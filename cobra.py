#!/usr/bin/env python
# Author: LinXing Chen, UC Berkeley

# cobra v1.2.3
# Contig Overlap Based Re-Assembly
# Modification date: Sep 3, 2023


import argparse
import concurrent.futures
import itertools
import os
from collections import defaultdict
from time import strftime
from typing import Callable, Iterable, Literal, NamedTuple, TextIO, TypeVar

import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import Seq

T = TypeVar("T")
try:
    from tqdm import tqdm
except ImportError:

    def tqdm(__object: Iterable[T], *nargs, **kwargs) -> Iterable[T]:
        return iter(__object)


try:
    from Bio.SeqUtils import gc_fraction as GC
except ImportError:
    from Bio.SeqUtils import GC  # type: ignore


def parse_args():
    parser = argparse.ArgumentParser(
        description="This script is used to get higher quality (including circular) virus genomes "
        "by joining assembled contigs based on their end overlaps."
    )
    requiredNamed = parser.add_argument_group("required named arguments")
    requiredNamed.add_argument(
        "-q",
        "--query",
        type=str,
        help="the query contigs file (fasta format), or the query name "
        "list (text file, one column).",
        required=True,
    )
    requiredNamed.add_argument(
        "-f",
        "--fasta",
        type=str,
        help="the whole set of assembled contigs (fasta format).",
        required=True,
    )
    requiredNamed.add_argument(
        "-a",
        "--assembler",
        type=str,
        choices=["idba", "megahit", "metaspades"],
        help="de novo assembler used, COBRA not tested for others.",
        required=True,
    )
    requiredNamed.add_argument(
        "-mink",
        "--mink",
        type=int,
        help="the min kmer size used in de novo assembly.",
        required=True,
    )
    requiredNamed.add_argument(
        "-maxk",
        "--maxk",
        type=int,
        help="the max kmer size used in de novo assembly.",
        required=True,
    )
    requiredNamed.add_argument(
        "-m",
        "--mapping",
        type=str,
        help="the reads mapping file in sam or bam format.",
        required=True,
    )
    requiredNamed.add_argument(
        "-c",
        "--coverage",
        type=str,
        help="the contig coverage file (two columns divided by tab).",
        required=True,
    )
    parser.add_argument(
        "-lm",
        "--linkage_mismatch",
        type=int,
        default=2,
        help="the max read mapping mismatches for "
        "determining if two contigs are spanned by "
        "paired reads. [2]",
    )
    parser.add_argument(
        "-tr",
        "--trim_readno",
        type=str,
        choices=["no", "trim", "auto"],
        default="no",
        help='whether trim the suffix `/1` and `/2` that distinguish paired reads or not ["no"]',
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="the name of output folder (default = '<query>_COBRA').",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=16,
        help="the number of threads for blastn. [16]",
    )
    parser.add_argument("-v", "--version", action="version", version="cobra v1.2.3")
    return parser.parse_args()


def get_log_info(steps: int, desc="COBRA Steps", log_file: TextIO | None = None):
    step = [0]
    step_format_str = "[{0:0>" + f"{len(str(steps))}" + "}/" + f"{steps}]"

    def steps_range():
        for i in tqdm(range(steps), desc=desc):
            if i:
                yield step_format_str.format(step[0])
                step[0] += 1
        while True:
            yield step_format_str.format(step[0])
            step[0] += 1

    stepin = steps_range()

    def _log_info(description: str, end="\n"):
        print(
            next(stepin),
            strftime("[%Y/%m/%d %H:%M:%S]"),
            description,
            end=end,
            file=log_file,
            flush=True,
        )

    print(desc, file=log_file, flush=True)
    return _log_info


# define functions for analyses
def compre_sets(p: int, li: set, lj: set):
    if p < 0:
        return ""
    if p == 0:
        return f"{len(li)}"
    i, j, ij = len(li), len(lj), len(li & lj)
    if i == j:
        return f"=="
    if i == ij:
        return ">"
    if j == ij:
        return "<"
    return f"{ij}"


def contig_name(end_name: str):
    """
    get contig name from end name
    """
    end_name_split = end_name.rsplit("_", 1)
    if end_name_split[-1] in ("L", "R", "Lrc", "Rrc"):
        return end_name_split[0]
    else:
        return end_name


def get_target(item: str):
    """
    to get the target for next run of joining
    """
    base_name, suffix = item.rsplit("_", 1)

    if suffix == "Lrc":
        return base_name + "_Rrc"
    elif suffix == "Rrc":
        return base_name + "_Lrc"
    elif suffix == "R":
        return base_name + "_L"
    else:
        return base_name + "_R"


def detect_self_circular(
    contig: str,
    one_path_end: set[str],
    two_paths_end: set[str],
    link_pair: dict[str, list[str]],
):
    """
    to determine the input sequences if they are self_circular genomes or not
    """
    end = contig + "_L"
    if end in one_path_end:
        if link_pair[end][0] == contig + "_R":
            # the other end could be joined with the current working end
            #           >contig
            #           |----->...|----->
            #           |- L ->
            #           |- R ->
            # |----->...|----->
            # >contig
            #           |----->.*.*.*
            #           >other contigs ends with `|----->` <
            #
            return "one_path_end"
    elif end in two_paths_end:
        if contig + "_R" in link_pair[end]:
            # the other end could be joined with the current working end
            #           >contig
            #           |----->...|----->
            #           |- L ->
            #           |- R ->
            # |----->...|----->
            # >contig
            # |----->...|----->
            # >the other contig ends with `|----->` <
            #
            if (
                link_pair[end][0].rsplit("_", 1)[0]
                == link_pair[end][1].rsplit("_", 1)[0]
            ):
                #           >contig
                #           |----->...|----->
                #           |- L ->
                #           |- R ->
                # |----->...|----->
                #                     contig< | FIXME: Why this happened?
                #           |~Rrc~>...|~Lrc~> |
                # |~Rrc~>...|~Lrc~>           |
                print(
                    "Warning: Unexpected self-multipile circular: ", contig, flush=True
                )
                print("       : from", end, "to", link_pair[end], flush=True)
                return ""
            return "two_paths_end"
    return ""


def check_self_circular_soft(sequence: Seq, min_over_len: int):
    """check overlap shorter than maxK, while not shorter than mink"""
    end_part = sequence[-min_over_len:]
    if sequence.count(end_part) > 1:
        expected_end = sequence.split(end_part, 1)[0] + end_part
        if sequence.endswith(expected_end):
            return len(expected_end)
    return -1


def get_contig2join(
    query_set: frozenset[str],
    orphan_end_query: frozenset[str],
    parsed_linkage: frozenset[tuple[str, str]],
    cov: dict[str, float],
    link_pair: dict[str, list[str]],
    one_path_end: frozenset[str],
    two_paths_end: frozenset[str],
    self_circular: frozenset[str],
):
    contig2join: dict[str, list[str]] = {}
    contig_checked: dict[str, list[str]] = {}
    contig2join_reason: dict[str, dict[str, str]] = {}
    for contig in query_set:
        contig2join[contig + "_L"] = []
        contig2join[contig + "_R"] = []
        contig_checked[contig + "_L"] = []
        contig_checked[contig + "_R"] = []
        contig2join_reason[contig] = {contig: "query"}

    path_circular_end: set[str] = set()

    def other_end_is_extendable(end: str):
        """True if the other end is extendable"""
        if contig_name(end) in self_circular:
            return False
        return len(link_pair[get_target(end)]) > 0

    def is_ok_to_add(end: str, contig: str):
        """
        if the other end of a potential joining contig cannot be extended,
        add only when it has a very similar coverage to that of the query contig
        """
        if contig_name(end) in self_circular:
            return False
        return (
            0.9 * cov[contig] <= cov[contig_name(end)]
            and cov[contig_name(end)] <= 1.11 * cov[contig]
        )

    def are_equal_paths(path1: str, path2: str):
        """
        True if two paths are equal

        .. code-block::
          >contig<
                 |-end-|
           .*.*.*|----->
          path1, path2      |           | exactly one contigs starts with |-R_2->    |
                 |- L ->... |           |-R_2->.*.*.*                                |
                 >contig2   | >contig2                                               |
                 |~Rrc~>*** | |-L_2->...|-R_2->                                      |
                   contig3< |             contig3<                                   |
                            | <-R_3-|   ...<-L_3-|                                   |
                            |              <-L_3-|*.*.*.                             |
                            |              | exactly one contigs starts with <-L_3-| |
          Expect <-L_3-| == reverse_complement(|-R_2->) inherent there are the same terminal of the same contig
          path1: |- L -> or |~Rrc~>
          path2: |- L -> or |~Rrc~>
          if |- L ->: just look for pair of |- R ->
          if |~Rrc~>: just look for pair of |~Lrc~>, equivalent to pair for |- L ->
        """
        path1_is_L = path1.rsplit("_", 1)[1].startswith("L")
        path2_is_L = path2.rsplit("_", 1)[1].startswith("L")
        path1_R_pair = link_pair[contig_name(path1) + ("_R" if path1_is_L else "_L")]
        path2_L_pair = link_pair[contig_name(path2) + ("_R" if path2_is_L else "_L")]
        return (
            len(path1_R_pair) == 1
            and len(path2_L_pair) == 1
            and contig_name(path1_R_pair[0]) == contig_name(path2_L_pair[0])
        )

    def the_dominant_one(path1: str, path2: str):
        """get the dominant path from two equal paths"""
        if cov[contig_name(path1)] >= cov[contig_name(path2)]:
            return path1
        else:
            return path2

    def the_better_one(paths: Iterable[str], contig: str):
        """Calculate the absolute differences in coverage"""
        return min(paths, key=lambda x: abs(cov[contig] - cov[contig_name(x)]))

    def could_circulate(
        point: str,
        contig: str,
        direction: Literal["L", "R"],
    ):
        """True if the path is circular with the current contig included"""
        contig_pair = ""
        other_direction = {"L": "_R", "R": "_L"}[direction]
        if len(link_pair[point]) == 2:  # point is the same as "target" in join_walker
            link_pair1, link_pair2 = link_pair[point]
            if contig_name(link_pair1) == contig:
                contig_pair = link_pair1
            elif contig_name(link_pair2) == contig:
                contig_pair = link_pair2
            if contig_pair and not cov[contig_name(point)] < 1.5 * cov[contig]:
                # 2 times is for repeat, but it is too risky, use 1.5 instead (same below)
                contig_pair = ""
        elif len(link_pair[point]) == 1:
            (link_pair1,) = link_pair[point]
            if contig_name(link_pair1) == contig:
                contig_pair = link_pair1
        return (
            contig_pair
            and contig_pair.endswith(other_direction)
            and (contig + other_direction, point) in parsed_linkage
        )

    def not_checked(end_list: Iterable[str], checked: list[str]):
        """True if all contig in end_list has been checked for adding or not"""
        return not any(contig_name(end) in checked for end in end_list)

    def join_walker(contig: str, direction: Literal["L", "R"]):
        """
        get potential joins for a given query
        """
        end = f"{contig}_{direction}"
        len_before_walk = len(contig2join[end])
        if len_before_walk == 0:
            contig_checked[end].append(contig)
            if end in one_path_end:
                link_pair1 = link_pair[end][0]
                if contig_name(link_pair1) != contig:
                    checked_reason = ""
                    # 1. link_pair1 is not point to a self-circulated contig
                    if other_end_is_extendable(link_pair1):
                        # 2.1. link_pair1 can be extend in the next step
                        checked_reason = "other_end_is_extendable"
                    elif is_ok_to_add(link_pair1, contig):
                        # 2.2. link_pair1 has similar coverage with contig
                        # however, it CANNOT be extend in the next step
                        checked_reason = "is_ok_to_add"
                    if checked_reason and (end, link_pair1) in parsed_linkage:
                        # 3. linkage between link_pair1 and end is supported by reads linkage
                        contig2join[end].append(link_pair1)
                        contig_checked[end].append(contig_name(link_pair1))
                        contig2join_reason[contig][
                            contig_name(link_pair1)
                        ] = checked_reason
            elif end in two_paths_end:
                # >contig<
                #        |-end-|
                #  .*.*.*|----->
                #        |----->.*.*.*
                #             >contig2<
                #        |----->.*.*.*
                #             >contig3<
                # no more contigs starts with `|----->`
                #
                link_pair1, link_pair2 = link_pair[end]
                if contig_name(link_pair1) != contig_name(link_pair2):
                    checked_reason = ""
                    if are_equal_paths(link_pair1, link_pair2):
                        #             >contig2<
                        #        |----->.*.   |=====>
                        # >contig<            |=====>.*.*.*
                        #  .*.*.*|----->            >contig4<
                        #        |-end-|
                        #        |----->   *.*|=====>
                        #             >contig3<
                        # no more contigs starts with `|----->` or ends with `|=====>`
                        #
                        link_pair_do = the_dominant_one(link_pair1, link_pair2)
                        checked_reason = "are_equal_paths"
                    elif (
                        cov[contig_name(link_pair1)] + cov[contig_name(link_pair2)]
                        >= cov[contig] * 0.5
                    ):
                        # 0.5 is ok, too big will get much fewer "Extended circular" ones.
                        # choise the one with more similar abundance
                        link_pair_do = the_better_one(link_pair[end], contig)
                        checked_reason = "the_better_one"
                    if (
                        checked_reason
                        and contig_name(link_pair_do) not in self_circular
                    ):
                        contig2join[end].append(link_pair_do)
                        contig_checked[end].append(contig_name(link_pair1))
                        contig_checked[end].append(contig_name(link_pair2))
                        contig2join_reason[contig][
                            contig_name(link_pair_do)
                        ] = checked_reason
                        # TODO: also record the discarded one
        else:
            target = get_target(contig2join[end][-1])
            if target in one_path_end:
                link_pair1 = link_pair[target][0]
                if checked := not_checked([link_pair1], contig_checked[end]):
                    # inherent contig_name(link_pair1) != contig
                    # the same as ablove
                    checked_reason = ""
                    if other_end_is_extendable(link_pair1):
                        checked_reason = "other_end_is_extendable"
                    elif is_ok_to_add(link_pair1, contig):
                        checked_reason = "is_ok_to_add"
                    if checked_reason and (target, link_pair1) in parsed_linkage:
                        contig2join[end].append(link_pair1)
                        contig_checked[end].append(contig_name(link_pair1))
                        contig2join_reason[contig][
                            contig_name(link_pair1)
                        ] = checked_reason
            elif target in two_paths_end:
                link_pair1, link_pair2 = link_pair[target]
                if checked := contig_name(link_pair1) != contig_name(link_pair2):
                    # >contig<
                    # .*.|----->
                    #    >contig of target<
                    #    |----->*.*.*.|----->
                    #                 |-----| target
                    #                 |----->.*.*.*
                    #                      >contig2<
                    #                 |----->.*.*.*
                    #                      >contig3<
                    # no more contigs starts with `|----->`
                    #
                    if cov[contig_name(target)] < 1.9 * cov[contig]:
                        # given a contig
                        #   |-query contig->[repeat region]|-contig2->[repeat region]|-contig3->
                        # where len(repeat region) > maxk, then the linkage looks like:
                        #   |-query contig->[repea>                           | 1* cov
                        #                   [repeat region]                   | 2* cov
                        #                           |egion]|-contig2->[repea> | 1* cov
                        #                           |egion]|-contig3->        | 1* cov
                        # where we can only confirm:
                        #   |-query contig->[repeat region]
                        # If we query the [repeat region], then never mind.
                        #
                        if not_checked([link_pair1, link_pair2], contig_checked[end]):
                            checked_reason = ""
                            if are_equal_paths(link_pair1, link_pair2):
                                link_pair_do = the_dominant_one(link_pair1, link_pair2)
                                checked_reason = "are_equal_paths"
                            elif (
                                cov[contig_name(link_pair1)]
                                + cov[contig_name(link_pair2)]
                                >= cov[contig] * 0.5
                            ):
                                # 0.5 is ok, too big will get much fewer "Extended circular" ones.
                                link_pair_do = the_better_one(link_pair[target], contig)
                                checked_reason = "the_better_one"
                            if (
                                checked_reason
                                and contig_name(link_pair_do) not in self_circular
                            ):
                                contig2join[end].append(link_pair_do)
                                contig_checked[end].append(contig_name(link_pair1))
                                contig_checked[end].append(contig_name(link_pair2))
                                contig2join_reason[contig][
                                    contig_name(link_pair_do)
                                ] = checked_reason
            else:
                checked = False
            if checked and could_circulate(target, contig, direction):
                # 1. target is extendable (the only link_pair here)
                path_circular_end.add(end)
                contig2join_reason[contig][contig_name(target)] = "could_circulate"
                # In the first walk, path will never add it self
                # In the following walks, paths will never be duplicated.
        return len_before_walk < len(contig2join[end])

    for contig in tqdm(
        query_set - (orphan_end_query | self_circular),
        desc="Detecting joins of contigs. ",
    ):
        # extend each contig from both directions
        while result_L := join_walker(contig, "L"):
            pass
        while result_R := join_walker(contig, "R"):
            pass
    return contig2join, contig2join_reason, path_circular_end


def get_contig2assembly(contig2join: dict[str, list[str]], path_circular_end: set[str]):
    contig2assembly: dict[str, set[str]] = {}
    for item in contig2join:
        if len(contig2join[item]) == 0:
            continue
        contig = contig_name(item)
        # detect extented query contigs
        if contig not in contig2assembly:
            contig2assembly[contig] = {contig}
        if (
            contig + "_L" in path_circular_end
            and contig + "_R" not in path_circular_end
        ):
            # here, only one path can be found to be circulated.
            contig2assembly[contig].update(
                (contig_name(i) for i in contig2join[contig + "_L"])
            )
        elif (
            contig + "_L" not in path_circular_end
            and contig + "_R" in path_circular_end
        ):
            contig2assembly[contig].update(
                (contig_name(i) for i in contig2join[contig + "_R"])
            )
        else:
            contig2assembly[contig].update((contig_name(i) for i in contig2join[item]))
    return contig2assembly


def join_seqs(
    contig: str,
    contig2join: dict[str, list[str]],
    path_circular_end: set[str],
):
    """
    get the join order of the sequences in a give path
    """
    left, right = contig + "_L", contig + "_R"

    order_all_ = [*reversed(contig2join[left]), contig]
    if left not in path_circular_end:
        # check if path is circular
        # then we can add right to it
        # TODO+FIXME: Compared with contig2assembly,
        #          1. what's the problem if path_circular_end in the middle of some path?
        #          2. what's the problem if some contig exist at both wind of a query?
        if contig2join[left] and contig2join[right]:
            # only when we add those in both contig2join[left] and contig2join[right]
            # and drop duplicates
            added_to_ = {contig_name(item) for item in order_all_}
            order_all_ += [
                item
                for item in contig2join[right]
                if contig_name(item) not in added_to_
            ]
        else:
            # either contig2join[left] or contig2join[right] is blank, just add one side is necessary
            order_all_ += contig2join[right]
    # elif left in path_circular:
    #   not extend right to order_all_ and just return
    return order_all_


def seqjoin2contig(seqjoin: Iterable[str]):
    used_queries = set()
    for item in seqjoin:
        contig = contig_name(item)
        if contig not in used_queries:
            yield item
            used_queries.add(contig)


def summary_fasta(
    fasta_file: str,
    length: int,
    cov: dict[str, float],
    self_circular: frozenset[str],
    self_circular_flexible_overlap: dict[str, int],
):
    """
    summary basic information of a fasta file
    """
    summary_file_headers = ["SeqID", "Length", "Coverage", "GC", "Ns"]
    if "self_circular" in fasta_file:
        summary_file_headers.append("DTR_length")
    with open(f"{fasta_file}.summary.tsv", "w") as so:
        print(*summary_file_headers, sep="\t", file=so)

        record: SeqIO.SeqRecord
        for record in SeqIO.parse(f"{fasta_file}", "fasta"):
            header = str(record.id).strip()
            seq: Seq = record.seq
            sequence_stats = [
                header,
                str(len(seq)),
                "",
                str(round(GC(seq), 3)),
                str(seq.count("N")),
            ]
            if "self_circular" in fasta_file:
                sequence_stats[2] = str(cov[header.split("_self")[0]])
                if header.split("_self")[0] in self_circular:
                    sequence_stats.append(str(length))
                else:
                    # header.split("_self")[0] in self_circular_flexible_overlap:
                    sequence_stats.append(
                        str(self_circular_flexible_overlap[header.split("_self")[0]])
                    )
            else:
                sequence_stats[2] = str(cov[header.split("_extended")[0]])
            print(*sequence_stats, sep="\t", file=so)


def get_direction(item: str):
    """
    get the direction in a joining path
    """
    return "reverse" if item.endswith("rc") else "forward"


def total_length(contig_list: list[str], header2len: dict[str, int]):
    """
    get the total length of all sequences in a joining path before overlap removing
    """
    total = sum(header2len[contig_name(item)] for item in contig_list)
    return total


def get_link_pair(assem_fa: str, maxk_length: int):
    header2seq: dict[str, Seq] = {}
    header2len: dict[str, int] = {}
    link_pair: dict[str, list[str]] = defaultdict(list)
    #
    d_L: dict[Seq, set[str]] = defaultdict(set)
    d_Lrc: dict[Seq, set[str]] = defaultdict(set)
    d_R: dict[Seq, set[str]] = defaultdict(set)
    d_Rrc: dict[Seq, set[str]] = defaultdict(set)
    #
    record: SeqIO.SeqRecord
    for record in tqdm(
        SeqIO.parse(assem_fa, "fasta"),
        desc="Reading contigs and getting the contig end sequences",
    ):
        header = str(record.id).strip()
        seq: Seq = record.seq
        header2seq[header] = seq
        header2len[header] = len(seq)
        # >contig
        # >>>>>>>...>>>>>>>
        # |- L ->   |- R ->
        # <~Lrc~|   <~Rrc~|
        # <<<<<<<***<<<<<<<
        #
        d_L[seq[:maxk_length]].add(header + "_L")
        d_Lrc[seq[:maxk_length].reverse_complement()].add(header + "_Lrc")
        d_R[seq[-maxk_length:]].add(header + "_R")
        d_Rrc[seq[-maxk_length:].reverse_complement()].add(header + "_Rrc")
    # get the shared seqs between direction pairs (L/Lrc, Lrc/L, L/R, R/L, R/Rrc, Rrc/R, Lrc/Rrc, Rrc/Lrc)
    #           >contig1
    #           >>>>>>>...>>>>>>>
    #           |- L ->
    #           |~Lrc~>
    # >>>>>>>***>>>>>>>
    #             contig2<
    # contig1_L - contig2_Lrc == contig1_Lrc - contig2_L
    d_L_d_Lrc_shared = set(d_L) & set(d_Lrc)
    #           >contig1
    #           >>>>>>>...>>>>>>>
    #           |- L ->
    #           |- R ->
    # >>>>>>>...>>>>>>>
    # >contig2
    # contig1_L - contig2_R == contig1_Lrc - contig2_Rrc
    d_L_d_R_shared = set(d_L) & set(d_R)
    # the d_R_d_L_shared will be included below
    # >contig1
    # >>>>>>>...>>>>>>>
    #           |- R ->
    #           |~Rrc~>
    #           >>>>>>>***>>>>>>>
    #                    contig2<
    # contig1_R - contig2_Rrc == contig1_Rrc - contig2_R
    d_R_d_Rrc_shared = set(d_R) & set(d_Rrc)
    # equals to                   |           >contig1
    #          contig1<           |           >>>>>>>...>>>>>>>
    # >>>>>>>***>>>>>>>           |           |- L ->
    #           |~Lrc~>           |           |- R ->
    #           |~Rrc~>           | >>>>>>>...>>>>>>>
    #           >>>>>>>***>>>>>>> | >contig2
    #                    contig2< |  the reverse_complement one
    d_Rrc_d_Lrc_shared = set(d_Rrc) & set(d_Lrc)
    ##
    # get link_pair between ends
    # link_pair = defaultdict(list)  # used to save all overlaps between ends
    for end in d_L_d_Lrc_shared:
        for left in d_L[end]:  # left is a seq name
            link_pair[left].extend(d_Lrc[end])
        for left_rc in d_Lrc[end]:
            link_pair[left_rc].extend(d_L[end])
    for end in d_L_d_R_shared:
        for left in d_L[end]:
            link_pair[left].extend(d_R[end])
        for right in d_R[end]:
            link_pair[right].extend(d_L[end])
    for end in d_R_d_Rrc_shared:
        for right in d_R[end]:
            link_pair[right].extend(d_Rrc[end])
        for right_rc in d_Rrc[end]:
            link_pair[right_rc].extend(d_R[end])
    for end in d_Rrc_d_Lrc_shared:
        for right_rc in d_Rrc[end]:
            link_pair[right_rc].extend(d_Lrc[end])
        for left_rc in d_Lrc[end]:
            link_pair[left_rc].extend(d_Rrc[end])
    return header2seq, header2len, link_pair


def check_Y_paths(
    link_pair: dict[str, list[str]],
    logfile: str,
):
    one_path_end: set[str] = set()  # the end of contigs with one potential join
    two_paths_end: set[str] = set()  # the end of contigs with two potential joins

    if logfile:
        p = open(logfile, "w")
    for item in tqdm(sorted(link_pair), desc="Collect one_path_end and two_path_end"):
        if len(link_pair[item]) == 1:  # and len(link_pair[link_pair[item][0]]) > 1:
            # | >contig3 or more contigs ends with `|----->` < |
            # |  .*.*.*|----->                                 |
            #   >contig1<
            #    .*.*.*|----->
            #          |----->.*.*.*
            #               >contig2<
            # |        no more contigs starts with `|----->` < |
            #
            # add one joining end to a list, its pair may have one or more joins
            one_path_end.add(item)
        elif (
            len(link_pair[item]) == 2
            and len(link_pair[link_pair[item][0]]) == 1
            and len(link_pair[link_pair[item][1]]) == 1
        ):
            # |          no more contigs ends with `|----->` < |
            #   >contig1<
            #    .*.*.*|----->
            #          |----->.*.*.*
            #               >contig2<
            #          |----->.*.*.*
            #               >contig3<
            # |        no more contigs starts with `|----->` < |
            #
            # add two joining end to a list, each of its pairs should only have one join
            two_paths_end.add(item)
        # print link pairs into a file for check if interested
        if logfile:
            print(item, *sorted(link_pair[item]), sep="\t", file=p)
    return one_path_end, two_paths_end


def get_cov(coverage_file: str):
    cov: dict[str, float] = {}  # the coverage of contigs
    with open(coverage_file) as coverage:
        # Sometimes will give a file with header, just ignore it once
        for line in coverage:
            header_cov, cov_value = line.strip().split("\t")[:2]
            try:
                cov[header_cov] = round(float(cov_value), 3)
            except ValueError:
                pass
            break
        for line in coverage:
            header_cov, cov_value = line.strip().split("\t")[:2]
            cov[header_cov] = round(float(cov_value), 3)
    return cov


def get_query_set(query_fa: str, uniset: dict[str, T]):
    query_set: set[str] = set()

    with open(query_fa) as f:
        first_line = next(f)
        if first_line.startswith(">"):
            # if the query file is in fasta format
            query_header = lambda x: (
                str(i.id).strip() for i in SeqIO.parse(x, "fasta")
            )
        else:
            # if the query file is in text format
            query_header = lambda x: (line.strip().split(" ")[0] for line in x)
    with open(query_fa) as query_file:
        for header in query_header(query_file):
            # some queries may not in the whole assembly, should be rare though.
            if header in uniset:
                query_set.add(header)
            else:
                print(
                    f"Query {header} is not in your whole contig fasta file, please check!",
                    flush=True,
                )
    return query_set


def get_parsed_linkage(
    mapping_file: str,
    trim_readno: Literal["no", "trim", "auto"],
    header2len: dict[str, int],
    orphan_end_query: frozenset[str],
    linkage_mismatch: int,
):
    """
    parsed_linkage: paired linkage supported by bam file
    parsed_contig_spanned_by_PE_reads: contig spanned by paired-end reads
    """
    # Initialize a defaultdict to store linked contigs
    parsed_linkage: set[tuple[str, str]] = set()
    linkage: dict[str, set[str]] = defaultdict(set)
    # Initialize a dictionary to store paired-end reads spanning contigs
    # Create a defaultdict(list) for each contig in header2seq, orphan_end first
    contig_spanned_by_PE_reads: dict[str, dict[str, list[int]]] = {
        contig: defaultdict(list) for contig in orphan_end_query
    }

    if trim_readno == "auto":
        with pysam.AlignmentFile(f"{mapping_file}", "rb") as map_file:
            for rmap in map_file:
                if rmap.query_name is not None:
                    if rmap.query_name[-2:] in ("/1", "/2"):
                        trim_readno = "trim"
                    else:
                        trim_readno = "no"
                    break
    parse_query_name: Callable[[str], str] = (
        (lambda x: x) if trim_readno == "no" else lambda x: x[:-2]
    )

    with pysam.AlignmentFile(f"{mapping_file}", "rb") as map_file:
        for rmap in tqdm(
            map_file,
            desc="Getting contig linkage based on sam/bam. Be patient, this may take long",
        ):
            if not rmap.is_unmapped and int(rmap.get_tag("NM")) <= linkage_mismatch:
                assert rmap.query_name is not None
                # mismatch should not be more than the defined threshold
                # Check if the read and its mate map to different contigs
                if rmap.reference_name != rmap.next_reference_name:
                    # get name just before use it to speed up
                    rmap_query_name = parse_query_name(rmap.query_name)
                    assert rmap.reference_name is not None
                    if header2len[rmap.reference_name] > 1000:
                        # determine if the read maps to the left or right end
                        # >contig1
                        # >>>>>>>>>>>>>>>>>>>>>>>>...>
                        # |           |           |
                        # ^- 0        ^- 500      ^- 1000
                        # +++++++ ... |              .       | linkage[read_i].add(contig1_L)
                        #         ... +++++++        .       | linkage[read_i].add(contig1_L)
                        #             |+++++++  ...  |       | linkage[read_i].add(contig1_R)
                        #                       ...  +++++++ | linkage[read_i].add(contig1_R)
                        # @read_i
                        #
                        linkage[rmap_query_name].add(
                            rmap.reference_name
                            + ("_L" if rmap.reference_start <= 500 else "_R")
                        )
                    else:
                        # >contig1
                        # >>>>>>>>>>>>>>>>>>>>...>
                        # |           |           |
                        # ^- 0        ^- 500      ^- 1000
                        # +++++++   ...           | linkage[read_i].add(contig1_L, contig1_R)
                        #           ...   +++++++ | linkage[read_i].add(contig1_L, contig1_R)
                        # @read_i
                        #
                        # add both the left and right ends to the linkage
                        linkage[rmap_query_name].add(rmap.reference_name + "_L")
                        linkage[rmap_query_name].add(rmap.reference_name + "_R")
                else:
                    # If the read and its mate map to the same contig, store the read mapped position (start)
                    if rmap.reference_name in orphan_end_query:
                        # >contig1, normally > 1000 bp
                        # >>>>>>>>>>>>>>......>>>>>>>>>>>>>
                        # |           | ...... |           |
                        # ^- 0        ^- 500   ^- -500     ^- -0
                        # +++++++ ... |        .             | contig_spanned_by_PE_reads[contig1][read_i].append(0)
                        #         ... +++++++  .             | contig_spanned_by_PE_reads[contig1][read_i].append(499)
                        #             |+++++++ |             |
                        #                      +++++++       |
                        #                       +++++++      | contig_spanned_by_PE_reads[contig1][read_i].append(-500)
                        # @read_i/1
                        #               +++++++
                        #               @read_i/2
                        #
                        # add both the left and right ends to the linkage
                        # only care about those in query
                        if (
                            rmap.reference_start <= 500
                            or header2len[rmap.reference_name] - rmap.reference_start
                            <= 500
                        ):
                            contig_spanned_by_PE_reads[rmap.reference_name][
                                parse_query_name(rmap.query_name)
                            ].append(rmap.reference_start)

    #
    for read in tqdm(linkage, desc="Parsing the linkage information"):
        # len(linkage[read]) in (1, 2, 3, 4)
        if len(linkage[read]) >= 2:  # Process only reads linked to at least two contigs
            for item, item_1 in itertools.combinations(linkage[read], 2):
                # Generate unique pairs of linked contigs for the current read using itertools.combinations
                if item.rsplit("_", 1)[1] != item_1.rsplit("_", 1)[1]:
                    # If the contigs have different ends (_L or _R), add the combinations to the parsed_linkage
                    parsed_linkage.add((item, item_1))
                    parsed_linkage.add((item_1, item))
                    parsed_linkage.add((item + "rc", item_1 + "rc"))
                    parsed_linkage.add((item_1 + "rc", item + "rc"))
                else:
                    # If the contigs have the same ends, add the combinations with reverse-complement (_rc) to the parsed_linkage
                    # Warning: contigs <= 1000 bp will be linked to it self: (contig1_L, contig1_Rrc) etc.
                    parsed_linkage.add((item, item_1 + "rc"))
                    parsed_linkage.add((item_1 + "rc", item))
                    parsed_linkage.add((item + "rc", item_1))
                    parsed_linkage.add((item_1, item + "rc"))

    # Initialize a set to store the contig spanned by paired-end reads
    parsed_contig_spanned_by_PE_reads: set[str] = set()
    for contig in orphan_end_query:
        # Check if the count is 0 and the contig has exactly two paired-end reads
        for PE in contig_spanned_by_PE_reads[contig]:
            if len(contig_spanned_by_PE_reads[contig][PE]) == 2:
                # Check if the absolute difference between the positions of the two paired-end reads is greater than or equal to
                # the length of contig minus 1000 bp
                if (
                    abs(
                        contig_spanned_by_PE_reads[contig][PE][0]
                        - contig_spanned_by_PE_reads[contig][PE][1]
                    )
                    >= header2len[contig] - 1000
                ):
                    # >contig1, normally > 1000 bp
                    # >>>>>>>>>>>>>>......>>>>>>>>>>>>>
                    # |           | ...... |           |
                    # ^- 0        ^- 500   ^- -500     ^- -0
                    # +++++++ ... |        |
                    #         ... +++++++  |
                    #   @read_i/1 .        |
                    #             .        +++++++
                    #             .        .     +++++++
                    #             .        .  @read_i/2
                    #             |--------| <- header2len[contig1] - 1000
                    #
                    parsed_contig_spanned_by_PE_reads.add(contig)

    return frozenset(parsed_linkage), frozenset(parsed_contig_spanned_by_PE_reads)


class QueryCovVar(NamedTuple):
    query_count: int
    cov_mean: float
    cov_std: float
    seq_len: float
    cov_adj_sum: float

    @classmethod
    def query(
        cls, contigs: Iterable[str], cov: dict[str, float], header2len: dict[str, int]
    ):
        """
        return stat of query coverage

        stat of queries:
        - query   count
        - cov     mean
        - cov     std
        - len     sum
        - len*cov sum
        """
        cov_lens = pd.Series({i: (cov[i], header2len[i]) for i in contigs})
        return cls(
            len(cov_lens),
            cov_lens.apply(lambda x: x[0]).mean(),
            cov_lens.apply(lambda x: x[0]).std(),
            cov_lens.apply(lambda x: x[1]).sum(),
            cov_lens.apply(lambda x: x[0] * x[1]).sum(),
        )


def extend_query(
    contigs: Iterable[str], header2seq: dict[str, Seq], maxk_length: int
) -> Seq:
    extend_seq = Seq("")
    last = Seq("")
    # print the sequences with their overlap removed
    for item in contigs:
        contig = contig_name(item)
        if item.endswith("rc"):
            seq = header2seq[contig].reverse_complement()
        else:
            # not a terminal, just the query itself
            if not (item.endswith("_R") or item.endswith("_L")):
                query = contig
                # assert query == contig
            seq = header2seq[contig]
        if last == "" or seq[:maxk_length] == last:
            extend_seq += seq[:-maxk_length]
            last = seq[-maxk_length:]
    return extend_seq + last


def write_and_run_blast(
    data_chunk: Iterable[tuple[str, Seq, int]],
    working_dir: str,
    prec_identity=70,
    retries=3,
):
    cobra_seq2lenblast: dict[str, tuple[int, list[list[str]]]] = {}
    os.makedirs(working_dir)
    outlog = f"{working_dir}/log"
    outfile = f"{working_dir}/blastdb_2.vs.1.tsv"
    for i in range(1, retries + 1):
        try:
            with (
                open(f"{working_dir}/blastdb_1.fa", "w") as blastdb_1,
                open(f"{working_dir}/blastdb_2.fa", "w") as blastdb_2,
            ):
                for header, seq, overlap in data_chunk:
                    real_seq_len = len(seq) - overlap
                    cobra_seq2lenblast[header] = real_seq_len, []
                    half = (len(seq) + 1) // 2
                    print(f">{header}_1\n{seq[:half]}", file=blastdb_1)
                    print(f">{header}_2\n{seq[half:real_seq_len]}", file=blastdb_2)
            os.system(f"ls -sh {blastdb_1.name} {blastdb_2.name} > {outlog} 2>&1")
            os.system(f"makeblastdb -in {blastdb_1.name} -dbtype nucl >> {outlog} 2>&1")
            os.system(
                "blastn"
                f" -task blastn"
                f" -db {blastdb_1.name}"
                f" -query {blastdb_2.name}"
                f" -out {outfile}"
                f" -outfmt 6 -evalue 1e-10 -perc_identity {prec_identity}"
                f" -num_threads 1 "
                f">> {outlog} 2>&1"
            )
            with open(outfile) as r:
                for line in r:
                    line_v = line.strip().split("\t")
                    contig_join = line_v[0].rsplit("_", 1)[0]
                    if (
                        contig_join == line_v[1].rsplit("_", 1)[0]
                        and line_v[0] != line_v[1]
                        and float(line_v[3]) >= 1000
                    ):
                        cobra_seq2lenblast[contig_join][1].append(line_v)
            return cobra_seq2lenblast
        except Exception:
            os.system(f"echo '' >> {outlog}")
            os.system(f"echo 'run {i} failed' >> {outlog}")
            os.system(f"echo '' >> {outlog}")
            if i == retries:
                raise
    raise NotImplementedError


def group_yield(data: Iterable[T], n=100):
    d = iter(data)
    while l := [i for i, _ in zip(d, range(n))]:
        yield l


def run_blast_half(
    query2compare: Iterable[tuple[str, Seq, int]], outfile: str, threads: int, n=100
):
    """
    Screening of contigs joined from closely related genomes using BLASTn comparison.
    https://www.nature.com/articles/s41564-023-01598-2/figures/8

    To prevent the joining of fragmented contigs from closely related (sub)populations,
      a BLASTn comparison is conducted
        between the first half and second half of each joined COBRA sequence.
    If the two parts share a region with a minimum length of 1000 bp
                                     and a minimum nucleotide similarity of 70%,
      all the query contigs involved in the join are labeled as "extended_failed"
        to indicate the failed extension.
    """
    with (
        concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor,
        open(outfile, "w") as results,
    ):
        os.makedirs(f"{outfile}_temp")
        futures = {
            executor.submit(write_and_run_blast, data_chunk, f"{outfile}_temp/tmp_{i}")
            for i, data_chunk in enumerate(group_yield(query2compare, n))
        }
        for future in tqdm(
            concurrent.futures.as_completed(futures),
            desc="Running blastn",
            total=len(futures),
        ):
            for query, (seq_len, blastns) in future.result().items():
                for blastn_v in blastns:
                    print(query, seq_len, *blastn_v[2:], sep="\t", file=results)
    return outfile


def get_failed_blast_half(outfile: str):
    """
    identify potential incorrect joins and remove them from corresponding category
    """
    contig_len_overlap: dict[str, float] = {}
    with open(outfile) as r:
        for line in r:
            line_v = line.strip().split("\t")
            if float(line_v[3]) >= 1000:
                if "_extended" in line_v[0]:
                    query = line_v[0].split("_extended")[0]
                else:
                    query = line_v[0]
                contig_len_overlap[query] = contig_len_overlap.get(query, 0) + float(
                    line_v[3]
                )

    return frozenset(
        contig for contig in contig_len_overlap if contig_len_overlap[contig] >= 1000
    )


def get_query2groups(contig2assembly: dict[str, set[str]], query_set: set[str]):
    """
    group all query paths if they shared paths

    .. code-block::
      k141_1: [k141_1, k141_2]                 |
      k141_2: [k141_1, k141_2, k141_3]         |
      k141_3: [k141_1, k141_2, k141_3]         | <- dup of k141_2 (query)
      k141_4: [k141_1, k141_2        , k141_4] |
      -----------------------------------------+
      [query]          ^- dup of k141_1 (contig in path)

    Yield: {
        k141_1: (frozenset(k141_1)                , frozenset()),
        k141_2: (frozenset(k141_1, k141_3)        , frozenset(k141_2, k141_3)),
        k141_4: (frozenset(k141_1        , k141_4), frozenset()),
    }
    """
    uassembly2contg = (
        # Series format of contig-assembly
        pd.Series(
            {k: frozenset(v) for k, v in contig2assembly.items()}, name="assembly"
        ).reset_index()
        # long format of contig-assembly
        .explode("assembly")
        # find contig (in path) exists in multipile query
        .groupby("assembly")["index"]
        # k141_1: [k141_1, k141_2]                 | [k141_1, k141_2, k141_4] :k141_1 |
        # k141_2: [k141_1, k141_2, k141_3]         | [k141_1, k141_2, k141_4] :       |
        #       : [k141_1, k141_2, k141_3]         | [        k141_2        ] :k141_3 |
        # k141_4: [k141_1, k141_2        , k141_4] | [                k141_4] :k141_4 |
        #
        .apply(frozenset)
    )
    query2groups_: dict[str, frozenset[str]] = {}
    queries: frozenset[str]
    # only call from query
    for query, queries in uassembly2contg.items():
        if query not in query_set:
            continue
        pan_queries: frozenset[str] = queries | {
            q for i in (queries & query2groups_.keys()) for q in query2groups_[i]
        }
        for query in pan_queries:
            query2groups_[query] = pan_queries
    # ucontig2assembly: dict[str, frozenset[str]] = (
    ucontig2assembly: dict[str, dict[Literal["index", "assembly"], frozenset[str]]] = (
        uassembly2contg.sort_index()
        .drop_duplicates()
        .reset_index()
        .explode("index")
        .groupby("index")["assembly"]
        .apply(frozenset)
        .reset_index()
        .groupby("assembly")["index"]
        .apply(lambda x: frozenset(x))
        .reset_index()
        .assign(Rep=lambda s: s["index"].apply(lambda x: next(iter(x))))
        .set_index("Rep")
        .T
    )
    for group in set(query2groups_.values()):
        yield {
            query: (v["assembly"], v["index"] if len(v["index"]) > 1 else frozenset())
            for query in group
            for v in [ucontig2assembly.get(query)]
            if v is not None
        }


def get_assembly2reason(
    query2groups: dict[
        int, dict[str, tuple[frozenset[str], frozenset[str] | frozenset]]
    ],
    contig2assembly: dict[str, set[str]],
    path_circular: frozenset[str],
    failed_join_potential: frozenset[str],
):
    """
    check query in each group, and report as a dict

    each item:
        query: (
            groupid,
            judgement,
            [sequences represented by it],
            frozenset(path with all the same query)
        )
    """
    check_assembly_reason: dict[
        str,
        tuple[
            int,
            Literal[
                "standalone",
                "conflict_query",
                "circular_in_sub",
                "circular_6_conflict",
                "circular_8_tight",
                "nolink_query",
                "longest",
            ],
            list[str],
            frozenset[str],
        ],
    ] = {}
    for groupi in query2groups:
        group = query2groups[groupi]
        # group: rep_query: {rep_query in the path}
        if len(group) == 1:
            (this,) = list(group)
            check_assembly_reason[this] = groupi, "standalone", [], group[this][1]
        else:
            subsets = sorted(group, key=lambda q: len(group[q][0]))
            this = subsets.pop(0)
            subset_trunks: dict[str, tuple[list[str], list[str]]] = {this: ([], [this])}
            subset_unextendable: dict[str, set[str]] = {}
            # assert all(i == v[-1] for i, v in subset_trunks.items())
            #        and values in v is sorted subset.
            while subsets:
                # start from the shortest path
                this = subsets.pop(0)
                this_feature: dict[
                    Literal[
                        "has_disjoint_found",
                        "is_subset_of",
                        "has_disjoint",
                        "sub_standalong",
                    ],
                    list[str],
                ] = {}
                for frag in subset_trunks | subset_unextendable.keys():
                    if group[frag][0] & group[this][0]:
                        if group[frag][0] < group[this][0]:
                            # k141_1: [k141_1]                 |
                            # k141_2: [k141_1, k141_3]         | *+ update longer
                            #
                            if frag in subset_unextendable:
                                this_feature.setdefault(
                                    "has_disjoint_found", []
                                ).append(frag)
                            else:
                                this_feature.setdefault("is_subset_of", []).append(frag)
                        elif group[frag][0] == group[this][0]:
                            assert False, "must be dereped previously"
                        elif frag in subset_trunks:
                            # k141_1: [k141_1]                 | keep
                            # k141_2: [k141_1, k141_3]         | - discard
                            # k141_4: [k141_1        , k141_4] | - discard
                            #
                            this_feature.setdefault("has_disjoint", []).append(frag)
                        elif frag in subset_unextendable:
                            this_feature.setdefault("has_disjoint_found", []).append(
                                frag
                            )
                        # else frag already has been captured by a subset_cannot_extend,
                        # or totally the same of something
                    else:
                        # k141_2: [k141_1, k141_3]         |
                        # k141_5: [              , k141_4] | + wait to create new item
                        #
                        this_feature.setdefault("sub_standalong", []).append(frag)
                if "has_disjoint" in this_feature:
                    check_assembly_reason[this] = (
                        groupi,
                        "conflict_query",
                        this_feature["has_disjoint"],
                        group[this][1],
                    )
                    for frag in this_feature["has_disjoint"] + this_feature.get(
                        "is_subset_of", []
                    ):
                        # k141_1: [k141_1, k141_3]                 |
                        # k141_3: [      , k141_3,       ]         | *- revert item
                        # k141_5: [      , k141_3, k141_4]         |
                        #
                        frags = subset_trunks[frag][1]
                        # find the not-conflict one
                        for fragi in range(len(frags)):
                            if (group[frags[fragi]][0] < group[this][0]) or not group[
                                frags[fragi]
                            ][0] & group[this][0]:
                                pass
                            else:
                                break
                        if fragi:
                            # common subset found
                            subset_trunks[frags[fragi - 1]] = (
                                subset_trunks.pop(frag)[0],
                                frags[:fragi],
                            )
                        else:
                            subset_trunks.pop(frag)
                            subset_unextendable[this] = {this}
                        # frozen it, or record in because this common subset cannot be bigger
                        subset_unextendable[frags[max(fragi - 1, 0)]] = set(
                            frags[fragi:]
                        )
                        check_assembly_reason[frags[-1]] = (
                            groupi,
                            f"conflict_query",
                            frags[fragi:],
                            group[frags[-1]][1],
                        )
                elif "has_disjoint_found" in this_feature:
                    check_assembly_reason[this] = (
                        groupi,
                        "conflict_query",
                        this_feature["has_disjoint_found"],
                        group[this][1],
                    )
                    for frag in this_feature["has_disjoint_found"]:
                        subset_unextendable[frag].add(this)
                elif "is_subset_of" in this_feature:
                    if len(frags := this_feature["is_subset_of"]) == 1:
                        # k141_1: [k141_1]                 |
                        # k141_6: [k141_1, k141_4]         | *+ update longer
                        #
                        subset_trunks[frags[0]][1].append(this)
                        subset_trunks[this] = subset_trunks.pop(frags[0])
                    else:
                        # k141_1: [k141_1]                 |
                        # k141_5: [      , k141_4]         |
                        # k141_6: [k141_1, k141_4]         | + select and create new item
                        #
                        standalong_subs = {x for x in frags for y in frags if x < y}
                        subset_trunks[this] = sorted(
                            i for i in frags if i not in standalong_subs
                        ), [this]
                else:  # if set(this_feature) == {"sub_standalong"}:
                    # k141_5: [k141_4]                 | + create new item
                    subset_trunks[this] = [], [this]
                # print(
                #    f"{this=}\n"
                #    f"{this_feature=}\n"
                #    f"{subset_trunks=}\n"
                #    f"{subset_unextendable=}\n"
                # )
            for subset in subset_trunks - subset_unextendable.keys():
                # subset          = [k141_1, k141_2, k141_3, k141_4]
                # subset_circular = [      , k141_2,       , k141_4]
                #  k141_1: 1 | not circular                              | * keep if k141_3 in k141_2
                #  k141_2: 0 | circular                                  | * keep if k141_3 not in k141_2
                #  k141_3: 6 | not circular but with a circular inside   |
                #  k141_4: 8 | circular and with another circular inside |
                #
                frags = list(subset_trunks[subset][1])
                if any(i for i in subset_trunks[subset][0] if i in path_circular):
                    # all the contig assemblies are based on another assembly,
                    # which is already circular and considered in another iter of this for-loop
                    check_assembly_reason[subset] = (
                        groupi,
                        "circular_in_sub",
                        frags,
                        group[subset][1],
                    )
                else:
                    subset_circular = [
                        i for i in subset_trunks[subset][1] if i in path_circular
                    ]
                    if subset_circular:
                        while (
                            subset_circular
                            and frags
                            # confirm potential `6` to speed up
                            and subset_circular[-1] != frags[-1]
                        ):
                            if (
                                contig2assembly[frags[-1]]
                                - contig2assembly[subset_circular[-1]]
                            ) & contig2assembly.keys():
                                # if contig < contig_1
                                # and contig_1 is not circular
                                # and contig_1 at the arm of `6`
                                #  >contig
                                # .|----->.
                                # ^       T
                                # |       |
                                # _       v
                                # .<-----|+|----->
                                #          >contig_1 (conflict and both failed)
                                #
                                # if so, remove the biggest `6`
                                frag = subset_circular[-1]
                                check_assembly_reason[frag] = (
                                    groupi,
                                    "circular_6_conflict",
                                    frags[frags.index(frag) :],
                                    group[frag][1],
                                )
                                # (next check for extensions on smaller circle)
                                frags = frags[: frags.index(frag)]
                            else:  # if subset_circular and subset_circular[0] != frags[-1]:
                                # check if any extension on the smallest `0`
                                # if so, any `6` will be considered as part of the biggest `8`
                                # and `8` is tightened as the smallest `0`
                                frag = frags[-1]
                                check_assembly_reason[frag] = (
                                    groupi,
                                    "circular_8_tight",
                                    frags[frags.index(frag) :],
                                    group[frag][1],
                                )
                                frags = frags[: frags.index(subset_circular[0]) + 1]
                                # assert subset_circular[0] != frags[-1], "expected break the while-loop"
                    if frags:
                        # remove query contigs that cannot extend itself
                        this = ""
                        for frag in frags:
                            if contig2assembly[frag] & failed_join_potential:
                                this = frag
                                break
                        if this:
                            check_assembly_reason[frags[-1]] = (
                                groupi,
                                "nolink_query",
                                frags[frags.index(frag) :],
                                group[frags[-1]][1],
                            )
                            frags = frags[: frags.index(frag)]
                        if frags:
                            check_assembly_reason[frags[-1]] = (
                                groupi,
                                "longest",
                                frags,
                                group[frags[-1]][1],
                            )
    return check_assembly_reason


def cobra(
    query_fa: str,
    assem_fa: str,
    mapping_file: str,
    coverage_file: str,
    maxk: int,
    mink: int,
    assembler: Literal["idba", "megahit", "metaspades"],
    trim_readno: Literal["no", "trim", "auto"],
    outdir: str = "",
    linkage_mismatch: int = 2,
    threads=1,
):
    ##
    # get information from the input files and parameters and save information
    # get the name of the whole contigs fasta file
    # folder of output
    if outdir:
        working_dir = f"{outdir}"
    else:
        # get the name of the query fasta file
        query_name = (
            f"{query_fa}".rsplit("/", 1)[1] if "/" in query_fa else f"{query_fa}"
        )
        working_dir = f"{query_name}_COBRA"

    # checking if output folder exists
    if os.path.exists(working_dir):
        raise FileExistsError(f"Output folder <{working_dir}> exists, please check.")
    else:
        os.mkdir(working_dir)

    # determine the length of overlap based on assembler and the largest kmer size
    maxk_length = maxk - 1 if assembler == "idba" else maxk

    # write input files information to log file
    log = open(f"{working_dir}/log", "w")  # log file
    print(
        "1. INPUT INFORMATION",
        "# Assembler: "
        + {"idba": "IDBA_UD", "metaspades": "metaSPAdes", "megahit": "MEGAHIT"}[
            assembler
        ],
        f"# Min-kmer: {mink}",
        f"# Max-kmer: {maxk}",
        f"# Overlap length: {maxk_length} bp",
        f"# Read mapping max mismatches for contig linkage: {linkage_mismatch}",
        f"# Query contigs: {os.path.abspath(query_fa)}",
        f"# Whole contig set: {os.path.abspath(assem_fa)}",
        f"# Mapping file: {os.path.abspath(mapping_file)}",
        f"# Coverage file: {os.path.abspath(coverage_file)}",
        f"# Output folder: {os.path.abspath(working_dir)}",
        sep="\n",
        file=log,
        flush=True,
    )

    print("\n", "2. PROCESSING STEPS", sep="\n", file=log, flush=True)
    ##
    # import the whole contigs and save their end sequences
    log_info_1 = get_log_info(5, "2.1. Loading assembly and mapping data", log)
    log_info_1("Reading contigs and getting the contig end pairs.", " ")
    header2seq, header2len, link_pair = get_link_pair(
        assem_fa=assem_fa, maxk_length=maxk_length
    )
    print(f"A total of {len(header2seq)} contigs were imported.", file=log, flush=True)

    ##
    # save all paired links to a file
    log_info_1("Joining contigs by contig end (maxK).", " ")
    one_path_end, two_paths_end = check_Y_paths(
        link_pair=link_pair, logfile=f"{working_dir}/COBRA_end_joining_pairs.tsv"
    )
    print(
        f"Among {len(link_pair)} link pairs, found"
        f" {len(one_path_end)} one path end,"
        f" {len(two_paths_end)} two paths end.",
        file=log,
        flush=True,
    )
    ##
    # read and save the coverage of all contigs
    log_info_1("Reading contig coverage information.")
    cov = get_cov(coverage_file)
    if len(cov) < len(header2seq):
        raise ValueError(
            "Some contigs do not have coverage information. Please check. COBRA exits."
        )

    ##
    # open the query file and save the information
    log_info_1("Getting query contig list.", " ")
    query_set = get_query_set(f"{query_fa}", uniset=header2seq)
    # distinguish orphan_end_query and non_orphan_end_query:
    orphan_end_query_1 = frozenset(
        header
        for header in query_set
        if header + "_L" not in link_pair.keys()
        and header + "_R" not in link_pair.keys()
    )
    print(
        f"A total of {len(query_set)} query contigs were imported, "
        f"with {len(orphan_end_query_1)} query with unique end (orphan).",
        file=log,
        flush=True,
    )

    ##
    # get the linkage of contigs based on paired-end reads mapping
    log_info_1("Getting linkage based on sam/bam. Be patient, this may take long.")
    parsed_linkage, parsed_contig_spanned_by_PE_reads = get_parsed_linkage(
        mapping_file,
        trim_readno=trim_readno,
        header2len=header2len,
        orphan_end_query=orphan_end_query_1,
        linkage_mismatch=linkage_mismatch,
    )

    ############################################################################
    ## Now all file are already loaded.                                       ##
    ############################################################################
    log_info_2 = get_log_info(
        8, "2.2. Analyzing assemblied paths and solving conflicts", log
    )
    debug = open(f"{working_dir}/debug.txt", "w")

    log_info_2("Detecting self_circular contigs, independent of mapping linkage.")
    self_circular = frozenset(
        contig
        for contig in tqdm(query_set, desc="Detecting self_circular contigs")
        if detect_self_circular(
            contig,
            one_path_end=one_path_end,
            two_paths_end=two_paths_end,
            link_pair=link_pair,
        )
    )

    # for debug
    print(f"# self_circular: {len(self_circular)}", file=debug, flush=True)
    print(sorted(self_circular), file=debug, flush=True)

    # orphan end queries info
    # determine potential self_circular contigs from contigs with orphan end
    min_over_len = mink - 1 if assembler == "idba" else mink
    # determine if there is DTR for those query with orphan ends, if yes, assign as self_circular as well
    self_circular_flexible_overlap = {
        contig: l
        for contig in parsed_contig_spanned_by_PE_reads
        if (l := check_self_circular_soft(header2seq[contig], min_over_len)) > 0
    }
    orphan_end_query = orphan_end_query_1 - self_circular_flexible_overlap.keys()

    # debug
    print(
        f"# self_circular_flexible_overlap: {len(self_circular_flexible_overlap)}",
        file=debug,
        flush=True,
    )
    print(sorted(self_circular_flexible_overlap), file=debug, flush=True)

    ##
    # walk the joins
    log_info_2("Detecting joins of contigs. ")
    contig2join, contig2join_reason, path_circular_end = get_contig2join(
        query_set=query_set,
        orphan_end_query=frozenset(orphan_end_query),
        parsed_linkage=parsed_linkage,
        cov=cov,
        link_pair=link_pair,
        one_path_end=one_path_end,
        two_paths_end=two_paths_end,
        self_circular=self_circular,
    )
    path_circular = frozenset(contig_name(i) for i in path_circular_end)
    print("# path_circular", file=debug, flush=True)
    print(sorted(path_circular), file=debug, flush=True)

    ##
    # save the potential joining paths
    log_info_2("Saving potential joining paths.")
    with open(f"{working_dir}/COBRA_potential_joining_paths.tsv", "w") as results:
        for item in sorted(contig2join):
            if contig_name(item) in self_circular:
                if item.endswith("_L"):
                    print(item, f"['{contig_name(item)}_R']", sep="\t", file=results)
                else:
                    print(item, f"['{contig_name(item)}_L']", sep="\t", file=results)
            elif contig2join[item]:
                print(item, contig2join[item], sep="\t", file=results)

    ##
    # get the fail_to_join contigs, but not due to orphan end
    # the simplest situation
    failed_join_potential = frozenset(
        contig
        for contig in query_set
        if (
            (contig + "_L" in link_pair or contig + "_R" in link_pair)
            and contig not in self_circular
            and len(contig2join[contig + "_L"]) == 0
            and len(contig2join[contig + "_R"]) == 0
        )
    )
    print("# failed_join_potential", file=debug, flush=True)
    print(sorted(failed_join_potential), file=debug, flush=True)

    ##
    # get the joining paths
    log_info_2("Getting the joining paths of contigs.")
    contig2assembly = get_contig2assembly(contig2join, path_circular_end)
    # for debug
    print("# contig2assembly", file=debug, flush=True)
    for k in sorted(contig2assembly):
        print(k, sorted(contig2assembly[k]), file=debug, flush=True)

    #
    # get the joining order of contigs, seems to be super of `contig2assembly`
    query2path = {
        query: join_seqs(
            query,
            contig2join=contig2join,
            path_circular_end=path_circular_end,
        )
        for query in contig2assembly
    }

    ##
    log_info_2("Getting joined seqeuences.")
    # query: seq, mark, len of overlap
    query2extension = {
        query: (
            extend_query(
                query2path[query], header2seq=header2seq, maxk_length=maxk_length
            ),
            "extended_circular" if query in path_circular else "extended_partial",
            maxk_length if query in path_circular else 0,
        )
        for query in tqdm(contig2assembly, desc="Extending path of query")
    }

    ##
    # Similar direct terminal repeats may lead to invalid joins
    log_info_2("Checking for invalid joining using BLASTn: close strains.")
    blast_half_file = run_blast_half(
        (
            i
            for j in (
                (
                    (contig, seq, overlap)
                    for contig, (seq, label, overlap) in query2extension.items()
                ),
                (
                    (contig, header2seq[contig], overlap)
                    for single_circle in (
                        self_circular,
                        self_circular_flexible_overlap,
                    )
                    for contig in single_circle
                    for overlap in (
                        self_circular_flexible_overlap.get(contig, maxk_length),
                    )
                ),
            )
            for i in j
        ),
        f"{working_dir}/blast_pairs.tsv",
        threads=threads,
        n=100,
    )
    # blast_half_file = f"{working_dir}/blast_pairs.tsv"
    failed_blast_half = get_failed_blast_half(blast_half_file)
    # for debug
    print("# failed_blast_half", file=debug, flush=True)
    print(sorted(failed_blast_half), file=debug, flush=True)
    debug.close()

    query2stat = {
        k: QueryCovVar.query(v, cov=cov, header2len=header2len)
        for k, v in tqdm(contig2assembly.items(), "stat query length and coverage")
    }

    log_info_2("Grouping paths by sharing queries to check for invalid queries.")
    query2groups = dict(
        enumerate(
            get_query2groups(contig2assembly=contig2assembly, query_set=query_set)
        )
    )
    check_assembly_reason = get_assembly2reason(
        query2groups=query2groups,
        contig2assembly=contig2assembly,
        path_circular=path_circular,
        failed_join_potential=failed_join_potential,
    )
    # for debug
    with open(f"{working_dir}/COBRA_check_assembly_reason.tsv", "w") as results:
        print("group", "reason", "query", "dups", sep="\t", file=results)
        for item in sorted(
            check_assembly_reason.items(), key=lambda x: (x[1][0], x[0])
        ):
            print(*item[1][:2], item[0], *item[1][2:], sep="\t", file=results)

    log_info_2("Filtering paths accoring to COBRA rules.")
    failed_joins: dict[Literal["conflict", "circular_6", "nolink"], dict[str, int]] = {
        l: {  # type: ignore [misc]
            j: v[0]
            for v in check_assembly_reason.values()
            if v[1] == v1
            for i in query2groups[v[0]]
            for j in contig2assembly[i]
            if j in query_set
        }
        for l, v1 in (
            ("conflict", "conflict_query"),
            ("circular_6", "circular_6_conflict"),
            ("nolink", "nolink_query"),
        )
    }
    failed_groups: dict[int, Literal["conflict", "circular_6", "nolink"]] = {
        i: q  # type: ignore [misc]
        for q in reversed(("conflict", "circular_6", "nolink"))
        for i in failed_joins[q].values()  # type: ignore [index]
    }
    checked_strict = {
        j: v[0]
        for k, v in check_assembly_reason.items()
        if v[1] in {"standalone", "longest"} and v[0] not in failed_groups
        for j in contig2assembly[k] & contig2assembly.keys()
    }
    failed_join_list = frozenset(i for f in failed_joins.values() for i in f.keys())

    assembly_rep = {
        max(v[3], key=lambda x: query2stat[x].seq_len) if len(v[3]) else i: v[0]
        for i, v in check_assembly_reason.items()
    }
    checked_strict_rep = {
        j: v[0]
        for k, v in check_assembly_reason.items()
        if v[1] in {"standalone", "longest"} and v[0] not in failed_groups
        for j in v[3] or {k}
        if j in assembly_rep and j not in failed_blast_half
    }

    redundant_circular_8 = {
        i: check_assembly_reason[i]
        for i in checked_strict.keys()
        if i in check_assembly_reason
        and check_assembly_reason[i][1] == "circular_8_tight"
    }
    print(
        *(
            "\t".join(
                [
                    *(compre_sets(pi - pj, li, lj) for pj, (j, lj) in enumerate(ll)),
                    f"{i}",
                ]
            )
            for ll in [
                [
                    ("query_set", query_set),
                    ("orphan_end_query", orphan_end_query),
                    (
                        "self_circular_flexible_overlap",
                        self_circular_flexible_overlap.keys(),
                    ),
                    ("self_circular", self_circular),
                    ("failed_join_potential", failed_join_potential),
                    ("failed_join_nolink", failed_joins["nolink"].keys()),
                    ("failed_join_conflict", failed_joins["conflict"].keys()),
                    ("failed_join_circular_6", failed_joins["circular_6"].keys()),
                    ("path_circular", path_circular),
                    ("checked_strict", checked_strict.keys()),
                    ("checked_strict_rep", checked_strict_rep.keys()),
                    ("assembly_rep", assembly_rep.keys()),
                    ("redundant_circular_8", redundant_circular_8.keys()),
                    ("failed_blast_half", failed_blast_half),
                ]
            ]
            for pi, (i, li) in enumerate(ll)
        ),
        sep="\n",
        file=log,
    )

    def check_check_assembly_reason():
        assert query_set - orphan_end_query_1 == (
            query_set | failed_join_potential | failed_join_list | checked_strict.keys()
        )
        assert not checked_strict.keys() & failed_join_potential
        assert not checked_strict.keys() & failed_joins["nolink"].keys()
        assert not checked_strict.keys() & failed_joins["conflict"].keys()
        assert not checked_strict.keys() & failed_joins["circular_6"].keys()
        assert not self_circular & self_circular_flexible_overlap.keys()
        assert not any(
            {
                query
                for query in checked_strict
                for contig in contig2assembly[query]
                if contig in failed_join_list
            }
        )

        def check_query(query: str):
            print(
                f"{query=}",
                f"{query in query_set=}",
                f"{contig2assembly.get(query)=}",
                f"{check_assembly_reason.get(query)=}",
                f"{query in path_circular=}",
                f"{query in failed_join_potential=}",
                f"{query in failed_joins['nolink']=}",
                f"{query in failed_joins['conflict']=}",
                f"{query in failed_joins['circular_6']=}",
                f"{query in checked_strict=}",
                sep="\n",
            )
            print(query2groups.get(check_assembly_reason.get(query, [-1])[0]))
            print({i: contig2assembly.get(i) for i in contig2assembly.get(query)})
            print({i: check_assembly_reason.get(i) for i in contig2assembly.get(query)})
            print(
                {
                    k: v
                    for k, v in query2groups.items()
                    if v.keys() & contig2assembly.get(query, {})
                }
            )

        check_query("k141_10804945")

    # make blastn database and run search if the database is not empty
    if (
        (not checked_strict_rep)
        and (not self_circular)
        and (not self_circular_flexible_overlap)
    ):
        # Of course we assume that sequence is not empty!
        print(
            "no query was extended by COBRA, exit! "
            "this is normal if you only provide few queries.",
            file=log,
            flush=True,
        )
        exit()

    ############################################################################
    ## Now output paths                                                       ##
    ############################################################################
    log_info_3 = get_log_info(5, "2.3. Output extended paths and circulated paths", log)
    ##
    # save the joining details information
    log_info_3(
        'Getting the joining details of unique "extended_circular" and "extended_partial" query contigs.'
    )
    with (
        open(
            f"{working_dir}/COBRA_category_ii-a_extended_circular_unique_joining_details.tsv",
            "w",
        ) as detail_circular,
        open(
            f"{working_dir}/COBRA_category_ii-b_extended_partial_unique_joining_details.tsv",
            "w",
        ) as detail_partial,
        open(f"{working_dir}/COBRA_joining_summary.tsv", "w") as assembly_summary,
    ):
        joining_detail_headers = [
            *("Final_Seq_ID", "Joined_Len", "Status"),
            *("Joined_Seq_ID", "Direction", "Joined_Seq_Len"),
            *("Start", "End", "Cov", "GC", "Joined_reason"),
        ]
        print(*joining_detail_headers, sep="\t", file=detail_circular)
        print(*joining_detail_headers, sep="\t", file=detail_partial)

        contig2join_details: dict[
            str,
            tuple[
                tuple[str, int, Literal["Circular", "Partial"]],
                list[
                    tuple[
                        str,
                        Literal["forward", "reverse"],
                        int,
                        int,
                        int,
                        float,
                        float,
                        str,
                    ]
                ],
            ],
        ] = {}
        for query in tqdm(
            checked_strict_rep,
            desc="Getting the joining details of extended query contigs.",
        ):
            site = 1
            query_extend_id = f"{query}_{query2extension[query][1]}"
            contig2join_details[query] = (
                query_extend_id,
                len(query2extension[query][0]) - query2extension[query][2],
                "Circular" if query in path_circular else "Partial",
            ), []
            for item in query2path[query]:
                contig = contig_name(item)
                if (direction := get_direction(item)) == "forward":
                    contig_start_end = (site, site + header2len[contig] - 1)
                else:
                    contig_start_end = (site, site + header2len[contig] - 1)
                contig2join_details[query][1].append(
                    (
                        contig,
                        direction,
                        header2len[contig],
                        *contig_start_end,
                        cov[contig],
                        round(GC(header2seq[contig]), 3),
                        contig2join_reason[query][contig],
                    )
                )
                site += header2len[contig] - maxk_length

        print(
            *("Query_Seq_ID", "Query_Seq_Len"),
            *("Final_Seq_ID", "Joined_Len", "Status", "Group_Id"),
            *("Total_Joined_Seqs", "Joined_seqs"),
            sep="\t",
            file=assembly_summary,
        )
        for query in sorted(contig2join_details):
            for detail_values in contig2join_details[query][1]:
                print(
                    *contig2join_details[query][0],
                    *detail_values,
                    sep="\t",
                    file=detail_circular if "circular" in query else detail_partial,
                )
            for contig in sorted(contig2assembly[query] & query_set):
                print(
                    *(contig, header2len[contig]),
                    *contig2join_details[query][0],
                    checked_strict_rep[query],
                    query2stat[query].query_count,
                    " ".join(seqjoin2contig(query2path[query])),
                    sep="\t",
                    file=assembly_summary,
                )

    log_info_3("Getting the joining details of failed query contigs.")
    with open(
        f"{working_dir}/COBRA_category_ii-c_extended_failed_paths.tsv", "w"
    ) as detail_failed:
        print(
            *("Group_Status", "Group_Id"),
            *("Failed_Seq_ID", "Failed_Status", "Failed_Len", "Failed_Path"),
            sep="\t",
            file=detail_failed,
        )
        reported_failed_groups: set[int] = set()
        for failed_reason in ("conflict", "circular_6", "nolink"):  # type: ignore[assignment]
            for groupi in sorted(
                set(failed_joins[failed_reason].values()) - reported_failed_groups  # type: ignore[index]
            ):
                for query in sorted(query2groups[groupi].keys()):
                    for contig in query2groups[groupi][query][1] or {query}:
                        print(
                            *(failed_reason, groupi, query),
                            check_assembly_reason.get(query, ("", "redundant"))[1],
                            len(query2extension[query][0]) - query2extension[query][2],
                            " ".join(query2path[query]),
                            sep="\t",
                            file=detail_failed,
                        )
                reported_failed_groups.add(groupi)
        for groupi, rep_query in sorted(
            (v[0], k)
            for k, v in check_assembly_reason.items()
            if v[1] in {"standalone", "longest"}
            and v[0] not in failed_groups.keys() | reported_failed_groups
            for j in v[3] or {k}
            if j in assembly_rep and j in failed_blast_half
        ):
            # now not failed by assembly, so just search "standalone" and "longest"
            for query in (
                contig
                for query in sorted(query2groups[groupi])
                for contig in sorted(query2groups[groupi][query][1] or {query})
            ):
                print(
                    *("blast_half", groupi, query),
                    check_assembly_reason[rep_query][1],
                    len(query2extension[query][0]) - query2extension[query][2],
                    " ".join(query2path[query]),
                    sep="\t",
                    file=detail_failed,
                )
            reported_failed_groups.add(groupi)
        for groupi, rep_query in sorted(
            (v[0], k)
            for k, v in redundant_circular_8.items()
            if v[0] not in failed_groups.keys() | reported_failed_groups
            and k in assembly_rep
        ):
            # now just search in path_circular
            # only report those not in the `longest`
            for query in sorted(
                contig2assembly[rep_query]
                & query_set
                - {
                    contig
                    for query in contig2assembly[rep_query] & checked_strict_rep.keys()
                    for contig in contig2assembly[query]
                }
            ):
                print(
                    *("circular_8", groupi, query),
                    check_assembly_reason[rep_query][1],
                    len(query2extension[query][0]) - query2extension[query][2],
                    " ".join(query2path[query]),
                    sep="\t",
                    file=detail_failed,
                )

    ##
    # save the joining summary information
    ##
    # save the joining status information of each query
    log_info_3("Saving joining status of all query contigs.")
    assembled_info = open(
        f"{working_dir}/COBRA_joining_status.tsv", "w"
    )  # shows the COBRA status of each query
    print(
        *("SeqID", "Length", "Coverage", "GC", "Group_Id", "Status", "Category"),
        sep="\t",
        file=assembled_info,
    )

    # for those could be extended to circular
    query_counts = {
        i: 0
        for i in (
            "extended_circular",
            "extended_partial",
            "extended_failed",
            "orphan_end",
            "self_circular",
        )
    }
    extended_fa_names = {}
    for ab, label in (("a", "extended_circular"), ("b", "extended_partial")):
        with open(
            f"{working_dir}/COBRA_category_ii-{ab}_{label}_unique.fa", "w"
        ) as extended_fasta:
            for query in sorted(checked_strict_rep):
                if query2extension[query][1] == label:
                    print(
                        f">{query}_{query2extension[query][1]}\n{query2extension[query][0]}",
                        file=(extended_fasta),
                        flush=True,
                    )
                    for contig in sorted(contig2assembly[query] & query_set):
                        query_counts[label] += 1
                        print(
                            contig,
                            header2len[contig],
                            cov[contig],
                            round(GC(header2seq[contig_name(contig)]), 3),
                            checked_strict_rep[query],
                            label,
                            f"category_ii-{ab}",
                            sep="\t",
                            file=assembled_info,
                        )
        extended_fa_names[label] = extended_fasta

    # for those cannot be extended
    with open(
        f"{working_dir}/COBRA_category_ii-c_extended_failed.fa", "w"
    ) as failed_join:
        label = "extended_failed"
        for query in sorted(failed_join_list):
            query_counts[label] += 1
            # assert query not in extended_circular_query
            # assert query not in extended_partial_query
            # assert query not in orphan_end_query
            print(f">{query}\n{header2seq[query]}", file=failed_join)
            print(
                query,
                header2len[query],
                cov[query],
                round(GC(header2seq[contig_name(query)]), 3),
                "",
                "extended_failed",
                "category_ii-c",
                sep="\t",
                file=assembled_info,
            )

    # for those due to orphan end
    with open(f"{working_dir}/COBRA_category_iii_orphan_end.fa", "w") as orphan_end:
        label = "orphan_end"
        for query in sorted(orphan_end_query):
            query_counts[label] += 1
            print(f">{query}\n{header2seq[query]}", file=orphan_end)
            print(
                query,
                header2len[query],
                cov[query],
                round(GC(header2seq[contig_name(query)]), 3),
                "",
                label,
                "category_iii",
                sep="\t",
                file=assembled_info,
            )

    # for self circular
    log_info_3("Saving self_circular contigs.")
    with open(
        f"{working_dir}/COBRA_category_i_self_circular.fa", "w"
    ) as circular_fasta:
        label = "self_circular"
        for query in sorted(self_circular):
            query_counts[label] += 1
            print(f">{query}_self_circular\n{header2seq[query]}", file=circular_fasta)
            print(
                query,
                header2len[query] - maxk_length,
                cov[query],
                round(GC(header2seq[contig_name(query)]), 3),
                "",
                label,
                "category_i",
                sep="\t",
                file=assembled_info,
            )
        for query in sorted(self_circular_flexible_overlap):
            query_counts[label] += 1
            print(f">{query}_self_circular\n{header2seq[query]}", file=circular_fasta)
            print(
                query,
                header2len[query] - self_circular_flexible_overlap[query],
                cov[query],
                round(GC(header2seq[contig_name(query)]), 3),
                "",
                label,
                "category_i",
                sep="\t",
                file=assembled_info,
            )
    assembled_info.close()

    for outio in (
        extended_fa_names["extended_circular"],
        extended_fa_names["extended_partial"],
        failed_join,
        orphan_end,
        circular_fasta,
    ):
        summary_fasta(
            outio.name,
            maxk_length,
            cov=cov,
            self_circular=self_circular,
            self_circular_flexible_overlap=self_circular_flexible_overlap,
        )

    ##
    # save new fasta file with all the others used in joining replaced by COBRA sequences excepting self_circular ones
    log_info_3("Saving the new fasta file.")
    query_extension_used = {
        contig: query
        for query in checked_strict_rep
        for contig in contig2assembly[query]
    }
    with open(f"{working_dir}/{assem_fa.rsplit('/', 1)[-1]}.new.fa", "w") as new:
        for header in sorted(header2seq.keys() - query_extension_used.keys()):
            print(f">{header}\n{header2seq[header]}", file=new)
    os.system(
        f"cat {extended_fa_names['extended_circular'].name} {extended_fa_names['extended_partial'].name} >> {new.name}"
    )
    ##
    # intermediate files

    ##
    # write the numbers to the log file
    print(
        "\n",
        "3. RESULTS SUMMARY",
        f"# Total queries: {len(query_set)}",
        f"# Category i   - self_circular: {query_counts['self_circular']}",
        f"# Category ii  - extended_circular: {query_counts['extended_circular']} (Unique: {len(checked_strict_rep.keys() & path_circular)})",
        f"# Category ii  - extended_partial: {query_counts['extended_partial']} (Unique: {len(checked_strict_rep.keys() - path_circular)})",
        f"# Category ii  - extended_failed (due to COBRA rules): {query_counts['extended_failed']}",
        f"# Category iii - orphan end: {query_counts['orphan_end']}",
        '# Check "COBRA_joining_status.tsv" for joining status of each query.',
        '# Check "COBRA_joining_summary.tsv" for joining details of "extended_circular" and "extended_partial" queries.',
        sep="\n",
        file=log,
        flush=True,
    )
    log.close()


def main():
    args = parse_args()
    cobra(
        query_fa=args.query,
        assem_fa=args.fasta,
        mapping_file=args.mapping,
        coverage_file=args.coverage,
        maxk=args.maxk,
        mink=args.mink,
        assembler=args.assembler,
        trim_readno=args.trim_readno,
        outdir=args.output,
        linkage_mismatch=args.linkage_mismatch,
        threads=args.threads,
    )


if __name__ == "__main__":
    main()
