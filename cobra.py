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
from pathlib import Path
import sys
from time import strftime
from typing import Callable, Iterable, KeysView, Literal, NamedTuple, TextIO, TypeVar

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


def parse_args(args=None):
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
        "--mapping-link-cache",
        type=str,
        nargs="*",
        help=" cache reads mapping to gzip info to speed up secondary run.",
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
        "--skip_joining",
        action="store_true",
        help="whether skip joining contigs output to disk",
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
    return parser.parse_args(args)


def get_log_info(steps: int, desc="COBRA Steps", log_file: TextIO | None = None):
    def steps_range():
        step = [0]
        for i in tqdm(range(steps), desc=desc):
            if i:
                yield step[0]
                step[0] += 1
        while True:
            yield step[0]
            step[0] += 1

    def _log_info(description: str, end="\n"):
        print(
            step_format_str.format(next(stepin)),
            strftime("[%Y/%m/%d %H:%M:%S]"),
            description,
            end=end,
            file=log_file,
            flush=True,
        )

    step_format_str = "[{0:0>" + f"{len(str(steps))}" + "}/" + f"{steps}]"
    stepin = steps_range()
    print(desc, file=log_file, flush=True)
    return _log_info


# define functions for analyses
def compre_sets(p: int, li: set, lj: set):
    if p < 0:
        return ""
    if p == 0:
        return f"{len(li)}"
    i, j, ij = len(li), len(lj), len(li & lj)
    if ij:
        if i == ij:
            if j == ij:
                return f"=="
            return ">"
        if j == ij:
            return "<"
    return f"{ij}"


def end2contig(end: str):
    """
    get contig name from end name
    """
    base_name, suffix = end.rsplit("_", 1)
    return base_name if suffix in ("L", "R", "Lrc", "Rrc") else end


def end2end2(end: str):
    """
    to get the target for next run of joining
    """
    base_name, suffix = end.rsplit("_", 1)

    if suffix == "Lrc":
        return base_name + "_Rrc"
    elif suffix == "Rrc":
        return base_name + "_Lrc"
    elif suffix == "R":
        return base_name + "_L"
    elif suffix == "L":
        return base_name + "_R"
    raise ValueError(f"Unexpected {suffix=}")


def detect_self_circular(
    contig: str,
    one_path_end: frozenset[str],
    two_paths_end: frozenset[str],
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
    orphan_query: frozenset[str],
    contig_pe_links: frozenset[tuple[str, str]],
    contig2cov: dict[str, float],
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
        return (
            end2contig(end) not in self_circular and len(link_pair[end2end2(end)]) > 0
        )

    def check_2path_to_add(end1: str, end2: str, contig: str):
        if link2_1contig(end1, end2):
            #             >contig2<
            #        |----->.*.   |=====>
            # >contig<            |=====>.*.*.*
            #  .*.*.*|----->            >contig4<
            #        |-end-|
            #        |----->   *.*|=====>
            #             >contig3<
            # no more contigs starts with `|----->` or ends with `|=====>`
            #
            checked_reason = "are_equal_paths"
            link_pair_do = the_dominant_one(end1, end2)
        elif (
            contig2cov[end2contig(end1)] + contig2cov[end2contig(end2)]
            >= contig2cov[contig] * 0.5
        ):
            # 0.5 is ok, too big will get much fewer "Extended circular" ones.
            # choise the one with more similar abundance
            checked_reason = "the_better_one"
            link_pair_do = the_better_one((end1, end2), contig)
        else:
            return "", ""
        if end2contig(link_pair_do) in self_circular:
            return "", ""
        return checked_reason, link_pair_do

    def check_1path_to_add(end: str, contig: str):
        """
        if the other end of a potential joining contig cannot be extended,
        add only when it has a very similar coverage to that of the query contig
        """
        # 2.1. link_pair1 can be extend in the next step
        if other_end_is_extendable(end):
            return "other_end_is_extendable"
        # however, it CANNOT be extend in the next step
        # 1. link_pair1 is not point to a self-circulated contig
        if end2contig(end) in self_circular:
            return ""
            # 2.2. link_pair1 has similar coverage with contig
        if (
            contig2cov[contig]  # don't div zero
            and 0.9 <= contig2cov[end2contig(end)] / contig2cov[contig] <= 1.11
        ):
            return "final_similar_cov"
        return ""

    def link2_1contig(end1: str, end2: str):
        """
        True if two paths both links to the only one same contig

        .. code-block::
          |-end1-| .*.*.*|-----> end1_2pair
          |-end2-| *.*.*.|-----> end2_2pair
                         |----->.*.*.*
                         * it means only one end matched by both pairs
        """
        end1_2pairs = link_pair[end2end2(end1)]
        end2_2pairs = link_pair[end2end2(end2)]
        return (
            len(end1_2pairs) == 1
            and len(end2_2pairs) == 1
            and end2contig(end1_2pairs[0]) == end2contig(end2_2pairs[0])
        )

    def the_dominant_one(path1: str, path2: str):
        """the paths with higher abundance"""
        return max(path1, path2, key=lambda x: contig2cov[end2contig(x)])

    def the_better_one(paths: Iterable[str], contig: str):
        """the path with more similar coverage"""
        return min(
            paths, key=lambda x: abs(contig2cov[contig] - contig2cov[end2contig(x)])
        )

    def could_circulate(
        end: str,
        contig: str,
        direction: Literal["L", "R"],
    ):
        """True if the path is circular with the current contig included"""
        contig_pair = ""
        other_direction = {"L": "_R", "R": "_L"}[direction]
        if end in two_paths_end:  # point is the same as "target" in join_walker
            link_pair1, link_pair2 = link_pair[end]
            if end2contig(link_pair1) == contig:
                contig_pair = link_pair1
            elif end2contig(link_pair2) == contig:
                contig_pair = link_pair2
            # 2 times is for repeat, but it is too risky, use 1.5 instead (same below)
            if contig_pair and 1.5 * contig2cov[contig] <= contig2cov[end2contig(end)]:
                contig_pair = ""
        elif end in one_path_end:
            (link_pair1,) = link_pair[end]
            if end2contig(link_pair1) == contig:
                contig_pair = link_pair1
        return (
            contig_pair
            and contig_pair.endswith(other_direction)
            and (contig + other_direction, end) in contig_pe_links
        )

    def not_checked(end_list: Iterable[str], checked: list[str]):
        """True if all contig in end_list has been checked for adding or not"""
        return not any(end2contig(end) in checked for end in end_list)

    def join_walker(contig: str, direction: Literal["L", "R"]):
        """
        get potential joins for a given query
        """
        end0 = f"{contig}_{direction}"
        len_before_walk = len(contig2join[end0])
        if len_before_walk == 0:
            contig_checked[end0].append(contig)
            if end0 in one_path_end:
                end1 = link_pair[end0][0]
                if end2contig(end1) != contig:
                    checked_reason = check_1path_to_add(end1, contig)
                    if checked_reason and (end0, end1) in contig_pe_links:
                        # 3. linkage between link_pair1 and end is supported by reads linkage
                        contig2join[end0].append(end1)
                        contig_checked[end0].append(end2contig(end1))
                        contig2join_reason[contig][end2contig(end1)] = checked_reason
            elif end0 in two_paths_end:
                # >contig<
                #        |-end-|
                #  .*.*.*|----->
                #        |----->.*.*.*
                #             >contig2<
                #        |----->.*.*.*
                #             >contig3<
                # no more contigs starts with `|----->`
                #
                end1, end2 = link_pair[end0]
                if end2contig(end1) != end2contig(end2):
                    checked_reason, end2add = check_2path_to_add(end1, end2, contig)
                    if checked_reason:
                        contig2join[end0].append(end2add)
                        contig_checked[end0].append(end2contig(end1))
                        contig_checked[end0].append(end2contig(end2))
                        contig2join_reason[contig][end2contig(end2add)] = checked_reason
                        # TODO: also record the discarded one
        else:
            end = end2end2(contig2join[end0][-1])
            if end in one_path_end:
                end1 = link_pair[end][0]
                if checked := not_checked([end1], contig_checked[end0]):
                    # inherent contig_name(link_pair1) != contig
                    # the same as ablove
                    checked_reason = check_1path_to_add(end1, contig)
                    if checked_reason and (end, end1) in contig_pe_links:
                        contig2join[end0].append(end1)
                        contig_checked[end0].append(end2contig(end1))
                        contig2join_reason[contig][end2contig(end1)] = checked_reason
            elif end in two_paths_end:
                end1, end2 = link_pair[end]
                if checked := end2contig(end1) != end2contig(end2):
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
                    if (
                        not_checked([end1, end2], contig_checked[end0])
                        and contig2cov[end2contig(end)] < 1.9 * contig2cov[contig]
                    ):
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
                        checked_reason, end2add = check_2path_to_add(end1, end2, contig)
                        if checked_reason:
                            contig2join[end0].append(end2add)
                            contig_checked[end0].append(end2contig(end1))
                            contig_checked[end0].append(end2contig(end2))
                            contig2join_reason[contig][
                                end2contig(end2add)
                            ] = checked_reason
            else:
                checked = False
            if checked and could_circulate(end, contig, direction):
                # 1. target is extendable (the only link_pair here)
                path_circular_end.add(end0)
                contig2join_reason[contig][end2contig(end)] = "could_circulate"
                # In the first walk, path will never add it self
                # In the following walks, paths will never be duplicated.
        return len_before_walk < len(contig2join[end0])

    for contig in tqdm(
        query_set - (orphan_query | self_circular),
        desc="Detecting joins of contigs. ",
    ):
        # extend each contig from both directions
        while join_walker(contig, "L"):
            pass
        while join_walker(contig, "R"):
            pass
    return contig2join, contig2join_reason, path_circular_end


def get_contig2assembly(contig2join: dict[str, list[str]], path_circular_end: set[str]):
    contig2assembly: dict[str, set[str]] = {}
    for item in contig2join:
        if len(contig2join[item]) == 0:
            continue
        contig = end2contig(item)
        # detect extented query contigs
        if contig not in contig2assembly:
            contig2assembly[contig] = {contig}
        if (
            contig + "_L" in path_circular_end
            and contig + "_R" not in path_circular_end
        ):
            # here, only one path can be found to be circulated.
            contig2assembly[contig].update(
                (end2contig(i) for i in contig2join[contig + "_L"])
            )
        elif (
            contig + "_L" not in path_circular_end
            and contig + "_R" in path_circular_end
        ):
            contig2assembly[contig].update(
                (end2contig(i) for i in contig2join[contig + "_R"])
            )
        else:
            contig2assembly[contig].update((end2contig(i) for i in contig2join[item]))
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
            added_to_ = {end2contig(item) for item in order_all_}
            order_all_ += [
                item for item in contig2join[right] if end2contig(item) not in added_to_
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
        contig = end2contig(item)
        if contig not in used_queries:
            yield item
            used_queries.add(contig)


def summary_fasta(
    fasta_file: str,
    length: int,
    contig2cov: dict[str, float],
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
                sequence_stats[2] = str(contig2cov[header.split("_self")[0]])
                if header.split("_self")[0] in self_circular:
                    sequence_stats.append(str(length))
                else:
                    # header.split("_self")[0] in self_circular_flexible_overlap:
                    sequence_stats.append(
                        str(self_circular_flexible_overlap[header.split("_self")[0]])
                    )
            else:
                sequence_stats[2] = str(contig2cov[header.split("_extended")[0]])
            print(*sequence_stats, sep="\t", file=so)


def get_direction(item: str):
    """
    get the direction in a joining path
    """
    return "reverse" if item.endswith("rc") else "forward"


def total_length(contig_list: list[str], contig2len: dict[str, int]):
    """
    get the total length of all sequences in a joining path before overlap removing
    """
    return sum(contig2len[end2contig(item)] for item in contig_list)


def get_link_pair(contig2seq: dict[str, Seq], maxk: int):
    d_L: dict[Seq, set[str]] = defaultdict(set)
    d_Lrc: dict[Seq, set[str]] = defaultdict(set)
    d_R: dict[Seq, set[str]] = defaultdict(set)
    d_Rrc: dict[Seq, set[str]] = defaultdict(set)
    #
    for header, seq in tqdm(
        contig2seq.items(), desc=f"Buiding graph based on {maxk}-mer end"
    ):
        # >contig
        # >>>>>>>...>>>>>>>
        # |- L ->   |- R ->
        # <~Lrc~|   <~Rrc~|
        # <<<<<<<***<<<<<<<
        d_L[seq[:maxk]].add(header + "_L")
        d_Lrc[seq[:maxk].reverse_complement()].add(header + "_Lrc")
        d_R[seq[-maxk:]].add(header + "_R")
        d_Rrc[seq[-maxk:].reverse_complement()].add(header + "_Rrc")
    # get the shared seqs between direction pairs (L/Lrc, Lrc/L, L/R, R/L, R/Rrc, Rrc/R, Lrc/Rrc, Rrc/Lrc)
    link_pair: dict[str, list[str]] = defaultdict(list)
    #           >contig1
    #           >>>>>>>...>>>>>>>
    #           |- L ->
    #           |~Lrc~>
    # >>>>>>>***>>>>>>>
    #       rc_contig2<
    # contig1_L - contig2_Lrc == contig1_Lrc - contig2_L
    for end in set(d_L) & set(d_Lrc):
        for left in d_L[end]:  # left is a seq name
            link_pair[left].extend(d_Lrc[end])
        for left_rc in d_Lrc[end]:
            link_pair[left_rc].extend(d_L[end])
    #           >contig1
    #           >>>>>>>...>>>>>>>
    #           |- L ->
    #           |- R ->
    # >>>>>>>...>>>>>>>
    # >contig2
    # contig1_L - contig2_R == contig1_Lrc - contig2_Rrc
    for end in set(d_L) & set(d_R):
        for left in d_L[end]:
            link_pair[left].extend(d_R[end])
        for right in d_R[end]:
            link_pair[right].extend(d_L[end])
    # the d_R_d_L_shared will be included below
    # >contig1
    # >>>>>>>...>>>>>>>
    #           |- R ->
    #           |~Rrc~>
    #           >>>>>>>***>>>>>>>
    #                 rc_contig2<
    # contig1_R - contig2_Rrc == contig1_Rrc - contig2_R
    for end in set(d_R) & set(d_Rrc):
        for right in d_R[end]:
            link_pair[right].extend(d_Rrc[end])
        for right_rc in d_Rrc[end]:
            link_pair[right_rc].extend(d_R[end])
    # equals to                   |           >contig1
    #       rc_contig1<           |           >>>>>>>...>>>>>>>
    # >>>>>>>***>>>>>>>           |           |- L ->
    #           |~Lrc~>           |           |- R ->
    #           |~Rrc~>           | >>>>>>>...>>>>>>>
    #           >>>>>>>***>>>>>>> | >contig2
    #                 rc_contig2< |  the reverse_complement one
    for end in set(d_Rrc) & set(d_Lrc):
        for right_rc in d_Rrc[end]:
            link_pair[right_rc].extend(d_Lrc[end])
        for left_rc in d_Lrc[end]:
            link_pair[left_rc].extend(d_Rrc[end])
    return link_pair


def check_Y_paths(
    link_pair: dict[str, list[str]],
    logfile: Path | None,
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
    return frozenset(one_path_end), frozenset(two_paths_end)


def get_cov(coverage_file: Path):
    contig2cov: dict[str, float] = {}  # the coverage of contigs
    with open(coverage_file) as coverage:
        # Sometimes will give a file with header, just ignore it once
        for line in coverage:
            header_cov, cov_value = line.strip().split("\t")[:2]
            try:
                contig2cov[header_cov] = float(cov_value)
            except ValueError:
                pass
            break
        for line in coverage:
            header_cov, cov_value = line.strip().split("\t")[:2]
            contig2cov[header_cov] = float(cov_value)
    return contig2cov


def get_query_set(query_fa: Path, uniset: dict[str, T]):
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
    return frozenset(query_set)


class MappingLinks:
    def __init__(self, mapping, mapping_link_cache):
        if mapping_link_cache is None:
            self.link_cache: Path | None = None
            self.mapping: Path | None = Path(mapping)
        elif len(mapping_link_cache) == 0:
            mapping_file = Path(mapping)
            if mapping_file.suffix == ".gz":
                self.link_cache = mapping_file
                self.mapping = None
            elif mapping_file.suffix in {".sam", ".bam"}:
                self.mapping = mapping_file
                self.link_cache = mapping_file.with_suffix(".link_cache.gz")
            else:
                self.mapping = mapping_file
                self.link_cache = mapping_file.with_suffix(
                    f"{mapping_file.suffix}.link_cache.gz"
                )
        elif len(mapping_link_cache) == 1:
            self.mapping = Path(mapping)
            self.link_cache = Path(mapping_link_cache[0])
        else:
            raise ValueError("Only one cache file is allowed.")

        if self.link_cache is not None and not self.link_cache.exists():
            self.link_cache.parent.mkdir(parents=True, exist_ok=True)

    def path(self):
        if self.cached():
            return self.link_cache
        return self.mapping

    def write_cache(self, contig_pe_links, orphan2pe_span):
        import gzip

        assert self.link_cache is not None
        with gzip.open(self.link_cache, "wb") as f:
            f.write(b">>>>>>>contig_pe_links\n")
            for i, j in contig_pe_links:
                f.write(f"{i}\t{j}\n".encode())
            f.write(b">>>>>>>orphan2pe_span\n")
            for i in orphan2pe_span:
                f.write(f"{i}\n".encode())

    def read_cache(self):
        import gzip

        assert self.link_cache is not None
        with gzip.open(self.link_cache, "rb") as f:
            contig_pe_links = set()
            orphan2pe_span = set()
            for line in f:
                if line.startswith(b">>>>>>>contig_pe_links"):
                    break
            for line in f:
                if line.startswith(b">>>>>>>orphan2pe_span"):
                    break
                contig_pe_links.add(tuple(line.decode().strip().split("\t")))
            for line in f:
                orphan2pe_span.add(line.decode().strip())
        return frozenset(contig_pe_links), frozenset(orphan2pe_span)

    def cached(self):
        return self.link_cache is not None and self.link_cache.exists()

    def get_link(
        self,
        trim_readno: Literal["no", "trim", "auto"],
        contig2len: dict[str, int],
        orphan_set: frozenset[str],
        linkage_mismatch: int = 2,
    ):
        """
        contig_pe_links: paired linkage supported by bam file
        orphan_pe_spanned:
            contig without kmer with other contigs, and the end may be spanned by paired-end reads
        """
        if self.cached():
            return self.read_cache()
        assert self.mapping is not None
        contig_pe_links, orphan_pe_spanned = self.read_link(
            self.mapping, trim_readno, contig2len, orphan_set, linkage_mismatch
        )
        if self.link_cache is not None:
            self.write_cache(contig_pe_links, orphan_pe_spanned)
        return contig_pe_links, orphan_pe_spanned

    @staticmethod
    def read_link(
        mapping_file: Path,
        trim_readno: Literal["no", "trim", "auto"],
        contig2len: dict[str, int],
        orphan_set: frozenset[str],
        linkage_mismatch: int = 2,
    ):
        """
        contig_pe_links: paired linkage supported by bam file
        orphan_pe_spanned:
            contig without kmer with other contigs, and the end may be spanned by paired-end reads
        """
        if trim_readno == "auto":
            with pysam.AlignmentFile(f"{mapping_file}", "rb") as map_file:
                for rmap in map_file:
                    if rmap.query_name is not None:
                        if rmap.query_name[-2:] in ("/1", "/2"):
                            trim_readno = "trim"
                        else:
                            trim_readno = "no"
                        break
        parse_readid: Callable[[str], str] = (
            (lambda x: x) if trim_readno == "no" else lambda x: x[:-2]
        )
        linkage: dict[str, set[str]] = defaultdict(set)
        orphan2pe_span: dict[str, dict[str, list[int]]] = {
            contig: defaultdict(list) for contig in orphan_set
        }
        with pysam.AlignmentFile(f"{mapping_file}", "rb") as map_file:
            for rmap in tqdm(
                map_file,
                desc="Getting contig linkage based on sam/bam. Be patient, this may take long",
            ):
                if rmap.is_unmapped:
                    continue
                assert rmap.query_name is not None
                assert rmap.reference_name is not None
                if linkage_mismatch < int(rmap.get_tag("NM")):
                    continue
                if rmap.reference_name != rmap.next_reference_name:
                    # Check if the read and its mate map to different contigs
                    # get name just before use it to speed up
                    rmap_readid = parse_readid(rmap.query_name)
                    if contig2len[rmap.reference_name] > 1000:
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
                        linkage[rmap_readid].add(
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
                        linkage[rmap_readid].add(rmap.reference_name + "_L")
                        linkage[rmap_readid].add(rmap.reference_name + "_R")
                else:
                    # If the read and its mate map to the same contig, store the read mapped position (start)
                    if rmap.reference_name in orphan_set:
                        # >contig1, normally > 1000 bp
                        # >>>>>>>>>>>>>>......>>>>>>>>>>>>>
                        # |           | ...... |           |
                        # ^- 0        ^- 500   ^- -500     ^- -0
                        # +++++++ ... |        .             | orphan2pe_span[contig1][read_i].append(0)
                        #         ... +++++++  .             | orphan2pe_span[contig1][read_i].append(499)
                        #             |+++++++ |             |
                        #                      +++++++       |
                        #                       +++++++      | orphan2pe_span[contig1][read_i].append(-500)
                        # @read_i/1
                        #               +++++++
                        #               @read_i/2
                        #
                        # add both the left and right ends to the linkage
                        # only care about those in query
                        if (
                            rmap.reference_start <= 500
                            or contig2len[rmap.reference_name] - rmap.reference_start
                            <= 500
                        ):
                            orphan2pe_span[rmap.reference_name][
                                parse_readid(rmap.query_name)
                            ].append(rmap.reference_start)

        contig_pe_links: set[tuple[str, str]] = set()
        for read in tqdm(linkage, desc="Parsing the linkage information"):
            # len(linkage[read]) in (1, 2, 3, 4)
            if (
                len(linkage[read]) >= 2
            ):  # Process only reads linked to at least two contigs
                for item, item_1 in itertools.combinations(linkage[read], 2):
                    # Generate unique pairs of linked contigs for the current read using itertools.combinations
                    if item.rsplit("_", 1)[1] != item_1.rsplit("_", 1)[1]:
                        # If the contigs have different ends (_L or _R), add the combinations to the parsed_linkage
                        contig_pe_links.add((item, item_1))
                        contig_pe_links.add((item_1, item))
                        contig_pe_links.add((item + "rc", item_1 + "rc"))
                        contig_pe_links.add((item_1 + "rc", item + "rc"))
                    else:
                        # If the contigs have the same ends, add the combinations with reverse-complement (_rc) to the parsed_linkage
                        # Warning: contigs <= 1000 bp will be linked to it self: (contig1_L, contig1_Rrc) etc.
                        contig_pe_links.add((item, item_1 + "rc"))
                        contig_pe_links.add((item_1 + "rc", item))
                        contig_pe_links.add((item + "rc", item_1))
                        contig_pe_links.add((item_1, item + "rc"))
        # Initialize a set to store the contig spanned by paired-end reads
        orphan_pe_spanned: set[str] = set()
        for contig in orphan_set:
            # Check if the count is 0 and the contig has exactly two paired-end reads
            for PE in orphan2pe_span[contig].values():
                if len(PE) == 2 and abs(PE[0] - PE[1]) >= contig2len[contig] - 1000:
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
                    #             |--------| <- contig2len[contig1] - 1000
                    #
                    orphan_pe_spanned.add(contig)
                    break

        return frozenset(contig_pe_links), frozenset(orphan_pe_spanned)


class QueryCovVar(NamedTuple):
    query_count: int
    cov_mean: float
    cov_std: float
    seq_len: float
    cov_adj_sum: float

    @classmethod
    def query(
        cls,
        contigs: Iterable[str],
        contig2: dict[str, float],
        contig2len: dict[str, int],
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
        cov_lens = pd.Series({i: (contig2[i], contig2len[i]) for i in contigs})
        return cls(
            len(cov_lens),
            cov_lens.apply(lambda x: x[0]).mean(),
            cov_lens.apply(lambda x: x[0]).std(),
            cov_lens.apply(lambda x: x[1]).sum(),
            cov_lens.apply(lambda x: x[0] * x[1]).sum(),
        )


def extend_query(contigs: Iterable[str], contig2seq: dict[str, Seq], maxk: int) -> Seq:
    extend_seq = Seq("")
    last = Seq("")
    # print the sequences with their overlap removed
    for end in contigs:
        contig = end2contig(end)
        if end.endswith("rc"):
            seq = contig2seq[contig].reverse_complement()
        else:
            seq = contig2seq[contig]
        if last == "" or seq[:maxk] == last:
            extend_seq += seq[:-maxk]
            last = seq[-maxk:]
    return extend_seq + last


def write_and_run_blast(
    data_chunk: Iterable[tuple[str, Seq, int]],
    working_dir: Path,
    prec_identity=70,
    retries=3,
):
    cobra_seq2lenblast: dict[str, tuple[int, list[list[str]]]] = {}
    os.makedirs(working_dir)
    outlog = working_dir / "log"
    outfile = working_dir / "blastdb_2.vs.1.tsv"
    for i in range(1, retries + 1):
        try:
            with (
                open(working_dir / "blastdb_1.fa", "w") as blastdb_1,
                open(working_dir / "blastdb_2.fa", "w") as blastdb_2,
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
    query2compare: Iterable[tuple[str, Seq, int]],
    outfile: Path,
    threads: int,
    chunk_size=100,
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
        temp_out = Path(f"{outfile}_temp")
        temp_out.mkdir()
        futures = {
            executor.submit(write_and_run_blast, data_chunk, temp_out / f"tmp_{i}")
            for i, data_chunk in enumerate(group_yield(query2compare, chunk_size))
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


class GroupAssemblyIndex(NamedTuple):
    special: frozenset[str]
    dup_queries: frozenset[str]


def query2groups(contig2assembly: dict[str, set[str]], query_set: frozenset[str]):
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
        # keys                     | special contigs                  | other query contigs
        k141_1: GroupAssemblyIndex(frozenset(k141_1)                , frozenset()),
        k141_2: GroupAssemblyIndex(frozenset(k141_1, k141_3)        , frozenset(k141_2, k141_3)),
        k141_4: GroupAssemblyIndex(frozenset(k141_1        , k141_4), frozenset()),
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
            query: GroupAssemblyIndex(
                v["assembly"], v["index"] if len(v["index"]) > 1 else frozenset()
            )
            for query in group
            for v in [ucontig2assembly.get(query)]
            if v is not None
        }


def group_this2feature(
    this: str,
    group: dict[str, GroupAssemblyIndex],
    subset_trunks: set[str] | KeysView[str],
    subset_unextendable: set[str] | KeysView[str],
):
    this_feature: dict[
        Literal["has_disjoint_found", "is_subset_of", "has_disjoint", "sub_standalong"],
        list[str],
    ] = {}
    ## first, compare this query to previous checked queries
    for frag in subset_trunks | subset_unextendable:
        if not group[frag].special & group[this].special:
            # there is no common sequence between this and frag, so they are disjoint
            # e.g.
            #    k141_2: [k141_1, k141_3]         |
            #    k141_5: [              , k141_4] | + wait to create new item
            this_feature.setdefault("sub_standalong", []).append(frag)
            # else we shouhd check their relationship
        elif frag in subset_unextendable:
            this_feature.setdefault("has_disjoint_found", []).append(frag)
        elif group[frag].special < group[this].special:
            # as sorted, group[frag].special never greater than group[this].special
            # e.g.
            #    k141_1: [k141_1]                 |
            #    k141_2: [k141_1, k141_3]         | *+ update longer
            this_feature.setdefault("is_subset_of", []).append(frag)
        else:
            assert (
                not group[frag].special >= group[this].special
            ), "must be uniq and smaller group first"
            # assert frag in subset_trunks
            # e.g.
            #    k141_1: [k141_1]                 | keep
            #    k141_2: [k141_1, k141_3]         | - discard
            #    k141_4: [k141_1        , k141_4] | - discard
            this_feature.setdefault("has_disjoint", []).append(frag)
    return this_feature


class AssemblyReason(NamedTuple):
    groupid: int
    judgement: Literal[
        "standalone",
        "conflict_query",
        "circular_in_sub",
        "circular_6_conflict",
        "circular_8_tight",
        "nolink_query",
        "longest",
    ]
    represent_seqs: list[str]
    dup_queries: frozenset[str]


class SubsetChunk(NamedTuple):
    standalong_subs: dict[str, "SubsetChunk"]
    frags: list[str]


def _get_subset_trunks(groupi: int, group: dict[str, GroupAssemblyIndex]):
    this, *subsets = sorted(group, key=lambda q: len(group[q].special))
    subset_trunks: dict[str, SubsetChunk] = {this: SubsetChunk({}, [this])}
    # the query in subset_trunks can not be extended,
    #  those contains query should be put in subset_unextendable[query]
    subset_unextendable: dict[str, set[str]] = {}
    check_assembly_reason: dict[str, AssemblyReason] = {}
    # assert all(i == v[-1] for i, v in subset_trunks.items())
    #        and values in v is sorted subset.
    while subsets:
        # start from the shortest path
        this = subsets.pop(0)
        this_feature = group_this2feature(
            this, group, subset_trunks.keys(), subset_unextendable.keys()
        )
        if "has_disjoint" in this_feature:
            # all things in this_feature["has_disjoint"], this_feature["is_subset_of"] and [this]
            #  cannot be extended, and should be reported as `conflict_query`
            # e.g.
            #    k141_1: [k141_1, k141_3]                 |
            #    k141_3: [      , k141_3,       ]         | *- revert item
            #    k141_5: [      , k141_3, k141_4]         |
            subset_unextendable[this] = {this}
            check_assembly_reason[this] = AssemblyReason(
                groupi,
                "conflict_query",
                this_feature["has_disjoint"],
                group[this].dup_queries,
            )
            disjoints = set(this_feature["has_disjoint"])
            while disjoints:
                frag = disjoints.pop()
                standalong_subs, frags = subset_trunks.pop(frag)
                # --- |       [  keep  ] if any
                # ----|--     [disjoint] frag
                #     |xxxxxx [disjoint] this
                # find the not-conflict one
                for fragi in range(len(frags)):
                    is_intersect = group[frags[fragi]].special & group[this].special
                    is_smaller = group[frags[fragi]].special < group[this].special
                    if is_intersect and not is_smaller:
                        break
                if fragi:
                    # Aha, as already assert len(frags) >= 1,
                    # here we assert len(frags) >= 2, and at least frags[0] is the common subset
                    # so we rsecue the common subset
                    subset_trunks[frags[fragi - 1]] = SubsetChunk(
                        standalong_subs, frags[:fragi]
                    )
                    subset_unextendable[frags[fragi - 1]] = set(frags[fragi:])
                elif standalong_subs:
                    # we should check those in .sdandalong_subs as well
                    for frag in standalong_subs:
                        if group[frag].special - group[this].special:
                            # disjoint happens here, dig into it
                            # if the common subset is not the first one, should keep it
                            disjoints.add(frag)
                        else:
                            # rsecue the common subset
                            subset_unextendable[frag] = {frag}
                        # will be moved in next loops if has disjoint
                        subset_trunks[frag] = standalong_subs[frag]
                else:
                    # all query in the extention failed
                    subset_unextendable[frag] = {frag}
                # frozen it, or record in because this common subset cannot be bigger
                check_assembly_reason[frags[-1]] = AssemblyReason(
                    groupi,
                    f"conflict_query",
                    frags[fragi:],
                    group[frags[-1]].dup_queries,
                )
        elif "has_disjoint_found" in this_feature:
            check_assembly_reason[this] = AssemblyReason(
                groupi,
                "conflict_query",
                this_feature["has_disjoint_found"],
                group[this].dup_queries,
            )
            for frag in this_feature["has_disjoint_found"]:
                subset_unextendable[frag].add(this)
        elif "is_subset_of" in this_feature:
            if len(frags := this_feature["is_subset_of"]) == 1:
                # e.g.
                #    k141_1: [k141_1]                 |
                #    k141_6: [k141_1, k141_4]         | *+ update longer
                subset_trunks[this] = subset_trunks.pop(frags[0])
                subset_trunks[this].frags.append(this)
            else:
                # e.g.
                #    k141_1: [k141_1]                 |
                #    k141_5: [      , k141_4]         |
                #    k141_6: [k141_1, k141_4]         | + select and create new item
                subset_trunks[this] = SubsetChunk(
                    {x: subset_trunks.pop(x) for x in frags}, [this]
                )
        else:  # if set(this_feature) == {"sub_standalong"}:
            # k141_5: [k141_4]                 | + create new item
            subset_trunks[this] = SubsetChunk({}, [this])
        # print(
        #    f"{this=}\n"
        #    f"{this_feature=}\n"
        #    f"{subset_trunks=}\n"
        #    f"{subset_unextendable=}\n"
        # )
    return (
        subset_trunks,
        {i for k in subset_unextendable.values() for i in k},
        check_assembly_reason,
    )


def get_assembly2reason(
    groupi: int,
    group: dict[str, GroupAssemblyIndex],
    contig2assembly: dict[str, set[str]],
    path_circular_potential: frozenset[str],
    contig_link_no_pe: frozenset[str],
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
    if len(group) == 1:
        for this in group:
            return {
                this: AssemblyReason(groupi, "standalone", [], group[this].dup_queries)
            }
    check_assembly_reason: dict[str, AssemblyReason] = {}
    subset_trunks, _unextendable, _failed_reason = _get_subset_trunks(groupi, group)
    check_assembly_reason.update(_failed_reason)
    subsets = subset_trunks.keys() - _unextendable
    while subsets:
        subset = subsets.pop()
        frags = subset_trunks[subset].frags[:]
        # subset          = [k141_1, k141_2, k141_3, k141_4]
        # subset_circular = [      , k141_2,       , k141_4]
        #  k141_1: 1 | not circular                              | * keep if k141_3 in k141_2
        #  k141_2: 0 | circular                                  | * keep if k141_3 not in k141_2
        #  k141_3: 6 | not circular but with a circular inside   |
        #  k141_4: 8 | circular and with another circular inside |
        #
        if any(
            i
            for i in subset_trunks[subset].standalong_subs
            if i in path_circular_potential
        ):
            # all the contig assemblies are based on another assembly,
            # which is already circular and considered in another iter of this for-loop
            check_assembly_reason[subset] = AssemblyReason(
                groupi,
                "circular_in_sub",
                frags,
                group[subset].dup_queries,
            )
            subsets.update(subset_trunks[subset].standalong_subs)
            subset_trunks.update(subset_trunks[subset].standalong_subs)
        else:
            subset_circular = [i for i in frags if i in path_circular_potential]
            while (
                subset_circular
                and frags
                # confirm potential `6` to speed up
                and subset_circular[0] != frags[-1]
            ):
                arm_contigs = (
                    contig2assembly[frags[-1]] - contig2assembly[subset_circular[-1]]
                )
                # if arm_contigs & contig_link_no_pe:
                #    # just drop and ignore it?
                #    frag = frags[-1]
                #    check_assembly_reason[frag] = AssemblyReason(
                #        groupi,
                #        "nolink_query",
                #        frags[frags.index(frag) :],
                #        group[frag].dup_queries,
                #    )
                #    frags = frags[: frags.index(frag)]
                # el...
                if (
                    arm_contigs & contig2assembly.keys()
                    or arm_contigs & contig_link_no_pe
                ):
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
                    check_assembly_reason[frag] = AssemblyReason(
                        groupi,
                        "circular_6_conflict",
                        frags[frags.index(frag) :],
                        group[frag].dup_queries,
                    )
                    # (next check for extensions on smaller circle)
                    frags = frags[: frags.index(frag)]
                    subset_circular = [i for i in frags if i in path_circular_potential]
                else:  # if subset_circular and subset_circular[0] != frags[-1]:
                    # check if any extension on the smallest `0`
                    # if so, any `6` will be considered as part of the biggest `8`
                    # and `8` is tightened as the smallest `0`
                    # it seems never be touched as subset_circular already handled as standalong_subs
                    frag = frags[-1]
                    check_assembly_reason[frag] = AssemblyReason(
                        groupi,
                        "circular_8_tight",
                        frags[frags.index(frag) :],
                        group[frag].dup_queries,
                    )
                    frags = frags[: frags.index(subset_circular[0]) + 1]
                    # assert subset_circular[0] != frags[-1], "expected break the while-loop"
                    # must end the while-loop
            if frags:
                # remove query contigs that cannot extend itself
                this = ""
                for frag in frags:
                    if contig2assembly[frag] & contig_link_no_pe:
                        this = frag
                        break
                if this:
                    check_assembly_reason[frags[-1]] = AssemblyReason(
                        groupi,
                        "nolink_query",
                        frags[frags.index(frag) :],
                        group[frags[-1]].dup_queries,
                    )
                    frags = frags[: frags.index(frag)]
                if frags:
                    check_assembly_reason[frags[-1]] = AssemblyReason(
                        groupi,
                        "longest",
                        frags,
                        group[frags[-1]].dup_queries,
                    )
    return check_assembly_reason


def get_assembly2reason2(
    groupi: int,
    group: dict[str, GroupAssemblyIndex],
    contig2assembly: dict[str, set[str]],
    path_circular_potential: frozenset[str],
    contig_link_no_pe: frozenset[str],
):
    if len(group) == 1:
        for this in group:
            return {
                this: AssemblyReason(groupi, "standalone", [], group[this].dup_queries)
            }
    largest = max(group, key=lambda q: len(group[q].special))
    if all(i.special <= group[largest].special for i in group.values()):
        return {
            largest: AssemblyReason(
                groupi, "longest", [largest], group[largest].dup_queries
            )
        }


def cobra(
    query_fa: Path,
    assem_fa: Path,
    mapping_links: MappingLinks,
    coverage_file: Path,
    maxk: int,
    mink: int,
    trim_readno: Literal["no", "trim", "auto"],
    working_dir: Path = Path(""),
    linkage_mismatch: int = 2,
    threads=1,
    logfile=sys.stdout,
    debugfile=sys.stderr,
    skip_joining=False,
):
    _log_info = lambda *n, **k: print(*n, **k, file=logfile, flush=True)
    _log_info(
        "1. INPUT INFORMATION",
        f"# Max-kmer: {maxk}",
        f"# Min-kmer: {mink}",
        f"# Overlap length: {maxk} bp",
        f"# Read mapping max mismatches for contig linkage: {linkage_mismatch}",
        f"# Query contigs: {query_fa.expanduser().absolute().resolve()}",
        f"# Whole contig set: {assem_fa.expanduser().absolute().resolve()}",
        f"# Mapping file: {mapping_links.path().expanduser().absolute().resolve()}",
        f"# Coverage file: {coverage_file.expanduser().absolute().resolve()}",
        f"# Output folder: {working_dir.expanduser().absolute().resolve()}",
        sep="\n",
    )
    _log_info("\n", "2. PROCESSING STEPS", sep="\n")
    # import the whole contigs and save their end sequences
    _log_info_1 = get_log_info(5, "2.1. Loading assembly and mapping data", logfile)
    _log_info_1("Reading contigs and getting the contig end pairs.", " ")
    contig2seq: dict[str, Seq] = {
        rec.id: rec.seq
        for rec in tqdm(SeqIO.parse(assem_fa, "fasta"), desc="Reading contigs")
    }
    _log_info(f"A total of {len(contig2seq)} contigs were imported.")

    # cache for faster
    contig2len = {i: len(seq) for i, seq in contig2seq.items()}
    link_pair = get_link_pair(contig2seq=contig2seq, maxk=maxk)
    # save all paired links to a file
    _log_info_1("Joining contigs by contig end (maxK).", " ")
    one_path_end, two_paths_end = check_Y_paths(
        link_pair=link_pair, logfile=working_dir / f"COBRA_end_joining_pairs.tsv"
    )
    _log_info(
        f"Among {len(link_pair)} link pairs, found"
        f" {len(one_path_end)} one path end,"
        f" {len(two_paths_end)} two paths end.",
    )
    orphan_pre_set = frozenset(
        header
        for header in contig2seq
        if f"{header}_L" not in link_pair and f"{header}_R" not in link_pair
    )
    # get the linkage of contigs based on paired-end reads mapping
    _log_info_1("Getting linkage based on sam/bam. Be patient, this may take long.")
    contig_pe_links, orphan_pe_spanned = mapping_links.get_link(
        trim_readno=trim_readno,
        contig2len=contig2len,
        orphan_set=orphan_pre_set,
        linkage_mismatch=linkage_mismatch,
    )

    # read and save the coverage of all contigs
    _log_info_1("Reading contig coverage information.")
    contig2cov = get_cov(coverage_file)
    assert not contig2seq.keys() - contig2cov.keys(), "Contigs missing coverage!"

    # open the query file and save the information
    _log_info_1("Getting query contig list.", " ")
    query_set = get_query_set(query_fa, uniset=contig2seq)
    # distinguish orphan_end_query and non_orphan_end_query:
    _log_info(
        f"A total of {len(query_set)} query contigs were imported, "
        f"with {len(orphan_pre_set & query_set)} query with unique end (orphan).",
    )
    ############################################################################
    ## Now all file are already loaded.                                       ##
    ############################################################################
    _log_info_2 = get_log_info(
        8, "2.2. Analyzing assemblied paths and solving conflicts", logfile
    )

    _log_info_2("Detecting self_circular contigs, independent of mapping linkage.")
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
    _log_debug = lambda *n, **k: print(*n, **k, file=debugfile, flush=True)
    _log_debug(f"# self_circular: {len(self_circular)}")
    _log_debug(sorted(self_circular))

    # orphan end queries info
    # determine if there is DTR for those query with orphan ends, if yes, assign as self_circular as well
    self_circular_flex = {
        contig: l
        for contig in orphan_pe_spanned
        if contig in query_set
        and (l := check_self_circular_soft(contig2seq[contig], mink)) > 0
    }
    orphan_query = orphan_pre_set - self_circular_flex.keys() & query_set
    _log_debug(f"# self_circular_flexible_overlap: {len(self_circular_flex)}")
    _log_debug(sorted(self_circular_flex))

    # walk the joins
    _log_info_2("Detecting joins of contigs. ")
    contig2join, contig2join_reason, path_circular_end = get_contig2join(
        query_set=query_set,
        orphan_query=orphan_query,
        contig_pe_links=contig_pe_links,
        contig2cov=contig2cov,
        link_pair=link_pair,
        one_path_end=one_path_end,
        two_paths_end=two_paths_end,
        self_circular=self_circular,
    )
    path_circular_potential = frozenset(end2contig(i) for i in path_circular_end)
    _log_debug("# path_circular")
    _log_debug(sorted(path_circular_potential))

    # save the potential joining paths
    _log_info_2("Saving potential joining paths.")
    with open(working_dir / f"COBRA_potential_joining_paths.tsv", "w") as results:
        for item in sorted(contig2join):
            if end2contig(item) in self_circular:
                print(item, [end2end2(item)], sep="\t", file=results)
            elif contig2join[item]:
                print(item, contig2join[item], sep="\t", file=results)

    # get the joining paths
    _log_info_2("Getting the joining paths of contigs.")
    contig2assembly = get_contig2assembly(contig2join, path_circular_end)
    # for debug
    _log_debug("# contig2assembly")
    for k in sorted(contig2assembly):
        _log_debug(k, sorted(contig2assembly[k]))

    # overlap found, but none end supported by reads linkage
    contig_link_no_pe = frozenset(
        contig
        for contig in query_set
        if (
            (contig + "_L" in link_pair or contig + "_R" in link_pair)
            and contig not in self_circular
            and len(contig2join[contig + "_L"]) == 0
            and len(contig2join[contig + "_R"]) == 0
        )
    )
    _log_debug("# contig_link_no_pe")
    _log_debug(sorted(contig_link_no_pe))

    # get the joining order of contigs, seems to be super of `contig2assembly`
    query2path = {
        query: join_seqs(
            query,
            contig2join=contig2join,
            path_circular_end=path_circular_end,
        )
        for _querys in (contig2assembly, contig_link_no_pe)
        for query in _querys
    }

    _log_info_2("Getting joined seqeuences.")
    # query: seq, mark, len of overlap
    query2extension = {
        query: (
            extend_query(query2path[query], contig2seq=contig2seq, maxk=maxk),
            *(
                ("extended_circular", maxk)
                if query in path_circular_potential
                else ("extended_partial", 0)
            ),
        )
        for query in tqdm(contig2assembly, desc="Extending path of query")
    }

    # Similar direct terminal repeats may lead to invalid joins
    _log_info_2("Checking for invalid joining using BLASTn: close strains.")
    blast_half_file = run_blast_half(
        (
            (contig, seq, overlap)
            for j in (
                ((contig, query2extension[contig]) for contig in contig2assembly),
                (
                    (contig, (contig2seq[contig], "", overlap))
                    for single_circle in (self_circular, self_circular_flex)
                    for contig in single_circle
                    for overlap in (self_circular_flex.get(contig, maxk),)
                ),
            )
            for contig, (seq, label, overlap) in j
        ),
        working_dir / "blast_pairs.tsv",
        threads=threads,
        chunk_size=100,
    )
    # blast_half_file = working_dir/f"blast_pairs.tsv"
    failed_blast_half = get_failed_blast_half(blast_half_file)
    # for debug
    _log_debug("# failed_blast_half")
    _log_debug(sorted(failed_blast_half))

    query2stat = {
        k: QueryCovVar.query(v, contig2=contig2cov, contig2len=contig2len)
        for k, v in tqdm(contig2assembly.items(), "stat query length and coverage")
    }

    _log_info_2("Grouping paths by sharing queries to check for invalid queries.")
    groups2ext_query = dict(
        enumerate(
            sorted(
                query2groups(contig2assembly=contig2assembly, query_set=query_set),
                key=lambda d: sorted(d),
            )
        )
    )
    check_assembly_reason = {
        k: v
        for groupi, group in groups2ext_query.items()
        for k, v in get_assembly2reason(
            groupi,
            group,
            contig2assembly=contig2assembly,
            path_circular_potential=path_circular_potential,
            contig_link_no_pe=contig_link_no_pe,
        ).items()
    }
    # for debug
    with open(working_dir / f"COBRA_check_assembly_reason.tsv", "w") as results:
        print("group", "reason", "query", "dups", sep="\t", file=results)
        for item in sorted(
            check_assembly_reason.items(), key=lambda x: (x[1][0], x[0])
        ):
            print(*item[1][:2], item[0], *item[1][2:], sep="\t", file=results)

    # region check assembly groups
    _log_info_2("Filtering paths accoring to COBRA rules.")
    failed_joins: dict[Literal["conflict", "circular_6", "nolink"], dict[str, int]] = {
        l: {  # type: ignore [misc]
            j: v[0]
            for v in check_assembly_reason.values()
            if v[1] == v1
            for i in groups2ext_query[v[0]]
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
        if j in assembly_rep
        and (
            any(end2contig(end) not in failed_blast_half for end in query2path[j])
            if any(end2contig(end) in path_circular_potential for end in query2path[j])
            else j not in failed_blast_half
        )
    }
    path_circular_rep = {
        q
        for q in checked_strict_rep
        if any(end2contig(end) in path_circular_potential for end in query2path[q])
    }

    redundant_circular_8 = {
        i: check_assembly_reason[i]
        for i in checked_strict.keys()
        if i in check_assembly_reason
        and check_assembly_reason[i][1] == "circular_8_tight"
    }
    _log_info(
        *(
            "\t".join(
                [*(compre_sets(pi - pj, li, lj) for pj, (j, lj) in enumerate(ll)), i]
            )
            for ll in (
                [
                    ("query_set", query_set),
                    ("orphan_end_query", orphan_query),
                    (
                        "self_circular_flexible_overlap",
                        self_circular_flex.keys(),
                    ),
                    ("self_circular", self_circular),
                    ("contig_link_no_pe", contig_link_no_pe),
                    ("failed_join [nolink]", failed_joins["nolink"].keys()),
                    ("failed_join [conflict]", failed_joins["conflict"].keys()),
                    ("failed_join [circular_6]", failed_joins["circular_6"].keys()),
                    ("path_circular_potential", path_circular_potential),
                    ("path_circular_rep", path_circular_rep),
                    ("checked_strict", checked_strict.keys()),
                    ("checked_strict_rep", checked_strict_rep.keys()),
                    ("assembly_rep", assembly_rep.keys()),
                    ("redundant_circular_8", redundant_circular_8.keys()),
                    ("failed_blast_half", failed_blast_half),
                ],
            )
            for pi, (i, li) in enumerate(ll)
        ),
        sep="\n",
    )
    # endregion check assembly groups

    def check_check_assembly_reason():
        assert query_set - orphan_pre_set == (
            query_set | contig_link_no_pe | failed_join_list | checked_strict.keys()
        )
        assert not checked_strict.keys() & contig_link_no_pe
        assert not checked_strict.keys() & failed_joins["nolink"].keys()
        assert not checked_strict.keys() & failed_joins["conflict"].keys()
        assert not checked_strict.keys() & failed_joins["circular_6"].keys()
        assert not self_circular & self_circular_flex.keys()
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
                f"{query                              =}",
                f"{query in query_set                 =}",
                f"{contig2assembly.get(query)         =}",
                f"{check_assembly_reason.get(query)   =}",
                f"{query in path_circular_potential             =}",
                f"{query in contig_link_no_pe     =}",
                f"{query in failed_joins['nolink']    =}",
                f"{query in failed_joins['conflict']  =}",
                f"{query in failed_joins['circular_6']=}",
                f"{query in checked_strict            =}",
                sep="\n",
            )
            print(groups2ext_query.get(check_assembly_reason.get(query, [-1])[0]))
            print({i: contig2assembly.get(i) for i in contig2assembly.get(query)})
            print({i: check_assembly_reason.get(i) for i in contig2assembly.get(query)})
            print(
                {
                    k: v
                    for k, v in groups2ext_query.items()
                    if v.keys() & contig2assembly.get(query, {})
                }
            )

        check_query("k141_10804945")

    # make blastn database and run search if the database is not empty
    ############################################################################
    ## Now output paths                                                       ##
    ############################################################################
    if not (checked_strict_rep or self_circular or self_circular_flex):
        # Of course we assume that sequence is not empty!
        _log_info(
            "no query was extended by COBRA, exit! "
            "this is normal if you only provide few queries."
        )
        # region ugly touch output
        with (
            open(working_dir / f"COBRA_joining_summary.tsv", "w") as assembly_summary,
            open(
                working_dir / f"COBRA_category_ii-c_extended_failed_paths.tsv", "w"
            ) as detail_failed,
            open(working_dir / f"COBRA_mark_failed.tsv", "w") as mark_failed,
        ):
            pass
        exit()

    _log_info_3 = get_log_info(
        5, "2.3. Output extended paths and circulated paths", logfile
    )
    # save the joining details information
    _log_info_3(
        'Getting the joining details of unique "extended_circular" and "extended_partial" query contigs.'
    )
    contig2join_detail = summary_join(
        working_dir=working_dir,
        checked_strict_rep=checked_strict_rep,
        query2extension=query2extension,
        path_circular_rep=path_circular_rep,
        query2path=query2path,
        contig2len=contig2len,
        contig2cov=contig2cov,
        contig2seq=contig2seq,
        contig2join_reason=contig2join_reason,
        contig2assembly=contig2assembly,
        query2stat=query2stat,
        query_set=query_set,
        maxk=maxk,
    )

    _log_info_3("Getting the joining details of failed query contigs.")
    with open(
        working_dir / f"COBRA_category_ii-c_extended_failed_paths.tsv", "w"
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
                for query in sorted(groups2ext_query[groupi].keys()):
                    for contig in groups2ext_query[groupi][query][1] or {query}:
                        if contig in contig_link_no_pe:
                            print(
                                *(failed_reason, groupi, query),
                                check_assembly_reason.get(query, ("", "redundant"))[1],
                                contig2len[query],
                                "",
                                sep="\t",
                                file=detail_failed,
                            )
                        else:
                            print(
                                *(failed_reason, groupi, query),
                                check_assembly_reason.get(query, ("", "redundant"))[1],
                                len(query2extension[query][0])
                                - maxk * (query in path_circular_rep),
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
                for query in sorted(groups2ext_query[groupi])
                for contig in sorted(groups2ext_query[groupi][query][1] or {query})
            ):
                if contig in contig_link_no_pe:
                    print(
                        *("blast_half", groupi, query),
                        check_assembly_reason[rep_query][1],
                        contig2len[query],
                        "",
                        sep="\t",
                        file=detail_failed,
                    )
                else:
                    print(
                        *("blast_half", groupi, query),
                        check_assembly_reason[rep_query][1],
                        len(query2extension[query][0])
                        - maxk * (query in path_circular_rep),
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
                if contig in contig_link_no_pe:
                    print(
                        *("circular_8", groupi, query),
                        check_assembly_reason[rep_query][1],
                        contig2len[query],
                        "",
                        sep="\t",
                        file=detail_failed,
                    )
                else:
                    print(
                        *("circular_8", groupi, query),
                        check_assembly_reason[rep_query][1],
                        len(query2extension[query][0])
                        - maxk * (query in path_circular_rep),
                        " ".join(query2path[query]),
                        sep="\t",
                        file=detail_failed,
                    )

    if skip_joining:
        _log_info(
            "\n",
            "3. RESULTS SUMMARY",
            f'# Check "{assembly_summary.name}" for successful joining queries.',
            f'# Check "{detail_failed.name}" for failed joining queries.',
            sep="\n",
        )
        return

    # save the joining summary information
    # save the joining status information of each query
    _log_info_3("Saving joining status of all query contigs.")
    assembled_info = open(working_dir / f"COBRA_joining_status.tsv", "w")
    # shows the COBRA status of each query
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
            working_dir / f"COBRA_category_ii-{ab}_{label}_unique.fa", "w"
        ) as extended_fasta:
            for query in sorted(checked_strict_rep):
                if label.endswith(contig2join_detail[query][0][2].lower()):
                    print(
                        f">{contig2join_detail[query][0][0]}\n{query2extension[query][0]}",
                        file=(extended_fasta),
                        flush=True,
                    )
                    for contig in sorted(contig2assembly[query] & query_set):
                        query_counts[label] += 1
                        print(
                            contig,
                            contig2len[contig],
                            contig2cov[contig],
                            round(GC(contig2seq[end2contig(contig)]), 3),
                            checked_strict_rep[query],
                            label,
                            f"category_ii-{ab}",
                            sep="\t",
                            file=assembled_info,
                        )
        extended_fa_names[label] = extended_fasta
    # for those cannot be extended
    with open(
        working_dir / f"COBRA_category_ii-c_extended_failed.fa", "w"
    ) as failed_join:
        label = "extended_failed"
        for query in sorted(failed_join_list):
            query_counts[label] += 1
            # assert query not in extended_circular_query
            # assert query not in extended_partial_query
            # assert query not in orphan_end_query
            print(f">{query}\n{contig2seq[query]}", file=failed_join)
            print(
                query,
                contig2len[query],
                contig2cov[query],
                round(GC(contig2seq[end2contig(query)]), 3),
                "",
                "extended_failed",
                "category_ii-c",
                sep="\t",
                file=assembled_info,
            )
    # for those due to orphan end
    with open(working_dir / f"COBRA_category_iii_orphan_end.fa", "w") as orphan_end:
        label = "orphan_end"
        for query in sorted(orphan_query):
            query_counts[label] += 1
            print(f">{query}\n{contig2seq[query]}", file=orphan_end)
            print(
                query,
                contig2len[query],
                contig2cov[query],
                round(GC(contig2seq[end2contig(query)]), 3),
                "",
                label,
                "category_iii",
                sep="\t",
                file=assembled_info,
            )

    # for self circular
    _log_info_3("Saving self_circular contigs.")
    with open(
        working_dir / f"COBRA_category_i_self_circular.fa", "w"
    ) as circular_fasta:
        label = "self_circular"
        for query in sorted(self_circular):
            query_counts[label] += 1
            print(f">{query}_self_circular\n{contig2seq[query]}", file=circular_fasta)
            print(
                query,
                contig2len[query] - maxk,
                contig2cov[query],
                round(GC(contig2seq[end2contig(query)]), 3),
                "",
                label,
                "category_i",
                sep="\t",
                file=assembled_info,
            )
        for query in sorted(self_circular_flex):
            query_counts[label] += 1
            print(f">{query}_self_circular\n{contig2seq[query]}", file=circular_fasta)
            print(
                query,
                contig2len[query] - self_circular_flex[query],
                contig2cov[query],
                round(GC(contig2seq[end2contig(query)]), 3),
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
            maxk,
            contig2cov=contig2cov,
            self_circular=self_circular,
            self_circular_flexible_overlap=self_circular_flex,
        )

    # save new fasta file with all the others used in joining replaced by COBRA sequences excepting self_circular ones
    _log_info_3("Saving the new fasta file.")
    query_extension_used = {
        contig: query
        for query in checked_strict_rep
        for contig in contig2assembly[query]
    }
    with open(working_dir / f"{assem_fa.name}.new.fa", "w") as new:
        for header in sorted(contig2seq.keys() - query_extension_used.keys()):
            print(f">{header}\n{contig2seq[header]}", file=new)
    os.system(
        f"cat {extended_fa_names['extended_circular'].name} {extended_fa_names['extended_partial'].name} >> {new.name}"
    )

    # write the numbers to the log file
    _log_info(
        "\n",
        "3. RESULTS SUMMARY",
        f"# Total queries: {len(query_set)}",
        f"# Category i   - self_circular: {query_counts['self_circular']}",
        f"# Category ii  - extended_circular: {query_counts['extended_circular']} (Unique: {len(checked_strict_rep.keys() & path_circular_potential)})",
        f"# Category ii  - extended_partial: {query_counts['extended_partial']} (Unique: {len(checked_strict_rep.keys() - path_circular_potential)})",
        f"# Category ii  - extended_failed (due to COBRA rules): {query_counts['extended_failed']}",
        f"# Category iii - orphan end: {query_counts['orphan_end']}",
        '# Check "COBRA_joining_status.tsv" for joining status of each query.',
        '# Check "COBRA_joining_summary.tsv" for joining details of "extended_circular" and "extended_partial" queries.',
        sep="\n",
    )


def summary_join(
    working_dir: Path,
    checked_strict_rep: dict[str, int],
    query2extension: dict[str, tuple[Seq, str, int]],
    path_circular_rep: set[str],
    query2path: dict[str, list[str]],
    contig2len: dict[str, int],
    contig2cov: dict[str, float],
    contig2seq: dict[str, Seq],
    contig2join_reason: dict[str, dict[str, str]],
    contig2assembly: dict[str, frozenset[str]],
    query2stat: dict[str, QueryCovVar],
    query_set: frozenset[str],
    maxk: int,
):
    class RepJoinDetail(NamedTuple):
        Final_Seq_ID: str
        Joined_Len: int
        Status: Literal["Circular", "Partial"]

    class ContigJoinDetail(NamedTuple):
        Query_Seq_ID: str
        Direction: Literal["forward", "reverse"]
        Joined_Seq_Len: int
        Start: int
        End: int
        Cov: float
        GC: float
        Joined_reason: str

    contig2join_detail: dict[str, tuple[RepJoinDetail, list[ContigJoinDetail]]] = {}
    for query in tqdm(
        checked_strict_rep,
        desc="Getting the joining details of extended query contigs.",
    ):
        site = 1
        join_status: Literal["Circular", "Partial"] = (
            "Circular" if query in path_circular_rep else "Partial"
        )
        query_extend_id = f"{query}_extended_{join_status.lower()}"
        contig2join_detail[query] = (
            RepJoinDetail(
                query_extend_id,
                len(query2extension[query][0]) - maxk * (query in path_circular_rep),
                join_status,
            ),
            [],
        )
        for item in query2path[query]:
            contig = end2contig(item)
            if (direction := get_direction(item)) == "forward":
                contig_start_end = (site, site + contig2len[contig] - 1)
            else:
                contig_start_end = (site, site + contig2len[contig] - 1)
            contig2join_detail[query][1].append(
                ContigJoinDetail(
                    contig,
                    direction,
                    contig2len[contig],
                    *contig_start_end,
                    contig2cov[contig],
                    round(GC(contig2seq[contig]), 3),
                    contig2join_reason[query][contig],
                )
            )
            site += contig2len[contig] - maxk
    with (
        open(
            working_dir
            / f"COBRA_category_ii-a_extended_circular_unique_joining_details.tsv",
            "w",
        ) as detail_circular,
        open(
            working_dir
            / f"COBRA_category_ii-b_extended_partial_unique_joining_details.tsv",
            "w",
        ) as detail_partial,
        open(working_dir / f"COBRA_joining_summary.tsv", "w") as assembly_summary,
    ):
        joining_detail_headers = [
            *RepJoinDetail._fields,
            *ContigJoinDetail._fields,
        ]
        print(*joining_detail_headers, sep="\t", file=detail_circular)
        print(*joining_detail_headers, sep="\t", file=detail_partial)
        print(
            *("Query_Seq_ID", "Query_Seq_Len"),  # contig, contig2len[contig]
            *RepJoinDetail._fields,
            *("Group_Id", "Total_Joined_Seqs", "Joined_seqs"),
            sep="\t",
            file=assembly_summary,
        )
        for query in sorted(contig2join_detail):
            for detail_values in contig2join_detail[query][1]:
                print(
                    *contig2join_detail[query][0],
                    *detail_values,
                    sep="\t",
                    file=(
                        detail_circular
                        if query in path_circular_rep
                        else detail_partial
                    ),
                )
            for contig in sorted(contig2assembly[query] & query_set):
                print(
                    *(contig, contig2len[contig]),
                    *contig2join_detail[query][0],
                    checked_strict_rep[query],
                    query2stat[query].query_count,
                    " ".join(seqjoin2contig(query2path[query])),
                    sep="\t",
                    file=assembly_summary,
                )
    return contig2join_detail


def main(args=None):
    args = parse_args(args)
    # get information from the input files and parameters and save information
    # get the name of the whole contigs fasta file
    # folder of output
    query_fa = Path(args.query)
    working_dir = Path(args.output) if args.output else Path(f"{query_fa.name}_COBRA")
    # check cache
    mapping_links = MappingLinks(args.mapping, args.mapping_link_cache)
    ##
    # get information from the input files and parameters and save information
    # get the name of the whole contigs fasta file
    # folder of output
    # checking if output folder exists
    if working_dir.exists():
        raise FileExistsError(f"Output folder <{working_dir}> exists, please check.")
    else:
        working_dir.mkdir(parents=True)
    # determine the length of overlap based on assembler and the largest kmer size
    maxk = args.maxk - 1 if args.assembler == "idba" else args.maxk
    # determine potential self_circular contigs from contigs with orphan end
    mink = args.mink - 1 if args.assembler == "idba" else args.mink

    with (
        open(working_dir / f"log", "w") as logfile,
        open(working_dir / f"debug.txt", "w") as debugfile,
    ):
        cobra(
            query_fa=query_fa,
            assem_fa=Path(args.fasta),
            mapping_links=mapping_links,
            coverage_file=Path(args.coverage),
            maxk=maxk,
            mink=mink,
            trim_readno=args.trim_readno,
            working_dir=working_dir,
            linkage_mismatch=args.linkage_mismatch,
            threads=args.threads,
            logfile=logfile,
            debugfile=debugfile,
            skip_joining=args.skip_joining,
        )


if __name__ == "__main__":
    main()
