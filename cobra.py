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
    from tqdm import tqdm as _tqdm

    def tqdm(__object: Iterable[T], *nargs, **kwargs) -> Iterable[T]:
        kwargs["disable"] = None
        return _tqdm(__object, *nargs, **kwargs)

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
        "--skip_new_assembly",
        action="store_true",
        help="whether skip new assembled contigs output to disk",
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


def group_yield(data: Iterable[T], n):
    d = iter(data)
    while l := [i for i, _ in zip(d, range(n))]:
        yield l


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

    def has_pe_link(end1: str, end2: str):
        return (end1, end2) in contig_pe_links

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
            and has_pe_link(contig + other_direction, end)
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
                    if checked_reason and has_pe_link(end0, end1):
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
                    # limited by mapping tools,
                    # we don't need has_pe_link(end0, end2add) here
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
                    # the same as above
                    checked_reason = check_1path_to_add(end1, contig)
                    if checked_reason and has_pe_link(end, end1):
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
                        if checked_reason:  # and has_pe_link(end, end2add):
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
    contig2cov: dict[str, float] = {}
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
    def __init__(self, mapping: str, mapping_link_cache: list[str] | None):
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
        contig2cov: dict[str, float],
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
        cov_lens = pd.Series({i: (contig2cov[i], contig2len[i]) for i in contigs})
        return cls(
            len(cov_lens),
            cov_lens.apply(lambda x: x[0]).mean(),
            cov_lens.apply(lambda x: x[0]).std(),
            cov_lens.apply(lambda x: x[1]).sum(),
            cov_lens.apply(lambda x: x[0] * x[1]).sum(),
        )


class ExtendQuery:
    def __init__(self, seq: Seq, label="", overlap_len=0):
        self.seq = seq
        self.overlap_len = overlap_len
        self.label = label
        self.len = len(seq) - overlap_len

    def __repr__(self):
        return f"<Query[{self.label}] {self.len} + {self.overlap_len}bp>"

    def to_fasta(self, name: str):
        return f">{name} len={self.len} overlap={self.overlap_len}\n{self.seq}"

    @classmethod
    def get_dict(
        cls,
        query2path: dict[str, list[str]],
        contig2seq: dict[str, Seq],
        maxk: int,
        path_circular_potential: frozenset[str],
        contig2assembly: Iterable[str],
    ):
        return {
            query: cls(
                cls.extend_query(query2path[query], contig2seq=contig2seq, maxk=maxk),
                (
                    "extended_circular"
                    if query in path_circular_potential
                    else "extended_partial"
                ),
                maxk * (query in path_circular_potential),
            )
            for query in tqdm(contig2assembly, desc="Extending path of query")
        }

    @staticmethod
    def extend_query(
        contigs: Iterable[str], contig2seq: dict[str, Seq], maxk: int
    ) -> Seq:
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


class BlastHalf:
    def __init__(self, outfile=Path("blast_pairs.tsv"), threads=10, chunk_size=200):
        self.outfile = outfile
        self.working_dir = outfile.parent
        self.threads = threads
        self.chunk_size = chunk_size

        self.min_overlap = 1000
        self.prec_identity = 70
        self.retries = 3

    def run(self, *packs: Iterable[tuple[str, ExtendQuery]]):
        self.run_blast_half(i for j in packs for i in j)
        return self.get_failed_blast_half()

    def write_and_run_blast(
        self, data_chunk: Iterable[tuple[str, ExtendQuery]], working_dir: Path
    ):
        cobra_seq2lenblast: dict[str, tuple[int, list[list[str]]]] = {}
        os.makedirs(working_dir)
        outlog = working_dir / "log"
        outfile = working_dir / "blastdb_2.vs.1.tsv"
        for i in range(1, self.retries + 1):
            try:
                with (
                    open(working_dir / "blastdb_1.fa", "w") as blastdb_1,
                    open(working_dir / "blastdb_2.fa", "w") as blastdb_2,
                ):
                    for header, eq in data_chunk:
                        cobra_seq2lenblast[header] = eq.len, []
                        half = eq.len // 2
                        print(f">{header}_1\n{eq.seq[:half]}", file=blastdb_1)
                        print(f">{header}_2\n{eq.seq[half:eq.len]}", file=blastdb_2)
                os.system(f"ls -sh {blastdb_1.name} {blastdb_2.name} > {outlog} 2>&1")
                os.system(
                    f"makeblastdb -in {blastdb_1.name} -dbtype nucl >> {outlog} 2>&1"
                )
                os.system(
                    "blastn"
                    f" -task blastn"
                    f" -db {blastdb_1.name}"
                    f" -query {blastdb_2.name}"
                    f" -out {outfile}"
                    f" -outfmt 6 -evalue 1e-10 -perc_identity {self.prec_identity}"
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
                if i == self.retries:
                    raise
        raise NotImplementedError

    def run_blast_half(self, query2compare: Iterable[tuple[str, ExtendQuery]]):
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
        temp_out = Path(f"{self.outfile}_temp")
        with (
            concurrent.futures.ThreadPoolExecutor(max_workers=self.threads) as executor,
            open(self.outfile, "w") as results,
        ):
            temp_out.mkdir()
            futures = {
                executor.submit(
                    self.write_and_run_blast, data_chunk, temp_out / f"tmp_{i}"
                )
                for i, data_chunk in enumerate(
                    group_yield(query2compare, self.chunk_size)
                )
            }
            for future in tqdm(
                concurrent.futures.as_completed(futures),
                desc="Running blastn",
                total=len(futures),
            ):
                for query, (seq_len, blastns) in future.result().items():
                    for blastn_v in blastns:
                        print(query, seq_len, *blastn_v[2:], sep="\t", file=results)

    def get_failed_blast_half(self):
        """
        identify potential incorrect joins and remove them from corresponding category
        """
        contig_len_overlap: dict[str, float] = {}
        with open(self.outfile) as r:
            for line in r:
                line_v = line.strip().split("\t")
                if float(line_v[3]) >= self.min_overlap:
                    if "_extended" in line_v[0]:
                        query = line_v[0].split("_extended")[0]
                    else:
                        query = line_v[0]
                    contig_len_overlap[query] = contig_len_overlap.get(
                        query, 0
                    ) + float(line_v[3])
        return frozenset(
            contig
            for contig in contig_len_overlap
            if contig_len_overlap[contig] >= self.min_overlap
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
        "complex_query",
        "longest",
    ]
    represent_seqs: list[str]
    dup_queries: frozenset[str]


class SubsetChunk(NamedTuple):
    standalong_subs: dict[str, "SubsetChunk"]
    frags: list[str]


def get_assembly2reason(
    groupi: int,
    group: dict[str, GroupAssemblyIndex],
    contig2assembly: dict[str, set[str]],
    path_circular_potential: frozenset[str],
    query_failed_join: frozenset[str],
):
    """
    for the group:
        1. check if the longest one is the superset of all others:
            True -> goto 2.
            False -> |
                check if there is sub standalone:
                    True -> cache the sub standalone, process it after 3.
                    False -> |
                        all those `not <= or >=` the longest one marked failed_wait
                        for each one in failed_wait:
                            add other `not <= or >=` to failed_wait
                goto 1.
        2. check if the longest one is avail, i.e. not the `6` or `8`
            `6` -> marked the longest and all the sub circular failed, and got 1.
            `8` -> reject the largest one and goto 1.
            False -> goto 3.
        3. check if any contig2assembly[frag] & query_failed_join
            True -> |
                Mark the frag as "complex_query",
                remove all paths containing this frag.
                goto 1.
            False -> |
                record the longest one
                if sth cached in 1.False.True:
                    True -> pop it and goto 1.
                    False -> return
    """
    if len(group) == 1:
        for this in group:
            judge: Literal["complex_query", "standalone"] = (
                "complex_query"
                if contig2assembly[this] & query_failed_join
                else "standalone"
            )
            return {this: AssemblyReason(groupi, judge, [], group[this].dup_queries)}
    subsets = sorted(
        group, key=lambda q: (len(group[q].special), len(group[q].dup_queries))
    )

    def clean_failed(*_failed_wait: str, reason="conflict_query"):
        nonlocal subsets
        failed_wait = set(_failed_wait)
        _this = _failed_wait[0]
        check_assembly_reason[_this] = AssemblyReason(
            groupi, reason, [], group[_this].dup_queries
        )
        while failed_wait:
            this = failed_wait.pop()
            for frag in subsets:
                conflict = group[frag].special & group[this].special
                if (
                    conflict
                    and (conflict < group[this].special)
                    and (conflict < group[frag].special)
                ):
                    failed_wait.add(frag)
            subsets = [i for i in subsets if i not in failed_wait and i != this]
            check_assembly_reason[_this].represent_seqs.append(this)

    check_assembly_reason: dict[str, AssemblyReason] = {}
    standalong_subs: list[list[str]] = []
    while subsets:
        # step 1
        this = subsets.pop(-1)
        conflicts = [
            frag for frag in subsets if not group[this].special >= group[frag].special
        ]
        if not conflicts:
            # the best one, meaning that this is representative to all others
            pass
        elif any(group[this].special & group[frag].special for frag in conflicts):
            clean_failed(this)
            continue
        else:
            # a full difference occurs
            subsets = [i for i in subsets if i not in conflicts]
            standalong_subs.append(conflicts)
        # step 2
        subset_circular = [i for i in subsets if i in path_circular_potential]
        if subset_circular:
            if this in path_circular_potential:
                clean_failed(this, reason="circular_8_tight")
            else:
                clean_failed(this, *subset_circular, reason="circular_6_conflict")
            continue
        # step 3
        no_pe = [i for i in subsets if contig2assembly[i] & query_failed_join]
        if no_pe or (contig2assembly[this] & query_failed_join):
            clean_failed(this, *no_pe, reason="complex_query")
            continue
        check_assembly_reason[this] = AssemblyReason(
            groupi,
            "longest",
            subsets,
            group[this].dup_queries,
        )
        if standalong_subs:
            subsets = standalong_subs.pop()
            continue
        break
    return check_assembly_reason


def log_group2graphviz(
    groups2ext_query: dict[int, dict[str, GroupAssemblyIndex]],
    contig2join: dict[str, list[str]],
    contig2cov: dict[str, float],
    contig_pe_links: frozenset[tuple[str, str]],
    file=sys.stdout,
):
    for groupi, group in groups2ext_query.items():
        if len(group) == 1:
            continue
        contigs: set[str] = set(group)
        contig_links: set[tuple[str, str]] = set()
        for contig in group:
            for end in (f"{contig}_L", f"{contig}_R"):
                last_end = end
                for next_end in contig2join.get(end, []):
                    if next_end.endswith("rc"):
                        next_end = next_end[:-2]
                    contig_links.add((last_end, next_end))
                    last_end = end2end2(next_end)
                    contigs.add(end2contig(next_end))
        gv_names: list[str] = []
        gv_covs: list[str] = []
        gv_links: list[str] = []
        for contig in contigs:
            gv_names.append(f'{contig}_L -> {contig}_R [label="{contig}"; color=blue]')
            gv_covs.append(
                f'{contig}_R -> {contig}_L [label="{contig2cov[contig]}"; color=blue]'
            )
        for pair in contig_links:
            color = ""
            if not (
                pair in contig_pe_links or (pair[0], pair[1] + "rc") in contig_pe_links
            ):
                color = "[color=gray]"
            gv_links.append(f"{pair[0]} -> {pair[1]} {color}")
        gv_str = (
            f"digraph group_{groupi} "
            + "{node[shape=box, style=rounded];"
            + (";".join(gv_names) + ";")
            + (";".join(gv_covs) + ";")
            + (";".join(gv_links) + ";")
            + "}"
        )
        print(gv_str, file=file)


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
    skip_new_assembly=False,
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
    orphan_query = (orphan_pre_set - self_circular_flex.keys()) & query_set
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
    # named failed_join_list in old versions
    query_failed_join = frozenset(
        contig
        for contig in query_set
        if (
            (contig + "_L" in link_pair or contig + "_R" in link_pair)
            and contig not in self_circular
            and len(contig2join[contig + "_L"]) == 0
            and len(contig2join[contig + "_R"]) == 0
        )
    )
    _log_debug("# query_failed_join (too complex to extend, so give up)")
    _log_debug(sorted(query_failed_join))

    # get the joining order of contigs, seems to be super of `contig2assembly`
    query2path = {
        query: join_seqs(
            query,
            contig2join=contig2join,
            path_circular_end=path_circular_end,
        )
        for _querys in (contig2assembly, query_failed_join)
        for query in _querys
    }

    _log_info_2("Getting joined seqeuences.")
    # query: seq, mark, len of overlap
    query2extension = (
        ExtendQuery.get_dict(
            query2path=query2path,
            contig2seq=contig2seq,
            maxk=maxk,
            path_circular_potential=path_circular_potential,
            contig2assembly=contig2assembly,
        )
        | {
            contig: ExtendQuery(contig2seq[contig], "overlap_maxk", maxk)
            for contig in self_circular
        }
        | {
            contig: ExtendQuery(contig2seq[contig], "overlap_flex", overlap)
            for contig, overlap in self_circular_flex.items()
        }
    )

    # Similar direct terminal repeats may lead to invalid joins
    _log_info_2("Checking for invalid joining using BLASTn: close strains.")
    failed_blast_half = BlastHalf(
        working_dir / "blast_pairs.tsv", threads=threads, chunk_size=100
    ).run(query2extension.items())
    # for debug
    _log_debug("# failed_blast_half")
    _log_debug(sorted(failed_blast_half))

    _log_info_2("Grouping paths by sharing queries to check for invalid queries.")
    groups2ext_query = dict(
        enumerate(
            sorted(
                query2groups(contig2assembly=contig2assembly, query_set=query_set),
                key=lambda d: sorted(d),
            )
        )
    )
    ext_query2group = {
        queryi: groupi
        for groupi, group in groups2ext_query.items()
        for queryi in group.keys()
        | {i for j in group.values() for i in j.special}
        | {i for j in group.values() for i in j.dup_queries}
    }
    with open(working_dir / f"COBRA_grouping_gv.tsv", "w") as results:
        log_group2graphviz(
            contig2join=contig2join,
            groups2ext_query=groups2ext_query,
            contig2cov=contig2cov,
            contig_pe_links=contig_pe_links,
            file=results,
        )
    check_assembly_reason = {
        k: v
        for groupi, group in groups2ext_query.items()
        for k, v in get_assembly2reason(
            groupi,
            group,
            contig2assembly=contig2assembly,
            path_circular_potential=path_circular_potential,
            query_failed_join=query_failed_join,
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
    failed_joins: dict[Literal["conflict", "circular_6", "complex"], dict[str, int]] = {
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
            ("complex", "complex_query"),
        )
    }
    failed_groups: dict[int, Literal["conflict", "circular_6", "complex"]] = {
        i: q  # type: ignore [misc]
        for q in reversed(("conflict", "circular_6", "complex"))
        for i in failed_joins[q].values()  # type: ignore [index]
    }
    checked_strict = {
        j: v[0]
        for k, v in check_assembly_reason.items()
        if v[1] in {"standalone", "longest"} and v[0] not in failed_groups
        for j in contig2assembly[k] & contig2assembly.keys()
    }

    assembly_rep = {
        max(v[3], key=lambda x: query2extension[x].len) if len(v[3]) else i: v[0]
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
    extended_status = (
        pd.Series(
            {
                "circular": path_circular_rep,
                "partial": checked_strict_rep.keys() - path_circular_rep,
                "failed_blast_half": path_circular_rep & failed_blast_half,
                "redundant_circular_8": redundant_circular_8.keys(),
            }
            | {k: v.keys() - query_failed_join for k, v in failed_joins.items()},
            name="RepQuery",
        )
        .rename_axis(index="StatusReason")
        .reset_index()
        .assign(
            Status=lambda s: s["StatusReason"].apply(
                lambda x: "extended_{}".format(
                    x if x in ("circular", "partial") else "failed"
                )
            )
        )
        .explode("RepQuery")
        .dropna(subset=["RepQuery"])
        .assign(
            Item=lambda df: df["RepQuery"].map(query2path),
            Contig=lambda df: df["Item"].apply(lambda x: [end2contig(i) for i in x]),
        )
    )
    extended_status_queries = (
        extended_status.assign(
            Query=lambda s: s["Contig"].apply(lambda x: set(x) & query_set)
        )
        .explode("Query")
        .groupby("Query")["Status"]
        .apply(lambda x: " ".join(sorted(set(x))))
        .reset_index()
        .groupby("Status")["Query"]
        .apply(frozenset)
    )
    _log_info(
        *(
            "\t".join(
                [*(compre_sets(pi - pj, li, lj) for pj, (j, lj) in enumerate(ll)), i]
            )
            for ll in (
                [
                    ("query_set", query_set),
                    ("orphan_end_query", orphan_query),
                    ("complex_end_query", query_failed_join),
                    ("self_circular [overlap_maxk]", self_circular),
                    ("self_circular [overlap_flex]", self_circular_flex.keys()),
                    ("assembly_rep", assembly_rep.keys()),
                    ("failed_join_queries", extended_status_queries["extended_failed"]),
                    ("failed_join [complex]", failed_joins["complex"].keys()),
                    ("failed_join [conflict]", failed_joins["conflict"].keys()),
                    ("failed_join [circular_6]", failed_joins["circular_6"].keys()),
                    ("failed_blast_half", failed_blast_half),
                    ("redundant_circular_8", redundant_circular_8.keys()),
                    ("checked_rep", checked_strict_rep.keys()),
                    (
                        "checked_partial_queries",
                        extended_status_queries["extended_partial"],
                    ),
                    (
                        "checked_circular_queries",
                        extended_status_queries["extended_circular"],
                    ),
                    ("path_circular_rep", path_circular_rep),
                ],
            )
            for pi, (i, li) in enumerate(ll)
        ),
        sep="\n",
    )
    # endregion check assembly groups

    # make blastn database and run search if the database is not empty
    summary_join_tsv = working_dir / f"COBRA_joining_summary.tsv"
    summary_fail_tsv = working_dir / f"COBRA_joining_failed_paths.tsv"
    if not (checked_strict_rep or self_circular or self_circular_flex):
        # Of course we assume that sequence is not empty!
        _log_info(
            "no query was extended by COBRA, exit! "
            "this is normal if you only provide few queries."
        )
        # region ugly touch output
        with (
            open(summary_join_tsv, "w"),
            open(summary_fail_tsv, "w"),
            open(working_dir / f"COBRA_mark_failed.tsv", "w"),
        ):
            return
        # endregion ugly touch output

    _log_info_3 = get_log_info(
        5, "2.3. Output extended paths and circulated paths", logfile
    )
    # save the joining details information
    _log_info_3("Getting the joining details of extended query contigs.")
    extended_status_detail = (
        extended_status.assign(
            FinalSeqID=lambda df: df.apply(
                lambda s: "_".join((s["RepQuery"], s["Status"])),
                axis=1,
            ),
            GroupID=lambda df: df["RepQuery"].apply(lambda x: ext_query2group[x]),
            JoinLen=lambda df: df["RepQuery"].map(query2extension).map(lambda x: x.len),
            StartOnJoin=lambda df: df.apply(
                lambda x: pd.Series(
                    [0]
                    + [contig2len[end2contig(contig)] - maxk for contig in x["Item"]][
                        :-1
                    ]
                )
                .cumsum()
                .values,
                axis=1,
            ),
            Direction=lambda df: df["Item"].apply(
                lambda x: ["reverse" if i.endswith("rc") else "forward" for i in x]
            ),
            Item=lambda df: df["Item"].apply(" ".join),
        )
        .sort_values(["Status", "StatusReason", "GroupID", "RepQuery"])
        .explode(["Direction", "Contig", "StartOnJoin"])
        .assign(
            JoinedReason=lambda df: df.apply(
                lambda x: contig2join_reason[x["RepQuery"]][x["Contig"]],
                axis=1,
            ),
            ContigLen=lambda df: df["Contig"].apply(lambda x: len(contig2seq[x])),
            Cov=lambda df: df["Contig"].map(contig2cov),
            GC=lambda df: df["Contig"].apply(lambda x: GC(contig2seq[x])),
            IsQuery=lambda df: df["Contig"].apply(lambda x: x in query_set),
        )
    )
    extended_status_detail.query("Status != 'extended_failed'").to_csv(
        summary_join_tsv, sep="\t", index=False
    )
    # extended_status_detail.query("Status != 'extended_failed'").query("IsQuery == True").groupby("Contig")["RepQuery"].count().max()

    _log_info_3("Getting the joining details of failed query contigs.")
    extended_failed_detail = (
        extended_status_detail.query("Status == 'extended_failed'")
        .query("IsQuery == True")[
            ["StatusReason", "GroupID", "Contig", "ContigLen", "Cov", "GC"]
        ]
        .drop_duplicates()
        .assign(
            FailedReason=lambda df: df["Contig"].apply(
                lambda x: check_assembly_reason.get(x, ["", "dup"])[1]
            ),
            Item=lambda df: df["Contig"].map(query2path).apply(" ".join),
        )
    )
    extended_failed_detail.to_csv(summary_fail_tsv, sep="\t", index=False)

    extended_self_circular = (
        pd.Series(
            {
                "overlap_maxk": self_circular,
                "overlap_flex": self_circular_flex,
            },
            name="Contig",
        )
        .rename_axis(index="StatusReason")
        .reset_index()
        .explode("Contig")
        .dropna(subset=["Contig"])
        .assign(
            Status="self_circular",
            Overlap=lambda df: df["Contig"].apply(
                lambda x: self_circular_flex.get(x, maxk)
            ),
        )
    )
    extended_ignored = (
        pd.Series(
            {"orphan_end": orphan_query, "complex_end": query_failed_join},
            name="Contig",
        )
        .rename_axis(index="StatusReason")
        .reset_index()
        .assign(Status="Ignored")
        .explode("Contig")
        .dropna(subset=["Contig"])
        .merge(
            extended_status_detail.groupby("Contig")[["RepQuery", "GroupID"]].agg(
                lambda x: " ".join(sorted({str(i) for i in x}))
            ),
            on="Contig",
            how="left",
        )
    )

    _log_info_3("Saving joining status of all query contigs.")
    query_status = concat_query_status(
        query_set=query_set,
        extended_self_circular=extended_self_circular,
        extended_status_detail=extended_status_detail,
        extended_failed_detail=extended_failed_detail,
        extended_ignored=extended_ignored,
        maxk=maxk,
        logfile=logfile,
        debugfile=debugfile,
    )
    query_status.to_csv(working_dir / f"COBRA_query_status.tsv", sep="\t", index=False)

    # save the joining status information of each query
    _log_info_3("Saving identified and modified contigs.")
    # for those could be extended
    extended_fa_names = {}
    for ab, label in (("a", "extended_circular"), ("b", "extended_partial")):
        with open(
            working_dir / f"COBRA_category_ii-{ab}_{label}_unique.fa", "w"
        ) as extended_fasta:
            extended_status_detail.query("IsQuery == True").query(
                f"Status == '{label}'"
            )[["RepQuery", "FinalSeqID"]].drop_duplicates().apply(
                lambda s: print(
                    query2extension[s["RepQuery"]].to_fasta(s["FinalSeqID"]),
                    file=extended_fasta,
                    flush=True,
                ),
                axis=1,
            )
        extended_fa_names[label] = extended_fasta

    _log_info_3(
        "Skip saving the new assembled contigs."
        if skip_new_assembly
        else "Summarise all query contigs and save the new assembled contigs."
    )
    if not skip_new_assembly:

        def contig2fasta(conig, *labels):
            return ">{} {}\n{}".format(conig, " ".join(labels), contig2seq[conig])

        # for self circular
        with open(
            working_dir / f"COBRA_category_i_self_circular.fa", "w"
        ) as circular_fasta:
            label = "self_circular"
            extended_self_circular.drop_duplicates().apply(
                lambda s: print(
                    contig2fasta(
                        s["Contig"],
                        label,
                        "DTR_length={}".format(
                            query2extension[s["Contig"]].overlap_len
                        ),
                    ),
                    file=circular_fasta,
                ),
                axis=1,
            )
        # for those cannot be extended
        with open(
            working_dir / f"COBRA_category_ii-c_extended_failed.fa", "w"
        ) as failed_join:
            label = "extended_failed"
            extended_failed_detail.drop_duplicates(subset=["Contig"]).apply(
                lambda s: print(contig2fasta(s["Contig"], label), file=failed_join),
                axis=1,
            )
        # for those due to orphan end or complex end
        for ab, label in (("a", "orphan_end"), ("b", "complex_end_rest")):
            with open(
                working_dir / f"COBRA_category_iii-{ab}_{label}.fa", "w"
            ) as extended_fasta:
                extended_ignored.query("RepQuery.isna()").query(
                    "StatusReason == '{}'".format(label.replace("_rest", ""))
                ).drop_duplicates().apply(
                    lambda s: print(
                        contig2fasta(s["Contig"], label), file=extended_fasta
                    ),
                    axis=1,
                )
            extended_fa_names[label] = extended_fasta

        for outio in (
            circular_fasta,
            extended_fa_names["extended_circular"],
            extended_fa_names["extended_partial"],
            failed_join,
            extended_fa_names["orphan_end"],
            extended_fa_names["complex_end_rest"],
        ):
            summary_fasta(
                outio.name,
                maxk,
                contig2cov=contig2cov,
                self_circular=self_circular,
                self_circular_flexible_overlap=self_circular_flex,
            )

        # save new fasta file with all the others used in joining replaced by COBRA sequences excepting self_circular ones
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
        "# Total queries: {}".format(len(query_set)),
        "# Category i   - self_circular: {}".format(
            sum(query_status["Status"] == "self_circular")
        ),
        "# Category ii  - extended_circular: {} (Unique: {})".format(
            sum(query_status["Status"] == "extended_circular"),
            len(
                query_status.query("Status == 'extended_circular'")[
                    "FinalSeqID"
                ].unique()
            ),
        ),
        "# Category ii  - extended_partial: {} (Unique: {})".format(
            sum(query_status["Status"] == "extended_partial"),
            len(
                query_status.query("Status == 'extended_partial'")[
                    "FinalSeqID"
                ].unique()
            ),
        ),
        "# Category ii  - extended_failed (due to COBRA rules): {}".format(
            sum(query_status["Status"] == "extended_failed")
        ),
        "# Category iii - orphan end: {}".format(
            sum(query_status["StatusReason"] == "orphan_end")
        ),
        "# Category iii - too complex to: {}".format(
            sum(query_status["StatusReason"] == "complex_end")
        ),
        '# Check "COBRA_query_status.tsv" for joining status of each query.',
        '# Check "COBRA_joining_summary.tsv" for joining details of "extended_circular" and "extended_partial" queries.',
        sep="\n",
    )


def concat_query_status(
    query_set: frozenset[str],
    extended_self_circular: pd.DataFrame,
    extended_status_detail: pd.DataFrame,
    extended_failed_detail: pd.DataFrame,
    extended_ignored: pd.DataFrame,
    maxk: int,
    logfile=sys.stdout,
    debugfile=sys.stderr,
):
    key_cols = ["Status", "StatusReason", "Contig"]

    query_df = pd.DataFrame({"Contig": sorted(query_set)})
    extend2concat = [
        extended_self_circular[[*key_cols, "Overlap"]],
        extended_status_detail.query("IsQuery == True")
        .query("Status != 'extended_failed'")
        .assign(Overlap=maxk, GroupID=lambda df: df["GroupID"].apply(str))[
            [*key_cols, "GroupID", "FinalSeqID", "Overlap"]
        ],
        extended_failed_detail.assign(
            Status="extended_failed",
            StatusReason=lambda df: df.merge(
                extended_ignored[["Contig", "StatusReason"]], on="Contig", how="left"
            ).apply(
                lambda s: " ".join(s[["StatusReason_x", "StatusReason_y"]].dropna()),
                axis=1,
            ),
        )[[*key_cols, "GroupID"]],
        extended_ignored.query("RepQuery.isna()")[key_cols],
    ]
    query_status = (
        query_df.merge(pd.concat(extend2concat))
        .groupby("Contig")
        .agg(lambda x: " | ".join(sorted({str(i) for i in x.dropna()})))
        .reset_index()
        .sort_values(key_cols)
    )
    abnormal_dup_queries = query_status.query("' | ' in StatusReason")
    if abnormal_dup_queries.shape[0]:
        print(
            f"\nWarning: a total of {abnormal_dup_queries.shape[0]} queries"
            f" were in duplicated catalog, "
            f"this may caused by bug in Cobra."
            f"\nPlease contact maintainer for help."
            f"\n",
            file=logfile,
        )
        print("# query contigs in duplicated catalog", file=debugfile)
        abnormal_dup_queries.to_csv(debugfile, sep="\t", index=False)
    abnormal_drop_queries = query_status.query("StatusReason == ''")
    if abnormal_drop_queries.shape[0]:
        print(
            f"\nWarning: a total of {abnormal_drop_queries.shape[0]} queries"
            f" were not in any catalog, "
            f"this may caused by bug in Cobra."
            f"\nPlease contact maintainer for help."
            f"\n",
            file=logfile,
        )
        print("# query contigs not in any catalog", file=debugfile)
        abnormal_drop_queries.to_csv(debugfile, sep="\t", index=False)
    return query_status


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
            skip_new_assembly=args.skip_new_assembly,
        )


if __name__ == "__main__":
    main()
