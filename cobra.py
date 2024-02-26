#!/usr/bin/env python
# Author: LinXing Chen, UC Berkeley

# cobra v1.2.3
# Contig Overlap Based Re-Assembly
# Modification date: Sep 3, 2023


import os
import Bio
from Bio import SeqIO
from Bio.Seq import reverse_complement
from collections import defaultdict
import argparse
import pysam
from time import strftime
import itertools

bio_version = Bio.__version__
if bio_version > '1.79':
    from Bio.SeqUtils import gc_fraction as GC
else:
    from Bio.SeqUtils import GC

parser = argparse.ArgumentParser(description="This script is used to get higher quality (including circular) virus genomes "
                                             "by joining assembled contigs based on their end overlaps.")
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("-q", "--query", type=str, help="the query contigs file (fasta format), or the query name "
                                                           "list (text file, one column).", required=True)
requiredNamed.add_argument("-f", "--fasta", type=str, help="the whole set of assembled contigs (fasta format).", required=True)
requiredNamed.add_argument("-a", "--assembler", type=str, choices=["idba", "megahit", "metaspades"],
                           help="de novo assembler used, COBRA not tested for others.", required=True)
requiredNamed.add_argument("-mink", "--mink", type=int, help="the min kmer size used in de novo assembly.", required=True)
requiredNamed.add_argument("-maxk", "--maxk", type=int, help="the max kmer size used in de novo assembly.", required=True)
requiredNamed.add_argument("-m", "--mapping", type=str, help="the reads mapping file in sam or bam format.", required=True)
requiredNamed.add_argument("-c", "--coverage", type=str, help="the contig coverage file (two columns divided by tab).",
                           required=True)
parser.add_argument("-lm", "--linkage_mismatch", type=int, default=2, help="the max read mapping mismatches for "
                                                                           "determining if two contigs are spanned by "
                                                                           "paired reads. [2]")
parser.add_argument("-o", "--output", type=str, help="the name of output folder (default = '<query>_COBRA').")
parser.add_argument("-t", "--threads", type=int, default=16, help="the number of threads for blastn. [16]")
parser.add_argument("-v", "--version", action='version', version='cobra v1.2.3')
args = parser.parse_args()

##
#
global one_path_end, link_pair, cov, parsed_linkage, self_circular, two_paths_end, contig2join, contig_checked, \
    contig2join_reason, path_circular, path_circular_end, order_all, added_to_contig, header2seq, \
    self_circular_non_expected_overlap, all_joined_query, extended_circular_query, extended_partial_query, \
    header2len, is_subset_of

#
cov = {}  # the coverage of contigs
header2seq = {}
header2len = {}
link_pair = {}  # used to save all overlaps between ends
parsed_linkage = set()  # Initialize an empty set to store the parsed linkage information
one_path_end = []  # the end of contigs with one potential join
two_paths_end = []  # the end of contigs with two potential joins
self_circular = set()
self_circular_non_expected_overlap = {}
contig2join = {}
contig_checked = {}
contig2join_reason = {}
path_circular_end = set()
path_circular = set()
is_subset_of = {}
extended_circular_query = set()
extended_partial_query = set()
order_all = {}
added_to_contig = {}
all_joined_query = set()


##
# define functions for analyses
def log_info(step, description, line_feed, log_file):
    localtime = ' [20' + strftime("%x").rsplit('/', 1)[1] + '/' + strftime("%x").rsplit('/', 1)[0] + ' ' + strftime("%X") + '] '
    print(step + localtime + description, end=line_feed, file=log_file, flush=True)


def determine_file_format(filename):
    """
    Determine the format of the input file.
    Returns either "fasta" or "txt" based on the content.
    """
    with open(filename, 'r') as f:
        first_line = f.readline().strip()
        if first_line.startswith(">"):
            return "fasta"
        else:
            return "txt"


def contig_name(end_name):
    """
    get contig name from end name
    """
    if end_name.endswith('_L') or end_name.endswith('_R') or end_name.endswith('_Lrc') or end_name.endswith('_Rrc'):
        return end_name.rsplit('_', 1)[0]
    else:
        return end_name


def get_target(item):
    """
    to get the target for next run of joining
    """
    suffix = item.rsplit('_', 1)[1]
    base_name = item.rsplit('_', 1)[0]

    if suffix == 'Lrc':
        return base_name + '_Rrc'
    elif suffix == 'Rrc':
        return base_name + '_Lrc'
    elif suffix == 'R':
        return base_name + '_L'
    else:
        return base_name + '_R'


def are_equal_paths(two_paths):
    """
    evaluate if two paths are equal, two_paths is a list including two ends, e.g., a_L, b_R
    """
    if two_paths[0].rsplit('_', 1)[1].startswith('L') and two_paths[1].rsplit('_', 1)[1].startswith('R'):
        return contig_name(two_paths[0]) + '_R' in one_path_end \
               and contig_name(two_paths[1]) + '_L' in one_path_end \
               and contig_name(link_pair[contig_name(two_paths[0]) + '_R'][0]) == \
               contig_name(link_pair[contig_name(two_paths[1]) + '_L'][0])
    elif two_paths[0].rsplit('_', 1)[1].startswith('R') and two_paths[1].rsplit('_', 1)[1].startswith('L'):
        return contig_name(two_paths[0]) + '_L' in one_path_end \
               and contig_name(two_paths[1]) + '_R' in one_path_end \
               and contig_name(link_pair[contig_name(two_paths[0]) + '_L'][0]) == \
               contig_name(link_pair[contig_name(two_paths[1]) + '_R'][0])
    elif two_paths[0].rsplit('_', 1)[1].startswith('R') and two_paths[1].rsplit('_', 1)[1].startswith('R'):
        return contig_name(two_paths[0]) + '_L' in one_path_end \
               and contig_name(two_paths[1]) + '_L' in one_path_end \
               and contig_name(link_pair[contig_name(two_paths[0]) + '_L'][0]) == \
               contig_name(link_pair[contig_name(two_paths[1]) + '_L'][0])
    elif two_paths[0].rsplit('_', 1)[1].startswith('L') and two_paths[1].rsplit('_', 1)[1].startswith('L'):
        return contig_name(two_paths[0]) + '_R' in one_path_end \
               and contig_name(two_paths[1]) + '_R' in one_path_end \
               and contig_name(link_pair[contig_name(two_paths[0]) + '_R'][0]) == \
               contig_name(link_pair[contig_name(two_paths[1]) + '_R'][0])
    else:
        return False


def the_dominant_one(two_paths):
    """
    get the dominant path from two equal paths
    """
    if cov[contig_name(two_paths[0])] >= cov[contig_name(two_paths[1])]:
        return two_paths[0]
    else:
        return two_paths[1]


def could_circulate(point, contig, direction):
    """
    check if the path is circular with the current contig included
    """
    if len(link_pair[point]) == 2:  # point is the same as "target" in join_walker
        if direction == 'L':
            if contig_name(link_pair[point][0]) == contig:
                if cov[contig_name(point)] < 1.5 * cov[contig]:
                    # 2 times is for repeat, but it is too risky, use 1.5 instead (same below)
                    return link_pair[point][0].endswith('_R') and (contig + '_R', point) in parsed_linkage
                else:
                    return False
            elif contig_name(link_pair[point][1]) == contig:
                if cov[contig_name(point)] < 1.5 * cov[contig]:
                    return link_pair[point][1].endswith('_R') and (contig + '_R', point) in parsed_linkage
                else:
                    return False
            else:
                return False
        elif direction == 'R':
            if contig_name(link_pair[point][0]) == contig:
                if cov[contig_name(point)] < 1.5 * cov[contig]:
                    return link_pair[point][0].endswith('_L') and (contig + '_L', point) in parsed_linkage
                else:
                    return False
            elif contig_name(link_pair[point][1]) == contig:
                if cov[contig_name(point)] < 1.5 * cov[contig]:
                    return link_pair[point][1].endswith('_L') and (contig + '_L', point) in parsed_linkage
                else:
                    return False
            else:
                return False
        else:
            return False

    elif len(link_pair[point]) == 1:
        if direction == 'L':
            if contig_name(link_pair[point][0]) == contig:
                return link_pair[point][0].endswith('_R') and (contig + '_R', point) in parsed_linkage
        elif direction == 'R':
            if contig_name(link_pair[point][0]) == contig:
                return link_pair[point][0].endswith('_L') and (contig + '_L', point) in parsed_linkage

    else:
        return False


def the_better_one(two_paths, contig):
    """
    Calculate the absolute differences in coverage
    """
    diff0 = abs(cov[contig] - cov[contig_name(two_paths[0])])
    diff1 = abs(cov[contig] - cov[contig_name(two_paths[1])])

    # If diff0 is 0, return two_paths[0]
    if diff0 == 0:
        return two_paths[0]
    else:
        # If diff1 is not 0, calculate the ratio and compare
        if diff1 != 0:
            ratio = diff0 / diff1

            if ratio > 1:
                return two_paths[1]
            elif ratio <= 1:
                return two_paths[0]
        # If diff1 is 0, return two_paths[1]
        else:
            return two_paths[1]


def other_end_is_extendable(end, contig):
    """
    check if the other end is extendable
    """
    if contig_name(end) not in self_circular:
        if get_target(end) in link_pair.keys():
            return len(link_pair[get_target(end)]) > 0
        else:
            return False
    else:
        return False


def is_ok_to_add(end, contig):
    """
    if the other end of a potential joining contig cannot be extended,
    add only when it has a very similar coverage to that of the query contig
    """
    if contig_name(end) not in self_circular:
        return 0.9 * cov[contig] <= cov[contig_name(end)] <= 1.11 * cov[contig]


def not_checked(end_list, checked):
    """
    to see if a potential contig has been checked for adding or not
    """
    total_checked = 0
    for end in end_list:
        if contig_name(end) in checked:
            total_checked += 1
    return total_checked == 0


def detect_self_circular(contig):
    """
    to determine the input sequences if they are self_circular genomes or not
    """
    end = contig + '_L'
    if end in one_path_end:
        if link_pair[end][0] == contig + '_R':
            # the other end could be joined with the current working end
            self_circular.add(contig)
        else:
            pass
    elif end in two_paths_end:
        if link_pair[end][0].rsplit('_', 1)[0] != link_pair[end][1].rsplit('_', 1)[0]:
            if contig + '_R' in link_pair[end]:
                self_circular.add(contig)
            else:
                pass
        else:
            pass
    else:
        pass


def join_walker(contig, direction):
    """
    get potential joins for a given query
    """
    global contig2join, contig_checked, contig2join_reason, path_circular, path_circular_end

    end = contig + '_' + direction
    a = len(contig2join[end])
    if a == 0:
        contig_checked[end].append(contig)
        if end in one_path_end:
            if contig_name(link_pair[end][0]) != contig:
                if other_end_is_extendable(link_pair[end][0], contig) and (end, link_pair[end][0]) in parsed_linkage \
                        and contig_name(link_pair[end][0]) not in self_circular:
                    contig2join[end].append(link_pair[end][0])
                    contig_checked[end].append(contig_name(link_pair[end][0]))
                    contig2join_reason[contig][contig_name(link_pair[end][0])] = 'other_end_is_extendable'
                elif is_ok_to_add(link_pair[end][0], contig) and (end, link_pair[end][0]) in parsed_linkage \
                        and contig_name(link_pair[end][0]) not in self_circular:
                    contig2join[end].append(link_pair[end][0])
                    contig_checked[end].append(contig_name(link_pair[end][0]))
                    contig2join_reason[contig][contig_name(link_pair[end][0])] = 'is_ok_to_add'
                else:
                    pass
            else:
                pass
        elif end in two_paths_end:
            if link_pair[end][0].rsplit('_', 1)[0] != link_pair[end][1].rsplit('_', 1)[0]:
                if are_equal_paths(link_pair[end]):
                    if contig_name(the_dominant_one(link_pair[end])) not in self_circular:
                        contig2join[end].append(the_dominant_one(link_pair[end]))
                        contig_checked[end].append(contig_name((link_pair[end][0])))
                        contig_checked[end].append(contig_name((link_pair[end][1])))
                        contig2join_reason[contig][contig_name(the_dominant_one(link_pair[end]))] = 'are_equal_paths'
                    else:
                        pass
                else:
                    if cov[contig_name(link_pair[end][0])] + cov[contig_name(link_pair[end][1])] >= cov[contig] * 0.5:
                        # 0.5 is ok, too big will get much fewer "Extended circular" ones.
                        if contig_name(the_better_one(link_pair[end], contig)) not in self_circular:
                            contig2join[end].append(the_better_one(link_pair[end], contig))
                            contig_checked[end].append(contig_name((link_pair[end][0])))
                            contig_checked[end].append(contig_name((link_pair[end][1])))
                            contig2join_reason[contig][contig_name(the_better_one(link_pair[end], contig))] = 'the_better_one'
                        else:
                            pass
                    else:
                        pass
            else:
                pass
        else:
            pass
    else:
        target = get_target(contig2join[end][-1])
        if target in one_path_end:
            if not_checked(link_pair[target], contig_checked[end]):
                if other_end_is_extendable(link_pair[target][0], contig) and (target, link_pair[target][0]) in parsed_linkage \
                        and contig_name(link_pair[target][0]) not in self_circular:
                    contig2join[end].append(link_pair[target][0])
                    contig_checked[end].append(contig_name(link_pair[target][0]))
                    contig2join_reason[contig][contig_name(link_pair[target][0])] = 'other_end_is_extendable'
                elif is_ok_to_add(link_pair[target][0], contig) and (target, link_pair[target][0]) in parsed_linkage \
                        and contig_name(link_pair[target][0]) not in self_circular:
                    contig2join[end].append(link_pair[target][0])
                    contig_checked[end].append(contig_name(link_pair[target][0]))
                    contig2join_reason[contig][contig_name(link_pair[target][0])] = 'is_ok_to_add'
                else:
                    pass
            else:
                if could_circulate(target, contig, direction):
                    path_circular.add(contig)
                    path_circular_end.add(contig + '_' + direction)
                    contig2join_reason[contig][contig_name(target)] = 'could_circulate'
        elif target in two_paths_end:
            if link_pair[target][0].rsplit('_', 1)[0] != link_pair[target][1].rsplit('_', 1)[0]:
                if cov[contig_name(target)] >= 1.9 * cov[contig]:
                    pass
                else:
                    if not_checked(link_pair[target], contig_checked[end]):
                        if are_equal_paths(link_pair[target]):
                            if contig_name(the_dominant_one(link_pair[target])) not in self_circular:
                                contig2join[end].append(the_dominant_one(link_pair[target]))
                                contig_checked[end].append(contig_name(link_pair[target][0]))
                                contig_checked[end].append(contig_name(link_pair[target][1]))
                                contig2join_reason[contig][contig_name(the_dominant_one(link_pair[target]))] = 'are_equal_paths'
                            else:
                                pass
                        else:
                            if cov[contig_name(link_pair[target][0])] + cov[contig_name(link_pair[target][1])] >= cov[contig] * 0.5:
                                # 0.5 is ok, too big will get much fewer "Extended circular" ones.
                                if contig_name(the_better_one(link_pair[target], contig)) not in self_circular:
                                    contig2join[end].append(the_better_one(link_pair[target], contig))
                                    contig_checked[end].append(contig_name((link_pair[target][0])))
                                    contig_checked[end].append(contig_name((link_pair[target][1])))
                                    contig2join_reason[contig][contig_name(the_better_one(link_pair[target], contig))] = 'the_better_one'
                                else:
                                    pass
                            else:
                                pass
                    else:
                        if could_circulate(target, contig, direction):
                            path_circular.add(contig)
                            path_circular_end.add(contig + '_' + direction)
                            contig2join_reason[contig][contig_name(target)] = 'could_circulate'
                        else:
                            pass
            else:
                pass
        else:
            pass

    return a < len(contig2join[end])


def join_seqs(contig):
    """
    get the join order of the sequences in a give path
    """

    global order_all, added_to_contig

    order_all[contig] = []
    left = contig + '_L'
    right = contig + '_R'
    added_to_contig[contig] = []

    if contig not in path_circular:
        if left in contig2join.keys() and right in contig2join.keys():
            for item in contig2join[left][::-1]:  # the order of contigs added into the left direction should be reversed
                order_all[contig].append(item)
                added_to_contig[contig].append(contig_name(item))

            order_all[contig].append(contig)  # the query contig itself should be added to the path as well

            for item in contig2join[right]:
                if contig_name(item) not in added_to_contig[contig]:
                    order_all[contig].append(item)

        elif left in contig2join.keys() and right not in contig2join.keys():
            for item in contig2join[left][::-1]:  # the order of contigs added into the left direction should be reversed
                order_all[contig].append(item)

            order_all[contig].append(contig)  # the query contig itself should be added to the path as well

        else:
            order_all[contig].append(contig)  # the query contig itself should be added to the path as well

            for item in contig2join[right]:
                order_all[contig].append(item)

    else:
        if left in path_circular_end and right in path_circular_end:
            for item in contig2join[left][::-1]:  # the order of contigs added into the left direction should be reversed
                order_all[contig].append(item)

            order_all[contig].append(contig)  # the query contig itself should be added to the path as well

        elif left in path_circular_end and right not in path_circular_end:
            for item in contig2join[left][::-1]:  # the order of contigs added into the left direction should be reversed
                order_all[contig].append(item)

            order_all[contig].append(contig)  # the query contig itself should be added to the path as well

        elif right in path_circular_end and left not in path_circular_end:
            order_all[contig].append(contig)  # the query contig itself should be added to the path as well

            for item in contig2join[right]:
                order_all[contig].append(item)


def retrieve(contig):
    """
    to retrieve and save all the contigs in the joining path of a query
    """
    if contig in order_all.keys() and len(order_all[contig]) > 0:
        out = open('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(contig), 'w')
        added = []
        for item in order_all[contig]:
            if contig_name(item) not in added:
                out.write('>' + contig_name(item) + '\n')
                out.write(header2seq[contig_name(item)] + '\n')
                added.append(contig_name(item))
            else:
                pass
        out.close()
    else:
        pass


def count_seq(fasta_file):
    """
    calculate the number of seqs in a fasta file
    """
    seq_num = 0
    a = open(fasta_file, 'r')
    for line in a.readlines():
        if line.startswith('>'):
            seq_num += 1
    a.close()
    return seq_num


def count_len(fasta_file):
    """
    calculate the length of sequences in a fasta file
    """
    seq_len = 0
    a = open(fasta_file, 'r')
    for line in a.readlines():
        if not line.startswith('>'):
            seq_len += len(line.strip())
    a.close()
    return seq_len


def summary_fasta(fasta_file):
    """
    summary basic information of a fasta file
    """
    summary_file = open('{0}.summary.txt'.format(fasta_file), 'w')
    if 'self_circular' in fasta_file:
        summary_file_headers = ['SeqID', 'Length', 'Coverage', 'GC', 'Ns', 'DTR_length', '\n']
        summary_file.write('\t'.join(summary_file_headers[:]))
        summary_file.flush()
    else:
        summary_file_headers = ['SeqID', 'Length', 'Coverage', 'GC', 'Ns', '\n']
        summary_file.write('\t'.join(summary_file_headers[:]))
        summary_file.flush()

    with open('{0}'.format(fasta_file), 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            header = str(record.id).strip()
            seq = str(record.seq)
            ns = seq.count('N')
            if header.split('_self')[0] in self_circular:
                length = args.maxk - 1 if args.assembler == "idba" else args.maxk
                sequence_stats = [header, str(len(seq)), str(cov[header.split('_self')[0]]), str(round(GC(seq), 3)), str(ns), str(length), '\n']
                summary_file.write('\t'.join(sequence_stats[:]))
            elif header.split('_self')[0] in self_circular_non_expected_overlap.keys():
                sequence_stats = [header, str(len(seq)), str(cov[header.split('_self')[0]]), str(round(GC(seq), 3)), str(ns),
                                  str(self_circular_non_expected_overlap[header.split('_self')[0]]), '\n']
                summary_file.write('\t'.join(sequence_stats[:]))
            else:
                sequence_stats = [header, str(len(seq)), str(cov[header.split('_extended')[0]]), str(round(GC(seq), 3)), str(ns), '\n']
                summary_file.write('\t'.join(sequence_stats[:]))
    f.close()


def get_joined_seqs(fasta_file):
    joined_seqs = []
    a = open(fasta_file, 'r')
    for line in a.readlines():
        if line.startswith('>'):
            joined_seqs.append(line.strip().split(' ')[0][1:])
        else:
            pass
    for item in joined_seqs:
        all_joined_query.add(item)
    return ','.join(joined_seqs[:])


def summarize(contig):
    """
    summary the retrieved contigs and joined information of each query
    """
    if contig not in is_subset_of.keys():
        b = count_seq('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(contig))  # number of retrieved contigs
        c = count_len('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(contig))  # total length of retrieved contigs
        d = count_len('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(contig))  # total length after joining
        e = get_joined_seqs('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(contig))
        if contig in extended_circular_query:
            return '\t'.join([str(b), e, str(c), str(d), str(d - header2len[contig]), 'Extended_circular'])
        elif contig in extended_partial_query:
            return '\t'.join([str(b), e, str(c), str(d), str(d - header2len[contig]), 'Extended_partial'])
    else:
        if is_subset_of[contig] in is_subset_of.keys():
            item = is_subset_of[is_subset_of[contig]]
            b = count_seq('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(item))  # number of retrieved contigs
            c = count_len('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(item))  # total length of retrieved contigs
            d = count_len('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(item))  # total length after joining
            e = get_joined_seqs('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(item))
            if contig in extended_circular_query:
                return '\t'.join([str(b), e, str(c), str(d), str(d - header2len[contig]), 'Extended_circular'])
            elif contig in extended_partial_query:
                return '\t'.join([str(b), e, str(c), str(d), str(d - header2len[contig]), 'Extended_partial'])
        else:
            item = is_subset_of[contig]
            b = count_seq('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(item))  # number of retrieved contigs
            c = count_len('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(item))  # total length of retrieved contigs
            d = count_len('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(item))  # total length after joining
            e = get_joined_seqs('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(item))
            if contig in extended_circular_query:
                return '\t'.join([str(b), e, str(c), str(d), str(d - header2len[contig]), 'Extended_circular'])
            elif contig in extended_partial_query:
                return '\t'.join([str(b), e, str(c), str(d), str(d - header2len[contig]), 'Extended_partial'])


def get_direction(item):
    """
    get the direction in a joining path
    """
    if item.endswith('rc'):
        return 'reverse'
    else:
        return 'forward'


def total_length(contig_list):
    """
    get the total length of all sequences in a joining path before overlap removing
    """
    total = 0
    for item in contig_list:
        total += header2len[contig_name(item)]
    return total


def main():
    ##
    # get information from the input files and parameters and save information
    # get the name of the query fasta file
    if '/' in args.query:
        query_name = '{0}'.format(args.query).rsplit('/', 1)[1]
    else:
        query_name = '{0}'.format(args.query)

    # get the name of the whole contigs fasta file
    if '/' in args.fasta:
        fasta_name = '{0}'.format(args.fasta).rsplit('/', 1)[1]
    else:
        fasta_name = '{0}'.format(args.fasta)

    # folder of output
    if not args.output:
        working_dir = '{0}_COBRA'.format(query_name)
    else:
        working_dir = '{0}'.format(args.output)

    # checking if output folder exists
    if os.path.exists('{0}'.format(working_dir)):
        print('Output folder <{0}> exists, please check.'.format(working_dir))
        exit()
    else:
        os.mkdir('{0}'.format(working_dir))

    # determine the length of overlap based on assembler and the largest kmer size
    if args.assembler == "idba":
        length = args.maxk - 1
    else:
        length = args.maxk

    # write input files information to log file
    log = open('{0}/log'.format(working_dir), 'w')  # log file
    log.write('1. INPUT INFORMATION' + '\n')
    log.flush()

    if args.assembler == 'idba':
        parameters = ['# Assembler: IDBA_UD',
                      '# Min-kmer: ' + str(args.mink).strip(),
                      '# Max-kmer: ' + str(args.maxk).strip(),
                      '# Overlap length: ' + str(length) + ' bp',
                      '# Read mapping max mismatches for contig linkage: ' + str(args.linkage_mismatch),
                      '# Query contigs: ' + os.path.abspath(args.query),
                      '# Whole contig set: ' + os.path.abspath(args.fasta),
                      '# Mapping file: ' + os.path.abspath(args.mapping),
                      '# Coverage file: ' + os.path.abspath(args.coverage),
                      '# Output folder: ' + os.path.abspath(working_dir), '\n']
    elif args.assembler == 'metaspades':
        parameters = ['# Assembler: metaSPAdes',
                      '# Min-kmer: ' + str(args.mink).strip(),
                      '# Max-kmer: ' + str(args.maxk).strip(),
                      '# Overlap length: ' + str(length) + ' bp',
                      '# Read mapping max mismatches for contig linkage: ' + str(args.linkage_mismatch),
                      '# Query contigs: ' + os.path.abspath(args.query),
                      '# Whole contig set: ' + os.path.abspath(args.fasta),
                      '# Mapping file: ' + os.path.abspath(args.mapping),
                      '# Coverage file: ' + os.path.abspath(args.coverage),
                      '# Output folder: ' + os.path.abspath(working_dir), '\n']
    else:
        parameters = ['# Assembler: MEGAHIT',
                      '# Min-kmer: ' + str(args.mink).strip(),
                      '# Max-kmer: ' + str(args.maxk).strip(),
                      '# Overlap length: ' + str(length) + ' bp',
                      '# Read mapping max mismatches for contig linkage: ' + str(args.linkage_mismatch),
                      '# Query contigs: ' + os.path.abspath(args.query),
                      '# Whole contig set: ' + os.path.abspath(args.fasta),
                      '# Mapping file: ' + os.path.abspath(args.mapping),
                      '# Coverage file: ' + os.path.abspath(args.coverage),
                      '# Output folder: ' + os.path.abspath(working_dir), '\n']

    log.write('\n'.join(parameters[:]))
    log.write('2. PROCESSING STEPS' + '\n')
    log.flush()

    ##
    # import the whole contigs and save their end sequences
    log_info('[01/23]', 'Reading contigs and getting the contig end sequences. ', '', log)
    # header2seq = {}
    # header2len = {}
    gc = {}
    L = {}
    R = {}
    Lrc = {}
    Rrc = {}

    with open('{0}'.format(args.fasta), 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            header = str(record.id).strip()
            seq = str(record.seq)
            header2seq[header] = seq
            gc[header] = str(round(GC(seq), 3))
            header2len[header] = len(seq)
            L[header + '_L'] = seq[0:length]  # the first x bp of left end
            Lrc[header + '_Lrc'] = reverse_complement(seq[0:length])  # the reverse sequence of first x bp of left end
            R[header + '_R'] = seq[-length:]  # the first x bp of right end
            Rrc[header + '_Rrc'] = reverse_complement(seq[-length:])  # the reverse sequence of first x bp of right end

    log.write('A total of {0} contigs were imported.'.format(len(header2seq.keys())) + '\n')

    ##
    # get potential joins
    log_info('[02/23]', 'Getting shared contig ends.', '\n', log)

    # link_pair = {}  # used to save all overlaps between ends

    d_L = defaultdict(set)
    d_Lrc = defaultdict(set)
    d_R = defaultdict(set)
    d_Rrc = defaultdict(set)

    for k, v in L.items():  # save header2seq in dictionary with seqs as keys
        d_L[v].add(k)
    for k, v in Lrc.items():
        d_Lrc[v].add(k)
    for k, v in R.items():
        d_R[v].add(k)
    for k, v in Rrc.items():
        d_Rrc[v].add(k)

    d_L_d_Lrc_shared = set(d_L.keys()).intersection(set(d_Lrc.keys()))
    # get the shared seqs between direction pairs (L/Lrc, Lrc/L, L/R, R/L, R/Rrc, Rrc/R, Lrc/Rrc, Rrc/Lrc)
    d_L_d_R_shared = set(d_L.keys()).intersection(set(d_R.keys()))
    # the d_R_d_L_shared will be included below
    d_R_d_Rrc_shared = set(d_R.keys()).intersection(set(d_Rrc.keys()))
    d_Rrc_d_Lrc_shared = set(d_Rrc.keys()).intersection(set(d_Lrc.keys()))

    ##
    # get link_pair between ends
    for end in d_L_d_Lrc_shared:
        for left in d_L[end]:  # left is a seq name
            for left_rc in d_Lrc[end]:  # left_rc is a seq name
                if left not in link_pair.keys():
                    link_pair[left] = [left_rc]
                else:
                    link_pair[left].append(left_rc)
        for left_rc in d_Lrc[end]:
            for left in d_L[end]:
                if left_rc not in link_pair.keys():
                    link_pair[left_rc] = [left]
                else:
                    link_pair[left_rc].append(left)

    for end in d_L_d_R_shared:
        for left in d_L[end]:
            for right in d_R[end]:
                if left not in link_pair.keys():
                    link_pair[left] = [right]
                else:
                    link_pair[left].append(right)
        for right in d_R[end]:
            for left in d_L[end]:
                if right not in link_pair.keys():
                    link_pair[right] = [left]
                else:
                    link_pair[right].append(left)

    for end in d_R_d_Rrc_shared:
        for right in d_R[end]:
            for right_rc in d_Rrc[end]:
                if right not in link_pair.keys():
                    link_pair[right] = [right_rc]
                else:
                    link_pair[right].append(right_rc)
        for right_rc in d_Rrc[end]:
            for right in d_R[end]:
                if right_rc not in link_pair.keys():
                    link_pair[right_rc] = [right]
                else:
                    link_pair[right_rc].append(right)

    for end in d_Rrc_d_Lrc_shared:
        for right_rc in d_Rrc[end]:
            for left_rc in d_Lrc[end]:
                if right_rc not in link_pair.keys():
                    link_pair[right_rc] = [left_rc]
                else:
                    link_pair[right_rc].append(left_rc)
        for left_rc in d_Lrc[end]:
            for right_rc in d_Rrc[end]:
                if left_rc not in link_pair.keys():
                    link_pair[left_rc] = [right_rc]
                else:
                    link_pair[left_rc].append(right_rc)

    ##
    # save all paired links to a file
    log_info('[03/23]', 'Writing contig end joining pairs.', '\n', log)

    p = open('{0}/COBRA_end_joining_pairs.txt'.format(working_dir), 'w')

    for item in link_pair.keys():
        for point in link_pair[item]:
            p.write(item + '\t' + point + '\n')  # print link pairs into a file for check if interested

        if len(link_pair[item]) == 1:  # and len(link_pair[link_pair[item][0]]) > 1:
            one_path_end.append(item)  # add one joining end to a list, its pair may have one or more joins
        elif len(link_pair[item]) == 2 and len(link_pair[link_pair[item][0]]) == 1 and len(
                link_pair[link_pair[item][1]]) == 1:
            two_paths_end.append(item)  # add two joining end to a list, each of its pairs should only have one join
        else:
            pass
    p.close()

    ##
    # read and save the coverage of all contigs
    log_info('[04/23]', 'Getting contig coverage information.', '\n', log)
    coverage = open('{0}'.format(args.coverage), 'r')
    for line in coverage.readlines():
        line = line.strip().split('\t')
        cov[line[0]] = round(float(line[1]), 3)
    coverage.close()

    if len(cov.keys()) < len(header2seq.keys()):
        print('Some contigs do not have coverage information. Please check. COBRA exits.')
        exit()
    else:
        pass

    ##
    # open the query file and save the information

    log_info('[05/23]', 'Getting query contig list. ', '', log)
    query_set = set()
    orphan_end_query = set()
    non_orphan_end_query = set()

    with open('{0}'.format(args.query), 'r') as query_file:
        if determine_file_format(args.query) == 'fasta':  # if the query file is in fasta format
            for record in SeqIO.parse(query_file, 'fasta'):
                header = str(record.id).strip()
                if header in header2seq.keys():  # some queries may not in the whole assembly, should be rare though.
                    query_set.add(header)
                else:
                    print('Query {0} is not in your whole contig fasta file, please check!'.format(header), flush=True)
        else:  # if the query file is in text format
            for line in query_file:
                header = line.strip().split(' ')[0]
                if header in header2seq.keys():  # some queries may not in the whole assembly, should be rare though.
                    query_set.add(header)
                else:
                    print('Query {0} is not in your whole contig fasta file, please check!'.format(header), flush=True)

    # distinguish orphan_end_query and non_orphan_end_query:

    for header in query_set:
        if header + '_L' not in link_pair.keys() and header + '_R' not in link_pair.keys():
            orphan_end_query.add(header)
        else:
            non_orphan_end_query.add(header)

    #
    log.write('A total of {0} query contigs were imported.'.format(len(query_set)) + '\n')
    log.flush()

    ##
    # get the linkage of contigs based on paired-end reads mapping
    log_info('[06/23]', 'Getting contig linkage based on sam/bam. Be patient, this may take long.', '\n', log)
    linkage = defaultdict(set)  # Initialize a defaultdict to store linked contigs
    contig_spanned_by_PE_reads = {}  # Initialize a dictionary to store paired-end reads spanning contigs

    for contig in orphan_end_query:
        contig_spanned_by_PE_reads[contig] = defaultdict(list)  # Create a defaultdict(list) for each contig in header2seq

    with pysam.AlignmentFile('{0}'.format(args.mapping), 'rb') as map_file:
        for line in map_file:
            if not line.is_unmapped and line.get_tag("NM") <= args.linkage_mismatch:
                # mismatch should not be more than the defined threshold
                if line.reference_name != line.next_reference_name:
                    # Check if the read and its mate map to different contigs
                    if header2len[line.reference_name] > 1000:
                        # If the contig length is greater than 1000, determine if the read maps to the left or right end
                        if line.reference_start <= 500:
                            linkage[line.query_name].add(line.reference_name + '_L')  # left end
                        else:
                            linkage[line.query_name].add(line.reference_name + '_R')  # right end
                    else:
                        # If the contig length is 1000 or less, add both the left and right ends to the linkage
                        linkage[line.query_name].add(line.reference_name + '_L')
                        linkage[line.query_name].add(line.reference_name + '_R')
                else:
                    # If the read and its mate map to the same contig, store the read mapped position (start)
                    if line.reference_name in orphan_end_query:
                        if line.reference_start <= 500 or header2len[line.reference_name] - line.reference_start <= 500:
                            contig_spanned_by_PE_reads[line.reference_name][line.query_name].append(
                                line.reference_start)
                        else:
                            pass
                    else:
                        pass
            else:
                pass

    #
    log_info('[07/23]', 'Parsing the linkage information.', '\n', log)

    for read in linkage.keys():
        if len(linkage[read]) >= 2:  # Process only reads linked to at least two contigs
            # for item in linkage[read]:
            for item, item_1 in itertools.combinations(linkage[read], 2):
                # Generate unique pairs of linked contigs for the current read using itertools.combinations
                if item.rsplit('_', 1)[1] != item_1.rsplit('_', 1)[1]:
                    # If the contigs have different ends (_L or _R), add the combinations to the parsed_linkage
                    parsed_linkage.add((item, item_1))
                    parsed_linkage.add((item_1, item))
                    parsed_linkage.add((item + 'rc', item_1 + 'rc'))
                    parsed_linkage.add((item_1 + 'rc', item + 'rc'))
                else:
                    # If the contigs have the same ends, add the combinations with reverse-complement (_rc) to the parsed_linkage
                    parsed_linkage.add((item, item_1 + 'rc'))
                    parsed_linkage.add((item_1 + 'rc', item))
                    parsed_linkage.add((item + 'rc', item_1))
                    parsed_linkage.add((item_1, item + 'rc'))
        else:
            pass

    linkage = None  # remove it as no longer used later

    #
    parsed_contig_spanned_by_PE_reads = set()  # Initialize a set to store the contig spanned by paired-end reads
    for contig in orphan_end_query:
        for PE in contig_spanned_by_PE_reads[contig].keys():  # Check if the count is 0 and the contig has exactly two paired-end reads
            if len(contig_spanned_by_PE_reads[contig][PE]) == 2:
                # Check if the absolute difference between the positions of the two paired-end reads is greater than or equal to
                # the length of contig minus 1000 bp
                if abs(contig_spanned_by_PE_reads[contig][PE][0] - contig_spanned_by_PE_reads[contig][PE][1]) >= \
                        header2len[contig] - 1000:
                    parsed_contig_spanned_by_PE_reads.add(contig)
                else:
                    pass
            else:
                pass

    ##
    #
    log_info('[08/23]', 'Detecting self_circular contigs. ', '\n', log)

    for contig in non_orphan_end_query:
        detect_self_circular(contig)

    debug = open('{0}/debug.txt'.format(working_dir), 'w')

    # for debug
    print('self_circular', file=debug, flush=True)
    print(self_circular, file=debug, flush=True)

    # orphan end queries info
    # determine potential self_circular contigs from contigs with orphan end

    min_over_len = 0
    if args.assembler == 'idba':
        min_over_len = args.mink - 1
    else:
        min_over_len = args.mink

    for contig in orphan_end_query:  # determine if there is DTR for those query with orphan ends, if yes, assign as self_circular as well
        if contig in parsed_contig_spanned_by_PE_reads:
            sequence = header2seq[contig]
            end_part = sequence[-min_over_len:]
            if sequence.count(end_part) == 2:
                expected_end = sequence.split(end_part)[0] + end_part
                if sequence.endswith(expected_end):
                    self_circular_non_expected_overlap[contig] = len(expected_end)
                else:
                    pass
            else:
                pass
        else:
            pass

    for contig in self_circular_non_expected_overlap.keys():
        orphan_end_query.remove(contig)

    # debug
    print(self_circular_non_expected_overlap, file=debug, flush=True)

    ##
    # walk the joins
    log_info('[09/23]', 'Detecting joins of contigs. ', '', log)

    for contig in query_set:
        contig2join[contig + '_L'] = []
        contig2join[contig + '_R'] = []
        contig_checked[contig + '_L'] = []
        contig_checked[contig + '_R'] = []
        contig2join_reason[contig] = {}
        contig2join_reason[contig][contig] = 'query'

    percentage = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    finished_join_walker = 0
    total_join_walker = len(query_set) - len(orphan_end_query) - len(self_circular)

    #
    for contig in query_set:
        if contig not in list(orphan_end_query) + list(self_circular):

            # extend each contig from both directions
            while True:
                result_L = join_walker(contig, 'L')
                if not result_L:
                    break

            while True:
                result_R = join_walker(contig, 'R')
                if not result_R:
                    break

            # calculate the percentage of query contigs have been checked for extension
            finished_join_walker += 1
            finished_percentage = int(finished_join_walker / total_join_walker * 100)

            if finished_percentage in percentage:
                log.write(str(finished_percentage) + '%, ')
                log.flush()
                percentage.remove(finished_percentage)
            else:
                pass

        else:
            pass

    ##
    # save the potential joining paths
    log.write('100% finished.' + '\n')
    log_info('[10/23]', 'Saving potential joining paths.', '\n', log)

    with open('{0}/COBRA_potential_joining_paths.txt'.format(working_dir), 'w') as results:
        for item in contig2join.keys():
            if contig_name(item) in self_circular:
                if item.endswith('_L'):
                    results.write(item + '\t' + "['" + contig_name(item) + '_R' + "']" + '\n')
                    results.flush()
                else:
                    results.write(item + '\t' + "['" + contig_name(item) + '_L' + "']" + '\n')
                    results.flush()
            else:
                results.write(item + '\t' + str(contig2join[item]) + '\n')
                results.flush()

    ##
    # get the fail_to_join contigs, but not due to orphan end
    failed_join_list = []
    for contig in query_set:
        if contig + '_L' in link_pair.keys() or contig + '_R' in link_pair.keys():
            if len(contig2join[contig + '_L']) == 0 and len(contig2join[contig + '_R']) == 0 \
                    and contig not in self_circular:
                failed_join_list.append(contig)
                print(contig, len(contig2join[contig + '_L']), len(contig2join[contig + '_R']), contig not in self_circular)
            else:
                pass
        else:
            pass

    ##
    # get the joining paths
    log_info('[11/23]', 'Checking for invalid joining: sharing queries.', '\n', log)
    contig2assembly = {}
    for item in contig2join.keys():
        contig = contig_name(item)
        if contig + '_L' in path_circular_end and contig + '_R' not in path_circular_end:
            if contig not in contig2assembly.keys():
                contig2assembly[contig] = set()
                contig2assembly[contig].add(contig)
                for point in contig2join[contig + '_L']:
                    contig2assembly[contig].add(contig_name(point))
            else:
                for point in contig2join[contig + '_L']:
                    contig2assembly[contig].add(contig_name(point))
        elif contig + '_L' not in path_circular_end and contig + '_R' in path_circular_end:
            if contig not in contig2assembly.keys():
                contig2assembly[contig] = set()
                contig2assembly[contig].add(contig)
                for point in contig2join[contig + '_R']:
                    contig2assembly[contig].add(contig_name(point))
            else:
                for point in contig2join[contig + '_R']:
                    contig2assembly[contig].add(contig_name(point))
        else:
            if contig not in contig2assembly.keys():
                contig2assembly[contig] = set()
                contig2assembly[contig].add(contig)
                for point in contig2join[item]:
                    contig2assembly[contig].add(contig_name(point))
            else:
                for point in contig2join[item]:
                    contig2assembly[contig].add(contig_name(point))

    # for debug
    print('path_circular', file=debug, flush=True)
    print(path_circular, file=debug, flush=True)
    print('1failed_join_list', file=debug, flush=True)
    print(failed_join_list, file=debug, flush=True)
    print('contig2assembly', file=debug, flush=True)
    print(contig2assembly, file=debug, flush=True)

    ##
    # find the redundant joining paths
    redundant = set()
    is_same_as = {}
    for contig in contig2assembly.keys():
        for contig_1 in contig2assembly.keys():
            if contig != contig_1 and contig2assembly[contig].issubset(contig2assembly[contig_1]):
                if contig2assembly[contig] != contig2assembly[contig_1]:
                    if contig in path_circular and contig_1 not in path_circular:
                        # in this rare case, should use contig, not contig_1, thus contig_1 is_subset_of contig
                        if contig_1 not in contig2assembly[contig]:
                            failed_join_list.append(contig)
                            failed_join_list.append(contig_1)
                            path_circular.remove(contig)
                        else:
                            redundant.add(contig_1)
                            path_circular.add(contig_1)
                            if contig not in is_subset_of.keys():
                                is_subset_of[contig_1] = contig
                            else:
                                is_subset_of[contig_1] = is_subset_of[contig]
                    elif contig not in path_circular and contig_1 in path_circular:
                        # in this case, should use contig_1, not contig, thus contig is_subset_of contig_1
                        path_circular.add(contig)
                        redundant.add(contig)
                        if contig_1 not in is_subset_of.keys():
                            is_subset_of[contig] = contig_1
                        else:
                            is_subset_of[contig] = is_subset_of[contig_1]
                    elif contig in path_circular and contig_1 in path_circular:
                        # in this rare case, should consider them as the same
                        contig2assembly[contig_1] = contig2assembly[contig]
                        if contig not in is_same_as.keys():
                            is_same_as[contig] = set()
                            is_same_as[contig].add(contig_1)
                        else:
                            is_same_as[contig].add(contig_1)
                    else:  # in this case, contig is_subset_of contig_1, which is different from the first case above.
                        redundant.add(contig)
                        if contig_1 not in is_subset_of.keys():
                            is_subset_of[contig] = contig_1
                        else:
                            is_subset_of[contig] = is_subset_of[contig_1]
                else:
                    if contig not in is_same_as.keys():
                        is_same_as[contig] = set()
                        is_same_as[contig].add(contig_1)
                    else:
                        is_same_as[contig].add(contig_1)
            else:
                pass

    # for debug
    print('2failed_join_list', file=debug, flush=True)
    print(failed_join_list, file=debug, flush=True)

    for contig in is_subset_of.keys():
        if is_subset_of[contig] in failed_join_list:
            failed_join_list.append(contig)

    # for debug
    print('3failed_join_list', file=debug, flush=True)
    print(failed_join_list, file=debug, flush=True)

    ##
    # remove the queries in multiple paths
    all = []  # all contigs in all the paths checked
    same_path = []  # two are the same path
    contig_shared_by_paths = set()  # the queries in multiple non-unique paths

    for contig in contig2assembly.keys():
        if contig not in redundant and contig not in failed_join_list:
            if contig2assembly[contig] not in same_path:
                same_path.append(contig2assembly[contig])
                for item in contig2assembly[contig]:
                    if item in query_set:
                        all.append(item)
                    else:
                        pass
            else:
                pass
        else:
            pass

    # for debug
    print('same_path', file=debug, flush=True)
    print(same_path, file=debug, flush=True)

    for contig in set(all):
        if all.count(contig) > 1:
            for contig_1 in contig2assembly.keys():
                if contig_1 not in redundant and contig not in failed_join_list:
                    if contig in contig2assembly[contig_1]:
                        contig_shared_by_paths.add(contig_1)
                    else:
                        pass
                else:
                    pass
        else:
            pass

    # for debug
    print('contig_shared_by_paths', file=debug, flush=True)
    print(contig_shared_by_paths, file=debug, flush=True)

    for contig in contig_shared_by_paths:
        if contig in contig2assembly.keys():
            del contig2assembly[contig]
            failed_join_list.append(contig)
        else:
            pass

    # for debug
    print('4failed_join_list', file=debug, flush=True)
    print(failed_join_list, file=debug, flush=True)

    for contig in contig_shared_by_paths:
        for item in is_subset_of.keys():
            if contig == is_subset_of[item]:
                del contig2assembly[item]
                failed_join_list.append(item)

                for item_1 in is_subset_of.keys():
                    if is_subset_of[item_1] == item:
                        del contig2assembly[item_1]
                        failed_join_list.append(item_1)
            else:
                pass

    # for debug
    print('5failed_join_list', file=debug, flush=True)
    print(failed_join_list, file=debug, flush=True)

    ##
    # determine the joining status of queries
    log_info('[12/23]', 'Getting initial joining status of each query contig.', '\n', log)

    for contig in contig2assembly.keys():
        if contig not in failed_join_list:
            if contig in redundant:
                if is_subset_of[contig] in path_circular:
                    extended_circular_query.add(contig)
                else:
                    extended_partial_query.add(contig)
            else:
                if contig in path_circular:
                    extended_circular_query.add(contig)
                else:
                    if contig not in failed_join_list and contig not in orphan_end_query and contig not in self_circular \
                            and contig not in self_circular_non_expected_overlap.keys():
                        extended_partial_query.add(contig)
                    else:
                        pass
        else:
            pass

    # for debug
    print('6failed_join_list', file=debug, flush=True)
    print(failed_join_list, file=debug, flush=True)

    ##
    # deal with cross-assignment queries
    log_info('[13/23]', 'Getting final joining status of each query contig.', '\n', log)
    for contig in failed_join_list:
        if contig in is_subset_of.keys() and is_subset_of[contig] in extended_circular_query:
            failed_join_list.remove(contig)
            extended_circular_query.add(contig)
        elif contig in is_subset_of.keys() and is_subset_of[contig] in extended_partial_query:
            failed_join_list.remove(contig)
            extended_partial_query.add(contig)
        elif contig in is_same_as.keys():
            for item in is_same_as[contig]:
                if item in extended_circular_query:
                    if contig in failed_join_list:
                        failed_join_list.remove(contig)
                        extended_circular_query.add(contig)
                elif item in extended_partial_query:
                    if contig in failed_join_list:
                        failed_join_list.remove(contig)
                        extended_partial_query.add(contig)
                else:
                    pass
        else:
            pass

    to_be_removed_from_epq = set()  # epq = extended_partial_query
    for contig in extended_partial_query:
        if contig in redundant:
            if is_subset_of[contig] in extended_circular_query:
                to_be_removed_from_epq.add(contig)
                extended_circular_query.add(contig)
            elif is_subset_of[contig] in failed_join_list:
                to_be_removed_from_epq.add(contig)
                failed_join_list.append(contig)
            else:
                pass
        elif contig in is_same_as.keys():
            for item in is_same_as[contig]:
                if item in extended_circular_query and contig in extended_partial_query:
                    to_be_removed_from_epq.add(contig)
                    extended_circular_query.add(contig)
                else:
                    pass
        else:
            pass

    for contig in to_be_removed_from_epq:
        extended_partial_query.remove(contig)

    # for is_same_as ones, get the redundant ones
    is_same_as_redundant = []
    for contig in is_same_as.keys():
        if contig not in is_same_as_redundant:
            for item in is_same_as[contig]:
                is_same_as_redundant.append(item)
            else:
                pass
        else:
            pass

    # for debug
    print('7failed_join_list', file=debug, flush=True)
    print(failed_join_list, file=debug, flush=True)

    # for debug
    print('path_circular', file=debug, flush=True)
    print(path_circular, file=debug, flush=True)
    print('redundant', file=debug, flush=True)
    print(redundant, file=debug, flush=True)
    print('is_subset_of', file=debug, flush=True)
    print(is_subset_of, file=debug, flush=True)
    print('is_same_as', file=debug, flush=True)
    print(is_same_as, file=debug, flush=True)

    #
    def all_contigs_in_the_path_are_good(contig_set):
        total = 0
        for contig in contig_set:
            if contig in failed_join_list:
                total += 1
            else:
                pass
        return total == 0

    ##
    # get the joining order of contigs
    log_info('[14/23]', 'Getting the joining order of contigs.', '\n', log)
    # order_all = {}
    # added_to_contig = {}

    for contig in contig2assembly.keys():
        # only those contigs left in contig2assembly after filtering
        # (see above "# remove the queries in multiple paths") will be
        # checked for join paths (join_seqs) to get order_all
        if len(contig2assembly[contig]) > 1 and contig not in failed_join_list and contig not in redundant:
            if all_contigs_in_the_path_are_good(contig2assembly[contig]):
                join_seqs(contig)
            else:
                for item in contig2assembly[contig]:
                    if item in extended_circular_query:
                        extended_circular_query.remove(item)
                        failed_join_list.append(item)
                    elif item in extended_partial_query:
                        extended_partial_query.remove(item)
                        failed_join_list.append(item)
                    else:
                        pass
        else:
            pass

    ##
    # get retrieved sequences
    log_info('[15/23]', 'Getting retrieved contigs.', '\n', log)
    os.chdir('{0}'.format(working_dir))
    os.mkdir('COBRA_retrieved_for_joining')
    retrieved = []
    for contig in order_all.keys():
        retrieve(contig)
        retrieved.append(contig)

    # for debug
    print('retrieved', file=debug, flush=True)
    print(retrieved, file=debug, flush=True)

    ##
    # writing joined sequences
    log_info('[16/23]', 'Saving joined seqeuences.', '\n', log)
    header2joined_seq = {}
    contig2extended_status = {}
    for contig in retrieved:
        a = open('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(contig), 'w')
        last = ''
        # print header regarding the joining status
        if contig in path_circular:
            a.write('>' + contig + '_extended_circular' + '\n')
            header2joined_seq[contig + '_extended_circular'] = ''
            contig2extended_status[contig] = contig + '_extended_circular'
        else:
            a.write('>' + contig + '_extended_partial' + '\n')
            header2joined_seq[contig + '_extended_partial'] = ''
            contig2extended_status[contig] = contig + '_extended_partial'

        # print the sequences with their overlap removed
        for item in order_all[contig][:-1]:
            if item.endswith('_R') or item.endswith('_L'):
                if last == '':
                    a.write(header2seq[item.rsplit('_', 1)[0]][:-length])
                    last = header2seq[item.rsplit('_', 1)[0]][-length:]
                    header2joined_seq[contig2extended_status[contig]] += header2seq[item.rsplit('_', 1)[0]][:-length]
                else:
                    if header2seq[item.rsplit('_', 1)[0]][:length] == last:
                        a.write(header2seq[item.rsplit('_', 1)[0]][:-length])
                        last = header2seq[item.rsplit('_', 1)[0]][-length:]
                        header2joined_seq[contig2extended_status[contig]] += header2seq[item.rsplit('_', 1)[0]][
                                                                             :-length]
                    else:
                        pass
            elif item.endswith('rc'):
                if last == '':
                    a.write(reverse_complement(header2seq[item.rsplit('_', 1)[0]])[:-length])
                    last = reverse_complement(header2seq[item.rsplit('_', 1)[0]])[-length:]
                    header2joined_seq[contig2extended_status[contig]] += reverse_complement(
                        header2seq[item.rsplit('_', 1)[0]])[:-length]
                else:
                    if reverse_complement(header2seq[item.rsplit('_', 1)[0]])[:length] == last:
                        a.write(reverse_complement(header2seq[item.rsplit('_', 1)[0]])[:-length])
                        last = reverse_complement(header2seq[item.rsplit('_', 1)[0]])[-length:]
                        header2joined_seq[contig2extended_status[contig]] += reverse_complement(
                            header2seq[item.rsplit('_', 1)[0]])[:-length]
                    else:
                        pass
            else:
                if last == '':
                    a.write(header2seq[contig][:-length])
                    last = header2seq[contig][-length:]
                    header2joined_seq[contig2extended_status[contig]] += header2seq[contig][:-length]
                else:
                    if header2seq[contig][:length] == last:
                        a.write(header2seq[contig][:-length])
                        last = header2seq[contig][-length:]
                        header2joined_seq[contig2extended_status[contig]] += header2seq[contig][:-length]
                    else:
                        pass

        if order_all[contig][-1].endswith('rc'):
            a.write(reverse_complement(header2seq[order_all[contig][-1].rsplit('_', 1)[0]]) + '\n')
            header2joined_seq[contig2extended_status[contig]] += reverse_complement(
                header2seq[order_all[contig][-1].rsplit('_', 1)[0]])
        elif order_all[contig][-1].endswith('_R') or order_all[contig][-1].endswith('_L'):
            a.write(header2seq[order_all[contig][-1].rsplit('_', 1)[0]] + '\n')
            header2joined_seq[contig2extended_status[contig]] += header2seq[order_all[contig][-1].rsplit('_', 1)[0]]
        else:
            a.write(header2seq[order_all[contig][-1]] + '\n')
            header2joined_seq[contig2extended_status[contig]] += header2seq[order_all[contig][-1]]

        a.close()

    print('contig2extended_status', file=debug, flush=True)
    print(contig2extended_status, file=debug, flush=True)
    print('header2joined_seq', file=debug, flush=True)
    print(header2joined_seq.keys(), file=debug, flush=True)

    ##
    # Similar direct terminal repeats may lead to invalid joins
    log_info('[17/23]', 'Checking for invalid joining using BLASTn: close strains.', '\n', log)
    blastdb_1 = open('blastdb_1.fa', 'w')
    blastdb_2 = open('blastdb_2.fa', 'w')
    cobraSeq2len = {}

    for contig in retrieved:
        a = open('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(contig), 'r')
        for record in SeqIO.parse(a, "fasta"):
            header = str(record.id).strip()
            seq = str(record.seq)
            cobraSeq2len[header.split('_extended', 1)[0]] = len(seq)

            if len(seq) % 2 == 0:
                half = int(len(seq) / 2)
            else:
                half = int((len(seq) + 1) / 2)

            blastdb_1.write('>' + header + '_1' + '\n')
            blastdb_1.write(seq[:half] + '\n')
            blastdb_2.write('>' + header + '_2' + '\n')
            blastdb_2.write(seq[half:] + '\n')

        a.close()

    for contig in self_circular:
        cobraSeq2len[contig] = header2len[contig] - length

        if header2len[contig] % 2 == 0:
            half = int(header2len[contig] / 2)
        else:
            half = int((header2len[contig] + 1) / 2)

        blastdb_1.write('>' + contig + '_1' + '\n')
        blastdb_1.write(header2seq[contig][:half] + '\n')
        blastdb_2.write('>' + contig + '_2' + '\n')
        blastdb_2.write(header2seq[contig][half:] + '\n')

    for contig in self_circular_non_expected_overlap.keys():
        cobraSeq2len[contig] = header2len[contig] - self_circular_non_expected_overlap[contig]

        if header2len[contig] % 2 == 0:
            half = int(header2len[contig] / 2)
        else:
            half = int((header2len[contig] + 1) / 2)

        blastdb_1.write('>' + contig + '_1' + '\n')
        blastdb_1.write(header2seq[contig][:half] + '\n')
        blastdb_2.write('>' + contig + '_2' + '\n')
        blastdb_2.write(header2seq[contig][half:] + '\n')

    blastdb_1.close()
    blastdb_2.close()

    # make blastn database and run search if the database is not empty
    if os.path.getsize('blastdb_1.fa') == 0:
        print('no query was extended, exit! this is normal if you only provide few queries.', file=log, flush=True)
        exit()
    else:
        os.system('makeblastdb -in blastdb_1.fa -dbtype nucl')
        os.system('blastn -task blastn -db blastdb_1.fa -query blastdb_2.fa -out blastdb_2.vs.blastdb_1 -evalue 1e-10 '
                  '-outfmt 6 -perc_identity 70 -num_threads {0}'.format(args.threads))

    # parse the blastn results
    contig2TotLen = {}
    r = open('blastdb_2.vs.blastdb_1', 'r')
    for line in r.readlines():
        line = line.strip().split('\t')
        if line[0].rsplit('_', 1)[0] == line[1].rsplit('_', 1)[0] and line[0] != line[1]:
            if float(line[3]) >= 1000:
                if '_extended' in line[0]:
                    if line[0].split('_extended')[0] not in contig2TotLen.keys():
                        contig2TotLen[line[0].split('_extended')[0]] = float(line[3])
                    else:
                        contig2TotLen[line[0].split('_extended')[0]] += float(line[3])
                else:
                    if line[0].rsplit('_', 1)[0] not in contig2TotLen.keys():
                        contig2TotLen[line[0].rsplit('_', 1)[0]] = float(line[3])
                    else:
                        contig2TotLen[line[0].rsplit('_', 1)[0]] += float(line[3])
            else:
                pass
        else:
            pass
    r.close()

    # identify potential incorrect joins and remove them from corresponding category
    for contig in contig2TotLen.keys():
        if contig2TotLen[contig] >= 1000:  # previously, if contig2TotLen[contig] / cobraSeq2len[contig] >= 0.05
            if os.path.exists('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(contig)):
                a = open('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(contig), 'r')
                for record in SeqIO.parse(a, "fasta"):
                    header = str(record.id).strip()
                    if header in extended_partial_query:
                        extended_partial_query.remove(header)
                        failed_join_list.append(header)
                    elif header in extended_circular_query:
                        extended_circular_query.remove(header)
                        failed_join_list.append(header)
                    else:
                        pass
                a.close()
            elif contig in self_circular:
                self_circular.remove(contig)
                failed_join_list.append(contig)
            elif contig in self_circular_non_expected_overlap.keys():
                del self_circular_non_expected_overlap[contig]
                failed_join_list.append(contig)
            else:
                pass
        else:
            pass

    # for debug
    print('extended_circular_query', file=debug, flush=True)
    print(extended_circular_query, file=debug, flush=True)
    print('extended_partial_query', file=debug, flush=True)
    print(extended_partial_query, file=debug, flush=True)
    print('contig2assembly', file=debug, flush=True)
    print(contig2assembly, file=debug, flush=True)

    ##
    # get the unique sequences of COBRA "Extended" query contigs for joining check
    log_info('[18/23]', 'Saving unique sequences of "Extended_circular" and "Extended_partial" for joining checking.',
             '\n', log)
    extended_circular_fasta = open('COBRA_category_ii-a_extended_circular_unique.fasta', 'w')
    extended_partial_fasta = open('COBRA_category_ii-b_extended_partial_unique.fasta', 'w')
    query2current = {}

    for contig in query_set:
        if contig in extended_circular_query and contig not in redundant and contig not in is_same_as_redundant:
            extended_circular_fasta.write('>' + contig + '_extended_circular' + '\n')
            extended_circular_fasta.write(header2joined_seq[contig2extended_status[contig]] + '\n')
            extended_circular_fasta.flush()
            for item in order_all[contig]:
                if contig_name(item) in query_set:
                    query2current[contig_name(item)] = contig + '_extended_circular'
                else:
                    pass
        elif contig in extended_partial_query and contig not in redundant and contig not in is_same_as_redundant:
            extended_partial_fasta.write('>' + contig + '_extended_partial' + '\n')
            extended_partial_fasta.write(header2joined_seq[contig2extended_status[contig]] + '\n')
            extended_partial_fasta.flush()
            for item in order_all[contig]:
                if contig_name(item) in query_set:
                    if contig_name(item) in extended_partial_query:
                        query2current[contig_name(item)] = contig + '_extended_partial'
                    elif contig_name(item) in failed_join_list:
                        query2current[contig_name(item)] = contig + '_extended_partial'
                        extended_partial_query.add(contig_name(item))
                        failed_join_list.remove(contig_name(item))
                else:
                    pass
        else:
            pass
    extended_partial_fasta.close()
    extended_circular_fasta.close()
    summary_fasta('COBRA_category_ii-a_extended_circular_unique.fasta')
    summary_fasta('COBRA_category_ii-b_extended_partial_unique.fasta')

    # for debug
    print('query2current', file=debug, flush=True)
    print(query2current, file=debug, flush=True)
    debug.close()

    ##
    # save the joining details information
    log_info('[19/23]',
             'Getting the joining details of unique "Extended_circular" and "Extended_partial" query contigs.', '\n',
             log)
    joining_detail_headers = ['Final_Seq_ID', 'Joined_Len', 'Status', 'Joined_Seq_ID', 'Direction', 'Joined_Seq_Len',
                              'Start', 'End', 'Joined_Seq_Cov', 'Joined_Seq_GC', 'Joined_reason']
    extended_circular_joining_details = open('COBRA_category_ii-a_extended_circular_unique_joining_details.txt', 'w')
    extended_circular_joining_details.write('\t'.join(joining_detail_headers[:]) + '\n')
    extended_partial_joining_details = open('COBRA_category_ii-b_extended_partial_unique_joining_details.txt', 'w')
    extended_partial_joining_details.write('\t'.join(joining_detail_headers[:]) + '\n')
    contig2join_details = {}

    for contig in order_all.keys():
        site = 1
        if contig in path_circular and contig not in redundant and contig not in is_same_as_redundant and contig not in failed_join_list:
            contig2join_details[contig + '_extended_circular'] = []
            for item in order_all[contig][:-1]:
                if get_direction(item) == 'forward':
                    contents = [contig + '_extended_circular',
                                str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                                'Circular', contig_name(item), 'forward', str(header2len[contig_name(item)]), str(site),
                                str(site + header2len[contig_name(item)] - 1),
                                str(cov[contig_name(item)]), str(gc[contig_name(item)]),
                                contig2join_reason[contig][contig_name(item)]]
                    contig2join_details[contig + '_extended_circular'].append('\t'.join(contents[:]))
                    site += header2len[contig_name(item)] - length
                else:
                    contents = [contig + '_extended_circular',
                                str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                                'Circular', contig_name(item) + '_rc', 'reverse', str(header2len[contig_name(item)]),
                                str(site + header2len[contig_name(item)] - 1),
                                str(site), str(cov[contig_name(item)]), str(gc[contig_name(item)]),
                                contig2join_reason[contig][contig_name(item)]]
                    contig2join_details[contig + '_extended_circular'].append('\t'.join(contents[:]))
                    site += header2len[contig_name(item)] - length

            last = order_all[contig][-1]
            if get_direction(last) == 'forward':
                contents = [contig + '_extended_circular',
                            str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                            'Circular', contig_name(last), 'forward', str(header2len[contig_name(last)]), str(site),
                            str(site + header2len[contig_name(last)] - 1),  # length - 1),
                            str(cov[contig_name(last)]), str(gc[contig_name(last)]),
                            contig2join_reason[contig][contig_name(last)], '\n']
                contig2join_details[contig + '_extended_circular'].append('\t'.join(contents[:]))

            else:
                contents = [contig + '_extended_circular',
                            str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                            'Circular', contig_name(last) + '_rc', 'reverse', str(header2len[contig_name(last)]),
                            str(site + header2len[contig_name(last)] - 1),  # - length - 1),
                            str(site), str(cov[contig_name(last)]), str(gc[contig_name(last)]),
                            contig2join_reason[contig][contig_name(last)], '\n']
                contig2join_details[contig + '_extended_circular'].append('\t'.join(contents[:]))
        elif contig in extended_partial_query and contig not in redundant and contig not in is_same_as_redundant:
            contig2join_details[contig + '_extended_partial'] = []
            for item in order_all[contig][:-1]:
                if get_direction(item) == 'forward':
                    contents = [contig + '_extended_partial',
                                str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                                'Partial', contig_name(item), 'forward', str(header2len[contig_name(item)]), str(site),
                                str(site + header2len[contig_name(item)] - 1),
                                str(cov[contig_name(item)]), str(gc[contig_name(item)]),
                                contig2join_reason[contig][contig_name(item)]]
                    contig2join_details[contig + '_extended_partial'].append('\t'.join(contents[:]))
                    site += header2len[contig_name(item)] - length
                else:
                    contents = [contig + '_extended_partial',
                                str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                                'Partial', contig_name(item) + '_rc', 'reverse', str(header2len[contig_name(item)]),
                                str(site + header2len[contig_name(item)] - 1),
                                str(site), str(cov[contig_name(item)]), str(gc[contig_name(item)]),
                                contig2join_reason[contig][contig_name(item)]]
                    contig2join_details[contig + '_extended_partial'].append('\t'.join(contents[:]))
                    site += header2len[contig_name(item)] - length

            last = order_all[contig][-1]  # for the last one in non-circular path, the end position should be different
            if get_direction(last) == 'forward':
                contents = [contig + '_extended_partial',
                            str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                            'Partial', contig_name(last), 'forward', str(header2len[contig_name(last)]), str(site),
                            str(site + header2len[contig_name(last)] - 1),
                            str(cov[contig_name(last)]), str(gc[contig_name(last)]),
                            contig2join_reason[contig][contig_name(last)], '\n']
                contig2join_details[contig + '_extended_partial'].append('\t'.join(contents[:]))
            else:
                contents = [contig + '_extended_partial',
                            str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                            'Partial', contig_name(last) + '_rc', 'reverse', str(header2len[contig_name(last)]),
                            str(site + header2len[contig_name(last)] - 1),
                            str(site), str(cov[contig_name(last)]), str(gc[contig_name(last)]),
                            contig2join_reason[contig][contig_name(last)], '\n']
                contig2join_details[contig + '_extended_partial'].append('\t'.join(contents[:]))
        else:
            pass

    for seq in contig2join_details.keys():
        if 'circular' in seq:
            extended_circular_joining_details.write('\n'.join(contig2join_details[seq][:]))
        else:
            extended_partial_joining_details.write('\n'.join(contig2join_details[seq][:]))
    extended_circular_joining_details.close()
    extended_partial_joining_details.close()

    ##
    # save the joining summary information
    log_info('[20/23]', 'Saving joining summary of "Extended_circular" and "Extended_partial" query contigs.', '\n',
             log)
    assembly_summary = open('COBRA_joining_summary.txt', 'w')
    assembly_summary_headers = ['Query_Seq_ID', 'Query_Seq_Len', 'Total_Joined_Seqs', 'Joined_seqs', 'Total_Joined_Len',
                                'Assembled_Len', 'Extended_Len', 'Status', 'Final_Seq_ID']
    assembly_summary.write('\t'.join(assembly_summary_headers[:]) + '\n')

    for contig in list(extended_circular_query) + list(extended_partial_query):
        if contig in query2current.keys():
            assembly_summary.write(
                '\t'.join([contig, str(header2len[contig]), summarize(contig), query2current[contig]]) + '\n')
        else:
            if contig in extended_circular_query:
                extended_circular_query.remove(contig)
                failed_join_list.append(contig)
            elif contig in extended_partial_query:
                extended_partial_query.remove(contig)
                failed_join_list.append(contig)
            else:
                pass
    assembly_summary.close()

    ##
    # save the joining status information of each query
    log_info('[21/23]', 'Saving joining status of all query contigs.', '\n', log)
    assembled_info = open('COBRA_joining_status.txt', 'w')  # shows the COBRA status of each query
    assembled_info.write(
        'SeqID' + '\t' + 'Length' + '\t' + 'Coverage' + '\t' + 'GC' + '\t' + 'Status' + '\t' + 'Category' + '\n')

    # for those could be extended to circular
    for contig in extended_circular_query:
        assembled_info.write(
            contig + '\t' + str(header2len[contig]) + '\t' + str(cov[contig]) + '\t' + gc[contig] + '\t'
            + 'Extended_circular' + '\t' + 'category_ii-a' + '\n')

    # for those could be extended ok
    for contig in extended_partial_query:
        assembled_info.write(
            contig + '\t' + str(header2len[contig]) + '\t' + str(cov[contig]) + '\t' + gc[contig] + '\t'
            + 'Extended_partial' + '\t' + 'category_ii-b' + '\n')

    # for those cannot be extended
    failed_join = open('COBRA_category_ii-c_extended_failed.fasta', 'w')
    for contig in set(failed_join_list):
        if contig not in extended_circular_query or contig not in extended_partial_query or contig not in orphan_end_query:
            failed_join.write('>' + contig + '\n')
            failed_join.write(header2seq[contig] + '\n')
            assembled_info.write(contig + '\t' + str(header2len[contig]) + '\t' + str(cov[contig]) + '\t' + gc[contig]
                                 + '\t' + 'Extended_failed' + '\t' + 'category_ii-c' + '\n')
        else:
            pass
    failed_join.close()
    summary_fasta('COBRA_category_ii-c_extended_failed.fasta')

    # for those due to orphan end
    orphan_end = open('COBRA_category_iii_orphan_end.fasta', 'w')
    for contig in orphan_end_query:
        orphan_end.write('>' + contig + '\n')
        orphan_end.write(header2seq[contig] + '\n')
        assembled_info.write(
            contig + '\t' + str(header2len[contig]) + '\t' + str(cov[contig]) + '\t' + gc[contig] + '\t'
            + 'Orphan_end' + '\t' + 'category_iii' + '\n')
    orphan_end.close()
    summary_fasta('COBRA_category_iii_orphan_end.fasta')

    # for self circular
    log_info('[22/23]', 'Saving self_circular contigs.', '\n', log)
    circular_fasta = open('COBRA_category_i_self_circular.fasta', 'w')

    for contig in self_circular:
        assembled_info.write(
            contig + '\t' + str(header2len[contig] - length) + '\t' + str(cov[contig]) + '\t' + gc[contig]
            + '\t' + 'Self_circular' + '\t' + 'category_i' + '\n')
        circular_fasta.write('>' + contig + '_self_circular' + '\n')
        circular_fasta.write(header2seq[contig] + '\n')

    for contig in self_circular_non_expected_overlap.keys():
        assembled_info.write(
            contig + '\t' + str(header2len[contig] - self_circular_non_expected_overlap[contig]) + '\t' +
            str(cov[contig]) + '\t' + gc[contig] + '\t' +
            'Self_circular' + '\t' + 'category_i' + '\n')
        circular_fasta.write('>' + contig + '_self_circular' + '\n')
        circular_fasta.write(header2seq[contig] + '\n')

    circular_fasta.close()
    assembled_info.close()
    summary_fasta('COBRA_category_i_self_circular.fasta')

    ##
    # save new fasta file with all the others used in joining replaced by COBRA sequences excepting self_circular ones
    log_info('[23/23]', 'Saving the new fasta file.', '\n', log)

    for contig in all_joined_query:
        del header2seq[contig]

    with open('{0}.new.fa'.format(fasta_name.rsplit('.', 1)[0]), 'w') as new:
        for header, sequence in header2seq.items():
            new.write(f">{header}\n{sequence}\n")

    os.system('cat {0}.new.fa COBRA_category_ii-a_extended_circular_unique.fasta '
              'COBRA_category_ii-b_extended_partial_unique.fasta '
              '>{0}.new.fa.'.format(fasta_name.rsplit('.', 1)[0]))
    os.system('mv {0}.new.fa. {0}.new.fa'.format(fasta_name.rsplit('.', 1)[0]))

    ##
    # intermediate files
    os.mkdir('intermediate.files')
    os.system(
        'mv COBRA_end_joining_pairs.txt COBRA_potential_joining_paths.txt COBRA_retrieved_for_joining intermediate.files')
    os.mkdir('intermediate.files/invalid.checking')
    os.system('mv blastdb_1.fa* blastdb_2.fa blastdb_2.vs.blastdb_1* intermediate.files/invalid.checking')

    ##
    # write the numbers to the log file
    log.write('\n')
    log.write('3. RESULTS SUMMARY' + '\n')
    log.write('# Total queries: ' + str(len(query_set)) + '\n' +
              '# Category i   - Self_circular: ' + str(count_seq('COBRA_category_i_self_circular.fasta')) + '\n' +
              '# Category ii  - Extended_circular: ' + str(len(extended_circular_query)) + ' (Unique: ' +
              str(count_seq('COBRA_category_ii-a_extended_circular_unique.fasta')) + ')\n' +
              '# Category ii  - Extended_partial: ' + str(len(extended_partial_query)) + ' (Unique: ' +
              str(count_seq('COBRA_category_ii-b_extended_partial_unique.fasta')) + ')\n' +
              '# Category ii  - Extended_failed (due to COBRA rules): ' + str(len(set(failed_join_list))) + '\n' +
              '# Category iii - Orphan end: ' + str(len(orphan_end_query)) + '\n' +
              '# Check "COBRA_joining_status.txt" for joining status of each query.' + '\n' +
              '# Check "COBRA_joining_summary.txt" for joining details of "Extended_circular" and "Extended_partial" queries.')
    log.flush()
    log.close()


if __name__ == '__main__':
    main()
