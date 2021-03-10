#!/home/linking/.pyenv/versions/3.8.2/bin/python
# Author: Lin-Xing Chen, UC Berkeley

# COBRA v1.1.19
# Contig Overlap Based Re-Assembly
# Release date: Feb 5, 2021

# Usage: check COBRA.py -h

import os
from Bio import SeqIO
from Bio.Seq import reverse_complement
from collections import defaultdict
from Bio.SeqUtils import GC
import argparse
import math
import pysam
from time import strftime

parser = argparse.ArgumentParser(description="This script is used to get higher quality (including circular) virus genomes by combining assembled contigs based on their end overlaps.")
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("-q", "--query", type=str, help="the query contigs file (fasta format).", required=True)
requiredNamed.add_argument("-qi",type=str, choices=["all", "some"], help="specify if the queries include <all> or only <some> of the ≥ 10 kbp virus contigs in the sample. "
                           "If 'some', then host contigs ≥ 10 kbp may be joined, thus the COBRA sequences should be checked with checkV to remove host fraction; "
                           "if 'all', for any contig that is ≥ 10 kbp, it will be joined only when it is in the <query>.", required=True)
requiredNamed.add_argument("-f", "--fasta", type=str, help="the whole set of assembled contigs (fasta format).", required=True)
requiredNamed.add_argument("-a", "--assembler", type=str, choices=["idba", "megahit", "metaspades"], help="de novo assembler used, COBRA not tested for others.", required=True)
requiredNamed.add_argument("-k", "--kmer", type=int, help="the max kmer size used in de novo assembly.", required=True)
requiredNamed.add_argument("-m", "--mapping", type=str, help="the reads mapping file in sam or bam format.", required=True)
parser.add_argument("-o", "--output", type=str, help="the name of output folder (default = <query>_COBRA').")
parser.add_argument("-c", "--coverage", type=str, help="the contig coverage file (two columns divided by tab), not mandatory if assembler is metaSPAdes.")
parser.add_argument("-v", "--version", action='version', version='COBRA v1.1.19')
args = parser.parse_args()


def log_info(description, line_feed, log_file):
    localtime = strftime("%c")
    print('[' + localtime + '] ' + description, end=line_feed, file=log_file, flush=True)


def contig_name(item):  # get contig name from end name
    if item.endswith('_L') or item.endswith('_R') or item.endswith('_Lrc') or item.endswith('_Rrc'):
        return item.rsplit('_', 1)[0]
    else:
        return item


def are_equal_paths(two_paths):  # evaluate if two paths are equal, two_paths is a list including two ends, e.g., a_L, b_R
    if two_paths[0].rsplit('_', 1)[1].startswith('L') and two_paths[1].rsplit('_', 1)[1].startswith('R'):
        return contig_name(two_paths[0]) + '_R' in one_join_end \
               and contig_name(two_paths[1]) + '_L' in one_join_end \
               and contig_name(link_pair[contig_name(two_paths[0]) + '_R'][0]) == \
               contig_name(link_pair[contig_name(two_paths[1]) + '_L'][0])
    elif two_paths[0].rsplit('_', 1)[1].startswith('R') and two_paths[1].rsplit('_', 1)[1].startswith('L'):
        return contig_name(two_paths[0]) + '_L' in one_join_end \
               and contig_name(two_paths[1]) + '_R' in one_join_end \
               and contig_name(link_pair[contig_name(two_paths[0]) + '_L'][0]) == \
               contig_name(link_pair[contig_name(two_paths[1]) + '_R'][0])
    elif two_paths[0].rsplit('_', 1)[1].startswith('R') and two_paths[1].rsplit('_', 1)[1].startswith('R'):
        return contig_name(two_paths[0]) + '_L' in one_join_end \
               and contig_name(two_paths[1]) + '_L' in one_join_end \
               and contig_name(link_pair[contig_name(two_paths[0]) + '_L'][0]) == \
               contig_name(link_pair[contig_name(two_paths[1]) + '_L'][0])
    elif two_paths[0].rsplit('_', 1)[1].startswith('L') and two_paths[1].rsplit('_', 1)[1].startswith('L'):
        return contig_name(two_paths[0]) + '_R' in one_join_end \
               and contig_name(two_paths[1]) + '_R' in one_join_end \
               and contig_name(link_pair[contig_name(two_paths[0]) + '_R'][0]) == \
               contig_name(link_pair[contig_name(two_paths[1]) + '_R'][0])
    else:
        return False


def could_circulate(point, contig, direction):  # check if the path is circular with the current contig included
    if len(link_pair[point]) == 2:  # point is the same as "target" in join_walker
        if not are_equal_paths(link_pair[point]):
            if direction == 'L':
                if contig_name(link_pair[point][0]) == contig:
                    # at least one of the ends of the point itself must not have > 1 join,
                    # otherwise it is a contig represents a intragenome repeat (same below)
                    if cov[contig_name(point)] < 1.5 * cov[contig]:
                        # 2 times is for repeat, but it is too risky, use 1.9 instead (same below)
                        return link_pair[point][0].endswith('_R') and is_linked(contig, point)
                    else:
                        contig2join[contig + '_L'].append(link_pair[point][1])
                        # if the point contig seems to a intragenome repeat one,
                        # the other join of it will be added to the path
                        # and the query will not be considered as circular here (same below)
                elif contig_name(link_pair[point][1]) == contig:
                    if cov[contig_name(point)] < 1.5 * cov[contig]:
                        return link_pair[point][1].endswith('_R') and is_linked(contig, point)
                    else:
                        contig2join[contig + '_L'].append(link_pair[point][0])
            elif direction == 'R':
                if contig_name(link_pair[point][0]) == contig:
                    if cov[contig_name(point)] < 1.5 * cov[contig]:
                        return link_pair[point][0].endswith('_L') and is_linked(contig, point)
                    else:
                        contig2join[contig + '_R'].append(link_pair[point][1])

                elif contig_name(link_pair[point][1]) == contig:
                    if cov[contig_name(point)] < 1.5 * cov[contig]:
                        return link_pair[point][1].endswith('_L') and is_linked(contig, point)
                    else:
                        contig2join[contig + '_R'].append(link_pair[point][0])
    elif len(link_pair[point]) == 1:
        if direction == 'L':
            if contig_name(link_pair[point][0]) == contig:
                return link_pair[point][0].endswith('_R')
        elif direction == 'R':
            if contig_name(link_pair[point][0]) == contig:
                return link_pair[point][0].endswith('_L')
    else:
        return False


def the_better_one(two_paths, contig):  # to choose the better path for join based on coverage
    if cov[contig]-cov[contig_name(two_paths[1])] != 0:
        if abs(cov[contig]-cov[contig_name(two_paths[0])])/abs(cov[contig]-cov[contig_name(two_paths[1])]) > 1:
            return two_paths[1]
        elif abs(cov[contig]-cov[contig_name(two_paths[0])])/abs(cov[contig]-cov[contig_name(two_paths[1])]) < 1:
            return two_paths[0]
        else:
            return two_paths[0]
    elif cov[contig]-cov[contig_name(two_paths[0])] == 0:
        return two_paths[0]
    else:
        return two_paths[1]


def the_dominant_one(two_paths):  # get the dominant path from two equal paths
    if cov[contig_name(two_paths[0])] >= cov[contig_name(two_paths[1])]:
        return two_paths[0]
    else:
        return two_paths[1]


def the_rare_one(two_paths):  # get the less dominant path from two equal paths
    if cov[contig_name(two_paths[0])] >= cov[contig_name(two_paths[1])]:
        return two_paths[1]
    else:
        return two_paths[0]


def other_end_is_extendable(end, contig):  # check if the other end is extendable
    if contig_name(end) not in self_circular:
        if 0.333 * cov[contig] <= cov[contig_name(end)] <= 2.2 * cov[contig]:
            # 0.333 is for within-population variations of three subpopulations
            # and 2.2 is for direct repeat at the end of linear phage genomes
            if end.rsplit('_', 1)[1].startswith('L'):
                if contig_name(end) + '_R' in link_pair.keys():
                    return len(link_pair[contig_name(end) + '_R']) > 0
                else:
                    return False
            elif end.rsplit('_', 1)[1].startswith('R'):
                if contig_name(end) + '_L' in link_pair.keys():
                    return len(link_pair[contig_name(end) + '_L']) > 0
                else:
                    return False
            else:
                return False
        else:
            return False
    else:
        return False


def is_ok_to_add(end, contig):
    # if the other end of a potential joining contig cant be extend,
    # add only when it is in the query list and has a good coverage
    if contig_name(end) not in self_circular:
        return 0.9 * cov[contig] <= cov[contig_name(end)] <= 1.1 * cov[contig] and contig_name(end) in query_list


def not_checked(end_list, checked):  # to see if a potential contig has been check for adding or not
    total_checked = 0
    for end in end_list:
        if contig_name(end) in checked:
            total_checked += 1
    return total_checked == 0


def get_target(item):  # to get the target for next run of joining
    if item.rsplit('_', 1)[1] == 'Lrc':
        return item.rsplit('_', 1)[0] + '_Rrc'
    elif item.rsplit('_', 1)[1] == 'Rrc':
        return item.rsplit('_', 1)[0] + '_Lrc'
    elif item.rsplit('_', 1)[1] == 'R':
        return item.rsplit('_', 1)[0] + '_L'
    else:
        return item.rsplit('_', 1)[0] + '_R'


def is_linked(a, b):  # to evaluate if two contigs are linked based on linkage information
    my_list = list()
    my_list.append(contig_name(a))
    my_list.append(contig_name(b))
    return my_list in linkage


def percentage_finished(finished, total, log_file):
    p_list = []
    p2number = {}
    for item in range(1, 101):
        p_list.append(int(total*item/100))
        p2number[int(total*item/100)] = item
    if len(finished)/2 in p_list:
        print('{0}% finished ... '.format(p2number[len(finished)/2]), end='', file=log_file, flush=True)
    else:
        pass


unequal_paths_end = []
def join_walker(contig, direction):
    # get potential joins for a given query
    end = contig + '_' + direction
    contig_checked[end].append(contig)
    a = len(contig2join[end])
    if a == 0:
        percentage_finished(finished_end, len(query_list)-len(DNA_break_query), log)
        if end in one_join_end:
            if contig_name(link_pair[end][0]) != contig:
                if other_end_is_extendable(link_pair[end][0], contig) and is_linked(contig, link_pair[end][0]):
                    if header2len[contig_name(link_pair[end][0])] >= 10000:
                        if args.qi == 'all':
                            if contig_name(link_pair[end][0]) in query_list:
                                contig2join[end].append(link_pair[end][0])
                                contig_checked[end].append(contig_name(link_pair[end][0]))
                            else:
                                pass
                        else:
                            contig2join[end].append(link_pair[end][0])
                            contig_checked[end].append(contig_name(link_pair[end][0]))
                    else:
                        contig2join[end].append(link_pair[end][0])
                        contig_checked[end].append(contig_name(link_pair[end][0]))
                elif is_ok_to_add(link_pair[end][0], contig) and is_linked(contig, link_pair[end][0]):
                    contig2join[end].append(link_pair[end][0])
                    contig_checked[end].append(contig_name(link_pair[end][0]))
                else:
                    pass
            else:
                if direction == 'L' and link_pair[end][0] == contig + '_R':
                    # the other end could be joined with the current working point
                    self_circular.append(contig_name(end))
                elif direction == 'R' and link_pair[end][0] == contig + '_L':
                    # the other end could be joined with the current working point
                    self_circular.append(contig_name(end))
        elif end in two_joins_end:
            if are_equal_paths(link_pair[end]):
                contig2join[end].append(the_dominant_one(link_pair[end]))
                contig_checked[end].append(contig_name((link_pair[end][0])))
                contig_checked[end].append(contig_name((link_pair[end][1])))
                if cov[contig] > cov[contig_name(the_dominant_one(link_pair[end]))]:
                    if contig not in contig2equalpath.keys():
                        contig2equalpath[contig] = []
                        contig2equalpath[contig].append(the_dominant_one(link_pair[end]))
                        contig2equalpath[contig].append(the_rare_one(link_pair[end]))
                    else:
                        contig2equalpath[contig].append(the_dominant_one(link_pair[end]))
                        contig2equalpath[contig].append(the_rare_one(link_pair[end]))
            else:
                unequal_paths_end.append(end)
                if direction == '_L' and contig_name(end) + '_R' in link_pair[end]:
                    if len(link_pair[contig_name(end) + '_R']) == 1:
                        # the other end must have only one join
                        self_circular.append(contig_name(end))
                elif direction == 'R' and contig_name(end) + '_L' in link_pair[end]:
                    if len(link_pair[contig_name(end) + '_L']) == 1:
                        # the other end must have only one join
                        self_circular.append(contig_name(end))
                elif cov[contig_name(link_pair[end][0])] + cov[contig_name(link_pair[end][1])] >= cov[contig] * 0.5:
                    # 0.5 is ok, too big will get much fewer "Extended circular" ones.
                    if header2len[contig_name(the_better_one(link_pair[end], contig))] >= 10000:
                        if args.qi == 'all':
                            if contig_name(the_better_one(link_pair[end], contig)) in query_list:
                                contig2join[end].append(the_better_one(link_pair[end], contig))
                                contig_checked[end].append(contig_name((link_pair[end][0])))
                                contig_checked[end].append(contig_name((link_pair[end][1])))
                            else:
                                pass
                        else:
                            contig2join[end].append(the_better_one(link_pair[end], contig))
                            contig_checked[end].append(contig_name((link_pair[end][0])))
                            contig_checked[end].append(contig_name((link_pair[end][1])))
                    else:
                        contig2join[end].append(the_better_one(link_pair[end], contig))
                        contig_checked[end].append(contig_name((link_pair[end][0])))
                        contig_checked[end].append(contig_name((link_pair[end][1])))
                else:
                    pass
        else:
            pass
    else:
        target = get_target(contig2join[end][-1])
        if target in one_join_end:
            if not_checked(link_pair[target], contig_checked[end]):
                if other_end_is_extendable(link_pair[target][0], contig) and is_linked(target, link_pair[target][0]):
                    if header2len[contig_name(link_pair[target][0])] >= 10000:
                        if args.qi == 'all':
                            if contig_name(link_pair[target][0]) in query_list:
                                contig2join[end].append(link_pair[target][0])
                                contig_checked[end].append(contig_name(link_pair[target][0]))
                            else:
                                pass
                        else:
                            contig2join[end].append(link_pair[target][0])
                            contig_checked[end].append(contig_name(link_pair[target][0]))
                    else:
                        contig2join[end].append(link_pair[target][0])
                        contig_checked[end].append(contig_name(link_pair[target][0]))
                elif is_ok_to_add(link_pair[target][0], contig) and is_linked(target, link_pair[target][0]):
                    if header2len[contig_name(link_pair[target][0])] >= 10000:
                        if args.qi == 'all':
                            if contig_name(link_pair[target][0]) in query_list:
                                contig2join[end].append(link_pair[target][0])
                                contig_checked[end].append(contig_name(link_pair[target][0]))
                            else:
                                pass
                        else:
                            contig2join[end].append(link_pair[target][0])
                            contig_checked[end].append(contig_name(link_pair[target][0]))
                    else:
                        contig2join[end].append(link_pair[target][0])
                        contig_checked[end].append(contig_name(link_pair[target][0]))
                else:
                    pass
            else:
                if could_circulate(target, contig, direction):
                    path_circular.append(contig)
        elif target in two_joins_end:
            if not_checked(link_pair[target], contig_checked[end]):
                if are_equal_paths(link_pair[target]):
                    contig2join[end].append(the_dominant_one(link_pair[target]))
                    contig_checked[end].append(contig_name(link_pair[target][0]))
                    contig_checked[end].append(contig_name(link_pair[target][1]))
                    if cov[contig] > cov[contig_name(the_dominant_one(link_pair[target]))]:
                        if contig not in contig2equalpath.keys():
                            contig2equalpath[contig] = []
                            contig2equalpath[contig].append(the_dominant_one(link_pair[target]))
                            contig2equalpath[contig].append(the_rare_one(link_pair[target]))
                        else:
                            contig2equalpath[contig].append(the_dominant_one(link_pair[target]))
                            contig2equalpath[contig].append(the_rare_one(link_pair[target]))
                else:
                    unequal_paths_end.append(target)
                    if cov[contig_name(link_pair[target][0])] + cov[contig_name(link_pair[target][1])] >= cov[contig] * 0.5:
                        # 0.5 is ok, too big will get much fewer "Extended circular" ones.
                        if header2len[contig_name(the_better_one(link_pair[target], contig))] >= 10000:
                            if args.qi == 'all':
                                if contig_name(the_better_one(link_pair[target], contig)) in query_list:
                                    contig2join[end].append(the_better_one(link_pair[target], contig))
                                    contig_checked[end].append(contig_name((link_pair[target][0])))
                                    contig_checked[end].append(contig_name((link_pair[target][1])))
                                else:
                                    pass
                            else:
                                contig2join[end].append(the_better_one(link_pair[target], contig))
                                contig_checked[end].append(contig_name((link_pair[target][0])))
                                contig_checked[end].append(contig_name((link_pair[target][1])))
                        else:
                            contig2join[end].append(the_better_one(link_pair[target], contig))
                            contig_checked[end].append(contig_name((link_pair[target][0])))
                            contig_checked[end].append(contig_name((link_pair[target][1])))
            else:
                if could_circulate(target, contig, direction):
                    path_circular.append(contig)
        else:
            pass
    if a == len(contig2join[end]):
        finished_end.add(end)
    return a < len(contig2join[end])


def retrieve(contig):  # to retrieve and save all the contigs in the joining path of a query
    if len(contig2join[contig + '_L']) > 0 or len(contig2join[contig + '_R']) > 0:
        out = open('{0}_retrieved.fa'.format(contig), 'w')
        added = []
        print('>' + contig, file=out)
        print(header2seq[contig], file=out)
        for item in contig2join[contig + '_L']:
            if contig_name(item) not in added:
                print('>' + contig_name(item), file=out)
                print(header2seq[contig_name(item)], file=out)
                added.append(contig_name(item))
            else:
                pass
        for item in contig2join[contig + '_R']:
            if contig_name(item) not in added:
                print('>' + contig_name(item), file=out)
                print(header2seq[contig_name(item)], file=out)
                added.append(contig_name(item))
            else:
                pass
    else:
        pass


def join_seqs(contig):  # get the join order of the sequences in a give path
    order_all[contig] = []
    left = contig + '_L'
    right = contig + '_R'
    added_to_contig[contig] = []
    if len(contig2join[left]) > 0:
        order_left[contig] = []
        for item in contig2join[left]:
            order_left[contig].append(item)
    if len(contig2join[right]) > 0:
        order_right[contig] = []
        for item in contig2join[right]:
           order_right[contig].append(item)

    if contig in order_left.keys() and len(order_left[contig]) > 0:
        for item in order_left[contig][::-1]:
            # the order of contigs added into the left direction should be reversed
            order_all[contig].append(item)
            added_to_contig[contig].append(contig_name(item))
    order_all[contig].append(contig)
    # the query contig itself should be added to the path as well
    if contig in order_right.keys() and len(order_right[contig]) > 0:
        for item in order_right[contig]:
            if contig_name(item) not in added_to_contig[contig]:
                # one contig should NOT be added twice into a single joining path
                order_all[contig].append(item)
            else:
                if cov[contig_name(item)] >= 1.9 * cov[contig]:
                    # if a contig thought to be a repeat, ok to add twice.
                    order_all[contig].append(item)


def calculate_seq_num(fasta_file):  # calculate the number of seqs in a fasta file
    seq_num = 0
    a = open(fasta_file, 'r')
    for line in a.readlines():
        if line.startswith('>'):
            seq_num += 1
    a.close()
    return seq_num


def calculate_seq_len(fasta_file):  # calculate the length of sequences in a fasta file
    seq_len = 0
    a = open(fasta_file, 'r')
    for line in a.readlines():
        if not line.startswith('>'):
            seq_len += len(line.strip())
    a.close()
    return seq_len


def summary_fasta(fasta_file):  # summary basic information of a fasta file
    summary_file = open('{0}.summary.txt'.format(fasta_file), 'w')
    if 'self_circular' in fasta_file:
        summary_file_headers = ['SeqID', 'Length', 'GC', 'Ns', 'DTR_length']
        print('\t'.join(summary_file_headers[:]), file=summary_file, flush=True)
    else:
        summary_file_headers = ['SeqID', 'Length', 'GC', 'Ns']
        print('\t'.join(summary_file_headers[:]), file=summary_file, flush=True)
    with open('{0}'.format(fasta_file), 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            header = str(record.id).strip()
            seq = str(record.seq)
            Ns = seq.count('N')
            if header.split('_self')[0] in self_circular:
                sequence_stats = [header, str(len(seq)), str(round(GC(seq), 3)), str(Ns), str(length)]
                print('\t'.join(sequence_stats[:]), file=summary_file, flush=True)
            elif header.split('_self')[0] in self_circular_non_expected_overlap.keys():
                sequence_stats = [header, str(len(seq)), str(round(GC(seq), 3)), str(Ns), str(self_circular_non_expected_overlap[header.split('_self')[0]])]
                print('\t'.join(sequence_stats[:]), file=summary_file, flush=True)
            else:
                sequence_stats = [header, str(len(seq)), str(round(GC(seq), 3)), str(Ns)]
                print('\t'.join(sequence_stats[:]), file=summary_file, flush=True)
    f.close()
    summary_file.close()


def summarize(contig):  # summary the retrieved contigs and joined information of a query
    if contig not in is_redundant_of.keys():
        b = calculate_seq_num('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(contig))  # number of retrieved contigs
        c = calculate_seq_len('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(contig))  # total length of retrieved contigs
        d = calculate_seq_num('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(contig))  # Total of sequences after joining
        e = calculate_seq_len('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(contig))  # total length after joining
        if contig in extended_circular_query:
            return str(b) + '\t' + str(c) + '\t' + str(d) + '\t' + str(e) + '\t' + 'Extended circular'
        elif contig in extended_10k_query:
            return str(b) + '\t' + str(c) + '\t' + str(d) + '\t' + str(e) + '\t' + 'Extended >= 10k'
        elif contig in extended_few_query:
            return str(b) + '\t' + str(c) + '\t' + str(d) + '\t' + str(e) + '\t' + 'Extended < 10k'
    else:
        if is_redundant_of[contig] in is_redundant_of.keys():
            item = is_redundant_of[is_redundant_of[contig]]
            b = calculate_seq_num('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(item))  # number of retrieved contigs
            c = calculate_seq_len('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(item))  # total length of retrieved contigs
            d = calculate_seq_num('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(item))  # Total of sequences after joining
            e = calculate_seq_len('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(item))  # total length after joining
            if item in extended_circular_query:
                return str(b) + '\t' + str(c) + '\t' + str(d) + '\t' + str(e) + '\t' + 'Extended circular'
            elif item in extended_10k_query:
                return str(b) + '\t' + str(c) + '\t' + str(d) + '\t' + str(e) + '\t' + 'Extended >= 10k'
            elif item in extended_few_query:
                return str(b) + '\t' + str(c) + '\t' + str(d) + '\t' + str(e) + '\t' + 'Extended < 10k'
        else:
            item = is_redundant_of[contig]
            b = calculate_seq_num('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(item))  # number of retrieved contigs
            c = calculate_seq_len('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(item))  # total length of retrieved contigs
            d = calculate_seq_num('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(item))  # Total of sequences after joining
            e = calculate_seq_len('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(item))  # total length after joining
            if item in extended_circular_query:
                return str(b) + '\t' + str(c) + '\t' + str(d) + '\t' + str(e) + '\t' + 'Extended circular'
            elif item in extended_10k_query:
                return str(b) + '\t' + str(c) + '\t' + str(d) + '\t' + str(e) + '\t' + 'Extended >= 10k'
            elif item in extended_few_query:
                return str(b) + '\t' + str(c) + '\t' + str(d) + '\t' + str(e) + '\t' + 'Extended < 10k'


def determine_status(contig):
    if contig not in failed_join_list:
        if contig in redundant:
            a = header2len[contig]  # length of Query contig
            b = calculate_seq_len('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(is_redundant_of[contig]))  # total length after joining
            if is_redundant_of[contig] in set(path_circular):
                extended_circular_query.append(contig)
            else:
                if b - a >= 10000:
                    extended_10k_query.append(contig)
                else:
                    extended_few_query.append(contig)
        else:
            a = header2len[contig]  # length of Query contig
            b = calculate_seq_len('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(contig))  # total length after joining
            if contig in set(path_circular):
                extended_circular_query.append(contig)
            else:
                if b - a >= 10000:
                    extended_10k_query.append(contig)
                else:
                    extended_few_query.append(contig)
    else:
        pass


def is_subset(set_b, set_a):  # to check if an joining path is part of another one
    not_in = 0
    for item in set_b:
        if not item in set_a:
            not_in += 1
    return not_in == 0


def get_direction(item):  # get the direction in an joining path
    if item.endswith('rc'):
        return 'reverse'
    else:
        return 'forward'


def total_length(contig_list):  # get the total length of all sequences in an joining path with overlap removing
    total = 0
    count = 0
    for item in contig_list:
        total += header2len[contig_name(item)]
    return total - count * length


def rename(text):  # rename sequence names for ggplot figures
    chars = ["-", '.']
    for char in chars:
        while char in text:
            text = text.replace(char, "_")
    return text


# check if the coverage file is provided
if args.assembler == 'idba' and not args.coverage:
    print('Please provide coverage information for {0} assembled contigs.'.format(args.assembler))
    exit()


# get the name of the query file
if '/' in args.query:
    query_name = '{0}'.format(args.query).rsplit('/')[1]
else:
    query_name = '{0}'.format(args.query)

if '/' in args.fasta:
    fasta_name = '{0}'.format(args.fasta).rsplit('/')[1]
else:
    fasta_name = '{0}'.format(args.fasta)


# folder of output
if not args.output:
    working_dir = '{0}_COBRA'.format(query_name)
else:
    working_dir = '{0}'.format(args.output)
if os.path.exists('{0}'.format(working_dir)):
    print('Output folder ({0}) exists, please check.'.format(working_dir))
    exit()
else:
    os.mkdir('{0}'.format(working_dir))


# log file
log = open('{0}/log'.format(working_dir), 'w')  # log file
star_num = int((150-len(' COBRA analyses for {0} '.format(query_name)))/2)
print('*' * star_num + ' COBRA analyses for {0} '.format(query_name) + '*' * star_num, file=log, flush=True)
print('', end='\n', file=log, flush=True)
if args.qi == 'some':
    parameters = ['Key parameters used: ', '** Assembler: ' + args.assembler, '** Max-kmer: ' + str(args.kmer),
                  '** Query file includes ONLY SOME predicted virus contigs (the obtained COBRA sequences should be checked for host fraction).']
else:
    parameters = ['Key parameters used: ', '** Assembler: ' + args.assembler, '** Max-kmer: ' + str(args.kmer), '** Query file includes ALL predicted virus contigs.']
print('\n'.join(parameters[:]), file=log, flush=True)
print('', end='\n', file=log, flush=True)


# determine the length of overlap based on assembler and the largest kmer size
if args.assembler == "idba":
    length = args.kmer - 1
else:
    length = args.kmer


# import the contigs in the pool and save the end sequence
log_info('[1/18] Reading contigs and getting contig ends ... ', '', log)
header2seq = {}
header2len = {}
gc = {}
L = {}
R = {}
Lrc = {}
Rrc = {}
f = open('{0}'.format(args.fasta), 'r')
for record in SeqIO.parse(f, "fasta"):
    header = str(record.id).strip()
    seq = str(record.seq)
    header2seq[header] = seq
    header2len[header] = len(seq)
    gc[header] = str(round(GC(seq), 3))
    L[header + '_L'] = seq[0:length]  # the first x bp of left end
    Lrc[header + '_Lrc'] = reverse_complement(seq[0:length]) # the reverse sequence of first x bp of left end
    R[header + '_R'] = seq[-length:]  # the first x bp of right end
    Rrc[header + '_Rrc'] = reverse_complement(seq[-length:]) # the reverse sequence of first x bp of right end
f.close()
print('A total of {0} contigs are imported.'.format(len(header2seq.keys())), file=log, flush=True)

# get potential joins
log_info('[2/18] Getting shared contig ends ...', '\n', log)

link_pair = {}  # used to save all overlaps between ends

d_L = defaultdict(set)
d_Lrc = defaultdict(set)
d_R = defaultdict(set)
d_Rrc = defaultdict(set)

for k,v in L.items():  # save header2seq in dictionary with seqs as keys
    d_L[v].add(k)
for k,v in Lrc.items():
    d_Lrc[v].add(k)
for k,v in R.items():
    d_R[v].add(k)
for k,v in Rrc.items():
    d_Rrc[v].add(k)

d_L_d_Lrc_shared = set(d_L.keys()).intersection(set(d_Lrc.keys()))
# get the shared seqs between direction pairs (L/Lrc, Lrc/L, L/R, R/L, R/Rrc, Rrc/R, Lrc/Rrc, Rrc/Lrc)
d_L_d_R_shared = set(d_L.keys()).intersection(set(d_R.keys()))
# the d_R_d_L_shared will be included below
d_R_d_Rrc_shared = set(d_R.keys()).intersection(set(d_Rrc.keys()))
d_Rrc_d_Lrc_shared = set(d_Rrc.keys()).intersection(set(d_Lrc.keys()))

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


# save all paired links to a file
log_info('[3/18] Writing contig end joining pairs ...', '\n', log)
one_join_end = []  # the end of contigs with one potential join
two_joins_end = []  # the end of contigs with two potential joins

p = open('{0}/COBRA_end_joining_pairs.txt'.format(working_dir, fasta_name), 'w')
for item in link_pair.keys():
    for point in link_pair[item]:
        print(item + '\t' + point, file=p)  # print link pairs into a file for check if interested

    if len(link_pair[item]) == 1:
        one_join_end.append(item)  # add one joining end to a list, its pair may have one or more joins
    elif len(link_pair[item]) == 2 and len(link_pair[link_pair[item][0]]) == 1 and len(link_pair[link_pair[item][1]]) == 1:
        two_joins_end.append(item)  # add two joining end to a list, each of its pairs should only have one join
    else:
        pass
p.close()


# coverage
log_info('[4/18] Getting contig coverage information ...', '\n', log)
cov = {}
if args.assembler == 'metaspades':
    if not args.coverage:
        for header in header2seq.keys():
            cov[header] = float(header.rsplit('_cov_', 1)[1]) # get the coverage from metaSPAdes contig headers
    else:
        coverage = open('{0}'.format(args.coverage), 'r')
        for line in coverage.readlines():
            line = line.strip().split('\t')
            cov[line[0]] = float(line[1])
        coverage.close()
else:
    coverage = open('{0}'.format(args.coverage), 'r')
    for line in coverage.readlines():
        line = line.strip().split('\t')
        cov[line[0]] = float(line[1])
    coverage.close()


# open query file and save the information
log_info('[5/18] Getting query contig list ...', '\n', log)
query_list = []
DNA_break_query = []
q = open('{0}'.format(args.query), 'r')
for record in SeqIO.parse(q, "fasta"):
    header = str(record.id).strip()
    query_list.append(header)
    if header + '_L' not in link_pair.keys() and header + '_R' not in link_pair.keys():
        DNA_break_query.append(header)
q.close()

#for debug
print(DNA_break_query, flush=True)

# get the linkage of contigs based on paired-end reads mapping
log_info('[6/18] Getting contig linkage based on sam/bam ... Be patient, this may take long ... ', '\n', log)
linkage = []
unexpected_break_query = []
map_file = pysam.AlignmentFile('{0}'.format(args.mapping))
for line in map_file:
    if line.reference_name != line.next_reference_name:
        linkage.append([line.reference_name, line.next_reference_name])
        if line.next_reference_name in DNA_break_query:
            contig = line.next_reference_name
            if L[contig + '_L'] in line.query_sequence[2:-1] \
                    or Lrc[contig + '_Lrc'] in line.query_sequence[2:-1] \
                    or R[contig + '_R'] in line.query_sequence[2:-1] \
                    or Rrc[contig + '_Rrc'] in line.query_sequence[2:-1]:
                DNA_break_query.remove(contig)
                unexpected_break_query.append(contig)
            else:
                pass
    else:
        if line.is_unmapped and line.reference_name in DNA_break_query:
            contig = line.reference_name
            if L[contig + '_L'] in line.query_sequence[2:-1] \
                    or Lrc[contig + '_Lrc'] in line.query_sequence[2:-1] \
                    or R[contig + '_R'] in line.query_sequence[2:-1] \
                    or Rrc[contig + '_Rrc'] in line.query_sequence[2:-1]:
                DNA_break_query.remove(contig)
                unexpected_break_query.append(contig)
            else:
                pass
        else:
            pass
map_file.close()


# check valid joins
contig2join = {}
contig_checked = {}
path_circular = []
self_circular = []

for contig in query_list:
    contig2join[contig + '_L'] = []
    contig2join[contig + '_R'] = []
    contig_checked[contig + '_L'] = []
    contig_checked[contig + '_R'] = []


# walk the joins
log_info('[7/18] Walking joins. ', '', log)
contig2equalpath = {}
finished_end = set()
for contig in query_list:
    if contig not in DNA_break_query:
        while join_walker(contig, 'L'):
            join_walker(contig, 'L')
        while join_walker(contig, 'R'):
            join_walker(contig, 'R')
    else:
        pass


# save the potential joining paths
print('100% finished ...', end='\n', file=log, flush=True)
log_info('[8/18] Saving potential joining paths ...', '\n', log)
results = open('{0}/COBRA_potential_joining_paths.txt'.format(working_dir), 'w')
for item in contig2join.keys():
    print(item, contig2join[item], file=results, flush=True)


# DNA break info
log_info('[9/18] Saving contig DNA break information ...', '\n', log)
# determine potential self_circular contigs from DNA break ones
self_circular_non_expected_overlap = {}
for contig in DNA_break_query:
    if header2seq[contig].count(header2seq[contig][-19:]) == 2:
        if header2seq[contig].endswith(header2seq[contig].split(header2seq[contig][-19:])[0] + header2seq[contig][-19:]):
            self_circular_non_expected_overlap[contig] = len(header2seq[contig].split(header2seq[contig][-19:])[0]) + 19
            DNA_break_query.remove(contig)
        else:
            pass
    else:
        pass


# get the failed to join contigs, but not due to DNA break
failed_join_list = []
for contig in query_list:
    if contig + '_L' in link_pair.keys() or contig + '_R' in link_pair.keys():
        if len(contig2join[contig + '_L']) == 0 and len(contig2join[contig + '_R']) == 0 and contig not in self_circular:
            failed_join_list.append(contig)
        else:
            pass
    else:
        pass


# get retrieved sequences
log_info('[10/18] Getting retrieved contigs ...', '\n', log)
os.chdir('{0}'.format(working_dir))
os.mkdir('COBRA_retrieved_for_joining')
for contig in query_list:
    if contig not in self_circular and contig not in failed_join_list and contig not in DNA_break_query:
        # self_circular ones do not need to retrieve
        retrieve(contig)
    else:
        pass
os.system('mv *_retrieved.fa COBRA_retrieved_for_joining')


# get the joining paths
log_info('[11/18] Analyzing for valid joining paths ...', '\n', log)
contig2assembly = {}
for item in contig2join.keys():
    contig = item.rsplit('_', 1)[0]
    if contig not in contig2assembly.keys():
        contig2assembly[contig] = set()
        contig2assembly[contig].add(contig)
        for point in contig2join[item]:
            contig2assembly[contig].add(point.rsplit('_', 1)[0])
    else:
        for point in contig2join[item]:
            contig2assembly[contig].add(point.rsplit('_', 1)[0])


# find the redundant joining paths
redundant = set()
is_redundant_of = {}
is_the_same_as = set()
for item in contig2assembly.keys():
    for item_1 in contig2assembly.keys():
        if item != item_1 and is_subset(contig2assembly[item], contig2assembly[item_1]):
            if contig2assembly[item] != contig2assembly[item_1]:
                if item in path_circular:
                    redundant.add(item_1)
                    if item not in is_redundant_of.keys():
                        is_redundant_of[item_1] = item
                    else:
                        is_redundant_of[item_1] = is_redundant_of[item]
                    if item_1 not in contig2assembly[item]:
                        failed_join_list.append(item_1)
                    else:
                        path_circular.append(item_1)
                else:
                    if item_1 in path_circular:
                        path_circular.append(item)
                        redundant.add(item)
                        if item_1 not in is_redundant_of.keys():
                            is_redundant_of[item] = item_1
                        else:
                            is_redundant_of[item] = is_redundant_of[item_1]
                    else:
                        redundant.add(item)
                        if item_1 not in is_redundant_of.keys():
                            is_redundant_of[item] = item_1
                        else:
                            is_redundant_of[item] = is_redundant_of[item_1]
            else:
                if item not in is_the_same_as:
                    is_the_same_as.add(item_1)
                else:
                    pass
        else:
            pass


# remove the queries in multiple paths
all = []  # all contigs in all the paths checked
same_path = []  # two are the same path
contigs_in_two_or_more_paths = set()  # contigs that appeared in multiple paths
paths_sharing_contig = set()  # the queries in multiple non-unique paths

for contig in contig2assembly.keys():
    if contig not in redundant and contig not in failed_join_list:
        if contig2assembly[contig] not in same_path:
            same_path.append(contig2assembly[contig])
            for item in contig2assembly[contig]:
                if item in query_list:
                    all.append(item)
                else:
                    pass
        else:
            pass
    else:
        pass

for contig in set(all):
    if all.count(contig) > 1:
        contigs_in_two_or_more_paths.add(contig)
    else:
        pass

for contig in contigs_in_two_or_more_paths:
    for contig_1 in contig2assembly.keys():
        if contig_1 not in redundant and contig not in failed_join_list:
            if contig in contig2assembly[contig_1]:
                paths_sharing_contig.add(contig_1)
            else:
                pass
        else:
            pass

for contig in paths_sharing_contig:
    if contig in contig2assembly.keys():
        del contig2assembly[contig]
        failed_join_list.append(contig)
        os.remove('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(contig))
    else:
        pass

    for item in is_redundant_of.keys():
        if contig == is_redundant_of[item]:
            del contig2assembly[item]
            failed_join_list.append(item)
            if os.path.exists('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(item)):
                os.remove('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(item))
            else:
                pass

            for item_1 in is_redundant_of.keys():
                if is_redundant_of[item_1] == item:
                    del contig2assembly[item_1]
                    failed_join_list.append(item_1)
                    if os.path.exists('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(item_1)):
                        os.remove('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(item_1))
                    else:
                        pass
        else:
            pass


# join the retrieved sequences of non-redundant paths
log_info('[12/18] Saving valid joining paths ...', '\n', log)

order_left = {}
order_right = {}
order_all = {}
added_to_contig = {}

for contig in contig2assembly.keys():
    # only those contigs left in contig2assembly after filtering
    # (see above "# remove the queries in multiple paths") will be
    # checked for join paths (join_seqs) to get order_all
    if len(contig2assembly[contig]) > 1:
        join_seqs(contig)


for contig in order_all.keys():
    a = open('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(contig), 'w')
    # print header regarding the joining status
    if contig in path_circular:
        print('>' + contig + '_extended_circular', file=a, flush=True)
    else:
        print('>' + contig + '_extended_partial', file=a, flush=True)
    # print the sequences with their overlap removed
    for item in order_all[contig][:-1]:
        if item.endswith('_R') or item.endswith('_L'):
            print(header2seq[item.rsplit('_', 1)[0]][:-length], end='', file=a, flush=True)
        elif item.endswith('rc'):
            print(reverse_complement(header2seq[item.rsplit('_', 1)[0]])[:-length], end='', file=a, flush=True)
        else:
            print(header2seq[contig][:-length], end='', file=a, flush=True)
    if contig in path_circular:
        if order_all[contig][-1].endswith('rc'):
            print(reverse_complement(header2seq[order_all[contig][-1].rsplit('_',1)[0]])[:-length], file=a, flush=True)
        elif order_all[contig][-1].endswith('_R') or order_all[contig][-1].endswith('_L'):
            print(header2seq[order_all[contig][-1].rsplit('_',1)[0]][:-length], file=a, flush=True)
        else:
            print(header2seq[order_all[contig][-1]][:-length], file=a, flush=True)
    else:
        if order_all[contig][-1].endswith('rc'):
            print(reverse_complement(header2seq[order_all[contig][-1].rsplit('_',1)[0]]), file=a, flush=True)
        elif order_all[contig][-1].endswith('_R') or order_all[contig][-1].endswith('_L'):
            print(header2seq[order_all[contig][-1].rsplit('_',1)[0]], file=a, flush=True)
        else:
            print(header2seq[order_all[contig][-1]], file=a, flush=True)
    a.close()


# determine the joining status of queries
log_info('[13/18] Getting initial joining status of each query contigs ...', '\n', log)
extended_circular_query = []
extended_10k_query = []
extended_few_query = []
for contig in query_list:
    if os.path.exists('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(contig)):
        determine_status(contig)
    else:
        pass


# deal with cross-assignment queries
for item in failed_join_list:
    if item in is_redundant_of.keys() and is_redundant_of[item] in extended_circular_query:
        failed_join_list.remove(item)
        if item not in extended_circular_query:
            extended_circular_query.append(item)
        else:
            pass
    elif item in is_redundant_of.keys() and is_redundant_of[item] in extended_10k_query:
        failed_join_list.remove(item)
        if item not in extended_10k_query:
            extended_10k_query.append(item)
        else:
            pass
    else:
        pass


os.mkdir('COBRA_extended_few')
for item in extended_few_query:
    if item in redundant:
        if is_redundant_of[item] in extended_circular_query:
            extended_few_query.remove(item)
            if item not in extended_circular_query:
                extended_circular_query.append(item)
            else:
                pass
        elif is_redundant_of[item] in extended_10k_query:
            extended_few_query.remove(item)
            if item not in extended_10k_query:
                extended_10k_query.append(item)
            else:
                pass
        else:
            pass
    else:
        if item not in extended_10k_query and item not in extended_circular_query and item not in is_the_same_as:
            os.system('cp COBRA_retrieved_for_joining/{0}_retrieved*fa COBRA_extended_few'.format(item))
            # save those queries extended but fewer than 10 kbp to a folder


# save the joining summary information
log_info('[14/18] Saving joining summary of retrieved contigs ...', '\n', log)
assembly_summary = open('COBRA_joining_summary.txt', 'w')
assembly_summary_headers = ['QuerySeqID', 'QuerySeqLen', 'TotRetSeqs', 'TotRetLen', 'AssembledSeqs', 'AssembledLen', 'Status']
print('\t'.join(assembly_summary_headers[:]), file=assembly_summary, flush=True)
for contig in extended_circular_query + extended_10k_query + extended_few_query:
    print(contig + '\t' + str(header2len[contig]) + '\t' + summarize(contig), file=assembly_summary, flush=True)
assembly_summary.close()


# get the unique sequences of COBRA "well-extended" joining and their summary
log_info('[15/18] Saving unique sequences of "Extended_circular" and "Extended_10k" ...', '\n', log)
os.mkdir('COBRA_extended_10k_unique')
os.mkdir('COBRA_extended_circular_dominant_unique')
for contig in query_list:
    if contig in extended_circular_query and contig not in redundant and contig not in is_the_same_as:
        os.system('cp COBRA_retrieved_for_joining/{0}_retrieved*fa COBRA_extended_circular_dominant_unique'.format(contig))
    elif contig in extended_10k_query and contig not in redundant and contig not in is_the_same_as:
        os.system('cp COBRA_retrieved_for_joining/{0}_retrieved*fa COBRA_extended_10k_unique'.format(contig))
    else:
        pass

os.system('cat COBRA_extended_10k_unique/*joined.fa >COBRA_extended_10k_unique.fasta')
os.system('cat COBRA_extended_circular_dominant_unique/*joined.fa >COBRA_extended_circular_dominant_unique.fasta')
os.system('cat COBRA_extended_few/*joined.fa >COBRA_extended_few.fasta')
summary_fasta('COBRA_extended_10k_unique.fasta')
summary_fasta('COBRA_extended_circular_dominant_unique.fasta')
summary_fasta('COBRA_extended_few.fasta')


# save the joining details information
log_info('[16/18] Saving joining details of unique "Extended_circular" and "Extended_10k" queries ...', '\n', log)
joining_detail_headers = ['AssembledSeqID', 'NameForFig', 'AssembledLen', 'Status', 'RetSeqID', 'Direction', 'RetSeqLen', 'Start', 'End', 'RetSeqCov', 'RetSeqGC']
joining_detail_circular = open('COBRA_extended_circular_dominant_unique_joining_details.txt', 'w')
print('\t'.join(joining_detail_headers[:]), file=joining_detail_circular, flush=True)
joining_detail_10k = open('COBRA_extended_10k_unique_joining_details.txt', 'w')
print('\t'.join(joining_detail_headers[:]), file=joining_detail_10k, flush=True)
joining_detail_few = open('COBRA_extended_few_joining_details.txt', 'w')
print('\t'.join(joining_detail_headers[:]), file=joining_detail_few, flush=True)


for contig in order_all.keys():
    site = 1
    if contig in path_circular and contig not in redundant and contig not in is_the_same_as:
        for item in order_all[contig]:
            if get_direction(item) == 'forward':
                contents = [contig + '_extended_circular', rename(contig + '_extended_circular'), str(total_length(order_all[contig]) - length * len(order_all[contig])),
                            'Circular', contig_name(item),  'forward', str(header2len[contig_name(item)]), str(site), str(site + header2len[contig_name(item)] - length - 1),
                            str(cov[contig_name(item)]), str(gc[contig_name(item)])]
                print('\t'.join(contents[:]), file=joining_detail_circular, flush=True)
                site += header2len[contig_name(item)] - length
            else:
                contents = [contig + '_extended_circular', rename(contig + '_extended_circular'), str(total_length(order_all[contig]) - length * len(order_all[contig])),
                            'Circular', contig_name(item) + '_rc', 'reverse', str(header2len[contig_name(item)]), str(site + header2len[contig_name(item)] - length - 1),
                            str(site), str(cov[contig_name(item)]), str(gc[contig_name(item)])]
                print('\t'.join(contents[:]), file=joining_detail_circular, flush=True)
                site += header2len[contig_name(item)] - length
    elif contig in extended_10k_query and contig not in redundant and contig not in is_the_same_as:
        for item in order_all[contig][:-1]:
            if get_direction(item) == 'forward':
                contents = [contig + '_extended_partial', rename(contig + '_extended_partial'), str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                            'Partial', contig_name(item), 'forward', str(header2len[contig_name(item)]), str(site), str(site + header2len[contig_name(item)] - length - 1),
                            str(cov[contig_name(item)]), str(gc[contig_name(item)])]
                print('\t'.join(contents[:]), file=joining_detail_10k, flush=True)
                site += header2len[contig_name(item)] - length
            else:
                contents = [contig + '_extended_partial', rename(contig + '_extended_partial'), str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                            'Partial', contig_name(item) + '_rc', 'reverse', str(header2len[contig_name(item)]), str(site + header2len[contig_name(item)] - length - 1),
                            str(site), str(cov[contig_name(item)]), str(gc[contig_name(item)])]
                print('\t'.join(contents[:]), file=joining_detail_10k, flush=True)
                site += header2len[contig_name(item)] - length

        last = order_all[contig][-1]  # for the last one in non-circular path, the end position should be different
        if get_direction(last) == 'forward':
            contents = [contig + '_extended_partial', rename(contig + '_extended_partial'), str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                        'Partial', contig_name(last), 'forward', str(header2len[contig_name(last)]), str(site), str(site + header2len[contig_name(last)] - 1),
                        str(cov[contig_name(last)]), str(gc[contig_name(last)])]
            print('\t'.join(contents[:]), file=joining_detail_10k, flush=True)
        else:
            contents = [contig + '_extended_partial', rename(contig + '_extended_partial'), str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                        'Partial', contig_name(last) + '_rc', 'reverse', str(header2len[contig_name(last)]), str(site + header2len[contig_name(last)] - 1),
                        str(site), str(cov[contig_name(last)]), str(gc[contig_name(last)])]
            print('\t'.join(contents[:]), file=joining_detail_10k, flush=True)
    elif contig in extended_few_query and contig not in redundant and contig not in is_the_same_as:
        for item in order_all[contig][:-1]:
            if get_direction(item) == 'forward':
                contents = [contig + '_extended_partial', rename(contig + '_extended_partial'), str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                            'Partial', contig_name(item), 'forward', str(header2len[contig_name(item)]), str(site), str(site + header2len[contig_name(item)] - length - 1),
                            str(cov[contig_name(item)]), str(gc[contig_name(item)])]
                print('\t'.join(contents[:]), file=joining_detail_few, flush=True)
                site += header2len[contig_name(item)] - length
            else:
                contents = [contig + '_extended_partial', rename(contig + '_extended_partial'), str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                            'Partial', contig_name(item) + '_rc', 'reverse', str(header2len[contig_name(item)]), str(site + header2len[contig_name(item)] - length - 1),
                            str(site), str(cov[contig_name(item)]), str(gc[contig_name(item)])]
                print('\t'.join(contents[:]), file=joining_detail_few, flush=True)
                site += header2len[contig_name(item)] - length

        last = order_all[contig][-1]  # for the last one in non-circular path, the end position should be different
        if get_direction(last) == 'forward':
            contents = [contig + '_extended_partial', rename(contig + '_extended_partial'), str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                        'Partial', contig_name(last), 'forward', str(header2len[contig_name(last)]), str(site), str(site + header2len[contig_name(last)] - 1),
                        str(cov[contig_name(last)]), str(gc[contig_name(last)])]
            print('\t'.join(contents[:]), file=joining_detail_few, flush=True)
        else:
            contents = [contig + '_extended_partial', rename(contig + '_extended_partial'), str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                        'Partial', contig_name(last) + '_rc', 'reverse', str(header2len[contig_name(last)]), str(site + header2len[contig_name(last)] - 1),
                        str(site), str(cov[contig_name(last)]), str(gc[contig_name(last)])]
            print('\t'.join(contents[:]), file=joining_detail_few, flush=True)
    else:
        pass
joining_detail_10k.close()
joining_detail_circular.close()
joining_detail_few.close()


## get the joining details for each unique extended_circular and extended_10k for gggenes figure
assembled2detail = {}
with open('COBRA_extended_circular_dominant_unique_joining_details.txt', 'r') as joining_detail_circular:
    for line in joining_detail_circular.readlines():
        if not line.startswith('AssembledSeqID'):
            line = line.strip().split('\t')
            if line[0] not in assembled2detail.keys():
                assembled2detail[line[0]] = []
                assembled2detail[line[0]].append(line)
            else:
                assembled2detail[line[0]].append(line)
joining_detail_circular.close()

with open('COBRA_extended_10k_unique_joining_details.txt', 'r') as joining_detail_10k:
    for line in joining_detail_10k.readlines():
        if not line.startswith('AssembledSeqID'):
            line = line.strip().split('\t')
            if line[0] not in assembled2detail.keys():
                assembled2detail[line[0]] = []
                assembled2detail[line[0]].append(line)
            else:
                assembled2detail[line[0]].append(line)
joining_detail_10k.close()

with open('COBRA_extended_few_joining_details.txt', 'r') as joining_detail_few:
    for line in joining_detail_few.readlines():
        if not line.startswith('AssembledSeqID'):
            line = line.strip().split('\t')
            if line[0] not in assembled2detail.keys():
                assembled2detail[line[0]] = []
                assembled2detail[line[0]].append(line)
            else:
                assembled2detail[line[0]].append(line)
joining_detail_few.close()

for contig in order_all.keys():
    if os.path.exists('COBRA_extended_circular_dominant_unique/{0}_retrieved.fa'.format(contig)):
        a = open('COBRA_extended_circular_dominant_unique/{0}_retrieved_joining_details.txt'.format(contig), 'w')
        print('\t'.join(joining_detail_headers[:]), file=a, flush=True)
        for item in assembled2detail[contig + '_extended_circular']:
            print('\t'.join(item[:]), file=a, flush=True)
        a.close()
    elif os.path.exists('COBRA_extended_10k_unique/{0}_retrieved.fa'.format(contig)):
        a = open('COBRA_extended_10k_unique/{0}_retrieved_joining_details.txt'.format(contig), 'w')
        print('\t'.join(joining_detail_headers[:]), file=a, flush=True)
        for item in assembled2detail[contig + '_extended_partial']:
            print('\t'.join(item[:]), file=a, flush=True)
        a.close()
    elif os.path.exists('COBRA_extended_few/{0}_retrieved.fa'.format(contig)):
        a = open('COBRA_extended_few/{0}_retrieved_joining_details.txt'.format(contig), 'w')
        print('\t'.join(joining_detail_headers[:]), file=a, flush=True)
        for item in assembled2detail[contig + '_extended_partial']:
            print('\t'.join(item[:]), file=a, flush=True)
        a.close()
    else:
        pass


# save the joining status information of each query
assembled_info = open('COBRA_joining_status.txt', 'w')  # shows the COBRA status of each query
print('SeqID' + '\t' + 'Length' + '\t' + 'Coverage' + '\t' + 'GC' + '\t' + 'Status', file=assembled_info, flush=True)
# for those could be extended to circular
for contig in extended_circular_query:
    print(contig + '\t' + str(header2len[contig]) + '\t' + str(cov[contig]) + '\t' + gc[contig] + '\t' + 'Extended_circular', file=assembled_info, flush=True)
# for those could be extended >= 10k
for contig in extended_10k_query:
    print(contig + '\t' + str(header2len[contig]) + '\t' + str(cov[contig]) + '\t' + gc[contig] + '\t' + 'Extended >= 10k', file=assembled_info, flush=True)
# for those could be extended < 10k
for contig in extended_few_query:
    print(contig + '\t' + str(header2len[contig]) + '\t' + str(cov[contig]) + '\t' + gc[contig] + '\t' + 'Extended < 10k', file=assembled_info, flush=True)
# for those cant be extended
failed_join = open('COBRA_extended_failed.fasta', 'w')
for contig in set(failed_join_list):
    if contig not in extended_circular_query or contig not in extended_10k_query or contig not in extended_few_query:
        print('>' + contig, file=failed_join, flush=True)
        print(header2seq[contig], file=failed_join, flush=True)
        print(contig + '\t' + str(header2len[contig]) + '\t' + str(cov[contig]) + '\t' + gc[contig] + '\t' + 'Extended_failed', file=assembled_info, flush=True)
    else:
        failed_join_list.remove(contig)
summary_fasta('COBRA_extended_failed.fasta')
# for those due to DNA break
for contig in DNA_break_query:
    print(contig + '\t' + str(header2len[contig]) + '\t' + str(cov[contig]) + '\t' + gc[contig] + '\t' + 'DNA_break_or_short_contig_missing', file=assembled_info, flush=True)
# for those due to unexpected break at the end
for contig in unexpected_break_query:
    print(contig + '\t' + str(header2len[contig]) + '\t' + str(cov[contig]) + '\t' + gc[contig] + '\t' + 'Unexpected_break_by_assembler', file=assembled_info, flush=True)
# for self circular
log_info('[17/18] Saving self_circular contigs ...', '\n', log)
circular_fasta = open('COBRA_self_circular_queries_trimmed.fasta', 'w')
for contig in set(self_circular):
    print(contig + '\t' + str(header2len[contig] - length) + '\t' + str(cov[contig]) + '\t' + gc[contig] + '\t' + 'Self_circular', file=assembled_info, flush=True)
    print('>' + contig + '_self_circular', file=circular_fasta, flush=True)
    print(header2seq[contig][length:], file=circular_fasta, flush=True)
for contig in self_circular_non_expected_overlap.keys():
    print(contig + '\t' + str(header2len[contig] - self_circular_non_expected_overlap[contig]) + '\t' + str(cov[contig]) + '\t' + gc[contig] + '\t' + 'Self_circular', file=assembled_info, flush=True)
    print('>' + contig + '_self_circular', file=circular_fasta, flush=True)
    print(header2seq[contig][self_circular_non_expected_overlap[contig]:], file=circular_fasta, flush=True)
circular_fasta.close()
summary_fasta('COBRA_self_circular_queries_trimmed.fasta')
assembled_info.close()


# save the equal path information
log_info('[18/18] Saving less dominant genomes of "joined_circular" ...', '\n', log)
file_rows = []
dominant2lessdominant = {}
for key in contig2equalpath.keys():
    # Get list of numbers
    value_list = contig2equalpath[key]
    # Make sure to coerce to strings for printing
    value_list = [str(x) for x in value_list]
    # How many rows will you need to print
    n_rows = math.ceil(len(value_list) / 2)
    printed = []  # some equal paths were saved twice (left and right directions)
    for row_idx in range(n_rows):
        start_idx = row_idx * 2
        end_idx = start_idx + 1
        dominant2lessdominant[value_list[start_idx]] = value_list[end_idx]
        if contig_name(value_list[start_idx]) not in printed:
            if key in extended_circular_query and key not in is_the_same_as:  # those equal paths in non-extended_circular will not be saved.
                file_row = [key, contig_name(value_list[start_idx]), str(header2len[contig_name(value_list[start_idx])]), str(cov[contig_name(value_list[start_idx])]),
                            contig_name(value_list[end_idx]), str(header2len[contig_name(value_list[end_idx])]), str(cov[contig_name(value_list[end_idx])]), 'Extended_circular']
                printed.append(contig_name(value_list[start_idx]))
                file_rows.append(file_row)
            elif key in extended_10k_query and key not in is_the_same_as:  # those equal paths in non-extended_circular will not be saved.
                file_row = [key, contig_name(value_list[start_idx]), str(header2len[contig_name(value_list[start_idx])]), str(cov[contig_name(value_list[start_idx])]),
                            contig_name(value_list[end_idx]), str(header2len[contig_name(value_list[end_idx])]), str(cov[contig_name(value_list[end_idx])]), 'Extended_10k']
                printed.append(contig_name(value_list[start_idx]))
                file_rows.append(file_row)
            elif key in extended_few_query and key not in is_the_same_as:  # those equal paths in non-extended_circular will not be saved.
                file_row = [key, contig_name(value_list[start_idx]), str(header2len[contig_name(value_list[start_idx])]), str(cov[contig_name(value_list[start_idx])]),
                            contig_name(value_list[end_idx]), str(header2len[contig_name(value_list[end_idx])]), str(cov[contig_name(value_list[end_idx])]), 'Extended_few']
                printed.append(contig_name(value_list[start_idx]))
                file_rows.append(file_row)
            else:
                pass
        else:
            pass


outfile = open('COBRA_extended_equal_paths_contigs_info.txt', 'w')
title = ['QuerySeqID', 'Dominant', 'DominantLen', 'DominantCov', 'LessDominant', 'LessDominantLen', 'LessDominantCov', 'Category']
print('\t'.join(title), file=outfile, flush=True)
for file_row in file_rows:
    print('\t'.join(file_row[:]), file=outfile, flush=True)
outfile.close()


# get the less dominant paths of extended_circular queries
dominant2lessdominant_keys = []
for key in dominant2lessdominant.keys():
    dominant2lessdominant_keys.append(key)

for item in dominant2lessdominant_keys: # some of the equal paths may have different directions
    if item.endswith('_Lrc') and dominant2lessdominant[item].endswith('_R'):
        dominant2lessdominant[contig_name(item) + '_L'] = contig_name(dominant2lessdominant[item]) + '_Rrc'
    elif item.endswith('_Rrc') and dominant2lessdominant[item].endswith('_L'):
        dominant2lessdominant[contig_name(item) + '_R'] = contig_name(dominant2lessdominant[item]) + '_Lrc'
    elif item.endswith('_Rrc') and dominant2lessdominant[item].endswith('_Rrc'):
        dominant2lessdominant[contig_name(item) + '_R'] = contig_name(dominant2lessdominant[item]) + '_R'
    elif item.endswith('_Lrc') and dominant2lessdominant[item].endswith('_Lrc'):
        dominant2lessdominant[contig_name(item) + '_L'] = contig_name(dominant2lessdominant[item]) + '_L'
    elif item.endswith('_R') and dominant2lessdominant[item].endswith('_R'):
        dominant2lessdominant[contig_name(item) + '_Rrc'] = contig_name(dominant2lessdominant[item]) + '_Rrc'
    elif item.endswith('_L') and dominant2lessdominant[item].endswith('_L'):
        dominant2lessdominant[contig_name(item) + '_Lrc'] = contig_name(dominant2lessdominant[item]) + '_Lrc'

os.mkdir('COBRA_extended_circular_rare_unique')
order_all_rare = {}
extended_circular_unique = [] # only the unique of those with equal path will be written below
rare_set = set()
with open('COBRA_extended_circular_dominant_unique.fasta', 'r') as f:
    for line in f.readlines():
        if line.startswith('>'):
            extended_circular_unique.append(line.strip()[1:].split('_retrieved')[0])
f.close()

for contig in extended_circular_query:
    if contig in order_all.keys() and os.path.exists('COBRA_extended_circular_dominant_unique/{0}_retrieved.fa'.format(contig)):
        order_all_rare[contig] = []
        for item in order_all[contig]:
            if item.endswith('_R') or item.endswith('_L') or item.endswith('_Lrc') or item.endswith('_Rrc'):
                if item not in dominant2lessdominant.keys():
                    order_all_rare[contig].append(item)
                else:
                    order_all_rare[contig].append(dominant2lessdominant[item])
                    rare_set.add(contig) # only those with equal path will be written below
            else:
                order_all_rare[contig].append(item)

for contig in rare_set:
    a = open('COBRA_extended_circular_rare_unique/{0}_rare_retrieved.fa'.format(contig), 'w')
    for item in order_all_rare[contig]:
        print('>' + contig_name(item), file=a, flush=True)
        print(header2seq[contig_name(item)], file=a, flush=True)
    a.close()

for contig in rare_set:
    a = open('COBRA_extended_circular_rare_unique/{0}_rare_retrieved_joined.fa'.format(contig), 'w')
    print('>' + contig + '_extended_circular_rare', file=a, flush=True)
    for item in order_all_rare[contig][:-1]:
        if item.endswith('_R') or item.endswith('_L'):
            print(header2seq[item.rsplit('_', 1)[0]][:-length], end='', file=a, flush=True)
        elif item.endswith('rc'):
            print(reverse_complement(header2seq[item.rsplit('_', 1)[0]])[:-length], end='', file=a, flush=True)
        else:
            print(header2seq[contig][:-length], end='', file=a, flush=True)
    if order_all_rare[contig][-1].endswith('rc'):
        print(reverse_complement(header2seq[order_all_rare[contig][-1].rsplit('_',1)[0]])[:-length], file=a, flush=True)
    elif order_all_rare[contig][-1].endswith('_R') or order_all_rare[contig][-1].endswith('_L'):
        print(header2seq[order_all_rare[contig][-1].rsplit('_',1)[0]][:-length], file=a, flush=True)
    else:
        print(header2seq[order_all_rare[contig][-1]][:-length], file=a, flush=True)
    a.close()

os.system('cat COBRA_extended_circular_rare_unique/*joined.fa >COBRA_extended_circular_rare_unique.fasta')
summary_fasta('COBRA_extended_circular_rare_unique.fasta')


# write the numbers to the log file
print('\n', end='', file=log, flush=True)
print('=' * 150, file=log, flush=True)
print('Final summary', file=log, flush=True)
print('Total queries: ' + str(calculate_seq_num('../{0}'.format(args.query))) + '\n' + 'Self circular sequences: ' +
      str(calculate_seq_num('COBRA_self_circular_queries_trimmed.fasta'))+ '\n' + 'Extended circular: ' + str(len(extended_circular_query))
      + ' (Unique: ' + str(calculate_seq_num('COBRA_extended_circular_dominant_unique.fasta')) + ')' + '\n' + 'Extended >= 10k: ' + str(len(extended_10k_query))
      + ' (Unique: ' + str(calculate_seq_num('COBRA_extended_10k_unique.fasta')) + ')' + '\n' + 'Extended < 10k: ' + str(len(extended_few_query)) + '\n' +
      'Failed to extended: ' + str(calculate_seq_num('COBRA_extended_failed.fasta')) + '\n' +
      'Failed due to DNA break or short contig missing: ' + str(len(DNA_break_query)) + '\n' +
      'Failed due to unexpected break by assembler: ' + str(len(unexpected_break_query)), file=log, flush=True)
print('=' * 150, file=log, flush=True)
log.close()

