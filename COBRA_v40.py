#!/home/linking/.pyenv/versions/3.8.2/bin/python
# Author: Lin-Xing Chen, UC Berkeley

# COBRA v1.1.33
# Contig Overlap Based Re-Assembly
# Modification date: Mar 14, 2021

# updated from v1.1.23
# Modification: less strict in determining a valid join at two aspects
# - 1. "is_ok_to_add": the candidate does not have to be a query.
# - 2. remove the "-qi' flag, given that no prediction tool can predict all virus contigs with a minimum length of 10k from a metagenome assembly

# updated from v1.1.24
# add the exact extended length to the summary file

# updated from v1.1.25
# a BLASTN section to identify mis-join of contigs from genomes with similar direct terminal repeats

# updated from v1.1.26
# remove "AssembledSeqs" for summary file

# updated from v1.1.27
# modify the BLASTn parameters to filter invalid joins

# updated from v1.1.32
# function of "other_end_is_extendable" range modified from 0.333-2.2 to 0.333-2

# updated from v1.1.33
# function of "other_end_is_extendable" range modified from 0.333-2 to 0.5-2

# updated from v1.1.34
# when determining if two contigs are spanned by paired reads, the mapping mismatch for both should NOT be larger than 2

# updated from v1.1.35
# coverage must be provided even using metaSPAdes for assembly
# make invalid check possible for self_circular ones as well
# fix problems in detecting extended_circular paths
# extended_ok (previously >= 10k + < 10k)

# Usage: check COBRA.py -h


import os
from Bio import SeqIO
from Bio.Seq import reverse_complement
from collections import defaultdict
from Bio.SeqUtils import GC
import argparse
import pysam
from time import strftime


parser = argparse.ArgumentParser(description="This script is used to get higher quality (including circular) virus genomes by combining assembled contigs based on their end overlaps.")
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("-q", "--query", type=str, help="the query contigs file (fasta format).", required=True)
requiredNamed.add_argument("-f", "--fasta", type=str, help="the whole set of assembled contigs (fasta format).", required=True)
requiredNamed.add_argument("-a", "--assembler", type=str, choices=["idba", "megahit", "metaspades"], help="de novo assembler used, COBRA not tested for others.", required=True)
requiredNamed.add_argument("-k", "--maxk", type=int, help="the max kmer size used in de novo assembly.", required=True)
requiredNamed.add_argument("-m", "--mapping", type=str, help="the reads mapping file in sam or bam format.", required=True)
requiredNamed.add_argument("-c", "--coverage", type=str, help="the contig coverage file (two columns divided by tab).")
parser.add_argument("-mm", "--mismatch", type=int, default=2, help="the max read mapping mismatches (default, 2) when determining if two contigs are spanned by paired reads.")
parser.add_argument("-o", "--output", type=str, help="the name of output folder (default = '<query>_COBRA').")
parser.add_argument("-v", "--version", action='version', version='COBRA v1.1.19')
args = parser.parse_args()


def log_info(step, description, line_feed, log_file):
    localtime = ' [' + strftime("%c") + '] '
    print(step + localtime + description, end = line_feed, file=log_file, flush=True)


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


def is_linked(a, b):  # to evaluate if two contigs are linked based on linkage information
    my_set = set()
    my_set.add(contig_name(a))
    my_set.add(contig_name(b))
    return my_set in linkage_parsed


def could_circulate(point, contig, direction):  # check if the path is circular with the current contig included
    if len(link_pair[point]) == 2:  # point is the same as "target" in join_walker
        if not are_equal_paths(link_pair[point]):
            if direction == 'L':
                if contig_name(link_pair[point][0]) == contig:
                    if cov[contig_name(point)] < 1.5 * cov[contig]:
                        # 2 times is for repeat, but it is too risky, use 1.5 instead (same below)
                        return link_pair[point][0].endswith('_R') and is_linked(contig, point)
                    else:
                        pass
                elif contig_name(link_pair[point][1]) == contig:
                    if cov[contig_name(point)] < 1.5 * cov[contig]:
                        return link_pair[point][1].endswith('_R') and is_linked(contig, point)
                    else:
                        pass
            elif direction == 'R':
                if contig_name(link_pair[point][0]) == contig:
                    if cov[contig_name(point)] < 1.5 * cov[contig]:
                        return link_pair[point][0].endswith('_L') and is_linked(contig, point)
                    else:
                        pass
                elif contig_name(link_pair[point][1]) == contig:
                    if cov[contig_name(point)] < 1.5 * cov[contig]:
                        return link_pair[point][1].endswith('_L') and is_linked(contig, point)
                    else:
                        pass
        else:
            pass
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


def get_target(item):  # to get the target for next run of joining
    if item.rsplit('_', 1)[1] == 'Lrc':
        return item.rsplit('_', 1)[0] + '_Rrc'
    elif item.rsplit('_', 1)[1] == 'Rrc':
        return item.rsplit('_', 1)[0] + '_Lrc'
    elif item.rsplit('_', 1)[1] == 'R':
        return item.rsplit('_', 1)[0] + '_L'
    else:
        return item.rsplit('_', 1)[0] + '_R'


def other_end_is_extendable(end, contig):  # check if the other end is extendable
    if contig_name(end) not in self_circular:
        if 0.5 * cov[contig] <= cov[contig_name(end)] <= 2 * cov[contig]:
            if get_target(end) in link_pair.keys():
                return len(link_pair[get_target(end)]) > 0
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
        return 0.9 * cov[contig] <= cov[contig_name(end)] <= 1.11


def not_checked(end_list, checked):  # to see if a potential contig has been check for adding or not
    total_checked = 0
    for end in end_list:
        if contig_name(end) in checked:
            total_checked += 1
    return total_checked == 0


def percentage_finished(finished, total, log_file):
    p_finished = int(len(finished)/2/total*100)
    if p_finished in range(1, 101) and p_finished % 5 == 0 and p_finished not in p_printed:
        print('{0}%, '.format(p_finished), end='', file=log_file, flush=True)
        p_printed.append(p_finished)
    else:
        pass


def join_walker(contig, direction):
    # get potential joins for a given query
    end = contig + '_' + direction
    a = len(contig2join[end])
    if a == 0:
        contig_checked[end].append(contig)
        percentage_finished(finished_end, len(query_list)-len(DNA_break_query), log)
        if end in one_join_end:
            if contig_name(link_pair[end][0]) != contig:
                if other_end_is_extendable(link_pair[end][0], contig) and is_linked(contig, link_pair[end][0]):
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
                    self_circular.add(contig)
                elif direction == 'R' and link_pair[end][0] == contig + '_L':
                    # the other end could be joined with the current working point
                    self_circular.add(contig)
        elif end in two_joins_end and link_pair[end][0].rsplit('_', 1)[0] != link_pair[end][1].rsplit('_', 1)[0]:
            if are_equal_paths(link_pair[end]):
                contig2join[end].append(the_dominant_one(link_pair[end]))
                contig_checked[end].append(contig_name((link_pair[end][0])))
                contig_checked[end].append(contig_name((link_pair[end][1])))
            else:
                if direction == '_L' and contig + '_R' in link_pair[end]:
                    self_circular.add(contig)
                elif direction == 'R' and contig + '_L' in link_pair[end]:
                    self_circular.add(contig)
                elif cov[contig_name(link_pair[end][0])] + cov[contig_name(link_pair[end][1])] >= cov[contig] * 0.5:
                    # 0.5 is ok, too big will get much fewer "Extended circular" ones.
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
                    contig2join[end].append(link_pair[target][0])
                    contig_checked[end].append(contig_name(link_pair[target][0]))
                elif is_ok_to_add(link_pair[target][0], contig) and is_linked(target, link_pair[target][0]):
                    contig2join[end].append(link_pair[target][0])
                    contig_checked[end].append(contig_name(link_pair[target][0]))
                else:
                    pass
            else:
                if could_circulate(target, contig, direction):
                    path_circular.add(contig)
                    path_circular_end.add(contig + '_' + direction)
        elif target in two_joins_end and link_pair[target][0].rsplit('_', 1)[0] != link_pair[target][1].rsplit('_', 1)[0]:
            if not_checked(link_pair[target], contig_checked[end]):
                if are_equal_paths(link_pair[target]):
                    contig2join[end].append(the_dominant_one(link_pair[target]))
                    contig_checked[end].append(contig_name(link_pair[target][0]))
                    contig_checked[end].append(contig_name(link_pair[target][1]))
                else:
                    if cov[contig_name(link_pair[target][0])] + cov[contig_name(link_pair[target][1])] >= cov[contig] * 0.5:
                        # 0.5 is ok, too big will get much fewer "Extended circular" ones.
                        contig2join[end].append(the_better_one(link_pair[target], contig))
                        contig_checked[end].append(contig_name((link_pair[target][0])))
                        contig_checked[end].append(contig_name((link_pair[target][1])))
            else:
                if could_circulate(target, contig, direction):
                    path_circular.add(contig)
                    path_circular_end.add(contig + '_' + direction)
        else:
            pass
    if a == len(contig2join[end]):
        finished_end.add(end)
    return a < len(contig2join[end])


def join_seqs(contig):  # get the join order of the sequences in a give path
    order_all[contig] = []
    left = contig + '_L'
    right = contig + '_R'
    added_to_contig[contig] = []

    if contig not in path_circular:
        if left in contig2join.keys() and right in contig2join.keys():
            for item in contig2join[left][::-1]: # the order of contigs added into the left direction should be reversed
                order_all[contig].append(item)
                added_to_contig[contig].append(contig_name(item))

            order_all[contig].append(contig) # the query contig itself should be added to the path as well

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
            #if get_contigs_in_path(contig2join[left]) == get_contigs_in_path(contig2join[right]):
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


def get_contigs_in_path(end_list):
    contigs = set()
    for end in end_list:
        contigs.add(end.rsplit('_',1)[0])
    return contigs


def retrieve(contig):  # to retrieve and save all the contigs in the joining path of a query
    if contig in order_all.keys() and len(order_all[contig]) > 0:
        out = open('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(contig), 'w')
        added = []
        for item in order_all[contig]:
            if contig_name(item) not in added:
                print('>' + contig_name(item), file=out)
                print(header2seq[contig_name(item)], file=out)
                added.append(contig_name(item))
            else:
                pass
    else:
        pass


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
        summary_file_headers = ['SeqID', 'Length', 'Coverage', 'GC', 'Ns', 'DTR_length']
        print('\t'.join(summary_file_headers[:]), file=summary_file, flush=True)
    else:
        summary_file_headers = ['SeqID', 'Length', 'Coverage', 'GC', 'Ns']
        print('\t'.join(summary_file_headers[:]), file=summary_file, flush=True)

    with open('{0}'.format(fasta_file), 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            header = str(record.id).strip()
            seq = str(record.seq)
            Ns = seq.count('N')
            if header.split('_self')[0] in self_circular:
                sequence_stats = [header, str(len(seq)), str(cov[header.split('_self')[0]]), str(round(GC(seq), 3)), str(Ns), str(length)]
                print('\t'.join(sequence_stats[:]), file=summary_file, flush=True)
            elif header.split('_self')[0] in self_circular_non_expected_overlap.keys():
                sequence_stats = [header, str(len(seq)), str(cov[header.split('_self')[0]]), str(round(GC(seq), 3)), str(Ns),
                                  str(self_circular_non_expected_overlap[header.split('_self')[0]])]
                print('\t'.join(sequence_stats[:]), file=summary_file, flush=True)
            else:
                sequence_stats = [header, str(len(seq)), str(cov[header.split('_extended')[0]]), str(round(GC(seq), 3)), str(Ns)]
                print('\t'.join(sequence_stats[:]), file=summary_file, flush=True)
    f.close()
    summary_file.close()


def summarize(contig):  # summary the retrieved contigs and joined information of a query
    if contig not in is_redundant_of.keys():
        b = calculate_seq_num('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(contig))  # number of retrieved contigs
        c = calculate_seq_len('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(contig))  # total length of retrieved contigs
        d = calculate_seq_len('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(contig))  # total length after joining
        if contig in extended_circular_query:
            return '\t'.join([str(b), str(c), str(d), str(d-header2len[contig]), 'Extended_circular'])
        elif contig in extended_ok_query:
            return '\t'.join([str(b), str(c), str(d), str(d-header2len[contig]), 'Extended_ok'])
    else:
        if is_redundant_of[contig] in is_redundant_of.keys():
            item = is_redundant_of[is_redundant_of[contig]]
            b = calculate_seq_num('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(item))  # number of retrieved contigs
            c = calculate_seq_len('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(item))  # total length of retrieved contigs
            d = calculate_seq_len('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(item))  # total length after joining
            if contig in extended_circular_query:
                return '\t'.join([str(b), str(c), str(d), str(d-header2len[contig]), 'Extended_circular'])
            elif contig in extended_ok_query:
                return '\t'.join([str(b), str(c), str(d), str(d-header2len[contig]), 'Extended_ok'])
        else:
            item = is_redundant_of[contig]
            b = calculate_seq_num('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(item))  # number of retrieved contigs
            c = calculate_seq_len('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(item))  # total length of retrieved contigs
            d = calculate_seq_len('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(item))  # total length after joining
            if contig in extended_circular_query:
                return '\t'.join([str(b), str(c), str(d), str(d-header2len[contig]), 'Extended_circular'])
            elif contig in extended_ok_query:
                return '\t'.join([str(b), str(c), str(d), str(d-header2len[contig]), 'Extended_ok'])


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
    print('Output folder <{0}> exists, please check.'.format(working_dir))
    exit()
else:
    os.mkdir('{0}'.format(working_dir))


# determine the length of overlap based on assembler and the largest kmer size
if args.assembler == "idba":
    length = args.maxk - 1
else:
    length = args.maxk

# determine the number of mismatches allowed for spanning paired reads
if not args.mismatch:
    mismatch = '2'
else:
    mismatch = str(args.mismatch)


# log file
log = open('{0}/log'.format(working_dir), 'w')  # log file
star_num = int((150-len(' COBRA analyses for {0} '.format(query_name)))/2)
log.write('*' * star_num + ' COBRA analyses for {0} '.format(query_name) + '*' * star_num + '\n' + '\n')

if args.assembler == 'idba':
    parameters = ['# Key parameters:', '# Assembler: IDBA_UD', '# Max-kmer: ' + str(args.maxk).strip(),
                  '# Overlap length: ' + str(length) + ' bp', '# Mismatch: ' + mismatch, '\n']
elif args.assembler == 'metaspades':
    parameters = ['# Key parameters:', '# Assembler: metaSPAdes', '# Max-kmer: ' + str(args.maxk).strip(),
                  '# Overlap length: ' + str(length) + ' bp', '# Mismatch: ' + mismatch, '\n']
else:
    parameters = ['# Key parameters:', '# Assembler: MEGAHIT', '# Max-kmer: ' + str(args.maxk).strip(),
                  '# Overlap length: ' + str(length) + ' bp', '# Mismatch: ' + mismatch, '\n']

log.write('\n'.join(parameters[:]))


# import the contigs in the pool and save the end sequence
log_info('[01/20]', 'Reading contigs and getting contig ends ... ', '', log)
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
    gc[header] = str(round(GC(seq), 3))
    header2len[header] = len(seq)
    L[header + '_L'] = seq[0:length]  # the first x bp of left end
    Lrc[header + '_Lrc'] = reverse_complement(seq[0:length])  # the reverse sequence of first x bp of left end
    R[header + '_R'] = seq[-length:]  # the first x bp of right end
    Rrc[header + '_Rrc'] = reverse_complement(seq[-length:])  # the reverse sequence of first x bp of right end
f.close()
log.write('A total of {0} contigs were imported.'.format(len(header2seq.keys())) + '\n')


# get potential joins
log_info('[02/20]', 'Getting shared contig ends ...', '\n', log)

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


# save all paired links to a file
log_info('[03/20]', 'Writing contig end joining pairs ...', '\n', log)
one_join_end = []  # the end of contigs with one potential join
two_joins_end = []  # the end of contigs with two potential joins

p = open('{0}/COBRA_end_joining_pairs.txt'.format(working_dir, fasta_name), 'w')
for item in link_pair.keys():
    for point in link_pair[item]:
        p.write(item + '\t' + point + '\n')  # print link pairs into a file for check if interested

    if len(link_pair[item]) == 1:
        one_join_end.append(item)  # add one joining end to a list, its pair may have one or more joins
    elif len(link_pair[item]) == 2 and len(link_pair[link_pair[item][0]]) == 1 and len(link_pair[link_pair[item][1]]) == 1:
        two_joins_end.append(item)  # add two joining end to a list, each of its pairs should only have one join
    else:
        pass
p.close()


# coverage
log_info('[04/20]', 'Getting contig coverage information ...', '\n', log)
cov = {}
coverage = open('{0}'.format(args.coverage), 'r')
for line in coverage.readlines():
    line = line.strip().split('\t')
    cov[line[0]] = float(line[1])
coverage.close()


# open query file and save the information
log_info('[05/20]', 'Getting query contig list ... ', '', log)
query_list = []
DNA_break_query = []
q = open('{0}'.format(args.query), 'r')
for record in SeqIO.parse(q, "fasta"):
    header = str(record.id).strip()
    if header in header2seq.keys():
        if header not in query_list:
            query_list.append(header)
            if header + '_L' not in link_pair.keys() and header + '_R' not in link_pair.keys():
                DNA_break_query.append(header)
            else:
                pass
        else:
            pass
    else:
        print('Query {0} is not in your whole contig fasta file, please check!'.format(header), flush=True)
q.close()

log.write('{0} query contigs were imported ...'.format(len(query_list)) + '\n')


# get the linkage of contigs based on paired-end reads mapping
log_info('[06/20]', 'Getting contig linkage based on sam/bam ... Be patient, this may take long ... ', '\n', log)
linkage = {}
map_file = pysam.AlignmentFile('{0}'.format(args.mapping))
for line in map_file:
    if line.reference_name != line.next_reference_name and line.get_tag("NM") <= args.mismatch:  # mismatch should not be more than 2
        if line.query_name not in linkage.keys():
            linkage[line.query_name] = []
            pair = set()
            pair.add(line.reference_name)
            pair.add(line.next_reference_name)
            linkage[line.query_name].append(pair)
        else:
            pair = set()
            pair.add(line.reference_name)
            pair.add(line.next_reference_name)
            linkage[line.query_name].append(pair)
    else:
        pass
map_file.close()

linkage_parsed = []
for read in linkage.keys():
    if len(linkage[read]) == 2:  # 2 means both the paired end that spanning two contigs have mismatches fewer than the number indicated above
        linkage_parsed.append(linkage[read][0])
        linkage_parsed.append(linkage[read][1])
    else:
        pass


# check valid joins
contig2join = {}
contig_checked = {}
path_circular = set()
self_circular = set()

for contig in query_list:
    contig2join[contig + '_L'] = []
    contig2join[contig + '_R'] = []
    contig_checked[contig + '_L'] = []
    contig_checked[contig + '_R'] = []

# for debug
print(contig2join, flush=True)

# walk the joins
log_info('[07/20]', 'Walking joins... ', '', log)
finished_end = set()
p_printed = []
path_circular_end = set()
for contig in query_list:
    if contig not in DNA_break_query:
        while join_walker(contig, 'L'):
            join_walker(contig, 'L')
        while join_walker(contig, 'R'):
            join_walker(contig, 'R')
    else:
        pass


# save the potential joining paths
log.write('100% finished ...' + '\n')
log_info('[08/20]', 'Saving potential joining paths ...', '\n', log)
results = open('{0}/COBRA_potential_joining_paths.txt'.format(working_dir), 'w')
for item in contig2join.keys():
    if contig_name(item) in self_circular:
        if item.endswith('_L'):
            results.write(item + '\t' + "['" + contig_name(item) + '_R' + "']" + '\n')
        else:
            results.write(item + '\t' + "['" + contig_name(item) + '_L' + "']" + '\n')
    else:
        results.write(item + '\t' + str(contig2join[item]) + '\n')


# DNA break info
log_info('[09/20]', 'Saving contig DNA break information ...', '\n', log)
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


# get the joining paths
log_info('[10/20]', 'Checking for invalid joining - sharing queries ...', '\n', log)
contig2assembly = {}
for item in contig2join.keys():
    contig = item.rsplit('_', 1)[0]
    if contig not in contig2assembly.keys():
        contig2assembly[contig] = set()
        contig2assembly[contig].add(contig)
        for point in contig2join[item]:
            contig2assembly[contig].add(contig_name(point))
    else:
        for point in contig2join[item]:
            contig2assembly[contig].add(contig_name(point))


# find the redundant joining paths
redundant = set()
is_redundant_of = {}
is_the_same_as = set()
for item in contig2assembly.keys():
    for item_1 in contig2assembly.keys():
        if item != item_1 and contig2assembly[item].issubset(contig2assembly[item_1]):
            if contig2assembly[item] != contig2assembly[item_1]:
                if item in path_circular and item_1 not in path_circular:
                    if item_1 not in contig2assembly[item]:
                        failed_join_list.append(item)
                        failed_join_list.append(item_1)
                        path_circular.remove(item)
                    else:
                        redundant.add(item_1)
                        path_circular.add(item_1)
                        if item not in is_redundant_of.keys():
                            is_redundant_of[item_1] = item
                        else:
                            is_redundant_of[item_1] = is_redundant_of[item]
                elif item_1 in path_circular and item not in path_circular:
                    path_circular.add(item)
                    redundant.add(item)
                    if item_1 not in is_redundant_of.keys():
                        is_redundant_of[item] = item_1
                    else:
                        is_redundant_of[item] = is_redundant_of[item_1]
                elif item in path_circular and item_1 in path_circular:
                    failed_join_list.append(item)
                    failed_join_list.append(item_1)
                    path_circular.remove(item)
                    path_circular.remove(item_1)
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

for scaffold in is_redundant_of.keys():
    if is_redundant_of[scaffold] in failed_join_list:
        failed_join_list.append(scaffold)

# for debug
print(redundant, flush=True)
print(is_redundant_of, flush=True)
print(is_the_same_as, flush=True)

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
    else:
        pass

    for item in is_redundant_of.keys():
        if contig == is_redundant_of[item]:
            del contig2assembly[item]
            failed_join_list.append(item)

            for item_1 in is_redundant_of.keys():
                if is_redundant_of[item_1] == item:
                    del contig2assembly[item_1]
                    failed_join_list.append(item_1)
        else:
            pass


# get the joining order of contigs
log_info('[11/20]', 'Getting the joining order of contigs ...', '\n', log)
order_all = {}
added_to_contig = {}

for contig in contig2assembly.keys():
    # only those contigs left in contig2assembly after filtering
    # (see above "# remove the queries in multiple paths") will be
    # checked for join paths (join_seqs) to get order_all
    if len(contig2assembly[contig]) > 1 and contig not in failed_join_list and contig not in redundant: # and contig not in is_the_same_as:
        join_seqs(contig)


# get retrieved sequences
log_info('[12/20]', 'Getting retrieved contigs ...', '\n', log)
os.chdir('{0}'.format(working_dir))
os.mkdir('COBRA_retrieved_for_joining')
retrieved = []
for contig in order_all.keys():
    retrieve(contig)
    retrieved.append(contig)
    # for debug


# for debug
print(order_all, flush=True)
print(path_circular, flush=True)
print(path_circular_end, flush=True)
print(retrieved, flush=True)


# writing joined sequences
log_info('[13/20]', 'Saving joined seqeuences ...', '\n', log)
for contig in retrieved:
    # for debug
    print(contig, flush=True)
    a = open('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(contig), 'w')
    last = ''
    # print header regarding the joining status
    if contig in path_circular:
        a.write('>' + contig + '_extended_circular' + '\n')
    else:
        a.write('>' + contig + '_extended_partial' + '\n')
    # print the sequences with their overlap removed
    for item in order_all[contig][:-1]:
        if item.endswith('_R') or item.endswith('_L'):
            if last == '':
                a.write(header2seq[item.rsplit('_', 1)[0]][:-length])
                last = header2seq[item.rsplit('_', 1)[0]][-length:]
            else:
                if header2seq[item.rsplit('_', 1)[0]][:length] == last:
                    a.write(header2seq[item.rsplit('_', 1)[0]][:-length])
                    last = header2seq[item.rsplit('_', 1)[0]][-length:]
                else:
                    pass
        elif item.endswith('rc'):
            if last == '':
                a.write(reverse_complement(header2seq[item.rsplit('_', 1)[0]])[:-length])
                last = reverse_complement(header2seq[item.rsplit('_', 1)[0]])[-length:]
            else:
                if reverse_complement(header2seq[item.rsplit('_', 1)[0]])[:length] == last:
                    a.write(reverse_complement(header2seq[item.rsplit('_', 1)[0]])[:-length])
                    last = reverse_complement(header2seq[item.rsplit('_', 1)[0]])[-length:]
                else:
                    pass
        else:
            if last == '':
                a.write(header2seq[contig][:-length])
                last = header2seq[contig][-length:]
            else:
                if header2seq[contig][:length] == last:
                    a.write(header2seq[contig][:-length])
                    last = header2seq[contig][-length:]
                else:
                    pass

    if contig in path_circular:
        if order_all[contig][-1].endswith('rc'):
            a.write(reverse_complement(header2seq[order_all[contig][-1].rsplit('_',1)[0]])[:-length] + '\n')
        elif order_all[contig][-1].endswith('_R') or order_all[contig][-1].endswith('_L'):
            a.write(header2seq[order_all[contig][-1].rsplit('_',1)[0]][:-length] + '\n')
        else:
            a.write(header2seq[order_all[contig][-1]][:-length] + '\n')
    else:
        if order_all[contig][-1].endswith('rc'):
            a.write(reverse_complement(header2seq[order_all[contig][-1].rsplit('_',1)[0]]) + '\n')
        elif order_all[contig][-1].endswith('_R') or order_all[contig][-1].endswith('_L'):
            a.write(header2seq[order_all[contig][-1].rsplit('_',1)[0]] + '\n')
        else:
            a.write(header2seq[order_all[contig][-1]] + '\n')

    a.close()


# determine the joining status of queries
log_info('[14/20]', 'Getting initial joining status of each query contig ...', '\n', log)
extended_circular_query = set()
extended_ok_query = set()

for contig in contig2assembly.keys():
    if contig not in failed_join_list:
        if contig in redundant:
            if is_redundant_of[contig] in path_circular:
                extended_circular_query.add(contig)
            else:
                extended_ok_query.add(contig)
        else:
            if contig in path_circular:
                extended_circular_query.add(contig)
            else:
                if contig not in failed_join_list and contig not in DNA_break_query and contig not in self_circular and contig not in self_circular_non_expected_overlap.keys():
                    extended_ok_query.add(contig)
    else:
        pass


# deal with cross-assignment queries
log_info('[15/20]', 'Getting final joining status of each query contig ...', '\n', log)
for contig in failed_join_list:
    if contig in is_redundant_of.keys() and is_redundant_of[contig] in extended_circular_query:
        failed_join_list.remove(contig)
        extended_circular_query.add(contig)
    elif contig in is_redundant_of.keys() and is_redundant_of[contig] in extended_ok_query:
        failed_join_list.remove(contig)
        extended_ok_query.add(contig)
    else:
        pass


for contig in extended_ok_query:
    if contig in redundant:
        if is_redundant_of[contig] in extended_circular_query:
            extended_ok_query.remove(contig)
            extended_circular_query.add(contig)
        elif is_redundant_of[contig] in failed_join_list:
            extended_ok_query.remove(contig)
            failed_join_list.add(contig)
        else:
            pass
    else:
        pass


# Similar direct terminal repeats may lead to invalid joins
log_info('[16/20]', 'Checking for invalid joining using BLASTn ...', '\n', log)
for_blastn = open('for.blastn.fasta', 'w')
cobraSeq2len = {}

for contig in retrieved:
    a = open('COBRA_retrieved_for_joining/{0}_retrieved_joined.fa'.format(contig), 'r')
    for record in SeqIO.parse(a, "fasta"):
        header = str(record.id).strip()
        seq = str(record.seq)
        cobraSeq2len[header.split('_extended',1)[0]] = len(seq)

        if len(seq) % 2 == 0:
            half = int(len(seq) / 2)
        else:
            half = int((len(seq) + 1) / 2)

        for_blastn.write('>' + header + '_1' + '\n')
        for_blastn.write(seq[:half] + '\n')
        for_blastn.write('>' + header + '_2' + '\n')
        for_blastn.write(seq[half:] + '\n')

    a.close()

for contig in self_circular:
    cobraSeq2len[contig] = header2len[contig]-length

    if header2len[contig] % 2 == 0:
        half = int(header2len[contig] / 2)
    else:
        half = int((header2len[contig] + 1) / 2)

    for_blastn.write('>' + contig + '_1' + '\n')
    for_blastn.write(header2seq[contig][:half] + '\n')
    for_blastn.write('>' + contig + '_2' + '\n')
    for_blastn.write(header2seq[contig][half:] + '\n')

for contig in self_circular_non_expected_overlap.keys():
    cobraSeq2len[contig] = header2len[contig]-self_circular_non_expected_overlap[contig]

    if header2len[contig] % 2 == 0:
        half = int(header2len[contig] / 2)
    else:
        half = int((header2len[contig] + 1) / 2)

    for_blastn.write('>' + contig + '_1' + '\n')
    for_blastn.write(header2seq[contig][:half] + '\n')
    for_blastn.write('>' + contig + '_2' + '\n')
    for_blastn.write(header2seq[contig][half:] + '\n')

for_blastn.close()


os.system('makeblastdb -in for.blastn.fasta -dbtype nucl')
os.system('blastn -task blastn -db for.blastn.fasta -query for.blastn.fasta -out for.blastn.fasta.blastn.self -evalue 1e-10 -outfmt 6 -perc_identity 70 -num_threads 1')

contig2TotLen = {}
r = open('for.blastn.fasta.blastn.self', 'r')
for line in r.readlines():
    line = line.strip().split('\t')
    if line[0].rsplit('_', 1)[0] == line[1].rsplit('_', 1)[0] and line[0] != line[1] and line[0].rsplit('_', 1)[1] == '1':
            if float(line[3]) >= 1000:
                if '_extended' in line[0]:
                    if line[0].split('_extended')[0] not in contig2TotLen.keys():
                        contig2TotLen[line[0].split('_extended')[0]] = float(line[3])
                    else:
                        contig2TotLen[line[0].split('_extended')[0]] += float(line[3])
                else:
                    if line[0].rsplit('_',1)[0] not in contig2TotLen.keys():
                        contig2TotLen[line[0].rsplit('_',1)[0]] = float(line[3])
                    else:
                        contig2TotLen[line[0].rsplit('_',1)[0]] += float(line[3])
            else:
                pass
    else:
        pass
r.close()


for contig in contig2TotLen.keys():
    if contig2TotLen[contig] / cobraSeq2len[contig] >= 0.05:  # previously, if contig2TotLen[contig] >= 1000:
        if os.path.exists('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(contig)):
            a = open('COBRA_retrieved_for_joining/{0}_retrieved.fa'.format(contig), 'r')
            for record in SeqIO.parse(a, "fasta"):
                header = str(record.id).strip()
                if header in extended_ok_query:
                    extended_ok_query.remove(header)
                    failed_join_list.append(header)
                    if os.path.exists('COBRA_category_ii_extended_ok/{0}_retrieved*fa'.format(header)):
                        os.remove('COBRA_category_ii_extended_ok/{0}_retrieved*fa'.format(header))
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
print(extended_circular_query, flush=True)
print(extended_ok_query, flush=True)


# save the joining summary information
log_info('[17/20]', 'Saving joining summary of retrieved contigs ...', '\n', log)
assembly_summary = open('COBRA_joining_summary.txt', 'w')
assembly_summary_headers = ['QuerySeqID', 'QuerySeqLen', 'TotRetSeqs', 'TotRetLen', 'AssembledLen', 'ExtendedLen', 'Status']
assembly_summary.write('\t'.join(assembly_summary_headers[:]) + '\n')
for contig in list(extended_circular_query) + list(extended_ok_query):
    assembly_summary.write('\t'.join([contig, str(header2len[contig]), summarize(contig)]) + '\n')
assembly_summary.close()


# get the unique sequences of COBRA "well-extended" joining and their summary
log_info('[18/20]', 'Saving unique sequences of "Extended_circular" and "Extended_ok" ...', '\n', log)
os.mkdir('COBRA_category_ii_extended_ok_unique')
os.mkdir('COBRA_category_ii_extended_circular_unique')
for contig in query_list:
    if contig in extended_circular_query and contig not in redundant and contig not in is_the_same_as:
        os.system('cp COBRA_retrieved_for_joining/{0}_retrieved*fa COBRA_category_ii_extended_circular_unique'.format(contig))
    elif contig in extended_ok_query and contig not in redundant and contig not in is_the_same_as:
        os.system('cp COBRA_retrieved_for_joining/{0}_retrieved*fa COBRA_category_ii_extended_ok_unique'.format(contig))
    else:
        pass

os.system('cat COBRA_category_ii_extended_ok_unique/*joined.fa >COBRA_category_ii_extended_ok_unique.fasta')
os.system('cat COBRA_category_ii_extended_circular_unique/*joined.fa >COBRA_category_ii_extended_circular_unique.fasta')
summary_fasta('COBRA_category_ii_extended_ok_unique.fasta')
summary_fasta('COBRA_category_ii_extended_circular_unique.fasta')


# save the joining details information
log_info('[19/20]', 'Saving joining details of unique "Extended_circular" and "Extended_ok" queries ...', '\n', log)
joining_detail_headers = ['AssembledSeqID', 'NameForFig', 'AssembledLen', 'Status', 'RetSeqID', 'Direction', 'RetSeqLen', 'Start', 'End', 'RetSeqCov', 'RetSeqGC']
joining_detail_circular = open('COBRA_category_ii_extended_circular_unique_joining_details.txt', 'w')
print('\t'.join(joining_detail_headers[:]), file=joining_detail_circular, flush=True)
joining_detail_ok = open('COBRA_category_ii_extended_ok_unique_joining_details.txt', 'w')
print('\t'.join(joining_detail_headers[:]), file=joining_detail_ok, flush=True)


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
    elif contig in extended_ok_query and contig not in redundant and contig not in is_the_same_as:
        for item in order_all[contig][:-1]:
            if get_direction(item) == 'forward':
                contents = [contig + '_extended_partial', rename(contig + '_extended_partial'), str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                            'Partial', contig_name(item), 'forward', str(header2len[contig_name(item)]), str(site), str(site + header2len[contig_name(item)] - length - 1),
                            str(cov[contig_name(item)]), str(gc[contig_name(item)])]
                print('\t'.join(contents[:]), file=joining_detail_ok, flush=True)
                site += header2len[contig_name(item)] - length
            else:
                contents = [contig + '_extended_partial', rename(contig + '_extended_partial'), str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                            'Partial', contig_name(item) + '_rc', 'reverse', str(header2len[contig_name(item)]), str(site + header2len[contig_name(item)] - length - 1),
                            str(site), str(cov[contig_name(item)]), str(gc[contig_name(item)])]
                print('\t'.join(contents[:]), file=joining_detail_ok, flush=True)
                site += header2len[contig_name(item)] - length

        last = order_all[contig][-1]  # for the last one in non-circular path, the end position should be different
        if get_direction(last) == 'forward':
            contents = [contig + '_extended_partial', rename(contig + '_extended_partial'), str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                        'Partial', contig_name(last), 'forward', str(header2len[contig_name(last)]), str(site), str(site + header2len[contig_name(last)] - 1),
                        str(cov[contig_name(last)]), str(gc[contig_name(last)])]
            print('\t'.join(contents[:]), file=joining_detail_ok, flush=True)
        else:
            contents = [contig + '_extended_partial', rename(contig + '_extended_partial'), str(total_length(order_all[contig]) - length * (len(order_all[contig]) - 1)),
                        'Partial', contig_name(last) + '_rc', 'reverse', str(header2len[contig_name(last)]), str(site + header2len[contig_name(last)] - 1),
                        str(site), str(cov[contig_name(last)]), str(gc[contig_name(last)])]
            print('\t'.join(contents[:]), file=joining_detail_ok, flush=True)
    else:
        pass
joining_detail_ok.close()
joining_detail_circular.close()


## get the joining details for each unique extended_circular and extended_ok for gggenes figure
assembled2detail = {}
with open('COBRA_category_ii_extended_circular_unique_joining_details.txt', 'r') as joining_detail_circular:
    for line in joining_detail_circular.readlines():
        if not line.startswith('AssembledSeqID'):
            line = line.strip().split('\t')
            if line[0] not in assembled2detail.keys():
                assembled2detail[line[0]] = []
                assembled2detail[line[0]].append(line)
            else:
                assembled2detail[line[0]].append(line)
joining_detail_circular.close()

with open('COBRA_category_ii_extended_ok_unique_joining_details.txt', 'r') as joining_detail_ok:
    for line in joining_detail_ok.readlines():
        if not line.startswith('AssembledSeqID'):
            line = line.strip().split('\t')
            if line[0] not in assembled2detail.keys():
                assembled2detail[line[0]] = []
                assembled2detail[line[0]].append(line)
            else:
                assembled2detail[line[0]].append(line)
joining_detail_ok.close()


for contig in order_all.keys():
    if os.path.exists('COBRA_category_ii_extended_circular_unique/{0}_retrieved.fa'.format(contig)):
        a = open('COBRA_category_ii_extended_circular_unique/{0}_retrieved_joining_details.txt'.format(contig), 'w')
        print('\t'.join(joining_detail_headers[:]), file=a, flush=True)
        for item in assembled2detail[contig + '_extended_circular']:
            print('\t'.join(item[:]), file=a, flush=True)
        a.close()
    elif os.path.exists('COBRA_category_ii_extended_ok_unique/{0}_retrieved.fa'.format(contig)):
        a = open('COBRA_category_ii_extended_ok_unique/{0}_retrieved_joining_details.txt'.format(contig), 'w')
        print('\t'.join(joining_detail_headers[:]), file=a, flush=True)
        for item in assembled2detail[contig + '_extended_partial']:
            print('\t'.join(item[:]), file=a, flush=True)
        a.close()
    else:
        pass


# save the joining status information of each query
assembled_info = open('COBRA_joining_status.txt', 'w')  # shows the COBRA status of each query
print('SeqID' + '\t' + 'Length' + '\t' + 'Coverage' + '\t' + 'GC' + '\t' + 'Status' + '\t' + 'Category', file=assembled_info, flush=True)

# for those could be extended to circular
for contig in extended_circular_query:
    print(contig + '\t' + str(header2len[contig]) + '\t' + str(cov[contig]) + '\t' + gc[contig] + '\t' + 'Extended_circular' + '\t' + 'category_ii', file=assembled_info, flush=True)

# for those could be extended ok
for contig in extended_ok_query:
    print(contig + '\t' + str(header2len[contig]) + '\t' + str(cov[contig]) + '\t' + gc[contig] + '\t' + 'Extended_ok' + '\t' + 'category_ii', file=assembled_info, flush=True)

# for those cant be extended
failed_join = open('COBRA_category_ii_extended_failed.fasta', 'w')
for contig in set(failed_join_list):
    if contig not in extended_circular_query or contig not in extended_ok_query or contig not in DNA_break_query:
        print('>' + contig, file=failed_join, flush=True)
        print(header2seq[contig], file=failed_join, flush=True)
        print(contig + '\t' + str(header2len[contig]) + '\t' + str(cov[contig]) + '\t' + gc[contig] + '\t' + 'Extended_failed' + '\t' + 'category_ii', file=assembled_info, flush=True)
    else:
        pass
failed_join.close()
summary_fasta('COBRA_category_ii_extended_failed.fasta')

# for those due to DNA break
DNA_break = open('COBRA_category_iii_DNA_break.fasta', 'w')
for contig in DNA_break_query:
    print('>' + contig, file=DNA_break, flush=True)
    print(header2seq[contig], file=DNA_break, flush=True)
    print(contig + '\t' + str(header2len[contig]) + '\t' + str(cov[contig]) + '\t' + gc[contig] + '\t' + 'DNA_break' + '\t' + 'category_iii', file=assembled_info, flush=True)
DNA_break.close()
summary_fasta('COBRA_category_iii_DNA_break.fasta')


# for self circular
log_info('[20/20]', 'Saving self_circular contigs ...', '\n', log)
circular_fasta = open('COBRA_category_i_self_circular_queries_trimmed.fasta', 'w')
for contig in self_circular:
    print(contig + '\t' + str(header2len[contig] - length) + '\t' + str(cov[contig]) + '\t' + gc[contig] + '\t' + 'Self_circular' + '\t' + 'category_i', file=assembled_info, flush=True)
    print('>' + contig + '_self_circular', file=circular_fasta, flush=True)
    print(header2seq[contig][length:], file=circular_fasta, flush=True)
for contig in self_circular_non_expected_overlap.keys():
    print(contig + '\t' + str(header2len[contig] - self_circular_non_expected_overlap[contig]) + '\t' + str(cov[contig]) + '\t' + gc[contig] + '\t' + 'Self_circular' + '\t' + 'category_i', file=assembled_info, flush=True)
    print('>' + contig + '_self_circular', file=circular_fasta, flush=True)
    print(header2seq[contig][self_circular_non_expected_overlap[contig]:], file=circular_fasta, flush=True)
circular_fasta.close()
summary_fasta('COBRA_category_i_self_circular_queries_trimmed.fasta')
assembled_info.close()


# intermediate files
os.mkdir('intermediate.files')
os.system('mv COBRA_end_joining_pairs.txt COBRA_potential_joining_paths.txt COBRA_retrieved_for_joining intermediate.files')
os.mkdir('intermediate.files/invalid.checking')
os.system('mv for.blastn.fasta* intermediate.files/invalid.checking')


# write the numbers to the log file
print('\n', end='', file=log, flush=True)
print('=' * 150, file=log, flush=True)
print('Final summary', file=log, flush=True)
print('Total queries: ' + str(len(query_list)) + '\n' + 'Category i - Self_circular: ' +
      str(calculate_seq_num('COBRA_category_i_self_circular_queries_trimmed.fasta'))+ '\n' + 'Category ii - Extended_circular: ' + str(len(extended_circular_query))
      + ' (Unique: ' + str(calculate_seq_num('COBRA_category_ii_extended_circular_unique.fasta')) + ')' + '\n' + 'Category ii - Extended_ok: ' + str(len(extended_ok_query))
      + ' (Unique: ' + str(calculate_seq_num('COBRA_category_ii_extended_ok_unique.fasta')) + ')' + '\n' +
      'Category ii - Failed due to COBRA rules: ' + str(len(set(failed_join_list))) + '\n' +
      'Category iii - Failed due to DNA break: ' + str(len(DNA_break_query)), file=log, flush=True)
print('=' * 150, file=log, flush=True)
log.close()

