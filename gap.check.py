#!/usr/bin/env python
# Author: Lin-Xing Chen, UC Berkeley

import os
from Bio import SeqIO
import argparse
import pysam
import itertools

parser = argparse.ArgumentParser(description="This script is used to identify the mapping gaps of assembled sequences.")
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("-f", "--fasta", type=str, help="the sequences to check in fasta format.", required=True)
parser.add_argument("-b", "--bam", type=str, help="the reads mapping file in bam format.")
parser.add_argument("-m", "--mismatch", type=int, default=2, help="the max read mapping mismatches (default=2).")
parser.add_argument("-1", "--read1", type=str, help="read 1 fastq file used for read mapping.")
parser.add_argument("-2", "--read2", type=str, help="read 2 fastq file used for read mapping.")
args = parser.parse_args()


#
def get_gaps(i):
    gaps = []
    for a, b in itertools.groupby(enumerate(i), lambda pair: pair[1] - pair[0]):
        b = list(b)
        gaps.append(str(b[0][1]) + '-' + str(b[-1][1]))
    return gaps


#
fasta_file = os.path.abspath(args.fasta)
fasta_name = fasta_file.strip().split('/')[-1]

# check for mapping file or read files for mapping
if args.bam:
    mapping_file = os.path.abspath(args.bam)
    mapping_name = mapping_file.strip().split('/')[-1]
    bamfile = pysam.AlignmentFile('{0}'.format(args.bam), 'rb')
    bamfile_mismatch = pysam.AlignmentFile("{0}.filter.m{1}".format(mapping_name, str(args.mismatch)), "wb", template=bamfile)
    for line in bamfile:
        if not line.is_unmapped:
            if line.get_tag("NM") <= args.mismatch:
                bamfile_mismatch.write(line)
            else:
                pass
        else:
            pass
    bamfile.close()
    bamfile_mismatch.close()
else:
    if args.read1 and args.read2:
        os.system('bowtie2-build {0} {0}'.format(fasta_name))
        os.system('bowtie2 -p 16 -X 2000 -x {0} -1 {1} -2 {2} | shrinksam -v | sambam >{0}.mapped.bam'.format(fasta_name, args.read1, args.read2))
        bamfile = pysam.AlignmentFile('{0}.mapped.bam'.format(fasta_name), 'rb')
        bamfile_mismatch = pysam.AlignmentFile("{0}.mapped.bam.filter.m{1}".format(fasta_name, str(args.mismatch)), "wb", template=bamfile)
        for line in bamfile:
            if not line.is_unmapped:
                if line.get_tag("NM") <= args.mismatch:
                    bamfile_mismatch.write(line)
                else:
                    pass
            else:
                pass
        bamfile.close()
        bamfile_mismatch.close()
    else:
        print('Please provide mapping file or read files for mapping!')
        exit()


# index the bam file using samtools
if args.bam:
    mapping_file = os.path.abspath(args.bam)
    mapping_name = mapping_file.strip().split('/')[-1]
    os.system('samtools index {0}.filter.m{1}'.format(mapping_name, str(args.mismatch)))
    os.system('samtools mpileup --reference {0} {1}.filter.m{2} >{1}.filter.m{2}.pileup.txt'.format(fasta_name, mapping_name, str(args.mismatch)))

    header2seq = {}
    header2position = {}
    fasta = open('{0}'.format(fasta_name), 'r')
    for record in SeqIO.parse(fasta, "fasta"):
        header = str(record.id).strip()
        seq = str(record.seq)
        header2seq[header] = seq
        header2position[header] = []
        for position in range(71, len(seq) - 69):
            header2position[header].append(position)
    fasta.close()

    pileup = open('{0}.filter.m{1}.pileup.txt'.format(mapping_name, str(args.mismatch)), 'r')
    for line in pileup.readlines():
        line = line.strip().split('\t')
        if line[2] in ['A', 'T', 'C', 'G'] and int(line[1]) in header2position[line[0]]:
            header2position[line[0]].remove(int(line[1]))
        else:
            pass
    pileup.close()

    # write results
    scaf2gaps = {}
    with open('{0}.filter.m{1}.gap.info.txt'.format(mapping_name, str(args.mismatch)), 'w') as info:
        info.write('Sequence' + '\t' + 'Length' + '\t' + 'Num_unmapped_bases' + '\t' + 'Ummpaed_bases' + '\t' + 'Num_of_gaps' + '\n')
        for header in header2position.keys():
            if len(header2position[header]) != 0:
                info.write(header + '\t' + str(len(header2seq[header])) + '\t' + str(len(header2position[header])) + '\t' + ','.join(get_gaps(header2position[header][:])) + '\t' +
                           str(len(get_gaps(header2position[header]))) + '\n')
                scaf2gaps[header] = ','.join(get_gaps(header2position[header][:]))
            else:
                info.write(header + '\t' + str(len(header2seq[header])) + '\t' + str(len(header2position[header])) + '\t' + 'NA' + '\t' + '0' + '\n')
    info.close()

else:
    os.system('samtools index {0}.mapped.bam.filter.m{1}'.format(fasta_name, str(args.mismatch)))
    os.system('samtools mpileup --reference {0} {0}.mapped.bam.filter.m{1} >{0}.mapped.bam.filter.m{1}.pileup.txt'.format(fasta_name, str(args.mismatch)))

    header2seq = {}
    header2position = {}
    fasta = open('{0}'.format(fasta_name), 'r')
    for record in SeqIO.parse(fasta, "fasta"):
        header = str(record.id).strip()
        seq = str(record.seq)
        header2seq[header] = seq
        header2position[header] = []
        for position in range(71, len(seq) - 69):
            header2position[header].append(position)
    fasta.close()

    pileup = open('{0}.mapped.bam.filter.m{1}.pileup.txt'.format(fasta_name, str(args.mismatch)), 'r')
    for line in pileup.readlines():
        line = line.strip().split('\t')
        if line[2] in ['A', 'T', 'C', 'G'] and int(line[1]) in header2position[line[0]]:
            header2position[line[0]].remove(int(line[1]))
        else:
            pass
    pileup.close()

    # write results
    scaf2gaps = {}
    with open('{0}.mapped.bam.filter.m{1}.gap.info.txt'.format(fasta_name, str(args.mismatch)), 'w') as info:
        info.write('Sequence' + '\t' + 'Length' + '\t' + 'Num_unmapped_bases' + '\t' + 'Unmapped_bases' + '\t' + 'Num_of_gaps' + '\n')
        for header in header2position.keys():
            if len(header2position[header]) != 0:
                info.write(header + '\t' + str(len(header2seq[header])) + '\t' + str(len(header2position[header])) + '\t' + ','.join(get_gaps(header2position[header][:])) + '\t' +
                           str(len(get_gaps(header2position[header]))) + '\n')
                scaf2gaps[header] = ','.join(get_gaps(header2position[header][:]))
            else:
                info.write(header + '\t' + str(len(header2seq[header])) + '\t' + str(len(header2position[header])) + '\t' + 'NA' + '\t' + '0' + '\n')
    info.close()


def get_ok_region(sequence):
    length = len(header2seq[sequence])
    region_start = [1]
    region_end = []
    ok_region = []

    for item in scaf2gaps[sequence].split(','):
        if item.split('-')[0] == '71':
            pass
        elif item.split('-')[1] == str(length - 70):
            pass
        else:
            region_end.append(int(item.split('-')[0]) - 1)
            region_start.append(int(item.split('-')[1]) + 1)

    region_end.append(length)

    if len(region_end) == len(region_end) == 1:
        ok_region.append(header2seq[sequence])
    else:
        for location in range(0,len(region_start)):
            ok_region.append(header2seq[sequence][region_start[location]-1:region_end[location]])

    return ok_region


# write the sequences with gaps by N to a file
out_file = open('{0}.mapped.bam.filter.m{1}.gaps.replaced.by.Ns.fasta'.format(fasta_name,str(args.mismatch)), 'w')
for sequence in header2seq.keys():
    if sequence not in scaf2gaps.keys():
        out_file.write('>' + sequence + '\n')
        out_file.write(header2seq[sequence] + '\n')
    else:
        if get_ok_region(sequence)[0] == header2seq[sequence]:
            out_file.write('>' + sequence + '\n')
            out_file.write(header2seq[sequence] + '\n')
        else:
            out_file.write('>' + sequence + '_gaps_replaced_by_Ns' + '\n')
            out_file.write('NNNNNNNNNN'.join(get_ok_region(sequence)[:]) + '\n')
out_file.close()
