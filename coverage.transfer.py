#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description="This script transfer the metabat depth file to a two-column file.")
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("-i", type=str, help="The input coverage file", required=True)
requiredNamed.add_argument("-o", type=str, help="The output coverage file", required=True)
args = parser.parse_args()

a = open('{0}'.format(args.i), 'r')
b = open('{0}'.format(args.o), 'w')

for line in a.readlines():
    if not line.startswith('contigName'):
        line = line.strip().split('\t')
        print(line[0] + '\t' + line[2], file=b)
    else:
        pass

b.close()
a.close()