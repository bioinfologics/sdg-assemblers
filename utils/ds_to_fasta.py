#!/usr/bin/env python3 -u
from ast import increment_lineno
import SDGpython as SDG
import argparse
from collections import Counter
import os

def print_step_banner(s):
    print('\n'+'*'*(len(s)+4))
    print(f'* {s} *')
    print('*'*(len(s)+4)+"\n")

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output_prefix", help="prefix for output files", type=str, required=True)
parser.add_argument("-p", "--paired_datastore", help="paired reads datastore", type=str, required=True)

args = parser.parse_args()



ws = SDG.WorkSpace()
peds = ws.add_paired_reads_datastore(args.paired_datastore)
print_step_banner("DUMPING FASTA FILES")
with open(f'{args.output_prefix}_1.fasta','w') as f1,open(f'{args.output_prefix}_2.fasta','w') as f2:
    for rid in range(1,peds.size(),2):
        f1.write(f'>{rid}\n{peds.get_read_sequence(rid)}\n')
        f2.write(f'>{rid}\n{peds.get_read_sequence(rid+1)}\n')