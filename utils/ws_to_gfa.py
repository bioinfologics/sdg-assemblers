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
parser.add_argument("-w", "--workspace", help="input workspace", type=str, required=True)
parser.add_argument("-u", "--unique_coverage", help="value for unique coverage on KCI", type=int, default=-1)
args = parser.parse_args()

ws=SDG.WorkSpace(f'{args.workspace}')
ws.sdg.write_to_gfa1(f'{args.output_prefix}.gfa')

if args.unique_coverage!=-1:
    kc=ws.get_kmer_counter("main")
    kc.set_kci_peak(args.unique_coverage)
    kc.update_graph_counts()
    kc.compute_all_kcis()
    print(ws.sdg.stats_by_kci())

    with open(f'{args.output_prefix}.csv','w') as of:
            of.write('Node,KCI,Colour\n')
            for x in ws.sdg.get_all_nodeviews():
                nid=x.node_id()
                kci=x.kci()
                if kci<.5: c='gray'
                elif kci<1.5: c='green'
                elif kci<2.5: c='blue'
                else: c='red'
                of.write(f'seq{nid},{kci :.2f},{c}\n')