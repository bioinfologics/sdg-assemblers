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
parser.add_argument("-w", "--workspace", help="input workspace", type=str, required=True)
parser.add_argument("-u", "--unique_coverage", help="value for unique coverage on KCI", type=int, required=True)
args = parser.parse_args()

ws=SDG.WorkSpace(f'{args.workspace}')
kc=ws.get_kmer_counter("main")
kc.set_kci_peak(args.unique_coverage)
kc.update_graph_counts()
kc.compute_all_kcis()
print(ws.sdg.stats_by_kci())