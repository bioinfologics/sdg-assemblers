#!/usr/bin/env -S python3 -u
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
parser.add_argument("-u", "--unique_coverage", help="value for unique coverage on KCI", type=int, required=True)
parser.add_argument("--tip_size", help="size of nodes for tip clipping", type=int, default=200)
parser.add_argument("--low_coverage", help="low coverage for node cleanup", type=int, default=5)
parser.add_argument("--low_coverage_size", help="size for low coverage node cleanup", type=int, default=200)
args = parser.parse_args()

ws=SDG.WorkSpace(f'{args.output_prefix}_01_dbg.sdgws')
peds=ws.get_paired_reads_datastore(ws.list_paired_reads_datastores()[0])

kc=ws.get_kmer_counter("main")
kc.set_kci_peak(args.unique_coverage)
kc.update_graph_counts()
print(ws.sdg.stats_by_kci())

print_step_banner("GRAPH CLEANING")

gc=SDG.GraphContigger(ws)
print("Tip clipping:")
gc.clip_tips(args.tip_size)
kc.update_graph_counts()
print(ws.sdg.stats_by_kci())

print("Low kcov node removal:")
gc.remove_low_kcov_nodes("main","pe",args.low_coverage,args.low_coverage_size)
kc.update_graph_counts()
print(ws.sdg.stats_by_kci())

print("Tip clipping:")
gc.clip_tips(args.tip_size)
kc.update_graph_counts()
print(ws.sdg.stats_by_kci())

peds.mapper.path_reads()

ws.dump(f'{args.output_prefix}_02_clean.sdgws')