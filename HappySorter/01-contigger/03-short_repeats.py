#!/usr/bin/env python3 -u
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
parser.add_argument("--pe_rep_min_support", help="paired read min support to expand canonical repeats", type=int, default=5)
parser.add_argument("--pe_rep_max_noise", help="paired read max noise to expand canonical repeats", type=int, default=5)
parser.add_argument("--pe_rep_snr", help="paired read snr to expand canonical repeats", type=int, default=10)
args = parser.parse_args()

ws=SDG.WorkSpace(f'{args.output_prefix}_02_clean.sdgws')
peds=ws.get_paired_reads_datastore(ws.list_paired_reads_datastores()[0])

kc=ws.get_kmer_counter("main")
kc.set_kci_peak(args.unique_coverage)
kc.update_graph_counts()
print(ws.sdg.stats_by_kci())

print_step_banner("CANONICAL REPEAT EXPANSION")

c=SDG.GraphContigger(ws)
c.solve_canonical_repeats(peds,min_support=args.pe_rep_min_support,max_noise=args.pe_rep_max_noise,snr=args.pe_rep_snr)
kc.update_graph_counts()
print(ws.sdg.stats_by_kci())
peds.mapper.path_reads()

##second round
c=SDG.GraphContigger(ws)
c.solve_canonical_repeats(peds,min_support=args.pe_rep_min_support,max_noise=args.pe_rep_max_noise,snr=args.pe_rep_snr)
kc.update_graph_counts()
print(ws.sdg.stats_by_kci())
peds.mapper.path_reads()


ws.dump(f'{args.output_prefix}_03_shortreps.sdgws')