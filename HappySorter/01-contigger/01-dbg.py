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
parser.add_argument("-k", "--k", help="k value for graph construction", type=int, default=63)
parser.add_argument("-c", "--min_coverage", help="min coverage for graph construction", type=int, default=3)
parser.add_argument("-b", "--disk_batches", help="disk batches for graph construction", type=int, default=1)
parser.add_argument("--kci_k", help="k value for KCI", type=int, default=31)
parser.add_argument("--load_unitigs", help="load pre-computed unitigs from fasta file", type=str, default="")
args = parser.parse_args()

print_step_banner("INITIAL GRAPH CONSTRUCTION")

ws = SDG.WorkSpace()
peds = ws.add_paired_reads_datastore(args.paired_datastore)
if args.load_unitigs:
    ws.sdg.load_from_bcalm(args.load_unitigs,args.k)
    if len(ws.sdg.get_all_nodeviews(include_disconnected=False))<.8*len(ws.sdg.get_all_nodeviews()):
        print(f"WARNING: {len(ws.sdg.get_all_nodeviews(include_disconnected=False))}/{len(ws.sdg.get_all_nodeviews())} loaded unitigs are disconnected, check k is set correctly")
else:
    SDG.GraphMaker(ws.sdg).new_graph_from_paired_datastore(peds,args.k,args.min_coverage,args.disk_batches)
    os.replace('small_K.freqs',f'{args.output_prefix}_small_K.freqs')
print(ws.sdg.simple_structure_stats())

kc=ws.add_kmer_counter("main",args.kci_k)
kc.add_count("pe",peds)
with open(f'{args.output_prefix}_pe_kc_spectra.csv','w') as of:
    for x,c in enumerate(kc.count_spectra("pe"),start=0):
        of.write(f'{x},{c}\n')

ws.dump(f'{args.output_prefix}_01_dbg.sdgws')