#!/usr/bin/env python3
import SDGpython as SDG
import argparse
from collections import Counter
import os

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output_prefix", help="prefix for output files", type=str, required=True)
parser.add_argument("-u", "--unique_coverage", help="value for unique coverage at 31-mers", type=int, required=True)
parser.add_argument("-m", "--min_coverage", help="min anchor coverage", type=int, required=True)
parser.add_argument("-M", "--max_coverage", help="max anchor coverage", type=int, required=True)
parser.add_argument("--anchor_k", help="k value for anchor reduction (-1: use kci_k from previosu step)", type=int, default=-1)
parser.add_argument( "--min_kci", help="min contig kmer coverage index", type=int, default=0)
parser.add_argument( "--max_kci", help="max_contig kmer coverage index", type=float, default=1.25)
parser.add_argument("-s", "--max_anchor_size", help="max anchor size (in 31-mers)", type=int, required=True)

args = parser.parse_args()

ws=SDG.WorkSpace(f'{args.output_prefix}_strided.sdgws')
ws2=SDG.WorkSpace()
mainkc=ws.get_kmer_counter('main')
oldpeds=ws.get_paired_reads_datastore(ws.list_paired_reads_datastores()[0])

if args.anchor_k==-1 or args.anchor_k==mainkc.k:
    kcname='main'
    print("using existing counter at k=%d"%mainkc.k)
    kc=ws.get_kmer_counter(kcname)
else:
    kcname='anchorkc'
    print("producing counter for new k value at k=%d"%args.anchor_k)
    kc=ws.add_kmer_counter(kcname,args.anchor_k)
    kc.add_count("pe",oldpeds)
    


kc.set_kci_peak(args.unique_coverage)
kc.compute_all_kcis()
print("producing anchors at k=%d"%kc.k)
SDG.GraphContigger(ws).contig_reduction_to_unique_kmers(ws2,kcname,"pe",args.min_coverage,args.max_coverage,args.max_anchor_size,args.min_kci,args.max_kci);

peds = ws2.add_paired_reads_datastore(ws.list_paired_reads_datastores()[0]+".prseq")
peds.mapper.path_reads(k=kc.k)
lords = ws2.add_long_reads_datastore(ws.list_long_reads_datastores()[0]+".loseq")
lrr=SDG.LongReadsRecruiter(ws2.sdg,lords,k=kc.k)
lrr.map()

ws2.dump(f'{args.output_prefix}_anchors.sdgws')
peds.mapper.dump_readpaths(f'{args.output_prefix}_anchors.pe_paths')
lrr.dump(f'{args.output_prefix}_anchors.lrr')
