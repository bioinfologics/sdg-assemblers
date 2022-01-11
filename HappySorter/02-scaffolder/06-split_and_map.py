#!/usr/bin/env python3 -u
import SDGpython as SDG
import argparse
from collections import Counter
import os
from math import ceil

def print_step_banner(s):
    print('\n'+'*'*(len(s)+4))
    print(f'* {s} *')
    print('*'*(len(s)+4)+"\n")

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_workspace", help="input workspace", type=str, required=True)
parser.add_argument("-o", "--output_prefix", help="prefix for output files", type=str, required=True)
parser.add_argument("-u", "--unique_coverage", help="value for unique coverage at 31-mers", type=int, required=True)
parser.add_argument("-m", "--min_anchormap_coverage", help="min read kcov value for anchormap (default: 0)", type=int, default=0)
parser.add_argument("-M", "--max_anchormap_coverage", help="max read kcov value for anchormap (default: 0, uses map without anchor filtering)", type=int, default=0)
parser.add_argument("-s", "--split_size", help="break contigs into chunks of roughly...", type=int, default=2000)
parser.add_argument("-l", "--long_datastore", help="long reads datastore (default: '' uses the first datastore in the input workspace)", type=str, default='')
args = parser.parse_args()

ws=SDG.WorkSpace(args.input_workspace)
peds=ws.get_paired_reads_datastore(ws.list_paired_reads_datastores()[0])
if args.long_datastore=='':
    lords=ws.get_long_reads_datastore(ws.list_long_reads_datastores()[0])
else:
    lords = ws.add_long_reads_datastore(args.long_datastore)

kc=ws.get_kmer_counter("main")
kc.set_kci_peak(args.unique_coverage)
kc.update_graph_counts()
print(ws.sdg.stats_by_kci())

print_step_banner("BREAK LONG CONTIGS")

def split_long_contigs(split_size):
    to_split=[x.node_id() for x in ws.sdg.get_all_nodeviews() if x.size()>split_size]
    print(f'{len(to_split)} nodes to split',flush=True)
    for nid in to_split:
        nv=ws.sdg.get_nodeview(nid)
        seq=nv.sequence()
        nns=int(ceil(nv.size()/ceil(nv.size()/split_size)))
        splits=list(range(0,len(seq),nns))
        last_nid=0
        for ci,cstart in enumerate(splits):
            if ci==0:
                new_nid=ws.sdg.add_node(seq[:splits[ci+1]+15],False)
                for lb in nv.prev():
                    ws.sdg.add_link(-lb.node().node_id(),new_nid,lb.distance())
            else:
                if ci==len(splits)-1: 
                    new_nid=ws.sdg.add_node(seq[cstart-15:],False)
                    for lf in nv.next():
                        ws.sdg.add_link(-new_nid,lf.node().node_id(),lf.distance())
                else:
                    new_nid=ws.sdg.add_node(seq[cstart-15:splits[ci+1]+15],False)
                ws.sdg.add_link(-last_nid,new_nid,-30)
            last_nid=new_nid
        ws.sdg.remove_node(nid)
                

split_long_contigs(args.split_size)
kc.update_graph_counts()
print(ws.sdg.stats_by_kci())


lrr=SDG.LongReadsRecruiter(ws.sdg,lords,k=31)

if args.max_anchormap_coverage:
    if kc.k==31:
        kc31=kc
        kc_name="main"
    else:
        kc31=ws.add_kmer_counter("main31",31)
        kc31.add_count("pe",peds)
        kc_name="main31"
    lrr.anchormap(kcname=kc_name, countname='pe', fmin=args.min_anchormap_coverage, fmax=args.max_anchormap_coverage, graph_fmin=1, graph_fmax=1)
else: lrr.map()

lrr.dump(f'{args.output_prefix}_06_split.lrr')
ws.dump(f'{args.output_prefix}_06_split.sdgws')
