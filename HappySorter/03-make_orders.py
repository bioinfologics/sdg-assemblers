#!/usr/bin/env python3
import SDGpython as SDG
import argparse
from collections import Counter
import os

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output_prefix", help="prefix for output files", type=str, required=True)
#parser.add_argument("-m", "--min_coverage", help="min anchor coverage", type=int, required=True)
#parser.add_argument("-M", "--max_coverage", help="max anchor coverage", type=int, required=True)
#parser.add_argument("-s", "--max_anchor_size", help="max anchor size (in 31-mers)", type=int, required=True)

args = parser.parse_args()

def write_orders_to_fasta(ts,filename):
    with open(filename,'w') as f:
        for sid,s in ts.sorters.items():
            f.write(f'>sorter{sid:06d}\n{"N".join([ws.sdg.get_nodeview(x).sequence() for x in s.order.as_signed_nodes()])}\n')


ws=SDG.WorkSpace(f'{args.output_prefix}_anchors.sdgws')
#kc=ws.get_kmer_counter('main')

peds = ws.get_paired_reads_datastore(ws.list_paired_reads_datastores()[0])
peds.mapper.load_readpaths(f'{args.output_prefix}_anchors.pe_paths')
lords = ws.get_long_reads_datastore(ws.list_long_reads_datastores()[0])
lrr=SDG.LongReadsRecruiter(ws.sdg,lords,k=31)
lrr.load(f'{args.output_prefix}_anchors.lrr')
lrr.simple_thread_reads()
rtg=lrr.rtg_from_threads(True)



