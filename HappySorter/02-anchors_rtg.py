#!/usr/bin/env python3
import SDGpython as SDG
import argparse
from collections import Counter
import os

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output_prefix", help="prefix for output files", type=str, required=True)
parser.add_argument("-m", "--min_coverage", help="min anchor coverage", type=int, required=True)
parser.add_argument("-M", "--max_coverage", help="max anchor coverage", type=int, required=True)
parser.add_argument("-s", "--max_anchor_size", help="max anchor size (in 31-mers)", type=int, required=True)

args = parser.parse_args()

ws=SDG.WorkSpace(f'{args.output_prefix}_strided.sdgws')
kc=ws.get_kmer_counter('main')

ws2=SDG.WorkSpace()
SDG.GraphContigger(ws).contig_reduction_to_unique_kmers(ws2,"main","pe",args.min_coverage,args.max_coverage,args.max_anchor_size);

peds = ws2.add_paired_reads_datastore(ws.list_paired_reads_datastores()[0]+".prseq")
peds.mapper.path_reads(k=31)
lords = ws2.add_long_reads_datastore(ws.list_long_reads_datastores()[0]+".loseq")
lrr=SDG.LongReadsRecruiter(ws2.sdg,lords,k=31)
lrr.map()

ws2.dump(f'{args.output_prefix}_anchors.sdgws')
peds.mapper.dump_readpaths(f'{args.output_prefix}_anchors.pe_paths')
lrr.dump(f'{args.output_prefix}_anchors.lrr')
