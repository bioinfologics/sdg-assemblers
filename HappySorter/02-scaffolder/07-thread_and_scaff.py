#!/usr/bin/env python3 -u
import SDGpython as SDG
import argparse
from collections import Counter
import os
from math import ceil
from statistics import median

def print_step_banner(s):
    print('\n'+'*'*(len(s)+4))
    print(f'* {s} *')
    print('*'*(len(s)+4)+"\n")

parser = argparse.ArgumentParser()

parser.add_argument("-o", "--output_prefix", help="prefix for output files", type=str, required=True)
parser.add_argument("-u", "--unique_coverage", help="value for unique coverage at 31-mers", type=int, required=True)
parser.add_argument("-s", "--min_size", help="min size to keep nodes", type=int, default=500)
parser.add_argument("--min_kci", help="min kci to keep nodes", type=float, default=.5)
parser.add_argument("--max_kci", help="max kci to keep nodes", type=float, default=1.5)
parser.add_argument("--min_hits", help="min hits to thread a node", type=int, default=3)
parser.add_argument("--min_links", help="threads linking two selected nodes", type=int, default=5)
args = parser.parse_args()

ws=SDG.WorkSpace(f'{args.output_prefix}_06_split.sdgws')
peds=ws.get_paired_reads_datastore(ws.list_paired_reads_datastores()[0])
lords=ws.get_long_reads_datastore(ws.list_long_reads_datastores()[0])

lrr=SDG.LongReadsRecruiter(ws.sdg,lords,k=31)
lrr.load(f'{args.output_prefix}_06_split.lrr')

kc=ws.get_kmer_counter("main")
kc.set_kci_peak(args.unique_coverage)
kc.update_graph_counts()
print(ws.sdg.stats_by_kci())

print_step_banner("REDUCED THREAD GRAPH")

lrr.simple_thread_reads(args.min_hits)
rtg=lrr.rtg_from_threads()
whitelisted_nvs=[ nv for nv in rtg.get_all_nodeviews() if nv.size()>=args.min_size and nv.kci()<=args.max_kci and nv.kci()>=args.min_kci]
rrtg=rtg.reduced_graph(set([x.node_id() for x in whitelisted_nvs]))

rrtg.dump(f'{args.output_prefix}_07_reduced.rtg')

def condense_rrtg(g, min_links=10, min_distance=-2000000):
    cg=SDG.DistanceGraph(g.sdg)
    for nv in rrtg.get_all_nodeviews(include_disconnected=False,both_directions=True):
        cc=Counter([lf.node().node_id() for lf in nv.next()])
        for dnid,c in cc.items():
            if abs(dnid)<abs(nv.node_id()) or c<min_links or c<len(nv.next())*.3: continue
            d=int(median([lf.distance() for lf in nv.next() if lf.node().node_id()==dnid]))
            if d>=min_distance: cg.add_link(-nv.node_id(),dnid,d)
    return cg

crtg=condense_rrtg(rrtg,min_links=args.min_links,min_distance=-200)
crtg.remove_transitive_links(5)

crtg.write_to_gfa1(f'{args.output_prefix}_07_crtg.gfa',selected_nodes=[x.node_id() for x in crtg.get_all_nodeviews(include_disconnected=False)])