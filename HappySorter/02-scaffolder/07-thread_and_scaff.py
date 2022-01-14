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
parser.add_argument("-u", "--unique_coverage", help="value for unique coverage on kci", type=int, required=True)
parser.add_argument("-s", "--min_size", help="min size to keep nodes", type=int, default=500)
parser.add_argument("--min_kci", help="min kci to keep nodes", type=float, default=.5)
parser.add_argument("--max_kci", help="max kci to keep nodes", type=float, default=1.5)
parser.add_argument("--min_hits", help="min hits to thread a node", type=int, default=3)
parser.add_argument("--min_links", help="threads linking two selected nodes", type=int, default=5)
parser.add_argument("--min_link_perc", help="link percentaje to join two selected nodes", type=float, default=.1)
parser.add_argument("--max_overlap", help="max overlap to join two selected nodes", type=int, default=200)
parser.add_argument("--max_thread_count", help="max threads to select a node (to avoid repeats)", type=int, default=200)
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

lrr.clean_sandwich_matches()
print_step_banner("REDUCED THREAD GRAPH")

lrr.simple_thread_reads(args.min_hits)
rtg=lrr.rtg_from_threads()
whitelisted_nvs=[ nv for nv in rtg.get_all_nodeviews() if nv.size()>=args.min_size and nv.kci()<=args.max_kci and nv.kci()>=args.min_kci and len(rtg.node_threads(nv.node_id()))<=args.max_thread_count and len(rtg.node_threads(nv.node_id()))>=args.min_links]
rrtg=rtg.reduced_graph(set([x.node_id() for x in whitelisted_nvs]))

rrtg.dump(f'{args.output_prefix}_07_reduced.rtg')

def condense_rrtg(g, min_links=10, min_distance=-2000000):
    cg=SDG.DistanceGraph(g.sdg)
    for nv in rrtg.get_all_nodeviews(include_disconnected=False,both_directions=True):
        cc=Counter([lf.node().node_id() for lf in nv.next()])
        for dnid,c in cc.items():
            if abs(dnid)<abs(nv.node_id()) or c<min_links or c<len(nv.next())*args.min_link_perc: continue
            d=int(median([lf.distance() for lf in nv.next() if lf.node().node_id()==dnid]))
            if d>=min_distance: cg.add_link(-nv.node_id(),dnid,d)
    return cg

crtg=condense_rrtg(rrtg,min_links=args.min_links,min_distance=-args.max_overlap)
crtg.remove_transitive_links(5)

print('Popping all tips')
to_remove=[]
for nv in crtg.get_all_nodeviews(include_disconnected=False):
    if nv.is_tip():
        to_remove.append(nv.node_id())
print(f'{len(to_remove)} tips found')
for nid in to_remove:
    rrtg.pop_node_from_all(nid)

crtg=condense_rrtg(rrtg,min_links=args.min_links,min_distance=-args.max_overlap)
crtg.remove_transitive_links(5)

def solvable_collapse(crtg,rtg,first_node,last_node,min_links=5,ratio=6):
    cgnvF=crtg.get_nodeview(first_node)
    cgnvL=crtg.get_nodeview(last_node)
    if len(cgnvF.prev())!=2 or len(cgnvL.next())!=2: return False
    pAtids=set(rtg.node_threads(cgnvF.prev()[0].node().node_id(),True))
    pBtids=set(rtg.node_threads(cgnvF.prev()[1].node().node_id(),True))
    nAtids=set(rtg.node_threads(cgnvL.next()[0].node().node_id(),True))
    nBtids=set(rtg.node_threads(cgnvL.next()[1].node().node_id(),True))
    AA=len(pAtids.intersection(nAtids))
    AB=len(pAtids.intersection(nBtids))
    BA=len(pBtids.intersection(nAtids))
    BB=len(pBtids.intersection(nBtids))
    if AA<AB:
        AA,AB,BA,BB=AB,AA,BB,BA
    if AA>=min_links and AB<min_links and AA>=ratio*AB and BB>=min_links and BA<min_links and BB>=ratio*BA: return True
    return False

collapsed_lines_to_remove=[l for l in crtg.get_all_lines(1) if len(l)<4 and solvable_collapse(crtg,rrtg,l[0],l[-1],args.min_links)]

while collapsed_lines_to_remove:
    print(f"removing {sum([len(x) for x in collapsed_lines_to_remove])} nodes from {len(collapsed_lines_to_remove)} collapsed lines",flush=True)
    for l in collapsed_lines_to_remove:
        for nid in map(abs,l):
            rrtg.pop_node_from_all(nid)
    crtg=condense_rrtg(rrtg,min_links=args.min_links,min_distance=-args.max_overlap)
    crtg.remove_transitive_links(5)
    collapsed_lines_to_remove=[l for l in crtg.get_all_lines(1) if len(l)<4 and solvable_collapse(crtg,rrtg,l[0],l[-1],args.min_links)]

crtg.write_to_gfa1(f'{args.output_prefix}_07_crtg.gfa',selected_nodes=[x.node_id() for x in crtg.get_all_nodeviews(include_disconnected=False)])
with open(f'{args.output_prefix}_07_crtg.csv','w') as of:
        of.write('Node,KCI,Colour\n')
        for x in rrtg.get_all_nodeviews(include_disconnected=False):
            nid=x.node_id()
            kci=x.kci()
            if kci<.5: c='gray'
            elif kci<1.5: c='green'
            elif kci<2.5: c='blue'
            else: c='red'
            of.write(f'seq{nid},{kci :.2f},{c}\n')
