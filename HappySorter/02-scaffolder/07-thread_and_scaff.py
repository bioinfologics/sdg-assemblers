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
parser.add_argument("--min_hits", help="min hits to thread a node", type=int, default=1)
parser.add_argument("--min_links", help="threads linking two selected nodes", type=int, default=5)
#parser.add_argument("--min_link_perc", help="link percentaje to join two selected nodes", type=float, default=.1)
#parser.add_argument("--max_overlap", help="max overlap to join two selected nodes", type=int, default=200)
parser.add_argument("--max_thread_count", help="max threads to select a node (to avoid repeats)", type=int, default=1000)
parser.add_argument("-l", "--long_datastore", help="long reads datastore (default: '' uses the first datastore in the input workspace)", type=str, default='')
args = parser.parse_args()

ws=SDG.WorkSpace(f'{args.output_prefix}_06_split.sdgws')
ws.dump(f'{args.output_prefix}_07_graphonly.sdgws',graph_only=True)
peds=ws.get_paired_reads_datastore(ws.list_paired_reads_datastores()[0])
if args.long_datastore=='':
    lords=ws.get_long_reads_datastore(ws.list_long_reads_datastores()[0])
else:
    lords = ws.get_long_reads_datastore(args.long_datastore)

lrr=SDG.LongReadsRecruiter(ws.sdg,lords,k=31)
lrr.load(f'{args.output_prefix}_06_split.lrr')

kc=ws.get_kmer_counter("main")
kc.set_kci_peak(args.unique_coverage)
kc.update_graph_counts()
print(ws.sdg.stats_by_kci())

print_step_banner("MERGING RTGs")
lrr.clean_sandwich_matches()
lrr.simple_thread_reads(4)
rtg=lrr.rtg_from_threads()
print(f"rtg with 4+ hits has {len(rtg.get_all_nodeviews(include_disconnected=False))} connected nodes")
for x in range(3,args.min_hits-1,-1):
    lrr.simple_thread_reads(x)
    rtg=rtg.merge(lrr.rtg_from_threads())
    print(f"after merging with {x}+ hits, rtg has {len(rtg.get_all_nodeviews(include_disconnected=False))} connected nodes")

rtg.dump(f'{args.output_prefix}_07_merged.rtg')

print_step_banner("REDUCED THREAD GRAPH")

rtg=lrr.rtg_from_threads()
whitelisted_nvs=[ nv for nv in rtg.get_all_nodeviews() if nv.size()>=args.min_size and nv.kci()<=args.max_kci and nv.kci()>=args.min_kci and len(rtg.node_threads(nv.node_id()))<=args.max_thread_count and len(rtg.node_threads(nv.node_id()))>=args.min_links]
rrtg=rtg.reduced_graph(set([x.node_id() for x in whitelisted_nvs]))
rrtg.dump(f'{args.output_prefix}_07_reduced.rtg')

print(f"reduced rtg has {len(rrtg.get_all_nodeviews(include_disconnected=False))} connected nodes")

print_step_banner("CLOSEST RELIABLE CONNECTIONS GRAPH")
def remove_all_transitive_links(rtg):
    to_remove=Counter()
    for nv in rtg.get_all_nodeviews(both_directions=True,include_disconnected=False):
        lns=nv.next()
        #if nv.node_id()==-2287973: print(f'starting links {lns}')
        if len(lns)<2: continue
        reached=set()
        for ln in lns:
            #if nv.node_id()==-2287973: print(f'evaluating {ln}, reached so far: {reached}')
            if ln.node().node_id() in reached:
                link=(min(-nv.node_id(),ln.node().node_id()),max(-nv.node_id(),ln.node().node_id()))
                to_remove[link]+=1
            reached.update([x.node().node_id() for x in ln.node().next()])
    for ln,c in to_remove.items():
        if c==2:
            rtg.remove_link(ln[0],ln[1])

rrcg=rrtg.closest_reliable_connections_graph(3,args.min_links)
remove_all_transitive_links(rrcg)
print(f"closest reliable connections graph has {len(rrcg.get_all_nodeviews(include_disconnected=False))} connected nodes")
rrcg.write_to_gfa1(f'{args.output_prefix}_07_crc.gfa',selected_nodes=[x.node_id() for x in rrcg.get_all_nodeviews(include_disconnected=False)])

# print('Popping all tips')
# to_remove=[]
# for nv in crtg.get_all_nodeviews(include_disconnected=False):
#     if nv.is_tip():
#         to_remove.append(nv.node_id())
# print(f'{len(to_remove)} tips found')
# for nid in to_remove:
    # rrtg.pop_node_from_all(nid)

def is_jumped_repeat(crtg,rtg,nid):
    cnv=crtg.get_nodeview(nid)
    nv=rtg.get_nodeview(nid)
    if len(cnv.prev())<2 or len(cnv.next())<2: return False
    through=0
    total=0
    for tid in rtg.node_threads(nid,True):
        total+=1
        tns=rtg.get_thread_nodes(tid)
        if tns[0]!=nid and tns[-1]!=nid: through+=1
    return through>=total*.8


def solvable_collapse(crtg,rtg,first_node,last_node,min_links=10,ratio=6):
    cgnvF=crtg.get_nodeview(first_node)
    cgnvL=crtg.get_nodeview(last_node)
    
    if first_node==last_node:
        if is_jumped_repeat(crtg,rtg,first_node): return True
#         if len(cgnvF.prev())>3 or len(cgnvL.next())>3: return True
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


collapsed_lines_to_remove=[l for l in rrcg.get_all_lines(1) if len(l)<10 and solvable_collapse(rrcg,rrtg,l[0],l[-1],20)]

while collapsed_lines_to_remove:
    print(f"removing {sum([len(x) for x in collapsed_lines_to_remove])} nodes from {len(collapsed_lines_to_remove)} collapsed lines",flush=True)
    for l in collapsed_lines_to_remove:
        for nid in map(abs,l):
            rrtg.pop_node_from_all(nid)
    rrcg=rrtg.closest_reliable_connections_graph(3,args.min_links)
    remove_all_transitive_links(rrcg)
    collapsed_lines_to_remove=[l for l in rrcg.get_all_lines(1) if len(l)<10 and solvable_collapse(rrcg,rrtg,l[0],l[-1],20)]
print(f"simplified closest reliable connections graph has {len(rrcg.get_all_nodeviews(include_disconnected=False))} connected nodes")    
rrcg.dump(f'{args.output_prefix}_07_crc_simp.dg')
rrcg.write_to_gfa1(f'{args.output_prefix}_07_crc_simp.gfa',selected_nodes=[x.node_id() for x in rrcg.get_all_nodeviews(include_disconnected=False)])

with open(f'{args.output_prefix}_07_crc_simp.csv','w') as of:
        of.write('Node,KCI,Colour\n')
        for x in rrtg.get_all_nodeviews(include_disconnected=False):
            nid=x.node_id()
            kci=x.kci()
            if kci<.5: c='gray'
            elif kci<1.5: c='green'
            elif kci<2.5: c='blue'
            else: c='red'
            of.write(f'seq{nid},{kci :.2f},{c}\n')
