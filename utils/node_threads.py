#!/usr/bin/env python3 -x
import SDGpython as SDG
import argparse
from collections import Counter
import os
from math import ceil
from statistics import median
from pylab import *

def print_step_banner(s):
    print('\n'+'*'*(len(s)+4))
    print(f'* {s} *')
    print('*'*(len(s)+4)+"\n")

parser = argparse.ArgumentParser()

parser.add_argument("-o", "--output_prefix", help="prefix for output files", type=str, required=True)
parser.add_argument("-w", "--workspace", help="input workspace", type=str, required=True)
parser.add_argument("-r","--input_rtg", help="input read threads graph", type=str, required=True)
parser.add_argument("-d","--input_dg", help="input distance graph (optional)", type=str, default="")
parser.add_argument("-u", "--unique_coverage", help="value for unique coverage on kci", type=int, required=True)
parser.add_argument("-m", "--min_threads", help="min_threads to include a node", type=int, default=1)
parser.add_argument("-n","--node", help="node to dump treads from", type=int, nargs='*', action='append', required=True)
args = parser.parse_args()


ws=SDG.WorkSpace(args.workspace)
kc=ws.get_kmer_counter("main")
kc.set_kci_peak(args.unique_coverage)
rtg=SDG.ReadThreadsGraph(ws.sdg)
rtg.load(args.input_rtg)
dg=SDG.DistanceGraph(ws.sdg)
if args.input_dg: dg.load(args.input_dg)

def node_order(node_lists):
    order=node_lists[0]
    for nl in node_lists:
        to_place=[]
        last_seen=-1
        for nid in nl:
            if nid in order:
                nindex=order.index(nid)
                if to_place:
                    if last_seen==nindex-1:
                        order=order[:last_seen+1]+to_place+order[nindex:]
                    to_place=[]
                last_seen=nindex
            else:
                to_place.append(nid)
        if to_place and last_seen==len(order)-1: order=order+to_place
    for nl in node_lists:
        to_place=[]
        last_seen=-1
        for nid in nl:
            if nid in order:
                nindex=order.index(nid)
                if to_place:
                    order=order[:last_seen+1]+to_place+order[last_seen+1:]
                    to_place=[]
                last_seen=nindex
            else:
                to_place.append(nid)
        if to_place: order=order+to_place
    return order

def node_label(nid):
    nv=ws.sdg.get_nodeview(nid)
    s=f'{nid} ( {nv.size()}bp, kci={nv.kci():.2f} ) '
    if dg.get_nodeview(nid).prev(): s+='<-'
    else: s+=' -'
    if dg.get_nodeview(nid).next(): s+='>'
    else: s+=' '
    return s
def thread_plot(rtg,tids):
    #each node gets an x position.
    y=0
    
    #place every node at an appropriate position the first time it is seen (may be suboptimal but not terribly)
    tnodes=[rtg.get_thread_nodes(tid) for tid in tids]
    tnc=Counter()
    for tn in tnodes:
        tnc.update(tn)
    tnodes=[[nid for nid in tn if tnc[nid]>=args.min_threads] for tn in tnodes]
    nord=node_order(tnodes)
    node_x={}
    for xp,nid in enumerate(nord):
        node_x[nid]=xp
    

    #reorder the tids by their first node
    tf=[]
    for tid in tids:
        for np in rtg.get_thread(tid):
            if np.node in node_x:
                tf.append([node_x[np.node],tid])
                break
    tids=[x[1] for x in sorted(tf)]
    yticks_y=[]
    yticks_l=[]
    for tid in tids:
        y-=1
        yticks_y.append(y)
        yticks_l.append("%d" % (tid))
        t_x=[]
        t_y=[]
        for np in rtg.get_thread(tid):
            if np.node in node_x: pos=node_x[np.node]
            elif -np.node in node_x: pos=-node_x[-np.node]
            else: continue    
            if len(t_x) and ( ((last_pos>0) != (pos>0)) or (last_pos>pos)):
                y-=1
            last_pos=pos
            t_x.append(abs(pos))
            t_y.append(y)
        #print(t_y,t_x)
        plot(t_x,t_y,'o-',linewidth=3)
    grid(color='black', axis='x',which='both', linestyle=':', linewidth=1)
    grid(color='black', axis='y',which='both', linestyle=':', linewidth=1, alpha=.5)
    xticks(range(len(node_x)),labels=[node_label(x[0]) for x in sorted(node_x.items(),key=lambda v:v[1])],rotation=90)
    yticks(yticks_y,labels=yticks_l)

tids=set()
for nid in args.node:
    for tid in rtg.node_threads(nid[0],True): tids.add(tid)
nodes=set()

figure(figsize=(40,80))
thread_plot(rtg,tids)
savefig(f'{args.output_prefix}_threadplot.png')


# tids=set()
# for nid in args.node:
#     for tid in rtg.node_threads(nid[0],True): tids.add(tid)

# threads={tid: rtg.get_thread_nodes(tid) for tid in tids}
# node_threads={}
# with open(f'{args.output_prefix}_thread_nodes.txt',"w") as tnf:
#     for tid,nl in threads.items():
#         tnf.write(f'{tid}: {nl}\n')
#         for nid in nl:
#             if abs(nid) not in node_threads:
#                 node_threads[abs(nid)]=[]
#             node_threads[abs(nid)].append(abs(tid))
#     tnf.write('seq'+',seq'.join(map(str,node_threads.keys()))+'\n')
#     if args.input_dg:
#         dg=SDG.DistanceGraph(ws.sdg)
#         dg.load(args.input_dg)
#         sel_nids=set(node_threads.keys()).intersection(set(nv.node_id() for nv in dg.get_all_nodeviews(include_disconnected=False)))
#         tnf.write('seq'+',seq'.join(map(str,sel_nids))+'\n')

# with open(f'{args.output_prefix}_thread_nodes.csv',"w") as tnf:
#     tnf.write(f'node,{",".join(str(x) for x in tids)}\n')
#     for node,ntids in node_threads.items():
#         s='seq'+str(node)
#         for tid in tids:
#             if tid in ntids: s+=',X'
#             else: s+=','
#         tnf.write(s+'\n')


    

# with open(f'{args.output_prefix}_joined.csv','w') as of:
#         of.write('Node,KCI,Colour\n')
#         for x in dg.get_all_nodeviews(include_disconnected=False):
#             nid=x.node_id()
#             kci=x.kci()
#             if kci<.5: c='gray'
#             elif kci<1.5: c='green'
#             elif kci<2.5: c='blue'
#             else: c='red'
#             of.write(f'seq{nid},{kci :.2f},{c}\n')
