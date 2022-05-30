#!/usr/bin/env python3 -x
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
parser.add_argument("--input_ws", help="input workspace", type=str, required=True)
parser.add_argument("--input_dg", help="input distance graph", type=str, required=True)
parser.add_argument("-u", "--unique_coverage", help="value for unique coverage on kci", type=int, required=True)
parser.add_argument("--input_rtg", help="input read threads", type=str, default="")
parser.add_argument("--input_lrr", help="input long read mappings", type=str, default="")
parser.add_argument("--join_gaps", help="keep gaps at negative distances",action='store_true',default=False)
args = parser.parse_args()

def find_reasonable_overlap(s1,s2,max_size,min_size=21):
    for x in range(max_size,min_size-1,-1):
        if s1[-x:]==s2[:x]: return x
    return 0

def join_graph(g,join_gaps=False,max_overlap=200,keep_old_nodes=False):
    node_groups=[]
    total_line_nodes=0
    for line in g.get_all_lines(2):
        node_groups+=[[line[0]]]
        total_line_nodes+=len(line)
        for i in range(1,len(line)):
            
            pnid=line[i-1]
            pnv=g.get_nodeview(pnid)
            nid=line[i]
            d=pnv.next()[0].distance()
            if d<0:
                d=max(63,int(d*1.5))
                nv=g.get_nodeview(nid)
                d=-find_reasonable_overlap(pnv.sequence(),nv.sequence(),d)
            if d==0: d=10 #put 10 Ns 
            if d<0 or join_gaps:
                node_groups[-1]+=[d,nid]
            else:
                node_groups+=[[nid]]
    #TODO: replace every node group for a single node representing them and migrate connections
    print(f"joining {total_line_nodes} nodes into {len(node_groups)} groups")
    for ng in node_groups:
        if len(ng)==1: continue
        pcons=[(-l.node().node_id(),l.distance()) for l in g.get_nodeview(ng[0]).prev()]
        ncons=[(l.node().node_id(),l.distance()) for l in g.get_nodeview(ng[-1]).next()]
        s=g.get_nodeview(ng[0]).sequence()
        for i in range(1,len(ng),2):
            d=ng[i]
            nid=ng[i+1]
            if d>=0:
                s+='N'*d+g.get_nodeview(nid).sequence()
            else:
                s+=g.get_nodeview(nid).sequence()[-d:]
        #print(pcons,len(s),ncons)
        new_nid=g.sdg.add_node(s)
        for pc in pcons:
            g.add_link(pc[0],new_nid,pc[1])
        for nc in ncons:
            g.add_link(-new_nid,nc[0],nc[1])
        for nid in ng[::2]:
            g.disconnect_node(nid)
            if not keep_old_nodes: g.sdg.remove_node(nid)

#open ws
ws=SDG.WorkSpace(args.input_ws)
dg=SDG.DistanceGraph(ws.sdg)

#open distance graph
dg.load(args.input_dg)
join_graph(dg,args.join_gaps)

dg.write_to_gfa1(f'{args.output_prefix}_joined.gfa',selected_nodes=[x.node_id() for x in dg.get_all_nodeviews(include_disconnected=False)])

kc=ws.get_kmer_counter("main")
kc.set_kci_peak(args.unique_coverage)
kc.update_graph_counts()
kc.compute_all_kcis()

with open(f'{args.output_prefix}_joined.csv','w') as of:
        of.write('Node,KCI,Colour\n')
        for x in dg.get_all_nodeviews(include_disconnected=False):
            nid=x.node_id()
            kci=x.kci()
            if kci<.5: c='gray'
            elif kci<1.5: c='green'
            elif kci<2.5: c='blue'
            else: c='red'
            of.write(f'seq{nid},{kci :.2f},{c}\n')
