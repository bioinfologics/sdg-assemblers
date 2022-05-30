#!/usr/bin/env python3 -u
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
parser.add_argument("-l", "--long_datastore", help="long reads datastore", type=str, required=True)
parser.add_argument("-u", "--unique_coverage", help="value for unique coverage at 31-mers", type=int, required=True)
parser.add_argument("--map_k", help="long read mapping k", type=int, default=25)
parser.add_argument("--lr_min_support", help="long read min support to expand canonical repeats", type=int, default=5)
parser.add_argument("--lr_snr", help="long read SNR to expand canonical repeats", type=int, default=10)
parser.add_argument("--lr_max_noise", help="long read max_noise to expand canonical repeats", type=int, default=10)
parser.add_argument("--rounds", help="rounds of repeat resolution", type=int, default=2)
parser.add_argument("--pe_rep_min_support", help="paired read min support to expand canonical repeats", type=int, default=5)
parser.add_argument("--pe_rep_max_noise", help="paired read max noise to expand canonical repeats", type=int, default=5)
parser.add_argument("--pe_rep_snr", help="paired read snr to expand canonical repeats", type=int, default=10)
parser.add_argument("--final_remap", help="remap short and long reads at the end (default: false, any mappings are invalid)", type=bool, default=False)
args = parser.parse_args()

ws=SDG.WorkSpace(f'{args.output_prefix}_04_strided.sdgws')
peds=ws.get_paired_reads_datastore(ws.list_paired_reads_datastores()[0])
lords = ws.add_long_reads_datastore(args.long_datastore)

kc=ws.get_kmer_counter("main")
kc.set_kci_peak(args.unique_coverage)
kc.update_graph_counts()
print(ws.sdg.stats_by_kci())


print_step_banner("LONG READ REPEAT RESOLUTION")
    
# def solve_with_llr2(nid,min_support=5,max_noise=2,snr=3,verbose=False):    
#     nA,nB=[x.node().node_id() for x in ws.sdg.get_nodeview(nid).next()]
#     pA,pB=[x.node().node_id() for x in ws.sdg.get_nodeview(nid).prev()]
#     apA,apB,anA,anB=abs(pA),abs(pB),abs(nA),abs(nB)
#     aset=set([apA,apB,anA,anB])
#     tAA=0
#     tAB=0
#     tBA=0
#     tBB=0
#     all_voting_tids=set(rtg.node_threads(pA).union(rtg.node_threads(pB))).intersection(set(rtg.node_threads(nA).union(rtg.node_threads(nB))))
#     for tid in all_voting_tids:
        
#         c=Counter([abs(x.node) for x in lrr.read_perfect_matches[tid] if abs(x.node) in aset])
#         if c[apA]>c[apB]*snr:
#             if c[anA]>c[anB]*snr: tAA+=1
#             elif c[anB]>c[anA]*snr: tAB+=1
#         elif c[apB]>c[apA]*snr:
#             if c[anA]>c[anB]*snr: tBA+=1
#             elif c[anB]>c[anA]*snr: tBB+=1
#     if verbose: print(nid,": ",tAA,tBB,tAB,tBA)
#     if min(tAA,tBB)>=min_support and max(tAB,tBA)<=max_noise and min(tAA,tBB)>=max(tAB,tBA)*snr:
#         return [(-pA,nA),(-pB,nB)]
#     if min(tAB,tBA)>=min_support and max(tAA,tBB)<=max_noise and min(tAB,tBA)>=max(tAA,tBB)*snr:
#         return [(-pA,nB),(-pB,nA)]
#     return []

def voting_rids(n1,n2,snr=3):
    rids1=set()
    rids2=set()
    if n1 not in node_rids or n2 not in node_rids: return (rids1,rids2)
    #decide where a rid goes between two competing nodes
    for rid in set(node_rids[n1].keys()).union(set(node_rids[n2].keys())):
        n1c=node_rids[n1][rid]
        n2c=node_rids[n2][rid]
        #print(n1c,n2c,end=' ')
        if n2c*snr<=n1c:
            rids1.add(rid)
            #print(' Node 1',end=' ')
        elif n1c*snr<=n2c:
            rids2.add(rid)
            #print(' Node 2',end=' ')
        #print()
    return (rids1,rids2)

def solve_with_voting_rids(nid,min_support=5,max_noise=2,snr=3):
    nA,nB=[x.node().node_id() for x in ws.sdg.get_nodeview(nid).next()]
    pA,pB=[x.node().node_id() for x in ws.sdg.get_nodeview(nid).prev()]
    rnA,rnB=voting_rids(nA,nB,snr)
    rpA,rpB=voting_rids(pA,pB,snr)
    vAA=len(rnA.intersection(rpA))
    vAB=len(rnA.intersection(rpB))
    vBA=len(rnB.intersection(rpA))
    vBB=len(rnB.intersection(rpB))
    
    #print(vAA,vAB,vBA,vBB)
    if min(vAA,vBB)>=min_support and max(vAB,vBA)<=max_noise and min(vAA,vBB)>=max(vAB,vBA)*snr:
        return [(-pA,nA),(-pB,nB)]
    if min(vAB,vBA)>=min_support and max(vAA,vBB)<=max_noise and min(vAB,vBA)>=max(vAA,vBB)*snr:
        return [(-pA,nB),(-pB,nA)]
    return []

lrr=SDG.LongReadsRecruiter(ws.sdg,lords,k=args.map_k)

for sp in range(1,args.rounds+1):
    lrr.map()
    node_rids={}
    for rid,rms in enumerate(lrr.read_perfect_matches):
        for m in rms:
            if m.node not in node_rids:
                node_rids[m.node]=Counter()
                node_rids[-m.node]=Counter()
            node_rids[m.node][rid]+=1
            node_rids[-m.node][-rid]+=1

    ge=SDG.GraphEditor(ws)

    total_cr=0
    solved_cr=0
    for x in ws.sdg.get_all_nodeviews():
        nid=x.node_id()
        if x.is_canonical_repeat():
            total_cr+=1
            sol=solve_with_voting_rids(nid,args.lr_min_support,args.lr_max_noise,args.lr_snr)
            if sol!=[]:
                ge.queue_node_expansion(nid,sol)
                solved_cr+=1
            if total_cr%100==0: print("%d / %d canonical repeats solved"%(solved_cr,total_cr),flush=True)

    print("%d / %d canonical repeats solved"%(solved_cr,total_cr))
    ge.apply_all()
    ws.sdg.join_all_unitigs()
    kc.update_graph_counts()
    kc.compute_all_kcis()
    print(ws.sdg.stats_by_kci())

print_step_banner("PE CANONICAL REPEAT RESOLUTION")
peds.mapper.path_reads()
c=SDG.GraphContigger(ws)
c.solve_canonical_repeats(peds,min_support=args.pe_rep_min_support,max_noise=args.pe_rep_max_noise,snr=args.pe_rep_snr)
kc.update_graph_counts()
print(ws.sdg.stats_by_kci())

ws.sdg.write_to_gfa1(args.output_prefix+"_05_long_repeats.gfa")
if args.final_remap:
    peds.mapper.path_reads()
    lrr.map()
    lrr.dump(f'{args.output_prefix}_05_long_repeats.lrr')
ws.dump(f'{args.output_prefix}_05_long_repeats.sdgws')
