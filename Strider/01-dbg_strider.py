#!/usr/bin/env python3
import SDGpython as SDG
import argparse
from collections import Counter
import os

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output_prefix", help="prefix for output files", type=str, required=True)
parser.add_argument("-p", "--paired_datastore", help="paired reads datastore", type=str, required=True)
parser.add_argument("-k", "--k", help="k value for graph construction", type=int, default=63)
parser.add_argument("-u", "--unique_coverage", help="value for unique coverage at 31-mers", type=int, required=True)
parser.add_argument("-c", "--min_coverage", help="min coverage for graph construction", type=int, default=3)
parser.add_argument("-b", "--disk_batches", help="disk batches for graph construction", type=int, default=1)
parser.add_argument("--low_node_coverage", help="low coverage for short node cleanup", type=int, default=5)
parser.add_argument("--low_bubble_coverage", help="low coverage for bubble cleanup", type=int, default=5)

args = parser.parse_args()


K=args.k
MIN_CVG=args.min_coverage
NUM_BATCHES=args.disk_batches


ws = SDG.WorkSpace()
peds = ws.add_paired_reads_datastore(args.paired_datastore)
SDG.GraphMaker(ws.sdg).new_graph_from_paired_datastore(peds,K,MIN_CVG,NUM_BATCHES)
os.replace('small_K.freqs',f'{args.output_prefix}_small_K.freqs')

gc=SDG.GraphContigger(ws)
gc.clip_tips(300)
kc=ws.add_kmer_counter("main")
kc.add_count("pe",peds)

LOW_CVG=args.low_node_coverage
to_delete=set()
for nv in ws.sdg.get_all_nodeviews():
    if nv.size()>125: continue
    nvlc=len([x for x in nv.kmer_coverage('main','pe') if x<=LOW_CVG])
    if nvlc:
        to_delete.add(nv.node_id())

print(len(to_delete))
for x in to_delete:
    ws.sdg.remove_node(x)

ws.sdg.join_all_unitigs()
kc.update_graph_counts()

BUBBLE_SIZE=125
LOW_CVG=args.low_bubble_coverage

to_delete=[]
for nv in ws.sdg.get_all_nodeviews():
    if nv.size()<=BUBBLE_SIZE and len(nv.parallels())==1:
        on=nv.parallels()[0]
        if on.size()<=BUBBLE_SIZE and abs(on.node_id())>nv.node_id():
            nvlc=len([x for x in nv.kmer_coverage('main','pe') if x<=LOW_CVG])
            onlc=len([x for x in on.kmer_coverage('main','pe') if x<=LOW_CVG])
            if nvlc and not onlc:
                #print('del',nv,nvlc,on,onlc)
                to_delete.append(nv.node_id())
            if onlc and not nvlc:
                #print('del',on,onlc,nv,nvlc)
                to_delete.append(abs(on.node_id()))
print(len(to_delete))
for x in to_delete:
    ws.sdg.remove_node(x)

ws.sdg.join_all_unitigs()

kc.update_graph_counts()
#for x in enumerate(kc.count_spectra("pe",100)[:-1]): print(x)
kc.set_kci_peak(args.unique_coverage)
print(ws.sdg.stats_by_kci())

ws.sdg.write_to_gfa1(args.output_prefix+"_dbg.gfa")

peds.mapper.path_reads()

s=SDG.Strider(ws)
s.add_datastore(peds)
s.stride_from_anchors(min_size=1)

def strider_run_from_cpp():
    
    linear_anchors=set()
    for nv in ws.sdg.get_all_nodeviews():
        nid=nv.node_id()
        if s.is_anchor[nid] and len(s.routes_fw[nid])>1 and len(s.routes_bw[nid])>1: linear_anchors.add(nid)

    print("%d anchors are linear" % len(linear_anchors))

    ge=SDG.GraphEditor(ws)
    used_nodes=[]
    paths=0
    accepted=0
    rejected=0
    nopaths=0
    for fnid in list(linear_anchors):
        for nid in [fnid,-fnid]:
            #print("\nNode",nid)
            if nid>0: fp=s.routes_fw[nid][1:]
            else: fp=s.routes_bw[-nid][1:]
            #print("FW",fp)
            p=[]
            for x in fp:
                if abs(x) in linear_anchors:
                    p=fp[:fp.index(x)+1]
                    #print(p)
                    if x<0: op= [-x for x in s.routes_fw[-x][1:][::-1]]
                    else: op= [-x for x in s.routes_bw[x][1:][::-1]]
                    #print(op)
                    #print(p[:-1],"==",op[op.index(nid):],"??")
                    if nid not in op or p[:-1]!=op[op.index(nid)+1:]: p=[]
                    break
            #print(p)
            if p and abs(p[-1])>abs(nid):
                try:
                    SDG.SequenceDistanceGraphPath(ws.sdg,[nid]+p).sequence()
                    #print([nid]+p)
                    paths+=1
                    if ge.queue_path_detachment([nid]+p,True):
                        accepted+=1
                    else:
                        rejected+=1
                except:
                    nopaths+=1
                for r in p[:-1]: used_nodes.append(abs(r))
    print("%d paths ( %d accepted, %d rejected) and %d nopaths"%(paths,accepted,rejected,nopaths))
    ge.apply_all()

    #delete "used up" tips
    for r in range(10):
        ge=SDG.GraphEditor(ws)
        for x in used_nodes:
            try:
                nv=ws.sdg.get_nodeview(x)
                if len(nv.prev())==0 or len(nv.next())==0:
                    ge.queue_node_deletion(x)
            except: pass

        ge.apply_all()

    ws.sdg.join_all_unitigs()
    ge.remove_small_components(10,1000,10000)

    ### Finishing-off the assembly: typical heuristics that can be more aggressive when things already look fine:

    #an incredibly crude loop resolver, only supporting one round across the loop:
    kc.update_graph_counts()
    peds.mapper.path_reads()
    print(len(ws.sdg.get_all_nodeviews()))
    ge=SDG.GraphEditor(ws)
    for nv in ws.sdg.get_all_nodeviews():
        if len(nv.prev())!=1 or len(nv.next())!=1 or nv.prev()[0].node().node_id()!=nv.next()[0].node().node_id(): continue
        r=nv.next()[0].node().node_id()
        orn=[x.node().node_id() for x in nv.next()[0].node().next() if x.node().node_id()!=nv.node_id()]
        if len(orn)!=1: continue
        orn=orn[0]
        orp=[x.node().node_id() for x in nv.next()[0].node().prev() if x.node().node_id()!=nv.node_id()]
        if len(orp)!=1: continue
        orp=orp[0]
        print("Loop on %d -> %d -> %d -> %d -> %d"%(orp,r,nv.node_id(),r,orn))

        c=Counter()
        nc=Counter()
        nid=nv.node_id()
        for rid in peds.mapper.paths_in_node[abs(nid)]:
            if nid<0: rid=-rid
            c.update([tuple(peds.mapper.read_paths[abs(rid)].path)])
            nc.update(peds.mapper.path_fw(rid,nid))
        #for x in c.most_common(): print(x)
        #print()
        #for x in nc.most_common(): print(x)
        if nid not in nc and orn in nc:
            print("loop goes around just once!")
            ge.queue_node_expansion(r,[[-orp,nid],[-nid, orn]])

    ge.apply_all()
    ws.sdg.join_all_unitigs()    

    print(len(ws.sdg.get_all_nodeviews()))

    ## Short, low information node removal
    from statistics import median
    kc.update_graph_counts()
    #peds.mapper.path_reads()
    removed=0
    for nv in ws.sdg.get_all_nodeviews():
        if nv.size()>1000: continue
        rc=nv.kmer_coverage("main","pe")
        gc=nv.kmer_coverage("main","sdg")
        ukci=[]
        skci=[]
        for i in range(len(rc)):
            if gc[i]==1: ukci.append(rc[i])
            else: skci.append(rc[i]/gc[i])
        if ukci==[]: ukci=[-1]
        if skci==[]: skci=[-1]
        ms=median(skci)
        mu=median(ukci)
        #print(nv.node_id(),)
        if mu>-1 and mu<=4 and ms>5*mu:
            #print("Node %d (%d bp), shared_coverage=%0.2f, unique_coverage=%0.2f"%(nv.node_id(),nv.size(),ms,mu))
            removed+=1
            ws.sdg.remove_node(nv.node_id())
    print("%d low coverage, low information nodes removed"%removed)
    ws.sdg.join_all_unitigs()

strider_run_from_cpp()

kc.update_graph_counts()
kc.compute_all_kcis()
print(ws.sdg.stats_by_kci())
ws.sdg.write_to_gfa1(args.output_prefix+"_strided.gfa")
ws.dump(args.output_prefix+"_strided.sdgws")
