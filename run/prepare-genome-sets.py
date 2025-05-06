#!/usr/bin/env python
from collections import defaultdict
import numpy as np
np.random.seed(666)
def main():
    #ds = ["95","96","97","98","99","99.5"]
    #ds = ["99.5","99","98"]
    #ks = [4,8,16,32]
    ds = ["99.99"]
    #ks = [4,8,16,32]
    ks = [8]
    for d in ds:
        for k in ks:
            path = f"genomes/SSU-{d}/hits.txt" 
            identities = defaultdict(list)
            with open(path) as f:
                for line in f:
                    fields = line[:-1].split("\t")
                    query_id, hit_id = fields[:2]
                    identity = float(fields[2])
                    qgenome_id, tgenome_id = query_id.split(":")[0], hit_id.split(":")[0]
                    identities[qgenome_id].append((tgenome_id,identity))
            ss = {}
            for qgenome_id in identities:
                genome_ids = []
                gsetname = f"{qgenome_id}-{d}-top-{k}"
                for tgenome_id,identity in identities[qgenome_id][:k]:
                    genome_ids.append(tgenome_id)
                    #print(d, qgenome_id,tgenome_id,identity,sep="\t")
                with open(f"genome-selection/{gsetname}.txt","w") as f:
                    print(qgenome_id,file=f)
                    for genome_id in genome_ids:
                        print(genome_id,file=f)
                #"""
                if len(identities[qgenome_id]) > 10:
                    genome_ids = identities[qgenome_id]
                    for j in range(3):
                        np.random.shuffle(genome_ids)
                        gsetname = f"{qgenome_id}-{d}-rand-{k}-{j}"
                        with open(f"genome-selection/{gsetname}.txt","w") as f:
                            print(qgenome_id,file=f)
                            for genome_id in genome_ids[:k]:
                                print(genome_id[0],file=f)                    
                #"""
                #s = ",".join(sorted(genome_ids))
                #print(s)
                #if s in ss:
                #    print(gsetname, ss[s])
                #    continue
                #ss[s] = gsetname
                #print(genome_ids)
                

if __name__ == "__main__":
    main()
