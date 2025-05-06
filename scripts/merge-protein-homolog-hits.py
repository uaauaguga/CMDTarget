#!/usr/bin/env python
import argparse
import sys
import os
import subprocess
from collections import defaultdict
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('merge')

def main():
    parser = argparse.ArgumentParser(description='Merge homolog hits of multiple query genome')
    parser.add_argument('--working-directory',  '-wd', type=str, required=True, help='The working directory')
    parser.add_argument('--query-ids',  '-qi', type=str, required=True, help='Query genome ids')
    parser.add_argument('--genome-set',  '-gs', type=str, required=True, help='Genome set')
    parser.add_argument('--output','-o', type=str, required=True , help="Merged table")
    parser.add_argument('--collapse','-c', type=str, required=True , help="Output collapsing table")
    args = parser.parse_args()
    fc = open(args.collapse,"w")
    fo = open(args.output,"w")
    profiles = {}
    logger.info("Load homolog hits ...")   
    query_ids = args.query_ids.split(",") 
    for query_id in query_ids: #open(args.query_ids).read().strip().split("\n"):
        path = f"{args.working_directory}/{query_id}/homologs/{args.genome_set}/best.hit.txt" 
        #asm_id = path.split("/")[1] 
        asm_id = query_id
        logger.info(f"Processing {path} ...")
        hits = defaultdict(list)
        with open(path) as f:
            for line in f:
                query_id, hit_id = line.strip().split("\t")
                hits[query_id].append(hit_id)
        for query_id in hits:
            hit_ids = sorted(hits[query_id])
            asm_ids = [hit_id.split(":")[0] for hit_id in hit_ids]
            asm_ids = set(asm_ids)
            if len(asm_ids) < 3:
                continue
            profile = ";".join(hit_ids)
            if profile in profiles:
                proxy_query_id = profiles[profile]
                #print(asm_id + ":" + query_id, proxy_query_id,profile)
                print(asm_id + ":" + query_id, proxy_query_id, sep="\t", file=fc)
            else:
                profiles[profile] = asm_id + ":" + query_id
                for hit_id in hit_ids:
                    print(asm_id + ":" + query_id, hit_id, sep="\t",file=fo)
    fc.close()
    fo.close()

if __name__ == "__main__":
    main()        
                               
