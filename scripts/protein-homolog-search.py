#!/usr/bin/env python
import argparse
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("protein homolog search")
import subprocess
import os
from collections import defaultdict
import sys

def main():
    parser = argparse.ArgumentParser(description='protein homolog search')
    parser.add_argument('--query', '-q',required=True,help="query proteins")
    parser.add_argument('--database','-db',required=True,help="target database")
    parser.add_argument('--proteins','-p',required=True,help="directions contain proteins")
    parser.add_argument('--output-directory','-od',required=True,help="output directory")
    parser.add_argument('--coverage','-c', type = float,default=0.9,help="alignment coverage required")
    parser.add_argument('--threads','-t',default=8,help="threads for mmseqs search")
    parser.add_argument('--rerun','-r',action="store_true",help="whether rerun homolog search even the results presenece")
    args = parser.parse_args()

    tdb = args.database
    if not os.path.exists(tdb + ".dbtype"):
        if args.proteins is None:
            logger.error("Target database not exists, sequences to build it should be provided.")
            sys.exit(1)   
        logger.info("Build target protein mmseqs database ...")
        cmd = ["mmseqs", "createdb", "--dbtype", "1", args.proteins, tdb ]
        subprocess.run(cmd)
    else:
        logger.info("Target database exists.")
                    
    qdb = os.path.join(args.output_directory,"proteins")
    if not os.path.exists(qdb + ".dbtype"):
        logger.info("build protein mmseqs database ...")
        cmd = ["mmseqs", "createdb", "--dbtype", "1", args.query, qdb ]
        subprocess.run(cmd)
    else:
        logger.info("protein db exists .")

    hdb = os.path.join(args.output_directory,"hits")
    if (not os.path.exists(hdb+".dbtype")) or args.rerun:
        logger.info("homolog search ...")
        cmd = ["mmseqs", "search", "--search-type", "1","-s", "7.5", "--max-seqs", "5000", qdb, tdb, hdb, hdb + ".tmp" ]
        subprocess.run(cmd)
    else:
        logger.info("hit db exists .")

    tsv = os.path.join(args.output_directory,"hits.tsv")
    if not os.path.exists(tsv) or args.rerun:
        logger.info("reformat homolog search hits ...")
        cmd = ["mmseqs", "convertalis", "--format-mode", "2", qdb, tdb, hdb, tsv ]
        subprocess.run(cmd)
    else:
        logger.info("hit file exists .")
            
    logger.info("extract best hits by genome ...")
    bitscore_by_genome = defaultdict(int)
    best_hits = {}
    with open(tsv) as f:
        for line in f:
            fields = line.strip().split("\t")
            hit_id, hstart, hend = fields[1], int(fields[8]), int(fields[9])
            query_id, qstart, qend = fields[0], int(fields[6]), int(fields[7])
            query_length, hit_length = int(fields[12]), int(fields[13])
            aligned_length = int(fields[3])
            coverage = min(aligned_length/query_length,aligned_length/hit_length)
            if coverage < args.coverage:
                continue     
            genome_id = hit_id[:hit_id.find(":")]
            bitscore = int(fields[11]) 
            if bitscore > bitscore_by_genome[(query_id, genome_id)]:
                best_hits[(query_id,genome_id)] = hit_id
                bitscore_by_genome[(query_id, genome_id)] = bitscore   

    logger.info("Save best hits ...")
    #best_hit = os.path.join(args.output_directory,f"best.hit.{args.coverage}.txt")
    best_hit = os.path.join(args.output_directory,f"best.hit.txt")
    fout = open(best_hit,"w")
    for query_id,genome_id in best_hits:
        hit_id = best_hits[(query_id,genome_id)]
        print(query_id, hit_id,file=fout,sep="\t")
    fout.close()

    logger.info("All done.")
              

if __name__ == "__main__":
    main()
