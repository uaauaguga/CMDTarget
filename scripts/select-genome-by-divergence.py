#!/usr/bin/env python
import os
import argparse
from collections import defaultdict
import subprocess
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('select genomes')
def main():
    parser = argparse.ArgumentParser(description='Select genomes for comparative analysis')
    parser.add_argument('--input',  '-i', type=str, required=True, help='Input rRNAs')
    parser.add_argument('--query-ids',  '-qs', type=str, required=True, help='Genome id of query genomes')
    parser.add_argument('--cutoff',  '-c', type=float, default=0.98, help='rRNA sequence identity cutoff')
    parser.add_argument('--output-directory','-od', type=str, required=True , help="Output results")
    args = parser.parse_args()

    query_ids = args.query_ids.split(",")
    qsequences = {}
    tsequences = {}
    logger.info("Load sequences ...")
    with open(args.input) as f:
        for header in f:
            sequence = next(f).strip()
            seq_id = header[1:].strip().split(" ")[0]
            genome_id = seq_id.split(":")[0]
            if genome_id in query_ids:
                qsequences[seq_id] = sequence
            else:
                tsequences[seq_id] = sequence
    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)

    qpath = os.path.join(args.output_directory,"query.fa")
    with open(qpath,"w") as fout:
        for seq_id in qsequences:
            sequence = qsequences[seq_id]
            print(f">{seq_id}",file=fout)
            print(sequence,file=fout)
    tpath = os.path.join(args.output_directory,"target.fa")
    with open(tpath,"w") as fout:
        for seq_id in tsequences:
            sequence = tsequences[seq_id]
            print(f">{seq_id}",file=fout)
            print(sequence,file=fout)
    flog = open(os.path.join(args.output_directory,"cd-hit-est.log"),"w")
    trpath = os.path.join(args.output_directory,"target.reduce.fa")
    cmd = ["cd-hit-est","-c",str(args.cutoff),"-i",tpath,"-o",trpath,"-r","0","-d","10000"]
    logger.info("Reduce target sequence identity ...") 
    subprocess.run(cmd,stdout=flog,stderr=flog)
    flog.close()

    logger.info("Search reduced sequences ...")
    flog = open(os.path.join(args.output_directory,"mmseqs.log"),"w")
    hits =  os.path.join(args.output_directory,"hits.txt")
    cmd = ["mmseqs","easy-search","--search-type","3",qpath,trpath,hits,"tmp"] 
    subprocess.run(cmd,stdout=flog,stderr=flog)
    flog.close()

    logger.info("Select top genomes ...")
    scores_by_hits = defaultdict(list)
    with open(hits) as f:
        for line in f:
            fields = line[:-1].split("\t") 
            query_id, target_id, score = fields[:3]
            score = float(score)
            scores_by_hits[target_id].append(score)
    records = []
    for target_id, scores in scores_by_hits.items():
        score = sum(scores)/len(scores)
        records.append((target_id, score))

    fscores = open(os.path.join(args.output_directory,"scores.txt"),"w")
    for record in sorted(records,key=lambda x:-x[1]):
        print(*record,sep="\t",file=fscores)
    fscores.close()
        
    logger.info(f"Please check: " + os.path.join(args.output_directory,"scores.txt")) 
    

     
if __name__ == "__main__":
    main()
