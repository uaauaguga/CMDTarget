#!/usr/bin/env python
import argparse
import os
import pandas as pd
from collections import defaultdict
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('target prediction')

def main():
    parser = argparse.ArgumentParser(description='performance evaluation')
    parser.add_argument('--input-directory',  '-id', type=str, required=True, help='input predicted sRNA scores')
    parser.add_argument('--known-targets',  '-kt', type=str, required=True, help='directory contains known targets')
    parser.add_argument('--output','-o', type=str , help="where to save performance")
    parser.add_argument('--validated','-v', action = "store_true", help="only use validated ones")
    parser.add_argument('--score-field','-sf', default="denoised score", help="field to use")
    parser.add_argument('--topk','-tk', type = int, default = 200, help="top predictions to consider")
    parser.add_argument('--min-homologs','-mh', type = int, default = 0, help="only consider targets with at least some homolog")
    parser.add_argument('--srna-ids','-si', required=True, help="sRNA ids to consider")
    args = parser.parse_args()

    logger.info("Load known interactions ...")

    known_targets = defaultdict(list)
    for txt in os.listdir(args.known_targets):
        path = os.path.join(args.known_targets,txt)
        sRNA_id = txt[:-4]        
        with open(path) as f:
            for line in f:                
                protein_id, evidence = line[:-1].split("\t")
                if args.validated and ("validated" not in evidence):
                    continue
                known_targets[sRNA_id].append(protein_id)

    logger.info("Evaluate predictions ...")

    sRNA_ids = open(args.srna_ids).read().strip().split("\n")
    sRNA_ids = set(sRNA_ids)


    fout = open(args.output,"w")
    n_detected = 0
    total_number = 0
    total_recall = 0
    for sRNA_id in sRNA_ids:
        if sRNA_id not in known_targets:
            #logger.info(f"{sRNA_id} not in known targets, skip it.")
            continue
        path = os.path.join(args.input_directory,sRNA_id+".txt")            
        if not os.path.exists(path):
            #logger.info(f"{sRNA_id} not in predicted targets, skip it.")
            continue
        table = pd.read_csv(path,sep="\t")
        if "homolog counts" in table.columns:
            table = table[table["homolog counts"] >= args.min_homologs]
        if "protein id" in table.columns:
            table.index = table["protein id"].map(lambda x:x.split(":")[1].split("-")[0])
        else:
            table = table.set_index(table.iloc[:,0])
        if args.score_field in table.columns:
            table = table.sort_values(by=args.score_field,ascending=False)
        else:
            table["key"] = table.apply(eval (args.score_field),axis=1)
            table = table.sort_values(by="key",ascending=False)
        recalls = table.index.isin(known_targets[sRNA_id]).cumsum()
        #total_recall += recalls[args.topk]/len(known_targets[sRNA_id])
        n_detected += recalls[args.topk]
        total_number += len(known_targets[sRNA_id])
        recall = round(recalls[args.topk-1]/len(known_targets[sRNA_id]),4)     
        print(sRNA_id, recall, recalls[args.topk-1], len(known_targets[sRNA_id]), sep="\t", file=fout)
        #print(sRNA_id, recall, sep="\t")
    fout.close()
    total_recall = n_detected/total_number
    logger.info(f"Total recall  is: {total_recall}.")
    logger.info("All done.")         



if __name__ == "__main__":
    main()
