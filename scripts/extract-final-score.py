#!/usr/bin/env python
from collections import defaultdict
import argparse
import os
import logging
import numpy as np
from copy import copy
import pandas as pd
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('final score')
import json

def main():
    parser = argparse.ArgumentParser(description='extract final scores')
    parser.add_argument('--comparative-scores','-cs',type=str, required=True, help="results of compatative scoring")
    parser.add_argument('--genome-scores','-gs',type=str, required=True, help="scores by genome")
    parser.add_argument('--collapse','-c',type=str, required=True, help='collapsing table of genes in query genomes')   
    parser.add_argument('--working-directory','-wd',type=str, required=True, help='working directory of the scoring') 
    parser.add_argument('--weights','-w',type=str, required=True, help='weights for scoring') 
    parser.add_argument('--srna-name','-sn',type=str, required=True, help='name of sRNA')
    parser.add_argument('--tag','-t',type=str, default = "default", help="tag for the prediction") 
    parser.add_argument('--representative-genome-ids','-rgi',type=str, required=True, help="representative genome ids") 
    parser.add_argument('--conservation-weight','-cw',type=float, default=0, help="weight of conservation score") 
    parser.add_argument('--hfq-label','-hl',type=str, required=True, help="hfq label") 
    args = parser.parse_args()
 
  
    logger.info("Load comparative scoring results ...")
    scores = pd.read_csv(args.comparative_scores,sep="\t")
    logger.info("Load collapsing table ...")
    expand_table = {}
    with open(args.collapse) as f:
        for line in f:
            q1, q2 = line.strip().split("\t")
            if q2 not in expand_table:
                expand_table[q2] = []
            expand_table[q2].append(q1)

    genomes_with_sRNA = set()
    genomes_with_query_proteins = defaultdict(set)
    logger.info("Extract conservation of query genes ...")
    representative_genome_ids = open(args.representative_genome_ids).read().strip().split("\n")
    representative_genome_ids = set(representative_genome_ids)
    print(representative_genome_ids)
    with open(args.genome_scores) as f:
        _ = next(f)
        for line in f:
            fields = line.strip().split("\t")
            query_id, target_id = fields[:2]
            genome_id = target_id.split(":")[0]
            if genome_id not in representative_genome_ids:  
                continue
            genomes_with_sRNA.add(genome_id)
            genomes_with_query_proteins[query_id].add(genome_id)
            if query_id in expand_table:
                for expand_id in expand_table[query_id]:
                    genomes_with_query_proteins[expand_id].add(genome_id)

    logger.info("Retrieve scores for collapsed genes ...")
    records = []
    for record in scores.to_records(index=False):
        record = list(record)
        query_id = record[0]
        records.append(record)
        record = copy(record) 
        if query_id in expand_table:
            for expand_id in expand_table[query_id]:
                record[0] = expand_id
                records.append(copy(record))                

    n_rep_genome_with_sRNA = len(genomes_with_sRNA)
    logger.info(str(n_rep_genome_with_sRNA) + " representative genomes contains ncRNA homolog.")
    for query_id in list(genomes_with_query_proteins.keys()):
        genomes_with_query_proteins[query_id] = genomes_with_query_proteins[query_id].intersection(genomes_with_sRNA) 

    columns = scores.columns
    scores = pd.DataFrame.from_records(records)
    scores.columns = columns
    scores["query genomes"] = scores["query id"].map(lambda x:x.split(":")[0])
    genome_ids = sorted(list(scores["query genomes"].unique()))
    scores["target genomes"] = scores["target id"].map(lambda x:x.split(":")[0])
    scores = scores[scores["query genomes"] == scores["target genomes"]]
    #scores = scores[scores["tree"]!="."]
    logger.info("Extract scores of query genomes ...")
    denoised_scores = scores.set_index("query id")["denoised"].to_dict()
    print(denoised_scores)
    logger.info("Insert denoised score to final results ...")
    srna_name = args.srna_name
    weights = json.load(open(args.weights))

    for genome_id in genome_ids:
        path = f"{args.working_directory}/{genome_id}/combined.{args.hfq_label}/{srna_name}.txt"
        if not os.path.exists(path):
            logger.info(f"{srna_name} not present in {genome_id}, skip it.")    
            continue
        logger.info(f"Process {genome_id} ...")            
        if not os.path.exists(f"{args.working_directory}/{genome_id}/combined.final"):
            os.makedirs(f"{args.working_directory}/{genome_id}/combined.final",exist_ok=True)
        if not os.path.exists(f"{args.working_directory}/{genome_id}/combined.final/{args.tag}"):
            os.makedirs(f"{args.working_directory}/{genome_id}/combined.final/{args.tag}",exist_ok=True)
        logger.info("Load scores ...")
        table = pd.read_csv(path,sep="\t")
        table["score"] = 0
        table["protein id"] = genome_id + ":" + table["protein id"]
        for key in weights:
            if weights[key] == 0:
                logger.info(f"{key} has zero weights, skip it.")
                continue
            table["score"] += table[key].fillna(0)*weights[key]
        #if args.normalize:
        #    mean = table["score"].mean()
        #    std = table["score"].std()
        #    table["score"] = (table["score"] - mean)/std
        table["homolog counts"] = table["protein id"].map(lambda x:len(genomes_with_query_proteins.get(x,set())))
        table["homolog scores"] = round(table["homolog counts"]/n_rep_genome_with_sRNA,4)
        table["score"] = table["score"].round(4)
        print(table["protein id"])
        table["denoised score"] = table["protein id"].map(lambda x:denoised_scores.get(x,np.nan)).round(4)        
        table["final score"] = table["denoised score"]
        if args.conservation_weight != 0:
            panelty = table["homolog scores"] - 0.4
            panelty[panelty>0] = 0 
            table["final score"] += panelty*args.conservation_weight
        index = table.index[table["final score"].isna()]
        table.loc[index,"final score"] = table.loc[index,"score"]
        n_missed = table["denoised score"].isna().sum()
        N = table.shape[0]
        logger.info(f"Denoising impact {N - n_missed} in {N} genes.")
        print(f"{args.working_directory}/{genome_id}/combined.final/{args.tag}/{srna_name}.txt") 
        table.to_csv(f"{args.working_directory}/{genome_id}/combined.final/{args.tag}/{srna_name}.txt",sep="\t",index=False)
    logger.info("All done.")         
    
    
if __name__ == "__main__":
    main()
    
            
