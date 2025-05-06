#!/usr/bin/env python
import numpy as np
import pandas as pd
import gseapy as gp
import re
import argparse
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("infer sRNA function")
import os

def get_gene_id(s):
    contig_id, start, end, gene_id = re.match(r"CDS.(.+):(\d+)-(\d+)\([+-]\)\.(.+)",s).groups()    
    return gene_id

def main():
    parser = argparse.ArgumentParser(description='infer sRNA function')
    parser.add_argument('--input', '-i', required=True, help="input count matrix")
    parser.add_argument('--output', '-o', required=True, help="output association")
    parser.add_argument('--annotation', '-a', required=True, help="eggnogmapper annotation")
    parser.add_argument('--metric', '-m', required=True, help="field use as metric")
    args = parser.parse_args()

    logger.info("Load scores ...")
    table = pd.read_csv(args.input,sep="\t",index_col=0)

    logger.info("Load annotations ...") 
    pathway_id2name = {}
    with open("reference/KEGG/pathway") as f:
        for line in f:
            pathway_id, pathway_name = line[:-1].split("\t")
            pathway_id2name[pathway_id] = pathway_name

    ko_id2name = {}
    with open("reference/KEGG/ko") as f:
        for line in f:
            ko_id, ko_name = line[:-1].split("\t")  
            ko_id2name[ko_id] = ko_name        

    ko_id2pathway_ids = {} 
    for pathway_id in pathway_id2name:
        with open(f"reference/KEGG/pathway2ko/{pathway_id}") as f:
            for line in f:
                line = line.strip()
                if len(line) == 0:
                    continue
                ko_id = line.split("\t")[1]
                if ko_id not in ko_id2pathway_ids:
                   ko_id2pathway_ids[ko_id] = set()
                ko_id2pathway_ids[ko_id].add(pathway_id)
    pathway_id2gene_ids = {}
    with open(args.annotation) as f:
        for line in f:
            if line.startswith("#"):
                continue 
            fields = line[:-1].split("\t")
            gene_id = fields[0] 
            ko_ids = fields[11].split(",")
            for ko_id in ko_ids:
                if ko_id in ko_id2pathway_ids:
                    for pathway_id in ko_id2pathway_ids[ko_id]:
                        if pathway_id not in pathway_id2gene_ids:
                            pathway_id2gene_ids[pathway_id] = set()
                        pathway_id2gene_ids[pathway_id].add(gene_id)
    metric = args.metric
    scores = table[metric]
    scores = scores[~scores.isna()]
    scores.index = scores.index.map(lambda x:x.split(":")[1])
    scores = scores.sort_values(ascending=False)        
    scores = scores + np.random.rand(scores.shape[0])*0.00001
    res = gp.prerank(rnk=scores,gene_sets=pathway_id2gene_ids,threads=4,min_size=5,max_size=10000,permutation_num=5000)
    res_table = res.res2d
    res_table["pathway"] = res_table["Term"].map(lambda x:pathway_id2name[x])        
    res_table.to_csv(args.output,sep="\t",index=False)
    print(res_table[res_table["NES"]>0].loc[:,["Name","ES","NES","NOM p-val","FDR q-val","pathway"]].sort_values(by="NOM p-val").head(50)) 

    

if __name__ == "__main__":
    main()
