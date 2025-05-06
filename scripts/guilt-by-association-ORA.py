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
from scipy.stats import fisher_exact


def main():
    parser = argparse.ArgumentParser(description='infer sRNA function')
    parser.add_argument('--input', '-i', required=True, help="input correlation matrix")
    parser.add_argument('--output-directory', '-od', required=True, help="output correlation")
    parser.add_argument('--annotation', '-a', required=True, help="GO annotation")
    args = parser.parse_args()

    logger.info("Load correlations ...")
    correlations = pd.read_csv(args.input,sep="\t",index_col=0)

    logger.info("Load annotations ...") 
    pathway_id2name = {}
    with open("reference/go-basic.txt") as f:
        for line in f:
            pathway_id, pathway_name = line[:-1].split("\t")
            pathway_id2name[pathway_id] = pathway_name

    pathway_id2gene_ids = {}
    with open(args.annotation) as f:
        for line in f:
            fields = line[:-1].split("\t")
            gene_id = fields[3] 
            go_ids = fields[4].split(",")
            for go_id in go_ids:
                if go_id not in pathway_id2gene_ids:
                    pathway_id2gene_ids[go_id] = set()
                pathway_id2gene_ids[go_id].add(gene_id)
    #for pathway_id in pathway_id2gene_ids:
    #    pathway_id2gene_ids[pathway_id] = list(pathway_id2gene_ids[pathway_id])
    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)
    k = 250
    for srna_id in correlations.columns:
        logger.info(f"processing {srna_id} ...")
        scores = correlations[srna_id]
        scores = scores[~scores.index.duplicated()]
        scores = scores + np.random.rand(scores.shape[0])*0.001
        scores = scores.sort_values(ascending=False)        
        gene_list=set(list(scores.index[:k]))
        background = list(scores.index)
        """
        for pathway_id in pathway_id2gene_ids:
            pgene_ids = pathway_id2gene_ids[pathway_id]
            n = np.intersect1d(gene_list,pgene_ids).shape[0]
            t = [[n,k],[len(pgene_ids),len(background)]]
            res = fisher_exact(t, alternative='greater')
            if res[1] < 0.01:
                print(pathway_id, pathway_id2name[pathway_id],*res, sep="\t")
                #print(pgene_ids)
        """
        res = gp.prerank(rnk=scores,gene_sets=pathway_id2gene_ids,threads=4,min_size=5,max_size=10000,permutation_num=5000)
        res_table = res.res2d
        res_table["pathway"] = res_table["Term"].map(lambda x:pathway_id2name.get(x,x))        
        res_table.to_csv(os.path.join(args.output_directory,srna_id + ".txt"),sep="\t",index=False)
        print(res_table.loc[:,["Name","ES","NES","NOM p-val","FDR q-val","pathway"]].head(50)) 

    

if __name__ == "__main__":
    main()
