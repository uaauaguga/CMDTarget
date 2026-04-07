#!/usr/bin/env python
import argparse
import pandas as pd
import os
from ete3 import Tree
import numpy as np
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('extract scores')

def score_to_quantile(xs):
    xs = np.array(xs)
    mask = ~np.isnan(xs)
    xs_not_nan = xs[mask]
    ranks = np.argsort(np.argsort(xs_not_nan))
    qs = np.full_like (xs,np.nan)
    qs[mask] = ranks/ranks.max()
    #qs[mask] = norm.ppf(ranks/ranks.max())
    return qs

def get_median_rate(tree):
    rates = []
    for n in tree.traverse(strategy="postorder"):
        if n.up is None:
            continue
        rate = ((float(n.score) - float(n.up.score))**2)/n.dist
        rates.append(rate)
    median_rate = np.median(rates)
    return median_rate

def get_scores(tree):
    genome_id2score = {}
    for node in tree.get_leaves():
        genome_id2score[node.name] = float(node.score)
    return genome_id2score, float(tree.score)

def main():
    parser = argparse.ArgumentParser(description='extract score from comparative denoising scores')
    parser.add_argument('--input', '-i', required=True, help="input comparative denoising results")
    parser.add_argument('--output-directory','-od', required=True, help="output directory")
    args = parser.parse_args()

    if not os.path.exists(args.output_directory):
        logger.info(f"{args.output_directory} not exists, create it .")
        os.mkdir(args.output_directory)

    logger.info("Load scores ...")
    median_rates = {}
    root_scores = {}
    records = []
    with open(args.input) as f:
        _ = next(f)
        for line in f:
            fields = line[:-1].split("\t")
            raw_score = float(fields[2]) 
            gene_id = fields[0]
            genome_id = fields[1].split(":")[0]
            records.append(("raw",gene_id,genome_id,raw_score))
            if fields[4] == ".":
                continue
            tree = Tree(fields[4].replace(")ROOT[",")["))
            scores, root_score = get_scores(tree)
            root_scores[gene_id] = root_score
            median_rate = get_median_rate(tree)
            median_rates[gene_id] = median_rate
            for genome_id in scores:
                records.append(("denoised",gene_id,genome_id,scores[genome_id])) 
    logger.info("Summarize scores by homolog .")
    denoised_score_by_homolog = pd.DataFrame.from_records(records)
    denoised_score_by_homolog.columns = ["data","gene id","genome id","score"]
    logger.info("Save raw scores ...")
    raw_score_by_homolog = denoised_score_by_homolog[denoised_score_by_homolog["data"] == "raw"]
    raw_score_by_homolog_matrix = raw_score_by_homolog.pivot(index="gene id",columns="genome id",values="score")
    raw_score_by_homolog_matrix.round(4).to_csv(args.output_directory + "/raw-zscore.txt",sep="\t")
    logger.info("Save raw quantiles ...")
    raw_quantile_by_homolog_matrix = raw_score_by_homolog_matrix.apply(lambda x:score_to_quantile(x.values),axis=0)
    raw_quantile_by_homolog_matrix.round(4).to_csv(args.output_directory + "/raw-quantile.txt",sep="\t")
    denoised_score_by_homolog = denoised_score_by_homolog[denoised_score_by_homolog["data"] == "denoised"]
    denoised_score_by_homolog_matrix = denoised_score_by_homolog.pivot(index="gene id",columns="genome id",values="score")    
    logger.info("Save denoised scores ...")
    denoised_score_by_homolog_matrix.round(4).to_csv(args.output_directory + "/denoised-zscore.txt",sep="\t")
    stds = denoised_score_by_homolog_matrix.std(axis=1).to_dict()#/denoised_score_by_homolog_matrix.mean(axis=1)).to_dict()
    denoised_quantile_by_homolog_matrix = denoised_score_by_homolog_matrix.apply(lambda x:score_to_quantile(x.values),axis=0)
    logger.info("Save denoised quantiles ...")
    denoised_quantile_by_homolog_matrix.round(4).to_csv(args.output_directory + "/denoised-quantile.txt",sep="\t")
    logger.info("Save rates ...")
    with open(args.output_directory + "/stats.txt","w") as fout:
        print("gene id","root score","std","median rate",sep="\t",file=fout)
        for gene_id in median_rates:
            print(gene_id, root_scores[gene_id], stds[gene_id], median_rates[gene_id], sep="\t", file=fout)
    fout.close()
    logger.info("All done.")
    
if __name__ == "__main__":
    main()
