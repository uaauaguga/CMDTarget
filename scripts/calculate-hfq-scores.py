#!/usr/bin/env python
import argparse
import subprocess
import os
import logging
import numpy as np
from collections import defaultdict
from scipy.stats import norm
import pandas as pd
from scipy.stats import gumbel_r
import io
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('calcualte hfq scores')


def clip(p):
    if p > 0.99999:
        p = 0.99999
    if p < 0.00001:
        p = 0.00001
    return p

def main():
    parser = argparse.ArgumentParser(description='prepare gene specific scores')
    parser.add_argument('--input','-i',type=str, required=True, help="score of hfq binding in bed format")
    parser.add_argument('--gene','-g',type=str, required=True, help='gene interval in bed format, strandness is required')    
    parser.add_argument('--output','-o',type=str, required=True, help="where to save the scores")
    parser.add_argument('--contig', '-c', type=str, required=True, help="contig length")
    parser.add_argument('--upstream',  type=int, default = 300, required=True, help="upstream distance")
    parser.add_argument('--downstream',  type=int, default = 200, required=True, help="downstream distance")
    args = parser.parse_args()
    
    upstream = args.upstream - 100
    downstream = args.downstream - 100

    lengths = {}
    logger.info("Load contig length ...")
    with open(args.contig) as f:
        for line in f:
            seq_id, length = line.strip().split("\t")[:2]
            length = int(length)
            lengths[seq_id] = length
            
    
    logger.info("Load hfq scores ...")
    fwd_scores = {}
    rev_scores = {}
    with open(args.input) as f:
        for line in f:
            fields = line.strip().split("\t")
            seq_id, start, end = fields[:3]
            start, end = int(start), int(end)
            strand = fields[5]
            if seq_id not in fwd_scores:
                fwd_scores[seq_id] = np.zeros(lengths[seq_id])
                rev_scores[seq_id] = np.zeros(lengths[seq_id])
            current_score = float(fields[4])                
            if strand == "+":
                previous_scores = fwd_scores[seq_id][start:end]
                previous_scores[previous_scores<current_score] = current_score
                fwd_scores[seq_id][start:end] = previous_scores
            else:
                previous_scores = rev_scores[seq_id][start:end]
                previous_scores[previous_scores<current_score] = current_score
                rev_scores[seq_id][start:end] = previous_scores
                
    logger.info("Fit hfq score background distribution ...")
    if len(lengths) > 1:
        seq_ids = list(lengths.keys())
        probs = np.array(list(lengths.values()))
        probs = probs/probs.sum()        
        sampled_seq_ids = np.random.choice(seq_ids,5000,p=probs)   
        sampled_seq_ids = list(sampled_seq_ids)
    else:
        sampled_seq_ids = [list(lengths.keys())[0]]*5000
    
    background_hfq_scores = []
    for seq_id in sampled_seq_ids:
        if lengths[seq_id] < 300:
            continue
        start = np.random.randint(lengths[seq_id]-300)
        end = start + 300
        if np.random.rand() > 0.5:
            p = fwd_scores[seq_id][start:end].max()
        else:
            p = rev_scores[seq_id][start:end].max()       
        p = clip(p) 
        logit = np.log(p/(1-p))
        background_hfq_scores.append(logit)
    loc, scale = gumbel_r.fit(background_hfq_scores)      
    logger.info(f"Backgound distribution parameter of the logits is loc={loc}, scale={scale}.")
 
    logger.info("Calculate Z score of each leader sequences ...")   
    zscores = {}
    logits = {}
    with open(args.gene) as f:
        for line in f:
            fields = line.strip().split("\t")
            seq_id, s, e = fields[:3]
            protein_id, strand = fields[3], fields[5]
            s, e = int(s), int(e)
            strand = fields[5] 
            if s > e:
                print(*fields)
                continue
            if strand == "+":
                start = max(0,s-upstream)
                end = min(lengths[seq_id],s+downstream)
                if start > end:
                    start, end = end-1, start
                p = fwd_scores[seq_id][start:end].max()
            else:
                start = max(0,e-downstream)
                end =  min(lengths[seq_id],e+upstream)
                if start > end:
                    start, end = end-1, start
                p = rev_scores[seq_id][start:end].max()  
            p = clip(p) 
            logit = np.log(p/(1-p))
            logits[protein_id] = logit
            zscores[protein_id] = norm.ppf(gumbel_r.cdf(logit, loc=loc,scale=scale))


    logger.info("Saving results ...")
    fout = open(args.output,"w")
    zscores2 = {}

    print("protein id","logit","bg zscore","leader zscore",sep="\t",file=fout)
    loc, scale = gumbel_r.fit(list(logits.values())) 
    for protein_id in logits:
        logit, zscore =  logits[protein_id], zscores[protein_id]
        zscore2 = norm.ppf(gumbel_r.cdf(logit, loc=loc,scale=scale)) 
        logit, zscore = round(logit,4), round(zscore,4)
        zscore2 = round(zscore2,4)
        print(protein_id, logit, zscore, zscore2, sep="\t", file=fout) 
    fout.close()
    logger.info("All done.")
    
if __name__ == "__main__":
    main()
    
            
