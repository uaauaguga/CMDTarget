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


epsilon = 1e-10
def clip(p):
    if p > 1-epsilon:
        p = 1-epsilon
    if p < epsilon:
        p = epsilon
    return p

def main():
    parser = argparse.ArgumentParser(description='prepare gene specific scores')
    parser.add_argument('--input','-i',type=str, required=True, help="score of hfq binding in bed format")
    parser.add_argument('--gene','-g',type=str, required=True, help='gene interval in bed format, strandness is required')    
    parser.add_argument('--output','-o',type=str, required=True, help="where to save the scores")
    parser.add_argument('--contig', '-c', type=str, required=True, help="contig length")
    args = parser.parse_args()
    
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
        sampled_seq_ids = np.random.choice(seq_ids,10000,p=probs)   
        sampled_seq_ids = list(sampled_seq_ids)
    else:
        sampled_seq_ids = [list(lengths.keys())[0]]*10000
    
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
                start = max(0,s-200)
                end = min(lengths[seq_id],s+100)
                if start > end:
                    start, end = end-1, start
                p = fwd_scores[seq_id][start:end].max()
            else:
                start = max(0,e-100)
                end =  min(lengths[seq_id],e+200)
                if start > end:
                    start, end = end-1, start
                p = rev_scores[seq_id][start:end].max()  
            p = clip(p) 
            logit = np.log(p/(1-p))
            logits[protein_id] = logit
    loc, scale = gumbel_r.fit(list(logits.values()))
    print(loc, scale) 


    logger.info("Saving results ...")
    fout = open(args.output,"w")

    print("protein id","logit","Z",sep="\t",file=fout)
    for protein_id in logits:
        logit =  logits[protein_id]
        Z = norm.ppf(gumbel_r.cdf(logit, loc=loc,scale=scale))
        print(protein_id, logit, Z, sep="\t", file=fout)

    idx = 0
    for logit in background_hfq_scores:       
        Z = norm.ppf(gumbel_r.cdf(logit, loc=loc,scale=scale))
        print("BG" + str(idx).zfill(5), logit, Z, sep="\t", file=fout)
        idx += 1
    
    fout.close()
    logger.info("All done.")
    
if __name__ == "__main__":
    main()
