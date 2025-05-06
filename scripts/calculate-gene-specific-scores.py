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
logger = logging.getLogger('calcualte gene scores')

def score_context(t):
    direction = t.loc["same"]
    distance = t.loc["distance"]
    if distance < -150:
        distance = -150
    if distance > 749:
        distance = 749
    if direction == 1:
        return context_scores.loc[distance,"upstream-conext-same-Z"]
    else:
        return context_scores.loc[distance,"upstream-conext-diff-Z"]

def main():
    parser = argparse.ArgumentParser(description='prepare gene specific scores')
    parser.add_argument('--hfq-score','-hs',type=str, required=True, help="score of hfq binding in bed format")
    parser.add_argument('--context-score','-cs',type=str, default="params/upstream-context-scores.txt", 
                        help="scoring rule of upstream context")
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
            
    logger.info("Load gene information ...")
    entries = defaultdict(list)
    with open(args.gene) as f:
        for line in f:
            fields = line.strip().split("\t")
            seq_id, start, end = fields[:3]
            name, strand = fields[3], fields[5]
            start, end = int(start), int(end)
            strand = fields[5]
            entries[seq_id].append((start,end,name,strand))
    records = []
    logger.info("Prepare gene table ...")
    for seq_id in entries:
        entries[seq_id] = sorted(entries[seq_id],key=lambda x:(x[0],x[1]))
        for i in range(1,len(entries[seq_id])-1):
            record = entries[seq_id][i]
            start, end = record[:2]
            length = end - start
            protein_id = record[2]
            strand = record[3]
            if strand == "+":             
                upstream_record = entries[seq_id][i-1]
                upstream_end = upstream_record[1]
                upstream_strand = upstream_record[3]
                same = int(strand == upstream_strand)
                distance = start - upstream_end
            else:
                upstream_record = entries[seq_id][i+1]
                upstream_start = upstream_record[0]
                upstream_strand = upstream_record[3]            
                same = int(strand == upstream_strand)
                distance = upstream_start - end
            records.append((seq_id, start, end, strand, protein_id, same, distance,length))
    gene_table = pd.DataFrame.from_records(records)
    gene_table.columns = ["seq id","start","end","strand","protein id","same","distance","length"]      
    
    logger.info("Load context scoring rules ...")
    global context_scores
    context_scores = pd.read_csv(args.context_score,sep="\t",index_col=0)
    gene_table["context Z score"] = gene_table.apply(score_context,axis=1)
    
    
    logger.info("Load hfq scores ...")
    fwd_scores = {}
    rev_scores = {}
    with open(args.hfq_score) as f:
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
        start = np.random.randint(lengths[seq_id]-300)
        end = start + 300
        if np.random.rand() > 0.5:
            p = fwd_scores[seq_id][start:end].max()
        else:
            p = rev_scores[seq_id][start:end].max()        
        score = np.log(p/(1-p))
        background_hfq_scores.append(score)
    loc, scale = gumbel_r.fit(background_hfq_scores)      
    logger.info(f"Backgound distribution parameter of the logits is loc={loc}, scale={scale}.")
    
    zscores = {}
    probabilities = {}
    logger.info("Calculate Z score of each leader sequences ...")
    for record in records:
        seq_id, s, e, strand, protein_id = record[:5]
        if strand == "+":
            start = max(0,s-200)
            end = min(lengths[seq_id],e+100)
            p = fwd_scores[seq_id][start:end].max()
        else:
            start = max(0,e-100)
            end =  min(lengths[seq_id],e+200)
            p = rev_scores[seq_id][start:end].max()     
        score = np.log(p/(1-p))
        zscore = norm.ppf(gumbel_r.cdf(score, loc=loc,scale=scale))
        zscores[protein_id] = zscore
        probabilities[protein_id] = p
    gene_table["hfq Z score"] = gene_table["protein id"].map(lambda x:zscores[x])  
    gene_table["hfq proba"] = gene_table["protein id"].map(lambda x:probabilities[x])
    
    logger.info("Saving results ...")
    
    gene_table.to_csv(args.output,sep="\t",index=False)
    
    logger.info("All done.")
    
if __name__ == "__main__":
    main()
    
            
