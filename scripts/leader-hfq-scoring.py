#!/usr/bin/env python
import os
import torch
from torch.functional import F
import numpy as np
from model import CNNClassifier
import argparse
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('scanning sequence')
from tqdm import tqdm
from scipy.signal import convolve
from scipy.stats import gumbel_r,norm
from pyfaidx import Fasta
rc_lut = {"A":"T","C":"G","G":"C","T":"A"}

epsilon = 1e-20
def clip(p):
    if p > 1 - epsilon:
        p = 1 - epsilon
    if p < epsilon:
        p = epsilon
    return p

def get_rc(sequence):
    rcs = []
    for c in sequence:
        c = rc_lut.get(c,'N')
        rcs.append(c)
    return "".join(rcs[::-1])

def collate_fn(examples):
    sequences = []
    sequence_labels = []
    group_labels = []
    for sequence, sequence_label, group_label in examples:
        sequences.append(sequence)
        sequence_labels.append(sequence_label)
        group_labels.append(group_label)
    sequence_labels = torch.tensor(sequence_labels)
    group_labels = torch.tensor(group_labels)
    # batch , channel , length
    sequences = torch.cat(sequences)
    return sequences, sequence_labels, group_labels

DNA_alp = dict(zip("ACGT",list(range(4))))
protein_alp = dict(zip("ACDEFGHIKLMNPQRSTVWY",list(range(20))))

def onehot(sequence,alp="DNA"):
    assert alp in ["DNA","protein"], "the alphabet should be either DNA or protein"
    if alp == "DNA":
        c2i = DNA_alp
        x = torch.zeros((4,len(sequence)))
    else:
        c2i = protein_alp
        x = torch.zeros((20,len(sequence)))
    tokens = []
    indices = []
    for i,c in enumerate(sequence):
        if c in c2i:
            tokens.append(c2i[c])
            indices.append(i)
    x[tokens,indices] = 1
    return x[None,:,:]

def main():
    parser = argparse.ArgumentParser(description='subsequence scoring')
    parser.add_argument('--fasta','-f',type=str, required=True, help="Input sequence")
    parser.add_argument('--bed','-b',type=str, required=True, help="CDS to consider")
    parser.add_argument('--batch-size','-bs',type=int,default=1024,help="Batch size for scanning")
    parser.add_argument('--device','-d',default="cuda:0",choices=["cuda:0","cuda:1","cuda:2","cuda:3","cpu"],help="Device to run the model")
    parser.add_argument('--model','-m', default = "models/hfq.scoring.pt", type=str,help="Where to load the model parameters")
    parser.add_argument('--output','-o',required=True,type=str,help="Where to save the predicted probabilities")
    parser.add_argument('--n-channels', '-nc', type=int, default=64, help="number of channels to use")
    parser.add_argument('--n-blocks', '-nb', type=int, default=10, help="number of blocks in the model")
    parser.add_argument('--kernel-size', '-k', type=int, default=5, help="convolution kernel size in the res-block")
    parser.add_argument('--length', '-l', type=int, default=100, help="length of an instance")
    parser.add_argument('--stride', '-s', type=int, default = 1, help="stride size for scanning")
    parser.add_argument('--cutoff', '-c', type=float, default = 0, help="cutoff")
    parser.add_argument('--offset', type=int, default = 0, help="offset relative to predicted ends")
    parser.add_argument('--upstream', type=int, default = 200, help="upstream distance")
    parser.add_argument('--downstream', type=int, default = 100, help="downstream distance")
    args = parser.parse_args()

    logger.info("load genome ...")
    fasta = Fasta(args.fasta)
    logger.info(f"load model weights from {args.model} ...")
    state_dict = torch.load(args.model, map_location = args.device)
    if 'daOut.weight' in state_dict:
        n_domains = state_dict['daOut.weight'].shape[0]
    else:
        n_domains = None
    model = CNNClassifier(n_domains = n_domains, n_channels=args.n_channels, kernel_size = args.kernel_size, n_blocks=args.n_blocks)
    model = model.to(args.device)
    model = model.eval()
    model.load_state_dict(state_dict) 
    model.daOut = None

    batched_instances = []
    batched_positions = []

    scores_by_gene = {}

    logger.info("run prediction ...")
    n = 0
    fbed = open(args.bed)
    for line in fbed:
        fields = line[:-1].split("\t")
        seq_id, start, end = fields[:3]
        strand = fields[5]
        start, end = int(start), int(end)
        if strand == "+":
            s, e = start - args.upstream, start + args.downstream
        else:
            s, e = end - args.downstream, end + args.upstream
        if s < 0:
            s = 0
        if e > len(fasta[seq_id]):
            e = len(fasta[seq_id])
        gene_id = fields[3]
        sequence = fasta[seq_id][s:e]
        if strand == "-":
            sequence = sequence.reverse.complement
        sequence = str(sequence)
        if (n+1)%10000 == 0:
            logger.info(f"{round(n/1000)} K sequence processed")
        L = len(sequence)
        n += 1
        if L < args.length:
            sequence = sequence + "N"*(args.length-len(sequence))
            L = args.length
        sequence = sequence.upper().replace("U","T")
        sequence = onehot(sequence)
        p = 0
        while p+args.length <= L:
            instance = sequence[:,:,p:p+args.length]
            batched_instances.append(instance)
            batched_positions.append((seq_id, p, gene_id))
            p += args.stride
            if len(batched_instances) == args.batch_size:
                X = torch.cat(batched_instances).to(args.device)
                y = model(X)                
                scores = torch.softmax(y,1)[:,1].detach().cpu().numpy()
                for i, (seq_id, position, gene_id) in enumerate(batched_positions):
                    score = scores[i]
                    if gene_id not in scores_by_gene:
                        scores_by_gene[gene_id] = np.zeros(args.upstream+args.downstream-99)
                    scores_by_gene[gene_id][position] = score
                batched_instances = []
                batched_positions = []
    if len(batched_instances) > 0:
        X = torch.cat(batched_instances).to(args.device)
        y = model(X)
        scores = torch.softmax(y,1)[:,1].detach().cpu().numpy()        
        for i, (seq_id, position, gene_id) in enumerate(batched_positions):
            score = scores[i]
            if gene_id not in scores_by_gene:
                scores_by_gene[gene_id] = np.zeros(args.upstream+args.downstream-99)
            scores_by_gene[gene_id][position] = score
    
    logger.info("processing the predictions ...")
    logger.info(f"results will be saved to {args.output} ...")

    max_scores = {}
    for gene_id in scores_by_gene:
        #print(scores_by_gene[gene_id])
        score = clip(scores_by_gene[gene_id].max())
        logit = np.log(score/(1-score))   
        max_scores[gene_id] = logit 
    fout = open(args.output,"w")
    print("gene id","score","Z",sep="\t",file=fout)
    loc, scale = gumbel_r.fit(list(max_scores.values()))
    for gene_id in max_scores:
        logit = max_scores[gene_id]
        Z = norm.ppf(gumbel_r.cdf(logit, loc=loc,scale=scale))    
        print(gene_id, logit, Z, file=fout, sep="\t")    

    fout.close()
    logger.info("all done .")

                
if __name__ == "__main__":
    main()
