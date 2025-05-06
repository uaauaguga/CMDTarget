#!/usr/bin/env python
import logging
import argparse
from scipy.stats import gumbel_r
import numpy as np
import pickle
from collections import defaultdict
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('background energy')
np.random.seed(666)

def count_frequency(sequence,k=2):
    counter = defaultdict(int)
    for i in range(len(sequence)-k):
        counter[sequence[i:i+k]] += 1
    frequency = np.zeros(4**k)
    idxlut = dict(zip(list("ACGT"),list(range(4))))
    for kmer in counter.keys():
        idx = 0
        skip = False
        for c in kmer:
            if c not in idxlut:
                skip = True
                break
            idx = idx*4 + idxlut[c]
        if not skip:
            frequency[idx] = counter[kmer]/10
    return list(frequency)

def load_fasta(path):
    sequences = {}
    with open(path) as f:
       for line in f:
           if line.startswith(">"):
               seq_id = line[1:].strip().split(" ")[0]
               sequences[seq_id] = ""
           else:
               sequences[seq_id] += line.strip()
    return sequences


def inference(X, params):
    X = X@params["linear_1.weight"].T + params["linear_1.bias"]
    X = np.maximum(X,0)
    X = X@params["linear_2.weight"].T + params["linear_2.bias"]
    X = np.maximum(X,0)
    X = X@params["linear_3.weight"].T + params["linear_3.bias"]
    X = np.maximum(X,0)
    X = X@params["linear_4.weight"].T + params["linear_4.bias"]
    X = np.maximum(X,0)
    return X@params["linear_5.weight"].T + params["linear_5.bias"]


def main():
    parser = argparse.ArgumentParser(description='predict RNA RNA interaction with intaRNA')
    parser.add_argument('--srnas',  '-rs', type=str, required=True, help='sRNA sequences')
    parser.add_argument('--targets',  '-ts', type=str, required=True, help='Target sequences')
    parser.add_argument('--output','-o', type=str , help="where to save output")
    parser.add_argument('--model', '-m', type=str, default = "models/20240404.model.pkl", help='model to use')
    parser.add_argument('--word-size', '-k', type=int, default = 2, help='word size to use')
    args = parser.parse_args()
    
    fout = open(args.output,"w")
    logger.info("load sequences ...")
    sRNAs = load_fasta(args.srnas)
    targets = load_fasta(args.targets) 
    logger.info("load weights for background distribution modeling ...")
    params = pickle.load(open(args.model,"rb"))


    print("sRNA id","target id","loc","scale",file=fout,sep="\t")    
    for sRNA_id in sRNAs:
        logger.info(f"Processing {sRNA_id} ...")
        frequencies = []
        target_ids = []
        for target_id in list(targets.keys()):
            sequence_1 = sRNAs[sRNA_id]
            sequence_2 = targets[target_id]
            frequency = count_frequency(sequence_1,args.word_size) + count_frequency(sequence_2,args.word_size) 
            frequencies.append(np.array(frequency))
            target_ids.append(target_id)
        X = np.array(frequencies)
        X = inference(X, params) 
        for i in range(X.shape[0]):
            loc, scale = X[i,0], X[i,1]
            print(sRNA_id, target_ids[i], loc, scale, sep="\t", file=fout)
    fout.close()
    logger.info("All done .")

if __name__ == "__main__":
    main()
