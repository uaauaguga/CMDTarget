#!/usr/bin/env python
import subprocess
import logging
import pandas as pd
from scipy.optimize import minimize
import json
import argparse
import numpy as np
import io
import os
import re
from collections import defaultdict
import pandas as pd
from ete3 import Tree
import subprocess
from multiprocessing import Pool
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('comparative scoring')

def prune(tree, names):
    nodes = []
    for node in tree.traverse():
        name = node.name
        if name in names:
            nodes.append(node)
    tree.prune(nodes)
    return tree

def scale_tree(mtree,height=None):
    root = mtree.get_tree_root()
    scaler_lut = {}
    leave_scalers = []
    for node in mtree.traverse(strategy="postorder"):
        if node.is_root():
            continue
        if node.is_leaf():
            # Set the distance for each leaf to the max depth
            scaler_lut[node] = node.get_distance(root)
            leave_scalers.append(scaler_lut[node])
        else:
            scalers = []
            for c in node.get_children():
                scalers.append(scaler_lut[c])
            scaler_lut[node] = np.exp(np.log(scalers).mean())#sum(scalers)/len(scalers)
    if height is None:
        height = np.median(leave_scalers)
    for node in mtree.traverse(strategy="postorder"):
        if node in scaler_lut:
            #print(height/scaler_lut[node])
            node.dist = height*node.dist/scaler_lut[node]

def check_height(mtree):
    heights = []
    root = mtree
    for node in mtree.traverse():
        if node.is_leaf():
            heights.append(node.get_distance(root))
    return np.median(heights),np.std(heights)


def init_scores(tree,scores,genome_id2target_id):
    tree_used = tree.copy()
    xs = []
    mask = []
    for clade in tree_used.traverse(strategy='postorder'):
        if clade.is_leaf():
            value = scores[genome_id2target_id[clade.name]]                        
            clade.score = value
            xs.append(clade.score)
            mask.append(True)
        else:
            lchild, rchild = clade.get_children()
            lbl = lchild.dist
            rbl = rchild.dist
            clade.score = ((lchild.score/lbl) + (rchild.score/rbl))/(1/lbl + 1/rbl)
            xs.append(clade.score)
            mask.append(False)
    return np.array(xs), np.array(mask)


def nll(xs,tree,scores,genome_id2target_id):
    values = []
    tree_used = tree.copy()
    i = 0
    cost = 0
    noise = args.noise_variance
    for clade in tree_used.traverse(strategy='postorder'):
        clade.score = xs[i]
        if clade.is_leaf():
            value = scores[genome_id2target_id[clade.name]]            
            cost = cost + (value - clade.score)**2/(2*noise)
        else:
            lchild, rchild = clade.get_children()
            lbl = lchild.dist
            rbl = rchild.dist
            lsignal = args.signal_rate*lbl
            rsignal = args.signal_rate*rbl
            cost = cost +  (clade.score - lchild.score)**2/(2*lsignal)            
            cost = cost +  (clade.score - rchild.score)**2/(2*rsignal)
            if clade.is_root():
                cost  = cost + clade.score**2/(2*args.root_variance)
        i += 1
    return cost


def denoise(query_id, pair_scores):
    logger.info(f"process {query_id} ...")
    used_tree = tree.copy()
    genome_ids_in_tree = [node.name for node in tree.traverse()]
    target_ids = pair_scores.keys()    
    genome_id2target_id = {}
    for target_id in target_ids:
        genome_id, protein_id = target_id.split(":")
        if genome_id not in genome_ids_in_tree:
            continue
        genome_id2target_id[genome_id] = target_id
    genome_ids = list(genome_id2target_id.keys())
    n_genomes = len(genome_ids)    
    used_tree = prune(used_tree, genome_ids)
    if args.reroot:
        # get midpoint as outgroup
        mroot = used_tree.get_midpoint_outgroup()
        # and set it as tree outgroup
        used_tree.set_outgroup(mroot)
    if args.scale_tree:
        height, std = check_height(used_tree)
        while std > 0.001:
            scale_tree(used_tree,height)
            _, std = check_height(used_tree)
    xs_init, mask = init_scores(used_tree,pair_scores,genome_id2target_id)
    result = minimize(lambda x:nll(x,used_tree,pair_scores,genome_id2target_id), xs_init, method='L-BFGS-B')
    i = 0
    records = []
    for clade in used_tree.traverse(strategy='postorder'):
        #clade.score = result.x[i]
        clade.add_features(score=result.x[i])
        i += 1
        if clade.is_leaf():
            genome_id = clade.name
            target_id = genome_id2target_id[genome_id]        
            raw_score = pair_scores[target_id]
            records.append((query_id, target_id, round(raw_score,4), round(result.x[i],4)))
    named_tree_with_score = used_tree.write(features=["score"],format=2)
    return records, named_tree_with_score

def main():
    parser = argparse.ArgumentParser(description='scoring sRNA targt interaction with comparative scoring')
    parser.add_argument('--input',  '-i', type=str, required=True, help='interaction scores')
    parser.add_argument('--output','-o', type=str , required=True, help="where to save output")
    parser.add_argument('--select','-s', type=str , help="which query to run")
    parser.add_argument('--root-variance', '-rv', type=float, default = 1,   help='root variance')
    parser.add_argument('--signal-rate',  '-sr', type=float, default = 0.08,   help='evolution rate')
    parser.add_argument('--noise-variance',  '-nv', type=float, default = 0.004,   help='noise variance')
    parser.add_argument('--tree','-t', type=str , required=True, help="the phylogenetic tree")
    parser.add_argument('--jobs','-j', type=int , default=16, help="number of process to run")
    parser.add_argument('--weights','-w',required=True, help="The weights for scoring")
    parser.add_argument('--normalize','-n', action = "store_true" , help="whether normalize the score")
    parser.add_argument('--reroot',action = "store_true" , help="whether perform midpoint rooting")
    parser.add_argument('--scale-tree',action = "store_true" , help="whether rescaling the branch length")
    global args
    args = parser.parse_args()


    weights = json.load(open(args.weights))
    scores = defaultdict(dict)
    logger.info("Load scores ...")
    table = pd.read_csv(args.input,sep="\t")
    
    table["score"] = 0
    for key in weights:
        if weights[key] == 0:
            logger.info(f"{key} has zero weight, skip it.")
            continue
        assert key in table.columns, f"{key} not present in input table"        
        #if key not in table.columns:
        #    logger.info(f"{key} not in input table. skip it.") 
        #    continue
        values = table[key].fillna(0)
        # clip outliers
        values[values>5] = 5
        values[values<-5] = -5       
        table["score"] += values*weights[key]
    if args.normalize:
        mean = table["score"].mean()
        std = table["score"].std()
        table["score"] = (table["score"] - mean)/std    
        logger.info(f"Score mean is {mean}, stderr is {std}.")
    for query_id, target_id, score in table.loc[:,["query id", "protein id", "score"]].to_records(index=False):
        scores[query_id][target_id] = score

    global tree
    logger.info("Load tree ...")
    tree = Tree(newick=args.tree,format=1,quoted_node_names=True)

    logger.info("Comparative scoring ...")

    pool = Pool(args.jobs)
    workers = []
    for query_id in scores:
        if args.select is not None:
            if query_id != args.select:
                continue
        pair_scores = scores[query_id]
        n_genomes = len(pair_scores)
        if n_genomes < 3:
            continue
        workers.append(pool.apply_async(func=denoise,args=(query_id, pair_scores)))
    fout = open(args.output,"w")
    print("query id", "target id", "raw", "denoised","tree",sep="\t",file=fout)
    for worker in workers:
        records, tree = worker.get()  
        for record in records:
            print(*record, tree, sep="\t", file=fout)
            tree = "." 
        fout.flush()
    fout.close()
    logger.info("All done.")
                

if __name__ == "__main__":
    main()
