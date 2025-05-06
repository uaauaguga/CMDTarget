#!/usr/bin/env python
import subprocess
import logging
import pandas as pd
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


def make_nex(path, scores):
    ntax = len(scores)
    fout = open(path,"w")
    print("#NEXUS", file=fout)
    print("BEGIN DATA;",file=fout)
    print(f"DIMENSIONS  NCHAR=1 NTAX={ntax};", file=fout)
    print("FORMAT datatype=Continuous missing=? gap=-;",file=fout)
    print("",file=fout)
    print("    MATRIX",file=fout)
    for target_id in scores:
        genome_id = target_id.split(":")[0]
        score = scores[target_id]
        if score == "nan":
            score = '?'
        print(f"    {genome_id} {score}", file=fout)
    print("",file=fout)
    print(";",file=fout)
    print("END;",file=fout)
    fout.close()

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


def mcmc(query_id, pair_scores):
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
    nex = os.path.join(args.temp_directory, query_id + ".nex")
    make_nex(nex, pair_scores)
    nwk = os.path.join(args.temp_directory, query_id + ".nwk")
    with open(nwk,"w") as f:
        f.write(used_tree.write(format=1))
    mcmc_path = os.path.join(args.temp_directory, query_id + ".txt")
    names = os.path.join(args.temp_directory, query_id + ".names.txt")
    cmd = ["rb", "scripts/bayesian-reconstruction.Rev", "--args", nex, nwk, mcmc_path, names, 
           str(args.root_variance),str(args.signal_rate_median), 
           str(args.signal_rate_dispersion), str(args.noise_variance_median), str(args.noise_variance_dispersion), 
           str(args.mcmc_runs), str(args.mcmc_iterations)] 
    print(" ".join(cmd)) 
    subprocess.run(cmd,stdout = subprocess.DEVNULL, stderr =  subprocess.DEVNULL)
    if not os.path.exists(names):
        logger.info(f"error in processing {nex}: no name file produced, skip it")
        logger.info(" ".join(cmd))
        return [], None
    #names: the named tree

    try:
        mcmc_table = pd.read_csv(mcmc_path,sep="\t")
    except:
        logger.info(f"error in processing {nex}: reading {mcmc_path} failed, skip it")
        logger.info(" ".join(cmd))
        return [], None
    named_tree = Tree(newick=names,format=3)
    means = mcmc_table.median(axis=0)
    stds = mcmc_table.std(axis=0)    
    records = []
    #for i in range(1,n_genomes+1):
    named_tree_with_score = None
    for nnode in named_tree.traverse():
        node_name = nnode.name
        p = node_name.find("[")
        genome_id, node_index = node_name[:p],node_name[p:]
        node_index = node_index[1:-1].split("=")[-1]
        key = f"signals[{node_index}]"        
        processed_score = means.loc[key]
        processed_score_std = round(stds.loc[key],4)
        if len(genome_id) > 0:
            target_id = genome_id2target_id[genome_id]        
            raw_score = pair_scores[target_id]
            rootSignal = round(means.loc["rootSignal"],4)
            sigma2 = round(means.loc["signalSigma2"],4)
            records.append((query_id, target_id, round(raw_score,4), round(processed_score,4), processed_score_std, rootSignal, sigma2))
        nnode.add_features(mean=processed_score, std=processed_score_std)
        nnode.name = genome_id
        named_tree_with_score = named_tree.write(features=["mean","std"],format=2)
        
    if not args.keep_temp:
        os.remove(nex)
        os.remove(names)
        os.remove(nwk)
        os.remove(mcmc_path) 
        for i in range(1,11):
            path = os.path.join(args.temp_directory, query_id + f"_run_{i}.txt")
            if os.path.exists(path):
                os.remove(path)
    return records, named_tree_with_score

def main():
    parser = argparse.ArgumentParser(description='scoring sRNA targt interaction with comparative scoring')
    parser.add_argument('--input',  '-i', type=str, required=True, help='interaction scores')
    parser.add_argument('--output','-o', type=str , required=True, help="where to save output")
    parser.add_argument('--select','-s', type=str , help="which query to run")
    parser.add_argument('--root-variance', '-rv', type=float, default = 1,   help='root variance')
    parser.add_argument('--signal-rate-median',  '-srm', type=float, default = 0.5,   help='median of the lognormally distributed evolution rate')
    parser.add_argument('--signal-rate-dispersion', '-srd', type=float, default = 0.5, help='dispersion  of the lognormally distributed evolution rate')
    parser.add_argument('--noise-variance-median',  '-nvm', type=float, default = 0.5,   help='mean of the lognormally distributed noise variance dispersion')
    parser.add_argument('--noise-variance-dispersion',  '-nvd', type=float, default = 0.5,   help='dipersion of lognormally distributed noise variance')
    parser.add_argument('--temp-directory','-td', type=str , help="where to save temp files")
    parser.add_argument('--tree','-t', type=str , required=True, help="the phylogenetic tree")
    parser.add_argument('--jobs','-j', type=int , default=16, help="number of process to run")
    parser.add_argument('--keep-temp','-kt', action = "store_true" , help="whether keep temp file")
    parser.add_argument('--weights','-w',required=True, help="The weights for scoring")
    parser.add_argument('--mcmc-runs',default=4, help="Number of mcmc runs")
    parser.add_argument('--mcmc-iterations', default=10000, help="Number of mcmc iterations per-run")
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
    if not os.path.exists(args.temp_directory):
        os.mkdir(args.temp_directory)

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
        workers.append(pool.apply_async(func=mcmc,args=(query_id, pair_scores)))
    fout = open(args.output,"w")
    print("query id", "target id", "raw", "denoised", "std","root", "sigma2","tree",sep="\t",file=fout)
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
