#!/usr/bin/env python
import subprocess
import logging
import argparse
import numpy as np
import re
import os
import io
from collections import defaultdict
from multiprocessing import Pool
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('target prediction')
from itertools import product
np.random.seed(666)


def prediction_intarna(sequence_1, sequence_2, seed=7):
    cmd = ["/BioII/lulab_b/jinyunfan/miniforge3/envs/IntaRNA-env/bin/IntaRNA","-q",sequence_1,"-t",sequence_2,"--outMode","C","--outNumber","1","--tAcc","C","--seedBP", str(seed)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr = subprocess.DEVNULL)
    f = io.TextIOWrapper(proc.stdout, encoding="unicode_escape")
    _ = next(f)
    records = []
    for line in f:
        if len(line.strip()) == 0:
            continue
        fields = line.strip().split(";")
        if len(fields) < 6:
            logger.info("Some thing wrong with " + " ".join(cmd))
            continue
        t, ts, te, q, qs, qe = fields[:6]
        energy = fields[-1]
        energy = float(energy)
        ts, te, qs, qe = int(ts), int(te), int(qs), int(qe)
        records.append((qs, qe,len(sequence_1), ts, te, len(sequence_2), energy))
    code =  proc.poll()     
    return records

def prediction_rnaup(sequence_1, sequence_2, seed=None):
    sequence = sequence_1 + "&" + sequence_2
    cmd = ["RNAup","--no_output_file","--interaction_pairwise","--include_both"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.DEVNULL)
    lines = proc.communicate(sequence.encode())[0].decode()
    line = lines.split("\n")[0]
    #((((((&))))))   6,11  :   6,11  (-8.20 = -8.23 + 0.04)
    p = line.rfind("(")
    energy = float(line[p+1:].split("=")[0].strip())
    fields = re.split(r"\s+",line[:p-1])
    iv_1, iv_2 = fields[1], fields[3]
    if len(sequence_1) > len(sequence_2):
        s1, e1 = iv_1.split(",")
        s2, e2 = iv_2.split(",")
    else:
        s1, e1  = iv_2.split(",")
        s2, e2 = iv_1.split(",")
    code =  proc.poll()
    return [(int(s1),int(e1),len(sequence_1),int(s2),int(e2),len(sequence_2),energy)]



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


def main():
    parser = argparse.ArgumentParser(description='predict RNA RNA interaction with IntaRNA')
    parser.add_argument('--srnas',  '-rs', type=str, required=True, help='sRNA sequences')
    parser.add_argument('--targets',  '-ts', type=str, required=True, help='Target sequences')
    parser.add_argument('--output','-o', type=str , required=True, help="where to save output")
    parser.add_argument('--method', '-m', type=str , default="IntaRNA", choices = ["IntaRNA","RNAup","RNAplex"], help="method for interaction enery calculation")
    parser.add_argument('--jobs','-j', type=int , default=128, help="number of process to run, not work for RNAplex")
    parser.add_argument('--seed','-s', type=int , default=7, help="seed length to use") 
    args = parser.parse_args()
       
    if args.method == "RNAplex":
        logger.info("Use RNA plex for energy scoring .")         
        lines0 = open(args.srnas).read() + open(args.targets).read()
        lines = ""
        for line in lines0.strip().split("\n"):
            if line.startswith(">"):
                line = ">" + line[line.find(":")+1:]
            lines += line
        wd = os.getcwd()
        profile_directory =  args.output + ".profile" 
        profile_directory = os.path.abspath(profile_directory) 
        if not os.path.exists(profile_directory):
            os.mkdir(profile_directory)
        os.chdir(profile_directory)              
        print(profile_directory)
        logger.info("Calculate accessibility with RNAplfold ...")
        #cmd = ["RNAplfold", "-W", "240", "-L", "160", "-u", "30", "--opening_energies"]        
        #proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.DEVNULL)
        #ret = proc.communicate(lines.encode())[0].decode()
        #proc.poll()
        os.chdir(wd)
        logger.info("Calculate interaction energy with RNAplex ...")
        cmd = ["RNAplex", "-l", "25","-q", args.srnas, "-t", args.targets,"-a", profile_directory] 
        tmp = args.output + ".tmp"
        ftmp = open(tmp,"w")
        subprocess.run(cmd, stdout=ftmp)
        ftmp.close()
        fout = open(args.output,"w") 
        print("sRNA id","sstart","send","target id","tstart","tend","energy",file=fout,sep="\t")    
        with open(tmp) as f:
            #((((((&)))))) 257,262 : 143,148 (-8.98 = -13.41 +  1.71 +  2.72) i:262,j:143 <-9.83>
            for line in f:
                target_id = line[1:-1]
                sRNA_id = next(f)[1:-1] 
                res = next(f)[:-1]
                fields = re.split(r"\s+",res)[1:]
                tstart, tend = fields[0].split(",")
                qstart, qend = fields[2].split(",")
                energy = float(fields[3][1:])
                print(sRNA_id, int(qstart)-1, int(qend), target_id, int(tstart)-1, int(tend), energy, sep="\t", file=fout)
        fout.close()
        os.remove(tmp)
    else: 
        logger.info(f"Use {args.method} for scoring ...")
        fout = open(args.output,"w")
        logger.info("load sequences ...")
        sRNAs = load_fasta(args.srnas)
        targets = load_fasta(args.targets) 
        logger.info(f"run {args.method} with {args.jobs} workers ...")
        pool = Pool(args.jobs)
        n_too_long = 0
        n_no_prediction = 0
        n_total = 0
        n_too_short = 0
        n_few_reads = 0
        if args.method == "IntaRNA":
            prediction = prediction_intarna 
        else:
            prediction = prediction_rnaup
        print("sRNA id","sstart","send","target id","tstart","tend","energy",file=fout,sep="\t")    
        for sRNA_id in sRNAs:
            workers = []
            for target_id in list(targets.keys()):
                sequence_1 = sRNAs[sRNA_id]
                sequence_2 = targets[target_id]
                workers.append((pool.apply_async(func=prediction, args=(sequence_1, sequence_2, args.seed)),sRNA_id, target_id))
            i = 0
            logger.info(f"{sRNA_id}: {len(workers)} interactions to process .")
            n_no_prediction = 0
            for worker, seq_id_1, seq_id_2 in workers:
                rs = worker.get()
                i += 1
                for r in rs:
                    qs, qe, lq, ts, te, lt, energy = r
                    print(seq_id_1, qs-1, qe, seq_id_2, ts-1, te, energy, sep="\t", file=fout)
                if len(rs) == 0:
                    n_no_prediction += 1
                    print(seq_id_1,"-1","-1", seq_id_2, "-1","-1", 0, sep="\t",file=fout)            
                fout.flush()
        logger.info(f"{n_no_prediction} have no prediction.")
        fout.close()
        logger.info("all done .")

if __name__ == "__main__":
    main()
