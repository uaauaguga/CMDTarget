#!/usr/bin/env python
from ushuffle import shuffle
from scipy.stats import gumbel_r
import argparse
from multiprocessing import Pool
import subprocess
import io
from collections import defaultdict
import pickle
import numpy as np
import argparse
import resource
soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
resource.setrlimit(resource.RLIMIT_NOFILE, (hard, hard))
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

def prediction(sequence_1, sequence_2, number=1, seed=7):
    cmd = ["/BioII/lulab_b/jinyunfan/miniforge3/envs/IntaRNA-env/bin/IntaRNA","-q",sequence_1,"-t",sequence_2,"--outMode","C","--outNumber",str(number),"--tAcc","C","--seedBP", str(seed)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr = subprocess.DEVNULL)
    f = io.TextIOWrapper(proc.stdout, encoding="unicode_escape")
    _ = next(f)
    records = []
    for line in f:
        if len(line.strip()) == 0:
            continue
        fields = line.strip().split(";")
        if len(fields) < 6:
            #logger.info("Some thing wrong with " + " ".join(cmd))
            continue
        t, ts, te, q, qs, qe = fields[:6]
        sequence, bp, energy = fields[-3:]
        energy = float(energy)
        ts, te, qs, qe = int(ts), int(te), int(qs), int(qe)
        records.append((qs, qe,len(sequence_1), ts, te, len(sequence_2), sequence, bp, energy))
    code =  proc.poll()
    return records

def inference5(X, params):
    X = X@params["linear_1.weight"].T + params["linear_1.bias"]
    X = np.maximum(X,0)
    X = X@params["linear_2.weight"].T + params["linear_2.bias"]
    X = np.maximum(X,0)
    X = X@params["linear_3.weight"].T + params["linear_3.bias"]
    X = np.maximum(X,0)
    X = X@params["linear_4.weight"].T + params["linear_4.bias"]
    X = np.maximum(X,0)
    return X@params["linear_5.weight"].T + params["linear_5.bias"]



def get_shuffled_scores(s1,s2,n=1000):    
    pool = Pool(8)
    workers = []
    for i in range(n):
        sequence_1 = shuffle(s1.encode(),2).decode()
        sequence_2 = shuffle(s2.encode(),2).decode()
        workers.append(pool.apply_async(func=prediction, args=(sequence_1, sequence_2)))
    scores = [] 
    for worker in workers:
        rs = worker.get()
        if len(rs)  == 0:
            energy = 0
        else:
            qs, qe, l1, ts, te, l2, sequence, bp, energy = rs[0]
        score = -energy/10
        scores.append(score)
    pool.close()
    pool.join()
    return scores #loc, scale

def get_predicted_scores(s1,s2):
    frequency = count_frequency(s1) + count_frequency(s2)
    frequency = np.array(frequency).reshape(1,-1)
    pred = inference5(frequency ,params5)
    locp, scalep = pred[0,0], pred[0,1]
    #meanp = np.euler_gamma*scalep + locp
    return locp, scalep


def main():
    parser = argparse.ArgumentParser(description='check performance of background model')
    parser.add_argument('--model', '-m', type=str, default = "models/20240404.model.pkl", help='model to use')
    parser.add_argument('--output', '-o', type=str, required = True, help='where to save output')
    args = parser.parse_args()
    fout = open(args.output,"w")
    global params5
    params5 = pickle.load(open(args.model,"rb"))

    print("approach","loc","scale","mean", "#", sep="\t", file=fout)
    for i in range(1000):
        L1 = np.random.randint(50,150)
        L2 = np.random.randint(250,350)
        s1 = "".join(["ACGT"[np.random.randint(4)] for i in range(L1)])
        s2 = "".join(["ACGT"[np.random.randint(4)] for i in range(L2)])
        #loc, scale = gumbel_r.fit(scores)
        for n in [2,4,6,8,10,20,30,40,50,60,100]:
            scores = get_shuffled_scores(s1,s2,n=n)
            np.random.shuffle(scores)            
            loc, scale = gumbel_r.fit(scores[:n])
            mean = np.euler_gamma*scale + loc
            print(i, f"shuffle-{i}",loc, scale, mean, n, sep="\t", file=fout)
        loc, scale = get_predicted_scores(s1,s2)
        mean = np.euler_gamma*scale + loc
        print(i, f"shuffle-{i}",loc, scale, mean, 0,sep="\t", file=fout)
        fout.flush()         
    fout.close()            


if __name__ == "__main__":
    main()
