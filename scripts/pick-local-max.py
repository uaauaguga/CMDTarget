#!/usr/bin/env python
import argparse
import numpy as np
import logging
from collections import defaultdict
import subprocess
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("select intervals")



def select_intervals(cached_scores, cached_ivs, cached_names, cached_attrs):
    cached_scores_by_strand = {"+":[],"-":[]}
    cached_iv_by_strand = {"+":[],"-":[]}
    cached_names_by_strand = {"+":[],"-":[]}
    cached_attrs_by_strand = {"+":[],"-":[]}
    for i in range(len(cached_scores)):
        seq_id, start, end, strand = cached_ivs[i]
        cached_scores_by_strand[strand].append(cached_scores[i])
        cached_iv_by_strand[strand].append((seq_id, start, end))
        cached_names_by_strand[strand].append(cached_names[i])
        cached_attrs_by_strand[strand].append(cached_attrs[i])
    records = []
    for strand in "+-":
        # current behavior is different from bedtools merge:
        # first merge intervals not considering strandness, then pick best interval from each strand
        # TODO: add an option to reproduce bedtools results
        if len(cached_scores_by_strand[strand]) == 0:
            continue
        i = np.argmax(cached_scores_by_strand[strand])
        seq_id, start, end = cached_iv_by_strand[strand][i]
        score = cached_scores_by_strand[strand][i]
        name = cached_names_by_strand[strand][i]
        attrs = cached_attrs_by_strand[strand][i]
        record = [seq_id, start, end, name, score, strand] + attrs
        records.append(record)
    return sorted(records,key=lambda x:x[1])


def main():
    parser = argparse.ArgumentParser(description='pick best interval from overlapped ones')
    parser.add_argument('--input', '-i',type=str,required=True,help="input intervals")
    parser.add_argument('--output','-o',type=str,required=True,help="output intervals")
    parser.add_argument('--window', '-w', type=int,  help="slop each interval up to a given window size brefore merging")
    args = parser.parse_args()

    logger.info("sort input intervals inplace ...")
    cmd = ["sort","-k1,1","-k2,2n","-o",args.input,args.input]
    subprocess.run(cmd)

    logger.info(f"load intervals from {args.input} ...")
    logger.info(f"picked intervals will be saved to {args.output} .")
    fin = open(args.input)
    seq_id = ""
    last_seq_id = ""
    last_used_start, last_used_end = -1, -1
    fout = open(args.output,"w")
    # overlapped intervals to select from
    cached_ivs = []
    cached_scores = []
    cached_attrs = []
    cached_names = []
    for line in fin:
        fields = line.strip().split("\t")
        seq_id, start, end, name, score, strand = fields[:6]
        start, end = int(start), int(end)
        if args.window is not None:
            center = int((start + end)/2)
            used_start, used_end = max(center - int(args.window/2),0), center + int(args.window/2)
            if used_start == used_end:
                used_end += 1
        else:
            used_start, used_end = int(start), int(end)
        score = float(score)
        attrs = fields[6:] if len(fields) > 6 else []
        if ((used_start > last_used_end) and (last_seq_id == seq_id)) or (last_seq_id != seq_id):
            if last_seq_id != seq_id:
                last_used_start, last_used_end = -1, -1
            # current entry does not overlap with previous one
            # save cached entries
            if len(cached_ivs) > 0:
                # group predictions by strand
                records = select_intervals(cached_scores, cached_ivs, cached_names, cached_attrs)
                for record in records:
                    print(*record, file=fout, sep="\t")
                # update the cache
                cached_ivs, cached_scores, cached_attrs, cached_names = [], [], [], []
        cached_ivs.append((seq_id, start, end, strand))
        cached_scores.append(score)
        cached_attrs.append(attrs)
        cached_names.append(name)
        last_seq_id = seq_id
        last_used_start = used_start
        last_used_end = max(last_used_end,used_end)

    if len(cached_ivs) > 0:
        records = select_intervals(cached_scores, cached_ivs, cached_names, cached_attrs)
        for record in records:
            print(*record, file=fout, sep="\t")

    fin.close()
    fout.close()

if __name__ == "__main__":
    main()

