#!/usr/bin/env python
import argparse
import re
import sys
import os
import subprocess
from pyfaidx import Fasta
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('extract marker sequence')
import re

def main():
    parser = argparse.ArgumentParser(description='Extract marker sequence')
    parser.add_argument('--genome','-g', type=str , required=True, help="Input genome sequence")
    parser.add_argument('--protein','-p', type=str , required=True, help="Input protein sequence")
    parser.add_argument('--bed','-b', type=str , required=True, help="Input CDS coordinate")
    parser.add_argument('--hmm',type=str , default = "models/rpoB.hmm",help="hmm model to use")  
    parser.add_argument('--output','-o', type=str, required=True , help="Output marker gene nucleotide sequences")
    args = parser.parse_args()

    
    logger.info("Search homolog protein sequences ...")
    tbl = args.output + ".tmp"
    cmd = ["hmmsearch","--tblout",tbl,"--noali", args.hmm, args.protein] 
    subprocess.run(cmd)
    
    hits = []
    with open(tbl) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = re.split(r"\s+",line[:-1])
            gene_id = fields[0] 
            hits.append(gene_id)
    
    if len(hits) > 0:
        logger.info("marker gene detected.")
        if len(hits) > 1:
            logger.info(str(len(hits)) + " marker gene detected, only use the best one.")
        gene_id = hits[0]
        logger.info("Extract gene coordinate ...")
        with open(args.bed) as f:
            for line in f:
                fields = line[:-1].split("\t")
                if fields[3] == gene_id:
                    seq_id, start, end, name, _, strand = fields[:6]
                    start, end = int(start), int(end)
                    break
        logger.info("Extract nucletide sequences ...")
        fasta = Fasta(args.genome)
        sequence = fasta[seq_id][start:end]
        if strand == "-":
            sequence = sequence.reverse.complement
        sequence = str(sequence)        
        loci_id = f"{seq_id}:{start}-{end}({strand})"
        with open(args.output,"w") as f:
            print(f">{loci_id}",file=f)
            print(sequence,file=f)
    else:
        logger.info("No marker gene detected.")


if __name__ == "__main__":
    main()
