#!/usr/bin/env python
import argparse
import gzip
import os
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('concatenate sequences')

def main():
    parser = argparse.ArgumentParser(description='concatenate fasta sequence')
    parser.add_argument('--input-directory', '-i', required=True, help="input directory")
    parser.add_argument('--output','-o', required=True, help="output file")
    parser.add_argument('--genome-ids','-gi',help="target genome ids to consider")
    args = parser.parse_args() 
         
    fastas = [ fa for fa in os.listdir(args.input_directory) if fa.endswith(".fa") ]
    logger.info(f"{len(fastas)} genomes presented in input directory")

    tgenome_ids = set()
    if args.genome_ids is not None:
        tgenome_ids = set(open(args.genome_ids).read().strip().split("\n"))

    logger.info(f"concatenate input sequences ...")
    fout = open(args.output,"w")
    
    n = 0
    for fasta in sorted(fastas):
        if not fasta.endswith(".fa"):
            continue
        path = os.path.join(args.input_directory, fasta)
        if path.endswith(".gz"):
            f = gzip.open(path)
            prefix = ".".join(fasta.split(".")[:-2])
        else:
            f = open(path)        
            prefix = ".".join(fasta.split(".")[:-1])
        if prefix not in tgenome_ids:
            continue
        n += 1
        for line in f:
            if path.endswith(".gz"):
                line = line.decode()
            if line.startswith(">"):
                line = ">" + prefix + ":" + line[1:]
            fout.write(line)
        f.close()
    logger.info(f"{n} genomes finally reserved.")
    fout.close()


if __name__ == "__main__":
    main()
