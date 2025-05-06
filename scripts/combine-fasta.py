#!/usr/bin/env python
import argparse
import gzip
import os
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('concatenate sequences')

def main():
    parser = argparse.ArgumentParser(description='chunkify fasta files ')
    parser.add_argument('--input-directory', '-i', required=True, help="input directory")
    parser.add_argument('--output','-o', required=True, help="output file")
    args = parser.parse_args() 
         
    fastas = [ fa for fa in os.listdir(args.input_directory) if not fa.endswith(".fai") ]
    logger.info(f"{len(fastas)} genomes presented in input directory")

    logger.info(f"concatenate input sequences ...")
    fout = open(args.output,"w")
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
        for line in f:
            if path.endswith(".gz"):
                line = line.decode()
            if line.startswith(">"):
                line = ">" + prefix + ":" + line[1:]
            fout.write(line)
        f.close()
    fout.close()


if __name__ == "__main__":
    main()
