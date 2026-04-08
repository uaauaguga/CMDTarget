#!/usr/bin/env python
import argparse
from collections import defaultdict
import os
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("group sRNAs by genome")

def main():
    parser = argparse.ArgumentParser(description='split sequence')
    parser.add_argument('--input-directory', '-id', type=str, required=True, help='input directory')
    parser.add_argument('--output-directory', '-od', type=str, required=True, help='output directory')
    args = parser.parse_args()
  
    logger.info("Load sRNAs ...")
    sequences = defaultdict(dict) 
    for fasta in os.listdir(args.input_directory):
        path = os.path.join(args.input_directory,fasta)
        sRNA_id = fasta[:fasta.rfind(".")]
        with open(path) as f:
            for header in f:
                genome_id = header[1:].strip().split(" ")[0]
                sequence = next(f).strip()
                sequences[genome_id][sRNA_id] = sequence
    
    logger.info("Group sRNA by genomes ...")
    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory) 
    for genome_id in sequences:
        os.mkdir(os.path.join(args.output_directory,genome_id))
        for sRNA_id in sequences[genome_id]:
            path = os.path.join(args.output_directory,genome_id,sRNA_id+".fa")
            sequence =  sequences[genome_id][sRNA_id]
            with open(path,"w") as f:
                f.write(f">{sRNA_id}\n")
                f.write(f"{sequence}\n")  
    
    logger.info("All done .")


if __name__ == "__main__":
    main()
