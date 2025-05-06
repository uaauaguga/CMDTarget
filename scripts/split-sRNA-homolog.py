#!/usr/bin/env python
import argparse
from collections import defaultdict
import os

def main():
    parser = argparse.ArgumentParser(description='split sequence')
    parser.add_argument('--input-directory', '-id', type=str, required=True, help='input directory')
    parser.add_argument('--output-directory', '-od', type=str, required=True, help='output directory')
    args = parser.parse_args()
  
    sequences = defaultdict(dict) 
    for fasta in os.listdir(args.input_directory):
        path = os.path.join(args.input_directory,fasta)
        sRNA_id = fasta[:fasta.rfind(".")]
        with open(path) as f:
            for header in f:
                genome_id = header[1:].strip().split(" ")[0]
                sequence = next(f).strip()
                sequences[genome_id][sRNA_id] = sequence
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


if __name__ == "__main__":
    main()
