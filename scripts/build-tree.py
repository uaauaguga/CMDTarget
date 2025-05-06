#!/usr/bin/env python
import argparse
import sys
import os
import subprocess
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('build tree')

def main():
    parser = argparse.ArgumentParser(description='Build tree')
    parser.add_argument('--input-directory',  '-id', type=str, required=True, help='Input directory contains marker sequences')
    parser.add_argument('--output-directory','-od', type=str, required=True , help="Output MSA and tree")
    parser.add_argument('--genome-ids','-gi', type=str, required=True , help="Genome ids to use")
    parser.add_argument('--marker','-m', type=str, default = "16S-rRNA" , help="Name of the marker gene")
    args = parser.parse_args()

    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory) 

    genome_ids = open(args.genome_ids).read().strip().split("\n")
    logger.info(f"Combine {args.marker} sequences ...")
    RNA = os.path.join(args.output_directory,f"{args.marker}.fa")
    fout = open(RNA,"w")
    #for fasta in os.listdir(args.input_directory):
    for genome_id in genome_ids:
        fasta = genome_id + ".fa"
        path = os.path.join(args.input_directory,fasta)
        asm_id = fasta[:fasta.rfind(".")]
        with open(path) as f:
            header = next(f).strip()
            sequence = next(f).strip()
        print(f">{asm_id}",file=fout)
        print(sequence,file=fout)
    fout.close()


    logger.info("Build MSA ...")
    cmd = ['mafft',"--maxiterate","1000","--localpair",RNA]
    RNA_MSA = os.path.join(args.output_directory,f"{args.marker}.afa")
    print(" ".join(cmd))
    fout = open(RNA_MSA,"w") 
    subprocess.run(cmd,stdout=fout)
    fout.close()


    logger.info("Build tree ...")
    cmd = ["FastTree","-gtr", "-nt",'-gamma', RNA_MSA]
    RNA_tree = os.path.join(args.output_directory,f"{args.marker}.nwk")
    fout = open(RNA_tree,"w")
    subprocess.run(cmd,stdout=fout)
    fout.close()
    
    logger.info("All done.")

if __name__ == "__main__":
    main()
     
            
                
        
