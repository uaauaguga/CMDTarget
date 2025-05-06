#!/usr/bin/env python
import argparse
import os
import subprocess
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('extract')

def main():
    parser = argparse.ArgumentParser(description='Prepare protein sequences')
    parser.add_argument('--fasta','-f', required=True, help="Input genome sequences")
    parser.add_argument('--gff','-g', required = True , help="Input genome annotation")
    parser.add_argument('--cds','-c', required=True , help="Output CDS coordinate")
    parser.add_argument('--protein','-p', required=True , help="Output protein sequences")
    args = parser.parse_args()

    logger.info("Extract CDS coordinates ...")

    cmd = ["scripts/gff2bed.py", "--gff", args.gff, "--bed", args.cds, "--feature", "CDS", "--name", "Name,locus_tag,gene"]
    subprocess.run(cmd)


    logger.info("Extract protein sequences ...")

    cmd = ["scripts/extract-protein-sequence-from-genome.py","-i",args.cds,"-g",args.fasta,'-o',args.protein] 
    subprocess.run(cmd)

    logger.info("All done.")


if __name__ == "__main__":
    main()

