#!/usr/bin/env python
import argparse
from Bio.Seq import Seq
from pyfaidx import Fasta
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('extract protein')

def main():
    parser = argparse.ArgumentParser(description='extract protein sequence given coordinate')
    parser.add_argument('--input', '-i', type=str, required=True, help='input bed file')
    parser.add_argument('--genome', '-g', type=str, required=True, help='input genome file')
    parser.add_argument('--output', '-o', type=str, required=True, help='output fasta file')
    args = parser.parse_args()
    
    fout = open(args.output,"w")
    genome = Fasta(args.genome)
    with open(args.input) as f:
        for line in f:
            seq_id, start, end, name, _, strand = line[:-1].split("\t")
            start, end = int(start), int(end)
            sequence = genome[seq_id][start:end]
            if len(sequence) % 3 != 0:
                logger.warning(f"Length {name} not a multiple of three, skip it.")
                continue
            if strand == "-":
                sequence = sequence.reverse.complement
            sequence = Seq(str(sequence)).translate()
            sequence = str(sequence)
            if "*" in sequence[:-1]:
                logger.warning(f"{name} has premature stop codon, skip it.")
                continue
            p = 0
            print(f">{name}",file=fout)          
            while p < len(sequence):
                print(sequence[p:p+70],file=fout)                
                p += 70
    fout.close()
    fout.close()

if __name__ == "__main__":
    main()
