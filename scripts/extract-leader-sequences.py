#!/usr/bin/env python
import argparse
from pyfaidx import Fasta

def main():
    parser = argparse.ArgumentParser(description='extract leader sequences')
    parser.add_argument('--input',  '-i', required=True, help="input CDS bed")
    parser.add_argument('--genome',  '-g', required=True, help="input genome sequences")
    parser.add_argument('--output', '-o', required=True, help="input leader fasta")
    parser.add_argument('--left',   '-l', type=int, default = 200,  help="left flanking length")
    parser.add_argument('--right',  '-r', type=int,default = 100,  help="right flanking length")    
    args = parser.parse_args()

    fasta  = Fasta(args.genome)
    fout = open(args.output,"w")
    with open(args.input) as fin:
        for line in fin:
            seq_id, start, end, name, score, strand = line.strip().split("\t")[:6]
            start, end = int(start), int(end)
            if strand == "+":
                s = max(0,start - args.left)
                e = min(start + args.right, len(fasta[seq_id]))
            else:
                s = max(0,end - args.right)
                e = min(end + args.left, len(fasta[seq_id]))
            sequence = fasta[seq_id][s:e]
            if strand == "-":
                sequence = sequence.reverse.complement
            sequence = str(sequence)
            #print(f">{name}:{seq_id}:{s}-{e}({strand})",file=fout)
            print(f">{name}",file=fout)
            print(sequence,file=fout)
    fout.close()


if __name__ == "__main__":
    main()
