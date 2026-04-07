#!/usr/bin/env python
import argparse
import sys
import os
import subprocess
from pyfaidx import Fasta
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('extract 16S rRNA')
import re

def main():
    parser = argparse.ArgumentParser(description='Extract 16S rRNA')
    parser.add_argument('--fasta','-f', type=str , required=True, help="Input genome sequence")
    parser.add_argument('--cm-model','-cm', type=str , default = "models/RF00177-SSU_rRNA_bacteria.cm",help="cm model of 16S rRNA")  
    parser.add_argument('--output','-o', type=str, required=True , help="Output rRNA sequences")
    args = parser.parse_args()

    logger.info("Run cmsearch ...")
    tmp = args.output + ".tmp"
    cmd = ["cmsearch","--tblout",tmp, "--noali", "--hmmonly", "--notrunc", "-E", "1e-5", args.cm_model, args.fasta ]
    print(" ".join(cmd))
    subprocess.run(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

    seq_id, start, end = "", -1, -1
    logger.info("Extract best hit ...")
    with open(tmp) as f:       
        for line in f:
            # cmsearch sort results by e value, simply use the fist hit
            if line.startswith("#"):
                continue
            fields = re.split(r"\s+",line.strip())
            seq_id = fields[0]
            start, end = fields[7], fields[8]
            strand = fields[9]
            start, end = int(start), int(end)
            if start > end:
                assert strand == "-"
                start, end = end, start
            start, end = start - 1, end
            if end - start < 1200:
                continue
            break
    os.remove(tmp)
    if end - start > 1200:
        print(start,end)
        logger.info("Extract 16S rRNA sequence ...")
        fasta = Fasta(args.fasta)
        sequence = fasta[seq_id][start:end]
        if strand == "-":
            sequence = sequence.reverse.complement
        sequence = str(sequence)        
        fout = open(args.output,"w")
        fout.write(f">{seq_id}:{start}-{end}({strand})\n")
        fout.write(sequence + "\n")
        fout.close()
        logger.info("All done.")
    else:
        logger.info("No 16S rRNA detected.")



if __name__ == "__main__":
    main()
