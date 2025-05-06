#!/bin/bash
u=200
d=100
e=32
outdir=genomes/hfq.${u}.${d}
mkdir -p $outdir
echo $outdir
for fasta in $(ls fasta.9.genomes | grep 'fa$');do
   genome_id=${fasta%.*}
   [ -s $outdir/${genome_id}.txt ] ||  scripts/leader-hfq-scoring.py --fasta fasta.9.genomes/$fasta --bed ../sRNATarget-revised/genomes/CDS/${genome_id}.bed --output $outdir/${genome_id}.txt --model models/hfq/${e}.pt --upstream $u --downstream $d -d cuda:1
done

