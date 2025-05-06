#!/bin/bash
u=200
d=100
for e in {0..19};do
  outdir=genomes/hfq.${u}.${d}.${e}
  mkdir -p $outdir
  for genome_id in $(cat ../FEMTarget/genomes/enterobacteria.genomes.9.txt);do
   [ -s $outdir/${genome_id}.txt ] ||  scripts/leader-hfq-scoring.py --fasta ../FEMTarget/genomes/fasta/${genome_id}.fa --bed ../FEMTarget/genomes/CDS/${genome_id}.bed --output $outdir/${genome_id}.txt --model models/hfq/${e}.pt --upstream $u --downstream $d -d cuda:1
  done
done

