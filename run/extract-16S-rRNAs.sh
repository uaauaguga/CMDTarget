#!/bin/bash
indir=genomes
mkdir -p $indir/SSU-rRNA/
#for genome_id in $(cat $indir/Enterobacteriaceae.ge4.txt | cut -f 1);do
for genome_id in $(cat genomes/genome-ids.txt);do
  #echo "processing $genome_id ..."
  #[ -s $indir/SSU-rRNA/${genome_id}.fa ] && echo $genome_id
  [ -s $indir/SSU-rRNA/${genome_id}.fa ] || scripts/extract-16S-rRNAs.py -f $indir/fasta/${genome_id}.fa -o $indir/SSU-rRNA/${genome_id}.fa > $indir/SSU-rRNA/${genome_id}.log 2>&1 
done
