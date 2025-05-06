#!/bin/bash
indir=genomes
mkdir -p $indir/rpoB/

for genome_id in $(ls genomes/fasta | grep '.fa$');do
  genome_id=${genome_id%.*}
  if [ -s genomes/CDS/${genome_id}.bed ];then
   echo "processing $genome_id ..."
   mkdir -p output/Enterobacteriaceae/$genome_id
   protein=output/Enterobacteriaceae/$genome_id/proteins.faa
   rpoB=genomes/rpoB/${genome_id}.fa
   [ -s $protein ] || scripts/extract-protein-sequence-from-genome.py -i genomes/CDS/${genome_id}.bed -g genomes/fasta/${genome_id}.fa -o $protein > output/Enterobacteriaceae/$genome_id/proteins.log 2>&1
   [ -s $rpoB ] || scripts/extract-protein-marker-nuc-sequences.py -g genomes/fasta/${genome_id}.fa -p $protein -b genomes/CDS/${genome_id}.bed --hmm models/rpoB.hmm -o $rpoB > genomes/rpoB/${genome_id}.log 2>&1
  fi
done
