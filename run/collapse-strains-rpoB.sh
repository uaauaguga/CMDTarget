#!/bin/bash
for species in Escherichia.coli Klebsiella.pneumoniae Salmonella.enterica;do
  mkdir -p genomes/rpoB.${species}
  for genome_id in $(cat genomes/genome-by-species/${species}.txt );do 
   [ -s genomes/rpoB/${genome_id}.fa ] && cat genomes/rpoB/${genome_id}.fa > genomes/rpoB.${species}/${genome_id}.fa
  done
  scripts/combine-fasta.py -i genomes/rpoB.${species} -o genomes/rpoB.${species}.fa
  cd-hit-est -i genomes/rpoB.${species}.fa -o genomes/rpoB.${species}.nodup.fa -c 0.999999999999999999999 -r 0 -d 100000
done
