#!/bin/bash
for gff in $(ls genomes/gff);do
  genome_id=${gff%.*}  
  [ -s genomes/CDS/${genome_id}.bed ] || echo $gff #scripts/gff2bed.py --gff genomes/gff/${genome_id}.gff --bed genomes/CDS/${genome_id}.bed --feature CDS --name Name,locus_tag,gene
done
