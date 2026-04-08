#!/bin/bash
mkdir -p genomes/CDS
for gff in $(ls genomes/gff);do
  genome_id=${gff%.*}  
  if [ ! -s genomes/CDS/${genome_id}.bed ];then
    echo "processing $genome_id ..."
    scripts/gff2bed.py --gff genomes/gff/${genome_id}.gff --bed genomes/CDS/${genome_id}.bed --feature CDS --name Name,locus_tag,gene > genomes/CDS/${genome_id}.log 2>&1
  else
    echo "$genome_id already present ."
  fi
done
