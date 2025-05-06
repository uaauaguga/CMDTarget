#!/bin/bash
u=200
d=100
e=32
outdir=genomes/hfq.${u}.${d}
mkdir -p $outdir
for fasta in $(ls genomes/fasta | grep 'fa$');do
   genome_id=${fasta%.*}
   if [ ! -s $outdir/${genome_id}.txt ];then 
     if [ -s  genomes/CDS/${genome_id}.bed ];then
        echo $genome_id
        scripts/leader-hfq-scoring.py --fasta genomes/fasta/$fasta --bed genomes/CDS/${genome_id}.bed --output $outdir/${genome_id}.txt --model models/hfq/${e}.pt --upstream $u --downstream $d -d cuda:1
      fi
   fi
done

