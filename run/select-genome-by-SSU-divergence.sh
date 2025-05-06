#!/bin/bash
scripts/select-genome-by-SSU-divergence.py -i genomes/SSU-rRNA.fa -q GCF_000005845.2,GCF_000210855.2,GCF_000742755.1 -od genomes/SSU-99.5 --cutoff 0.995
cat genomes/SSU-99.5/query.fa genomes/SSU-99.5/target.reduce.fa > genomes/SSU-99.5/combined.reduce.fa

scripts/select-genome-by-SSU-divergence.py -i genomes/SSU-99.5/combined.reduce.fa -q GCF_000005845.2,GCF_000210855.2,GCF_000742755.1 -od genomes/SSU-99 --cutoff 0.99
cat genomes/SSU-99/query.fa genomes/SSU-99/target.reduce.fa > genomes/SSU-99/combined.reduce.fa

scripts/select-genome-by-SSU-divergence.py -i genomes/SSU-99/combined.reduce.fa -q GCF_000005845.2,GCF_000210855.2,GCF_000742755.1 -od genomes/SSU-98 --cutoff 0.98
cat genomes/SSU-98/query.fa genomes/SSU-98/target.reduce.fa > genomes/SSU-98/combined.reduce.fa

scripts/select-genome-by-SSU-divergence.py -i genomes/SSU-98/combined.reduce.fa -q GCF_000005845.2,GCF_000210855.2,GCF_000742755.1 -od genomes/SSU-97 --cutoff 0.97
cat genomes/SSU-97/query.fa genomes/SSU-97/target.reduce.fa > genomes/SSU-97/combined.reduce.fa

scripts/select-genome-by-SSU-divergence.py -i genomes/SSU-97/combined.reduce.fa -q GCF_000005845.2,GCF_000210855.2,GCF_000742755.1 -od genomes/SSU-96 --cutoff 0.96
cat genomes/SSU-96/query.fa genomes/SSU-96/target.reduce.fa > genomes/SSU-96/combined.reduce.fa

scripts/select-genome-by-SSU-divergence.py -i genomes/SSU-96/combined.reduce.fa -q GCF_000005845.2,GCF_000210855.2,GCF_000742755.1 -od genomes/SSU-95 --cutoff 0.95
cat genomes/SSU-95/query.fa genomes/SSU-95/target.reduce.fa > genomes/SSU-95/combined.reduce.fa


