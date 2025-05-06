#!/bin/bash
for json in $(ls  genome-selection-config | grep  '99.99' | grep top);do
#for json in $(ls  genome-selection-config | awk '(NR>=10)&&(NR<20){print}');do
  ls  genome-selection-config/$json
  genomeset=${json%.*}
  snakemake --rerun-incomplete --jobs 40 --configfile genome-selection-config/$json > log/${genomeset}.genome.selection.log 2>&1
done
