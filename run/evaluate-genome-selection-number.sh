#!/bin/bash
for json in $(ls  genome-selection-config | grep 'GCF_000005845.2');do
  genomeset=${json%.*}
  snakemake --rerun-incomplete --jobs 40 --configfile genome-selection-config/$json > log/${genomeset}.genome.selection.log 2>&1
done
