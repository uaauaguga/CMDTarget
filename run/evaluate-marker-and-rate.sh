#!/bin/bash
#for json in $(ls   grids.20250425 | grep 'rpoB' | grep '0.7' | grep 'pair');do # finished
for json in $(ls   grids.20250425 | grep 'rpoB' | grep '0.7' | grep 'sRNA');do # node45 gpu07, to run
#for json in $(ls   grids.20250425 | grep 'rpoB' | grep '1.0' | grep 'pair');do #node47 node48, to run
#for json in $(ls   grids.20250425 | grep 'rpoB' | grep '1.0' | grep 'sRNA');do # node47 node47, to run
#for json in $(ls   grids.20250425 | grep '16S' | grep '0.7' | grep 'pair');do  #node48 done
#for json in $(ls   grids.20250425 | grep '16S' | grep '0.7' | grep 'sRNA');do #node48 node45, to run
#for json in $(ls   grids.20250425 | grep '16S' | grep '1.0' | grep 'pair');do # node46, to run, hub
#for json in $(ls   grids.20250425 | grep '16S' | grep '1.0' | grep 'sRNA');do # node46, to run
#for json in $(ls   grids.20250425 );do
  genomeset=${json%.*}
  echo $genomeset 
  snakemake  --rerun-triggers mtime --snakefile CMDTarget.snakefile --rerun-incomplete --jobs 32 --configfile  grids.20250425/$json --rerun-incomplete  > log/${genomeset}.genome.selection.log 2>&1
done
