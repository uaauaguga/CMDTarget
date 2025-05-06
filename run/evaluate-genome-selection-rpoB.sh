#!/bin/bash
#for json in $(ls  genome-selection-rpoB-config | awk '(NR>0)&&(NR<=50){print}' );do #gpu07 done
#for json in $(ls  genome-selection-rpoB-config | awk '(NR>50)&&(NR<=70){print}' );do #node45 done
#for json in $(ls  genome-selection-rpoB-config | awk '(NR>70)&&(NR<=90){print}' );do #node46
#for json in $(ls  genome-selection-rpoB-config | awk '(NR>90)&&(NR<=110){print}' );do #node47
#for json in $(ls  genome-selection-rpoB-config | awk '(NR>110)&&(NR<=130){print}' );do #node48 done
#for json in $(ls  genome-selection-rpoB-config | awk '(NR>130)&&(NR<=150){print}' );do #cnode
#for json in $(ls  genome-selection-rpoB-config | awk '(NR>150)&&(NR<=170){print}' );do #hub
for json in $(ls  genome-selection-rpoB-config);do
  genomeset=${json%.*}
  #snakemake --rerun-triggers mtime --snakefile CMDTarget.snakefile --rerun-incomplete --jobs 16 --configfile genome-selection-rpoB-config/$json -np
  #echo 
  snakemake --rerun-triggers mtime --snakefile CMDTarget.snakefile --rerun-incomplete --jobs 30 --configfile genome-selection-rpoB-config/$json  #> log/${genomeset}.genome.selection.log 2>&1
done
