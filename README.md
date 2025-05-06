# CMDTarget

## What is CMDTarget?

- CMDTarget (Comparative Modeling based Denosing for sRNA Target prediction) is a computational tool for sRNA target prediction and evolutionary analysis.
- While conceptually similar to CopraRNA, CMDTarget has several distinct features:
  - CMDTarget explicitly leverage information of evolutionary conservation, while unlike CopraRNA, it does not force homolog sRNA-mRNA pairs across different species to have the same interaction score by explicitly reconstruct interaction scores of both exant and ancestral sRNA-mRNA pairs.
  - CMDTarget include a module for Hfq binding prediction, and could optionally incorporate predicted Hfq binding scores for sRNA target prediction.
  - If you have multiple sRNAs, or a large genome set for comparative analysis, comparative genomics based sRNA target prediction can be computationally intensive. CMDTarget was built on snakemake, CMDTarget is friendly for HPC environments.


## Installation and dependency

### Dependecy
- IntaRNA: for interaction energy calculation
- snakemake: for pipeline management
- ete3: for tree operation
- scipy: implement L-BFGS-B for optimize scores in tips and internal nodes
- cd-hit: for sequence clustering
- cmsearch: for 16S rRNA sequence detection
- mafft: for MSA construction 
- FastTree: for 16S rRNA tree building
- mmseqs: for homolog search of sRNA and protein
- pytorch: for Hfq binding score prediction


## Using CMDTarget


### Prepare genomes for comparative analysis

- First we have to set up data for comparative analysis. Several steps need some manual inspection.
- You have multiple ways to curate your genome set for comparative analysis. Here we provide a simple strategy.
- Suppose you a genome of interest (called query genome), you have to select several phylogenetically related genome for comparative analysis. 
- You may start from GTDB genomes of the clade that contains these genomes. You may only keep refseq genome, or genome with completeness greater than certain cutoff. You get the assembly ids, starts with "GCF". Saving these assembly ids in a text file.

```{bash}
mkdir -p genomes
scripts/select-GTDB-genomes.py --clade f__Enterobacteriaceae --output genomes/Enterobacteriaceae.txt
```

- If your query genome(s) does not present in this list, please manually add it


- Then fetch sequences and annotations of these genomes. 

```{bash}
#First download ncbi assembly summary.
wget -P genomes https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
# fetch genomes of interests
mkdir -p genomes/gff genomes/fasta
scripts/fetch-refseq-genomes.py --genome-ids genomes/Enterobacteriaceae.txt --fasta-directory genomes/fasta --gff-directory genomes/gff --summary genomes/assembly_summary_refseq.txt
```

- Determine which genomes contain sRNA(s) for target prediction. 

To do so, first put your sRNAs in fasta file. 
```{bash}
#cat sRNAs.fa
#each sRNA sequence should be named as "genome_id:sRNA_id"
>GCF_000005845.2:CpxQ
TTTTCCTTGCCATAGACACCATCCCTGTCTTCCCCCACATGCTGTGGGGGTTTTTTTT
>GCF_000005845.2:FnrS
GCAGGTGAATGCAACGTCAAGCGATGGGCGTTGCGCTCCATATTGTCTTACTTCCTTTTTTGAATTACTGCATAGCACAATTGATTCGTACGACGCCGACTTTGATGAGTCGGCTTTTTTTT
>GCF_000005845.2:GcvB
ACTTCCTGAGCCGGAACGAAAAGTTTTATCGGAATGCGTGTTCTGGTGAACTTTTGGCTTACGGTTGTGATGTTGTGTTGTTGTGTTTGCAATTGGTCTGCGATTCAGACCATGGTAGCAAAGCTACCTTTTTTCACTTCCTGTACATTTACCCTGTCTGTCCATAGTGATTAATGTAGCACCGCCTAATTGCGGTGCTTTTTTTT
>GCF_000005845.2:MicA
GAAAGACGCGCATTTGTTATCATCATCCCTGAATTCAGAGATGAAATTTTGGCCACTCACGAGTGGCCTTTTT
>GCF_000005845.2:OmrA
CCCAGAGGTATTGATTGGTGAGATTATTCGGTACGCTCTTCGTACCCTGTCTCTTGCACCAACCTGCGCGGATGCGCAGGTTTTTTTT
>GCF_000005845.2:RybB
GCCACTGCTTTTCTTTGATGTCCCCATTTTGTGGAGCCCATCAACCCCGCCATTTCGGTTCAAGGTTGATGGGTTTTTT
>GCF_000005845.2:RyhB
GCGATCAGGAAGACCCTCGCGGAGAACCTGAAAGCACGACATTGCTCACATTGCTTCCAGTATTACTTAGCCAGCCGGGTGCTGGCTTTTTTTTT
>GCF_000005845.2:SdsR
GGCAAGGCAACTAAGCCTGCATTAATGCCAACTTTTAGCGCACGGCTCTCTCCCAAGAGCCATTTCCCTGGACCGAATACAGGAATCGTGTTCGGTCTCTTTTT
>GCF_000005845.2:SgrS
GATGAAGCAAGGGGGTGCCCCATGCGTCAGTTTTATCAGCACTATTTTACCGCGACAGCGAAGTTGTGCTGGTTGCGTTGGTTAAGCGTCCCACAACGATTAACCATGCTTGAAGGACTGATGCAGTGGGATGACCGCAATTCTGAAAGTTGACTTGCCTGCATCATGTGTGACTGAGTATTGGTGTAAAATCACCCGCCAGCAGATTATACCTGCTGGTTTTTTTT
>GCF_000005845.2:Spot42
GTAGGGTACAGAGGTAAGATGTTCTATCTTTCAGACCTTTTACTTCACGTAATCGGATTTGGCTGAATATTTTAGCCGCCCCAGTCAGTAATGACTGGGGCGTTTTTTATT
>GCF_000005845.2:CyaR
GCTGAAAAACATAACCCATAAAATGCTAGCTGTACCAGGAACCACCTCCTTAGCCTGTGTAATCTCCCTTACACGGGCTTATTTTTT
>GCF_000005845.2:DsrA
AACACATCAGATTTCCTGGTGTAACGAATTTTTTAAGTGCTTCTTGCTTAAGCAAGTTTCATCCCGACCCCCTCAGGGTCGGGATTTTTTT
>GCF_000005845.2:OmrB
CCCAGAGGTATTGATAGGTGAAGTCAACTTCGGGTTGAGCACATGAATTACACCAGCCTGCGCAGATGCGCAGGTTTTTTTT
>GCF_000005845.2:OxyS
GAAACGGAGCGGCACCTCTTTTAACCCTTGAAGTCACTGCCCGTTTCGAGAGTTTCTCAACTCGAATAACTAAAGCCAACGTGAACTTTTGCGGATCTCCAGGATCCGCTTTTTTTT
>GCF_000005845.2:RprA
ACGGTTATAAATCAACATATTGATTTATAAGCATGGAAATCCCCTGAGTGAAACAACGAATTGCTGTGTGTAGTCTTTGCCCATCTCCCACGATGGGCTTTTTTTT
>GCF_000005845.2:RybA
GTGCTATATCTGTATGTAATGCAATCATCCCTCAAGGATCGACGGGATTAGCAAGTCAGGAGGTCTTATGAATGAGTTCAAGAGGTGTATGCGCGTGTTTAGTCATTCTCCCTTTAAAGTACGGTTAATGCTGCTCTCTATGTTGTGCGATATGGTCAACAACAAACCGCAGCAAGATAAACCTTCCGATAAATAGCGGCGTCGCGGTACGCCGCTTCACTCCTGCTTTCATGCAGGCATAACGCGTTTTGGTCTGAAAAACCCCACTTTTTGTCGGATTTGCAATCCCCTTCGCAAAAGATTTGTTCGTCAGTAGTTGACCTGAACGGCGGCTCGCTCT
```

```{bash}
cat sRNAs.fa | grep '>' | cut -f 2 -d ':' > sRNA-ids.txt
```

 
Perform homolog search, and extract hitted sequences:
```{bash}
# perform homolog search, homolog sequence saved in genomes/sRNA-hits/hits.fa
scripts/sRNA-homolog-search.py -q sRNAs.fa -gd genomes/fasta -od genomes/sRNA-hits --threads 12
# group sRNAs by homolog
scripts/group-sRNAs-by-homolog.py -i genomes/sRNA-hits/hits.fa -od genomes/sRNA-hits/hits
# further group sRNA by genomes, make each file contain a single sequence
scripts/split-sRNA-homolog.py -id genomes/sRNA-hits/hits -od genomes/sRNA-hits/hits.groupped
```

You may only consider genomes with at least a subset of query sRNAs:
```{bash}
cat  genomes/sRNA-hits/counts.by.genome.txt | awk '$2>=4{print}' > genomes/Enterobacteriaceae.ge4.txt 
```

- Extract 16S rRNA
```{bash}
run/extract-16S-rRNAs.sh
scripts/combine-fasta.py  --input-directory genomes/SSU-rRNA --output genomes/SSU-rRNA.fa
```

- Select sequence based on 16 rRNA diversity
```{bash}
scripts/select-genome-by-SSU-divergence.py -i genomes/SSU-rRNA.fa -q GCF_000005845.2,GCF_000210855.2,GCF_000742755.1 -od genomes/SSU-96 --cutoff 0.96
cut -f 1 -d ':'  genomes/SSU-96/scores.txt | head -n 10 > genomes/SSU-96-10-genomes.txt
```


### Predict Hfq binding scores (Optional)
- We provide a model for Hfq binding prediction. 
- If your species of interest belong to Gammaproteobacteria, incorporating such information could improve sRNA target prediction performance
- If your species of interest does not harbor a homolog of Hfq protein, please skip this step as the prediction does not make sense.

```{bash}
scripts/leader-hfq-scoring.py --fasta genomes/fasta/$fasta --bed genomes/CDS/${genome_id}.bed --output $outdir/${genome_id}.txt --model models/hfq/${e}.pt --upstream $u --downstream $d -d cuda:1
```


### Prepare config file

- Here is a example config file in json format:

```{json}
{ "genome-set-name": "SSU-96-10-genomes", # name of genome set, can be arbitrary
  "genome-ids": "genomes/SSU-96-10-genomes.txt", # genome used for comparative analysis
  "query-ids": ["GCF_000005845.2"], # genome for sRNA target prediction, multiple genome can be provided
  "sRNA-ids": "sRNA-ids.txt", # sRNA ids to consider
  "indir": "genomes", # input information
  "outdir": "output/Enterobacteriaceae", # where to save the output
  "normalize": 1, # whether renormalize the energy Z score, useful when using "pair zscore" mode
  "hfq": "", # hfq scores to use, optional.
  "weights": {"pair zscore": 1.0, "hfq leader zscore": 0.0}, # weight of energy score and hfq score
  "denoise": {"srm": 10, "nvm": 10, "rescale": 0, "reroot": 1}, # parameter for comparative analysis
  "tag": "202504" # a tag in output , arbitrary
} 
```

### Run the comparative sRNA target prediction pipeline

- Once the data and config file is ready, the comparative analysis is fully automatic

```{bash}
snakemake --configfile example.json --jobs 16
```

### Interpretation of the results

Raw binding energy
- output/Enterobacteriaceae/GCF_000005845.2/energies

## Citation
