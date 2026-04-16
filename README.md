# CMDTarget

## What is CMDTarget?

- CMDTarget (Comparative Modeling based Denosing for sRNA Target prediction) is a computational tool for sRNA target prediction and evolutionary analysis.
- While conceptually similar to CopraRNA, CMDTarget has several distinct features:
  - CMDTarget explicitly leverage information of evolutionary conservation, while unlike CopraRNA, it does not force homolog sRNA-mRNA pairs across different species to have the same interaction score by explicitly reconstruct interaction scores of both exant and ancestral sRNA-mRNA pairs.
  - CMDTarget include a module for Hfq binding prediction, and could optionally incorporate predicted Hfq binding scores for sRNA target prediction.
  - If you have multiple sRNAs, or a large genome set for comparative analysis, comparative genomics based sRNA target prediction can be computationally intensive. CMDTarget was built on snakemake, CMDTarget is friendly for HPC environments.

## Installation and dependency

### Dependecy

- We recommand using conda to automatically resolve for dependency
```{bash}
conda env create -f environment.yml 
#alternatively, the following command should work
#conda install python snakemake ete3 intarna mafft cd-hit FastTree mmseqs2 hmmer infernal pyfaidx biopython tqdm ypytorch-cuda pytorch -c pytorch -c nvidia -c conda-forge -c bioconda
```

- IntaRNA: for interaction energy calculation
- snakemake: for pipeline management
- ete3: for tree operation
- scipy: implement L-BFGS-B for optimizing scores in tips and internal nodes, also required by ete3
- cd-hit: for sequence clustering
- mafft: for MSA construction 
- FastTree: for 16S rRNA tree building
- mmseqs: for homolog search of sRNA and protein
- pyfaidx: for genome sequence processing
- hmmer: for detection of rpoB marker gene
- infernal: cmsearch for 16S rRNA sequence detection, required is using 16S rRNAs as phylogenetic marker
- pytorch: for Hfq binding score prediction


## Using CMDTarget


### Input of CMDTarget

- The required input of CMDTarget include:
  - Genome sequences
  - Genome annotations in bed format
  - Homologuous sRNA sequences
- These inputs should be organized as follow:
  - `indir`: input directory, should contains the following subdirectory
    - `{indir}/fasta`: Genome sequence in fasta format, named as '{genome_id}.fa' 
    - `{indir}/gff`: genome annotation in gff format, named as '{genome_id}.gff' 
    - `{indir}/{hits}/hits.groupped`
      - organized like '{indir}/{hits}/hits.groupped/{genome_id}/{sRNA_id}.fa', each file contains 1 sRNA sequence
      - homologuous sRNA sequences share the same id
- Easy to use script for preparing input data was provided, see below.

### Prepare input for CMDTarget
- You have multiple ways to curate your genome set for comparative analysis. Here we provide a simple strategy.
- For testing, you can simply use 10 Escherichia genome sequences and annotations in `genomes.tar.gz`, including E.coli (GCF_000005845.2) 
- Suppose you a genome of interest (called query genome), you have to select several phylogenetically related genome for comparative analysis. 
- You may start from GTDB genomes of the clade that contains these genomes. You may only keep refseq genome, or genome with completeness greater than certain cutoff. You get the assembly ids, starts with "GCF". Saving these assembly ids in a text file.

```{bash}
python scripts/select-GTDB-genomes.py --clade g__Escherichia --output genomes/Escherichia.txt
```

- If your query genome(s) does not present in this list, please manually add it

- Then fetch sequences and annotations of these genomes. 

```{bash}
#First download ncbi assembly summary.
wget -P genomes https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
# fetch genomes of interests
mkdir -p genomes/gff genomes/fasta
scripts/fetch-refseq-genomes.py --genome-ids genomes/Escherichia.txt --fasta-directory genomes/fasta --gff-directory genomes/gff --summary genomes/assembly_summary_refseq.txt
```

- Reformat CDS coordinate from gff to bed
```{bash}
bash run/gff2bed.sh
```

- Determine which genomes contain sRNA(s) for target prediction. 

To do so, first put your sRNAs in a fasta file. 
```{bash}
#cat sRNAs.fa
>GCF_000005845.2:RyhB
GCGATCAGGAAGACCCTCGCGGAGAACCTGAAAGCACGACATTGCTCACATTGCTTCCAGTATTACTTAGCCAGCCGGGTGCTGGCTTTTTTTTT
>GCF_000005845.2:Spot42
GTAGGGTACAGAGGTAAGATGTTCTATCTTTCAGACCTTTTACTTCACGTAATCGGATTTGGCTGAATATTTTAGCCGCCCCAGTCAGTAATGACTGGGGCGTTTTTTATT
>GCF_000005845.2:CyaR
GCTGAAAAACATAACCCATAAAATGCTAGCTGTACCAGGAACCACCTCCTTAGCCTGTGTAATCTCCCTTACACGGGCTTATTTTTT
```

Perform the homolog search, and extract hit sequences:
```{bash}
# perform homolog search, homolog sequence saved in genomes/sRNA-hits/hits.fa
scripts/sRNA-homolog-search.py -q sRNAs.fa -gd genomes/fasta -od genomes/sRNA-hits --threads 12 -gi genomes/Escherichia.txt
# group sRNAs by homolog
scripts/group-sRNAs-by-homolog.py -i genomes/sRNA-hits/hits.fa -od genomes/sRNA-hits/hits
# further group sRNA by genomes, make each file contain a single sequence
scripts/split-sRNA-homolog.py -id genomes/sRNA-hits/hits -od genomes/sRNA-hits/hits.groupped
```

### Prepare config file

- Here is an example config file in json format:

```{json}
#example.json
{ "indir": "genomes",
  "outdir": "output/Enterobacteriaceae",
  "hits": "sRNA-hits",
  "genome-set-name": "Escherichia",
  "genome-ids": "genomes/Escherichia.txt",
  "query-ids": ["GCF_000005845.2"],
  "sRNA-ids": "sRNA-ids.txt",
  "marker": "rpoB",
  "weights": {"pair zscore": 0.7, "hfq zscore": 0.3},
  "denoise": {"srm": 10, "nvm": 10, "rescale": 0, "reroot": 1},
  "steps": ["marker detection","single species scoring","comparative denoising"]}
```

- `genome-ids`:
  - A text file specify the genomes used in comparative analysis.
  - One genome id per line
  - For all genome_id with '{indir}/fasta/{genome_id}.fa', '{indir}/gff/{genome_id}.gff' present, only those specified in this file are used
- `query-ids`:
  - A list specify the query genome(s). The genome should present in `genome-ids`
- `sRNA-ids`:
  - A text file specify sRNAs to consider
  - For all sRNAs with `{indir}/{hits}/hits.groupped/{genome_id}/{sRNA_id}.fa` present, only those present in this file are considered
- `marker`:
  - Phylogenetic marker for comparative analysis
- `weights`: weights of energy scores and hfq scores
- `denoise`: parameter for denoising
- `steps`: steps to run. one or more steps can be specified. 
  - "marker detection": detect phyloegentic marker
  - "single species scoring": scoring binding using single species mode
  - "comparative denoising": run comparative denoising. as this step depends on "marker detection" and "single species scoring", if only specify this step, the other two steps will be called

### Run the comparative sRNA target prediction pipeline

- Once the data and config file are ready, the comparative analysis is fully automatic

1. We recommand first run the "marker detection" step to make sure the marker gene present in genomes that will be used for downstream analysis

```{json}
#example.marker.detection.json
{ "indir": "genomes",
  "outdir": "output/Enterobacteriaceae",
  "hits": "sRNA-hits",
  "genome-set-name": "Escherichia",
  "genome-ids": "genomes/Escherichia.txt",
  "query-ids": ["GCF_000005845.2"],
  "sRNA-ids": "sRNA-ids.txt",
  "marker": "rpoB",
  "weights": {"pair zscore": 0.7, "hfq zscore": 0.3},
  "denoise": {"srm": 10, "nvm": 10, "rescale": 0, "reroot": 1},
  "steps": ["marker detection"]}
```

```{bash}
snakemake --configfile example.marker.detection.json --jobs 16 
```
- Marker genes were extracted to `{indir}/{marker}/{genome_id}.fa` 

2. Define genomes set used for downstream analysis

- You may keep all genomes with marker genes
```{bash}
ls genomes/rpoB/ | grep '.fa$' | sed 's/.fa$//' > genome-ids.txt
```
- When you have many genomes, up to 8 genomes withs <0.95 sequence identity is generally sufficient

```{bash}
# combine marker sequences
scripts/combine-fasta.py -i genomes/rpoB -o genomes/Escherichia.rpoB.fa -gi  genomes/Escherichia.txt
# get sequence identity to query marker, the result present in Escherichia.sim2query/scores.txt
scripts/select-genome-by-divergence.py -i genomes/Escherichia.rpoB.fa --query-ids GCF_000005845.2 --cutoff 0.95 -od Escherichia.sim2query/
```

3. Refine your genome set by modify "genome-ids" in config file, run single genome scoring and comparative analysis

```{json}
{ "indir": "genomes",
  "outdir": "output/Enterobacteriaceae",
  "hits": "sRNA-hits",
  "genome-set-name": "ent9", # an alias genome set used for comparative analysis
  "genome-ids": "genome-ids.refined.txt", 
  "query-ids": ["GCF_000005845.2"], # genome for sRNA target prediction, multiple genomes can be provided
  "sRNA-ids": "sRNA-ids.txt",
  "marker": "rpoB",
  "weights": {"pair zscore": 0.7, "hfq zscore": 0.3},
  "denoise": {"srm": 10, "nvm": 10, "rescale": 0, "reroot": 1},
  "steps": ["single species scoring","comparative denoising"]}
```

```{bash}
#--resources gpu=1 restrict jobs for runing Hfq score prediction
snakemake --configfile example.marker.detection.json --jobs 16 --resources gpu=1
```

### Interpretation of the results

- Raw hfq scores: `{indir}/hfq`
- Raw binding energy: `{outdir}/{asm_id}/energies`
- Single species combined scores by genome: `{outdir}/{asm_id}/combined.wt.hfq/{sRNA_id}.txt`
- Single species combined scores by homolog: `{outdir}/scores-by-homolog-pair/wt.hfq--{genome-set-name}/{sRNA_id}.txt`
- The phyloegentic tree: `{outdir}/phylogeny/{genome-set-name}/{marker}.nwk`

### FAQs


## Citation
