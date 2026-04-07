#!/usr/bin/env python
import pandas as pd
import os
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('filter GTDB genomes')
import argparse
def get_taxonomy(s, l = "c"):
    t = "."
    for c in s.split(";"):
        if c.startswith(l):
            t = c[3:] 
    return t

def main():
    parser = argparse.ArgumentParser(description='extract genomes')
    parser.add_argument('--gtdb-metadata', '-gm', type=str, default="genomes/bac120_metadata_r207.rep.tsv")
    parser.add_argument('--clade', '-c', type=str, required=True, help='clade to extract')
    parser.add_argument('--output', '-o', type=str, required=True, help='output assemblies')
    args = parser.parse_args()
    level, clade = args.clade.split("__")
    table = pd.read_csv(args.gtdb_metadata,sep="\t",index_col=0)
    table[level] = table["gtdb_taxonomy"].map(lambda x:get_taxonomy(x, level))    
    mask = table[level] == clade
    N = mask.sum()
    logger.info(f"{N} genomes belong to {args.clade} .")
    mask = mask & table.index.map(lambda x:x.startswith("RS"))
    N = mask.sum()
    logger.info(f"{N} refseq genomes.")
    mask = mask & (table["checkm_completeness"] > 97)
    N = mask.sum()
    logger.info(f"{N} with completeness > 97.")
    mask = mask & (table["ssu_count"] > 0)
    N = mask.sum()
    logger.info(f"{N} with 16S rRNA.")
    genome_ids = list(table[mask].index.map(lambda x:x[x.find("_")+1:]))
    logger.info(f"Saving selected genome ids to {args.output} .")
    with open(args.output,"w") as f:
        f.write("\n".join(genome_ids) + "\n")
    logger.info("All done.")
                
if __name__ == "__main__":
    main()
