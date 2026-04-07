#!/usr/bin/env python
import argparse
import os
import subprocess 
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('download')

def main():
    parser = argparse.ArgumentParser(description='Download genomes from ncbi')
    parser.add_argument('--genome-ids', '-gi', type=str, required=True, help='Input refseq genome ids')
    parser.add_argument('--fasta-directory','-fd', default = "genomes/fasta", help="Directory for fasta file")
    parser.add_argument('--gff-directory','-gd', default = "genomes/gff" , help="Directory for gff file")
    parser.add_argument('--summary','-s', default = "genomes/assembly_summary_refseq.txt" , help="Assembly summary")
    args = parser.parse_args()

    if not os.path.exists(args.fasta_directory):
        os.makedirs(args.fasta_directory)
    if not os.path.exists(args.gff_directory):
        os.makedirs(args.gff_directory)
    logger.info("Load genome ids ...")
    genome_ids = open(args.genome_ids).read().strip().split("\n")
    genome_ids = set(genome_ids)
    N = len(genome_ids)
    logger.info(f"{N} genomes specified.")

     
    asm_id2url = {}
    logger.info("Extract urls ...")    
    with open(args.summary) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            asm_id = fields[0]
            if asm_id not in genome_ids:
                continue
            url = fields[19]
            asm_id2url[asm_id] = url

    n = len(asm_id2url)
    if n != N:
        genomes_not_included = [ asm_id for asm_id in genome_ids if asm_id not in asm_id2url]
        logger.warning(f"{N - n} genomes not present in assembly_summary_refseq.txt:")
        logger.warning(",".join(genomes_not_included))

    for asm_id in asm_id2url:
        logger.info(f"Processing {asm_id} ...")
        url = asm_id2url[asm_id]
        asm_id_long = url.split("/")[-1]         
        fna = os.path.join(args.fasta_directory,asm_id + ".fa")    
        if os.path.exists(fna):
            logger.info(f"{asm_id} genome already retrieved.")
        else:
            logger.info(f"Downloading {asm_id} genome ...")            
            fna_url = url + "/" + asm_id_long + "_genomic.fna.gz"
            fnaz = os.path.join(args.fasta_directory, asm_id + ".fa.gz")  
            subprocess.run(["wget", "-O", fnaz, fna_url])
            subprocess.run(["gunzip",fnaz])
        gff = os.path.join(args.gff_directory,asm_id + ".gff")
        if  os.path.exists(gff): 
            logger.info(f"{asm_id} annotation already retrieved.")
        else:        
            logger.info(f"Downloading {asm_id} annotation ...")
            gff_url = url + "/" + asm_id_long + "_genomic.gff.gz" 
            gffz = os.path.join(args.gff_directory, asm_id + ".gff.gz")
            subprocess.run(["wget", "-O", gffz, gff_url])
            subprocess.run(["gunzip",gffz])

if __name__ == "__main__":
    main()
