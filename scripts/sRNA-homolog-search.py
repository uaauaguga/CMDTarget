#!/usr/bin/env python
import argparse
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("pick hits")
import subprocess
import os
from collections import defaultdict
from pyfaidx import Fasta

def main():
    parser = argparse.ArgumentParser(description='sRNA homolog search')
    parser.add_argument('--query', '-q',required=True,help="input query sRNAs")
    parser.add_argument('--genome-directory','-gd',required=True,help="target genomes")
    parser.add_argument('--output-directory','-od',required=True,help="output directory")
    parser.add_argument('--coverage','-c', type=float, default=0.8,help="alignment coverage required")
    parser.add_argument('--threads','-t',default=8,help="threads for mmseqs search")
    args = parser.parse_args()
    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)
    genome_ids = []
    genome_fasta = os.path.join(args.output_directory,"genomes.fa")
    if not os.path.exists(genome_fasta):  
        logger.info("prepare genome sequences ...")
        fout = open(genome_fasta,"w")
        for fasta in os.listdir(args.genome_directory):
            if not (fasta.endswith(".fa") or fasta.endswith(".fna")):
                continue
            path = os.path.join(args.genome_directory,fasta)
            genome_id = fasta[:fasta.rfind(".")]
            genome_ids.append(genome_id)
            with open(path) as f:
                for line in f:
                    if line.startswith(">"):
                        line = line[:1] + genome_id + ":" + line[1:]
                    fout.write(line)
        fout.close()


    tdb = os.path.join(args.output_directory,"genomes")
    if not os.path.exists(tdb + ".dbtype"):
        logger.info("build genome mmseqs database ...")
        cmd = ["mmseqs", "createdb", "--dbtype", "2", genome_fasta, tdb ]
        subprocess.run(cmd)
    else:
        logger.info("genome db exists .")

    qdb = os.path.join(args.output_directory,"sRNAs")
    if not os.path.exists(qdb + ".dbtype"):
        logger.info("build sRNA mmseqs database ...")
        cmd = ["mmseqs", "createdb", "--dbtype", "2", args.query, qdb ]
        subprocess.run(cmd)
    else:
        logger.info("sRNA db exists .")

    hdb = os.path.join(args.output_directory,"hits")
    if not os.path.exists(hdb+".dbtype"):
        logger.info("homolog search ...")
        cmd = ["mmseqs", "search", "--search-type", "3", "-k", "9", "-s", "7.5", "--max-seqs", "5000", qdb, tdb, hdb, "tmp" ]
        subprocess.run(cmd)
    else:
        logger.info("hit db exists .")


    tsv = os.path.join(args.output_directory,"hits.tsv")
    if not os.path.exists(tsv):
        logger.info("reformat homolog search hits ...")
        cmd = ["mmseqs", "convertalis", "--format-mode", "2", qdb, tdb, hdb, tsv ]
        subprocess.run(cmd)
    else:
        logger.info("hit file exists .")
            

    logger.info("extract best hits by genome ...")
    bitscore_by_genome = defaultdict(int)
    best_hits = {}
    with open(tsv) as f:
        for line in f:
            fields = line.strip().split("\t")
            contig_id, tstart, tend = fields[1], int(fields[8]), int(fields[9])
            qstart, qend = int(fields[6]), int(fields[7])
            assert tstart < tend
            query_length, target_length = int(fields[12]), int(fields[13])
            aligned_length = int(fields[3])
            coverage = aligned_length/query_length
            if coverage < args.coverage:
                continue     
            genome_id = contig_id[:contig_id.find(":")]
            bitscore = int(fields[11]) 
            if qstart < qend:
                strand = "+"
            else:
                strand = "-"
                qstart, qend = qend, qstart 
            query_id = fields[0]
            qstart -= 1
            tstart -= 1
            if bitscore > bitscore_by_genome[(query_id, genome_id)]:
                left_offset, right_offset = qstart, query_length - qend
                if strand == "+":
                    start, end = tstart - left_offset, tend + right_offset
                else:
                    start, end = tstart - right_offset, tend + left_offset
                if (start < 0) or (end > target_length):
                    continue
                best_hits[(query_id,genome_id)] = (contig_id, start, end, strand)
                bitscore_by_genome[(query_id, genome_id)] = bitscore
   
    fout = open(os.path.join(args.output_directory,"hits.fa"),"w")
    fasta = Fasta(genome_fasta)
    logger.info("save homolog sequences ...")    

    counts_by_genome = defaultdict(list)
    for query_id,genome_id in sorted(list(best_hits.keys())):
        counts_by_genome[genome_id].append(query_id)
        contig_id, start, end, strand = best_hits[(query_id,genome_id)]
        sequence = fasta[contig_id][start:end]
        if strand == "-":
            sequence = sequence.reverse.complement
        sequence = str(sequence)
        print(f">{query_id}--{contig_id}:{start}-{end}({strand})",file=fout)
        print(sequence,file=fout)
    fout.close() 
    fstat = open(os.path.join(args.output_directory,"counts.by.genome.txt"),"w")
    for genome_id in genome_ids:
        print(genome_id, len(counts_by_genome[genome_id]), ",".join(counts_by_genome[genome_id]), sep="\t",file=fstat)
    fstat.close()
    logger.info("All done.") 

              

if __name__ == "__main__":
    main()
