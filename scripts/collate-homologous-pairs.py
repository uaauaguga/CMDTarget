#!/usr/bin/env python
import argparse
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('combine scores')

def main():
    parser = argparse.ArgumentParser(description='Extract related scores')
    parser.add_argument('--protein-group',  '-pg', type=str, required=True, help='Protein group')
    parser.add_argument('--srna',  '-s', type=str, required=True, help='sRNA sequence path')
    parser.add_argument('--srna-name',  '-sn', type=str, help='sRNA name')
    parser.add_argument('--input-directory',  '-id', type=str, required=True, help='Input directory')
    parser.add_argument('--output','-o', type=str, required=True , help="Output scores")
    parser.add_argument('--hfq-label','-hl', type=str, help="Hfq label")
    parser.add_argument('--genome-ids','-gi', help="genome ids to consider")
    args = parser.parse_args()
    """
    for an sRNA, for genomes with its homolog
    extract interaction scores with hitted proteins of query proteins
    """ 
    
    candidate_genomes = open(args.genome_ids).read().strip().split("\n")
    logger.info("Load homolog lookup table ...")
    homolog2query = {} 
    with open(args.protein_group) as f:
        for line in f:
            query_id, homolog_id = line.strip().split("\t")
            if homolog_id not in homolog2query:
                homolog2query[homolog_id] = []
            homolog2query[homolog_id].append(query_id)
    if args.srna_name is None:
        srna_name = args.srna.split("/")[-1]
        srna_name = srna_name[:srna_name.rfind(".")]
        
    logger.info("Gather genomes contain this sRNA ...")
    genome_ids = []
    with open(args.srna)  as f:
        for line in f:
            if line.startswith(">"):
                genome_id = line[1:].strip()
                if genome_id not in candidate_genomes:
                    continue
                genome_ids.append(genome_id)
    entries_by_query_id = {}
    for genome_id in genome_ids:
        path = f"{args.input_directory}/{genome_id}/combined.{args.hfq_label}/{srna_name}.txt"
        with open(path) as f:
            header = next(f)
            for line in f:
                fields = line[:-1].split("\t")
                protein_id = genome_id + ":" + fields[0]
                fields[0] = protein_id
                if protein_id not in homolog2query:
                    continue
                for query_id in homolog2query[protein_id]:
                    if query_id not in entries_by_query_id:
                        entries_by_query_id[query_id] = []
                    entries_by_query_id[query_id].append(fields)
    logger.info("Saving gathered scores ...")
    fout = open(args.output,"w")
    fout.write("query id\t" + header) 
    for query_id in entries_by_query_id:
        for fields in entries_by_query_id[query_id]:
            print(query_id, *fields,sep="\t", file=fout)
    fout.close()
    logger.info("All done.") 
                    
                 

    
    


if __name__ == "__main__":
    main()
