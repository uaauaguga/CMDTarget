#!/usr/bin/env python
import argparse
from collections import defaultdict
import os
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("group sRNA homologs")


def main():
    parser = argparse.ArgumentParser(description='split sequence')
    parser.add_argument('--input', '-i', type=str, required=True, help='input sequences')
    parser.add_argument('--output-directory', '-od', type=str, required=True, help='output directory')
    args = parser.parse_args()
    

    logger.info("Load hits ...")
    sRNA_ids = {}
    skipped_pairs = set()
    lines = defaultdict(list)
    with open(args.input) as f:
        for header in f:
            sequence = next(f).strip()
            query_id, hit_id = header[1:].strip().split("--")
            p = query_id.find(":")
            query_genome_id, sRNA_id  = query_id[:p], query_id[p+1:]
            hit_genome_id, contig_id = hit_id.split(":")[:2]
            if (sRNA_id in sRNA_ids) and (sRNA_ids[sRNA_id] != query_genome_id):
                if (sRNA_id, query_genome_id) not in skipped_pairs:
                    skipped_pairs.add((sRNA_id, query_genome_id))
                    print(f"{sRNA_id}: {query_genome_id} => {sRNA_ids[sRNA_id]}")
                continue
            sRNA_ids[sRNA_id] = query_genome_id
            contig_id = contig_id.split(".")[0]
            lines[sRNA_id].append(f">{hit_genome_id}")
            lines[sRNA_id].append(sequence) 

    logger.info("Save hits ...")
    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)
    for sRNA_id in lines:
        fout = open(os.path.join(args.output_directory,sRNA_id + ".fa"),"w")
        print("\n".join(lines[sRNA_id]),file=fout)
        fout.close()

    logger.info("All done .")
                    


if __name__ == "__main__":
    main()
