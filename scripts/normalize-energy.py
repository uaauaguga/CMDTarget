#!/usr/bin/env python
import logging
import argparse
from scipy.stats import gumbel_r, norm
import numpy as np
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('normalization')


def main():
    parser = argparse.ArgumentParser(description='Normalize RNA-RNA interaction score to make them comparable across different species')
    parser.add_argument('--input',  '-i', type=str, required=True, help='Input pairwise scoring')
    parser.add_argument('--output', '-o', type=str, required=True, help='Output energy')
    parser.add_argument('--params', '-p', type=str, help='Gumbel parameter')
    args = parser.parse_args()

    logger.info("Load interaction scores ...")
    entries_by_sRNA = {}
    energies_by_sRNA = {}
    with open(args.input) as f:
        header = next(f)
        for line in f:
            fields = line.strip().split("\t")
            #AgrA    3       15      NP_414544.1-b0003-thrB:NC_000913.32800-3733(+)  77      92      -7.3296
            sRNA_id = fields[0]
            if sRNA_id not in entries_by_sRNA:
                entries_by_sRNA[sRNA_id] = []
                energies_by_sRNA[sRNA_id] = []
            entries_by_sRNA[sRNA_id].append(fields[1:])
            energy = float(fields[-1])            
            if np.isnan(energy):
                energy = 0
            energies_by_sRNA[sRNA_id].append(energy)
    
    if args.params is not None:
        pair_to_params = {}
        with open(args.params) as f:
            _ = next(f)
            for line in f: 
                fields = line.strip().split("\t")
                sRNA_id, target_id, loc, scale = fields
                loc, scale = float(loc), float(scale)
                pair_to_params[(sRNA_id, target_id)]  = loc, scale
    sRNA_to_params = {}
    for sRNA_id in energies_by_sRNA:
        loc, scale = gumbel_r.fit(-np.array(energies_by_sRNA[sRNA_id])/10)
        sRNA_to_params[sRNA_id] = loc, scale

    fout = open(args.output,"w")
    header = header.strip() + "\tsRNA pvalue\tsRNA zscore"
    if args.params is not None:
        header = header + "\tpair pvalue\tpair zscore"
    header += "\n"
    fout.write(header)
    for sRNA_id in entries_by_sRNA:
        sloc, sscale = sRNA_to_params[sRNA_id]
        logger.info(f"{sRNA_id}: loc={sloc}, scale={scale}.")
        for fields in entries_by_sRNA[sRNA_id]:
            target_id = fields[2]            
            energy = float(fields[-1])
            if np.isnan(energy):
                energy = 0
            pvalue = 1 - gumbel_r.cdf(-energy/10, loc=sloc,scale=sscale)
            zscore = norm.ppf(1-pvalue)
            fields += [round(pvalue,4),round(zscore,4)]
            if args.params is not None:
                loc, scale = pair_to_params[(sRNA_id, target_id)]
                pvalue = 1 - gumbel_r.cdf(-energy/10, loc=loc,scale=scale)
                zscore = norm.ppf(1-pvalue)
                fields += [round(pvalue,4),round(zscore,4)]
            print(sRNA_id, *fields, sep="\t", file=fout)
    fout.close()

if __name__ == "__main__":
    main()
