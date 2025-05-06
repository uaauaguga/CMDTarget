#!/usr/bin/env python
import argparse
import pandas as pd
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('combine scores')

def main():
    parser = argparse.ArgumentParser(description='Incorporate gene specific score to sRNA-target interaction')
    parser.add_argument('--energy-score',  '-es', type=str, required=True, help='Energy scores')
    parser.add_argument('--hfq-score',  '-hs', type=str, help='Hfq scores')
    parser.add_argument('--hfq-label',  '-hl', type=str, help='Hfq label')
    parser.add_argument('--output','-o', type=str, required=True , help="Output combined energy score")
    args = parser.parse_args()

    logger.info("Load energy scores ...")
    energy = pd.read_csv(args.energy_score,sep="\t",index_col=0)
    energy["protein id"] = energy["target id"].map(lambda x:x.split(":")[0])
    energy = energy.set_index("protein id")

    if args.hfq_label != "wo.hfq":
        logger.info("Load hfq scores ...")
        hfq = pd.read_csv(args.hfq_score,sep="\t",index_col=0)
        hfq.columns = "hfq " + hfq.columns
        energy = energy.join(hfq)        

    logger.info("Saving results ...")
    energy.to_csv(args.output,sep="\t")
    logger.info("All done.")

if __name__ == "__main__":
    main()
     
            
                
        
