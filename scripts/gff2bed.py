#!/usr/bin/env python
import argparse
import sys 
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('gff2bed')

def parseAttr(s):
    info = {}
    s = s.strip()
    for data in s.split(";"):
        data = data.strip()
        if len(data) == 0:
            continue
        if "=" in data:
            #key, value = data.split("=")
            p = data.find("=")
            key, value = data[:p], data[p+1:]
        else:
            #key,value = data.split(" ")
            p = data.find(" ")
            key, value = data[:p], data[p+1:]
        key = key.replace('"','')
        value = value.replace('"','')
        info[key] = value
    return info



def main():
    parser = argparse.ArgumentParser(description='Convert gff3 format to bed format')
    parser.add_argument('--gff', '-g', type=str, required=True, help='Input gff3 file')
    parser.add_argument('--bed','-b',type=str, required=True, help='Output bed file')
    parser.add_argument('--feature','-f',type=str, help="Keep records in gff3 file where column 3 equal to this value")
    parser.add_argument('--name','-n', type=str , help="Use this field in gff column 9 as feature name in bed file, can be multiple value separated by coma")
    parser.add_argument('--value','-v', type=str , help="Use this field in gff column 9 as feature value in bed file")
    parser.add_argument('--keep-feature','-kf', action="store_true" , help="Keep feature in the name field")
    args = parser.parse_args()

    logger.info(f"Processing {args.gff} ...")
    logger.info(f"Results will be saved to {args.bed} ...")
    fin = open(args.gff)
    fout = open(args.bed,"w")

    encountered_names = {}
    for line in fin:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 3:
            continue
        if (args.feature is not None) and ( fields[2] != args.feature):
            continue
        name = []
        if args.keep_feature:
            name.append(fields[2])
        attrs = parseAttr(fields[8])
        skip = False 
        if args.name is not None:
            for k in args.name.split(","):            
                name.append(attrs.get(k,attrs['locus_tag']))
        name = "-".join(name)
        if name in encountered_names:            
            encountered_names[name] += 1
            suffix = encountered_names[name]
            logger.info(f"{name} encountered, add a suffix -{suffix} to avoid duplicate")
            name = name + "-" + str(encountered_names[name])            
            logger.info(line[:-1] + "\t" + name)
        else:
            encountered_names[name] = 0
        value = []
        skip = False
        if args.value is not None:
            for k in args.value.split(","):
                if k not in attrs:
                    skip = True                     
                else:                
                    value.append(attrs[k])
        if skip:
            continue
        value = ",".join(value)
        chrom, start, end, strand = fields[0], int(fields[3]) - 1, int(fields[4]), fields[6]
        print(f"{chrom}\t{start}\t{end}\t{name}\t{value}\t{strand}",file=fout) 
    fout.close()
    logger.info("All done.")


if __name__ == "__main__":
    main()
