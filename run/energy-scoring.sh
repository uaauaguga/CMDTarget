#nohup snakemake --configfile configs/CyaR.json --jobs 16 > log/CyaR.log 2>&1 & #cnode all
#nohup snakemake --configfile configs/DsrA.json --jobs 16 > log/DsrA.log 2>&1 & #hub
#nohup snakemake --configfile configs/OmrB.json --jobs 16 > log/OmrB.log 2>&1 & #45 all
#nohup snakemake --configfile configs/OxyS.json --jobs 20 > log/OxyS.log 2>&1 & #46 all
#nohup snakemake --configfile configs/RprA.json --jobs 20 > log/RprA.log 2>&1 & #47 all
#nohup snakemake --configfile configs/RybA.json --jobs 20 > log/RybA.log 2>&1 & #gpu all
