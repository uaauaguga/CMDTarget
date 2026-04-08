import json
import os

#load all genome ids
genome_ids = open(config["genome-ids"]).read().strip().split("\n")
#get query genome ids
query_ids = config["query-ids"]
for query_id in query_ids:
    if query_id not in genome_ids:
        #print(f"{query_id} not in genome set, insert it.")
        genome_ids.append(query_id)    

print("number of genome:",len(genome_ids))
# define input and output directory
indir = config["indir"]
outdir = config["outdir"]

sRNA_ids_considered0 = None
if "sRNA-ids" in config:
    if isinstance(config["sRNA-ids"],list):
        sRNA_ids_considered0 = config["sRNA-ids"]
    else:
        sRNA_ids_considered0 = open(config["sRNA-ids"]).read().strip().split("\n")    

genome_ids_by_srna = {}
sRNA_ids_considered = []

hits = config["hits"]
for genome_id in genome_ids:
    for fasta in os.listdir(f"{indir}/{hits}/hits.groupped/{genome_id}"):
        sRNA_id = fasta[:fasta.rfind(".")]
        if (sRNA_ids_considered0 is not None) and (sRNA_id not in sRNA_ids_considered0):
            continue
        if sRNA_id not in genome_ids_by_srna:
            genome_ids_by_srna[sRNA_id] = []
        genome_ids_by_srna[sRNA_id].append(genome_id) 

for sRNA_id in genome_ids_by_srna:
    if len(genome_ids_by_srna[sRNA_id]) >= 4:
        sRNA_ids_considered.append(sRNA_id)  

weights = config["weights"]
hfqlabel = ["wo.hfq","wt.hfq"][int(weights.get("hfq zscore",0) > 0)]

params_string = []
weights_string = []

for k,v in weights.items():
    k = k.replace(" ",".")
    params_string.append(f"{k}-{v}")
    weights_string.append(f"{k}-{v}")
weights_string = "_".join(weights_string)

if "normalize" in config:
    k = "normalize"
    v = config[k]
    params_string.append(f"{k}-{v}")

k = "marker"
v = config.get(k,"rpoB")
marker = v
params_string.append(f"{k}-{v}")

denoise = {}
denoise = config["denoise"]
for k,v in denoise.items():    
    params_string.append(f"{k}-{v}")

params_string = "_".join(params_string)

genome_set = config["genome-set-name"]
def get_output(wildcards):
    paths = []
    #["single species scoring","marker detection","comparative denoising"]
    for genome_id in genome_ids:
        sRNA_ids = []
        if not os.path.exists(f"{indir}/{hits}/hits.groupped/{genome_id}/"):
            continue
        for fa in os.listdir(f"{indir}/{hits}/hits.groupped/{genome_id}/"):
            if fa[:-3] not in sRNA_ids_considered:
                continue
            sRNA_ids.append(fa[:-3])
        if "single species scoring" in config["steps"]:
            paths += expand(outdir + f"/{genome_id}/energies/"+"{sRNA_id}.txt",sRNA_id = sRNA_ids)
        if "marker detection" in config["steps"]:
            paths += expand(indir + f"/{marker}/" + "{genome_id}.done",genome_id=genome_id)
    if "comparative denoising" in config["steps"]:
        paths += expand(outdir + "/{query_id}/homologs/{genome_set}/hits.dbtype",query_id=query_ids, genome_set = genome_set)
        paths += [outdir + f"/proteins/{genome_set}/table.txt"]     
        paths += expand(outdir + f"/comparative-analysis-denoised-ml/{hfqlabel}--{genome_set}/" + params_string + "/{sRNA_id}.txt", sRNA_id = sRNA_ids_considered)
        paths += expand(outdir + f"/checkpoints/{hfqlabel}.{genome_set}/" + params_string + "/{sRNA_id}.txt",sRNA_id = sRNA_ids_considered)
    return paths
                
rule all:
    input:
        score = get_output

rule prepare_bed:
    input:
        gff =  indir + "/gff/{genome_id}.gff"
    output:
        bed = indir + "/CDS/{genome_id}.bed"
    shell:
        """
        scripts/gff2bed.py --gff {input.gff} --bed {output.bed} --feature CDS --name Name,locus_tag,gene 
        """

rule prepare_protein:
    input:
        fasta = indir + "/fasta/{genome_id}.fa",
        bed = indir + "/CDS/{genome_id}.bed"
    output:
        protein = indir + "/proteins/{genome_id}.faa",
    log:
        log = indir + "/proteins/{genome_id}.log"
    shell:
        """
        scripts/extract-protein-sequence-from-genome.py -i {input.bed} -g {input.fasta} -o {output.protein} > {log.log} 2>&1
        """        

rule extract_leader:
    input:
        fasta = indir +"/fasta/{genome_id}.fa",
        CDS = indir + "/CDS/{genome_id}.bed"
    output:
        leader = outdir + "/{genome_id}/leader.fa",
        fai = indir +"/fasta/{genome_id}.fa.fai"
    log:
        log = outdir + "/{genome_id}/leader.log"
    shell:
        """
        scripts/extract-leader-sequences.py -i {input.CDS} -g {input.fasta} -o {output.leader}  > {log.log} 2>&1
        """ 

rule energy_scoring:
    input:
        sRNA = indir + "/" + hits + "/hits.groupped/{genome_id}/{sRNA_id}.fa",
        leaders = outdir + "/{genome_id}/leader.fa"
    output:
        energy = outdir + "/{genome_id}/energies/{sRNA_id}.txt"
    log:
        log = outdir + "/{genome_id}/energies/{sRNA_id}.log"
    threads: 4
    shell:
        """
        scripts/calculate-energy-score.py -rs "{input.sRNA}" -ts "{input.leaders}" -o "{output.energy}" -j 8 > {log.log} 2>&1
        """

rule predict_dinucbg_params:
    input:
        sRNA = indir + "/" + hits + "/hits.groupped/{genome_id}/{sRNA_id}.fa",
        leaders = outdir + "/{genome_id}/leader.fa"
    output:
        gumbel = outdir + "/{genome_id}/background.params/{sRNA_id}.txt"
    log:
        log = outdir + "/{genome_id}/background.params/{sRNA_id}.log"
    shell:
        """
        scripts/calculate-dinucleotide-background-params.py -rs "{input.sRNA}" -ts "{input.leaders}" -o "{output.gumbel}" > {log.log} 2>&1
        """


rule energy_normalization:
    input:
        energy = outdir + "/{genome_id}/energies/{sRNA_id}.txt",
        params = outdir + "/{genome_id}/background.params/{sRNA_id}.txt"
    output:
        energy = outdir + "/{genome_id}/energies.normalized/{sRNA_id}.txt"
    log:
        log = outdir + "/{genome_id}/energies.normalized/{sRNA_id}.log"
    shell:
        """
        scripts/normalize-energy.py -i "{input.energy}" -o "{output.energy}" -p "{input.params}" > {log.log} 2>&1
        """

rule hfq_binding_prediction:
    input:
        fasta = indir +"/fasta/{genome_id}.fa",
        CDS = indir + "/CDS/{genome_id}.bed"
    output:
        hfq = indir + "/hfq/{genome_id}.txt" if hfqlabel != "wo.hfq" else config["genome-ids"] # do nothing if not use hfq
    params:
        upstream = 200,
        downstream = 100,
        device = 'cuda:0'
    log:
        log =  indir + "/hfq/{genome_id}.log" 
    shell:
        """
        scripts/leader-hfq-scoring.py --fasta {input.fasta} --bed {input.CDS} --output {output.hfq} \
        --upstream {params.upstream} --downstream {params.downstream} -d {params.device} > {log.log} 2>&1
        """


rule combine_hfq_scores:
    input:
        energy = outdir + "/{genome_id}/energies.normalized/{sRNA_id}.txt",
        hfq = indir + "/hfq/{genome_id}.txt" if hfqlabel != "wo.hfq" else config["genome-ids"] # when not use hfq, the label is arbitrary
    output:
        score = outdir + "/{genome_id}/combined." + hfqlabel + "/{sRNA_id}.txt"
    log:
        log = outdir + "/{genome_id}/combined." + hfqlabel + "/{sRNA_id}.log"
    params:
        hfqlabel = hfqlabel
    shell:
        """   
        scripts/combine-scores.py -hl {params.hfqlabel} -es "{input.energy}" -hs {input.hfq} -o "{output.score}" > {log.log} 2>&1
        """


rule combine_proteins:
    input:
        proteins = expand(indir + "/proteins/{genome_id}.faa",genome_id=genome_ids)
    output:
        proteins = outdir + "/proteins/{genome_set}/db.faa"
    params:
        od = outdir + "/proteins"
    run:
        import os
        print("combine proteins ...")
        fout = open(output.proteins,"w")
        for path in input.proteins:
            genome_id = path.split("/")[-1]
            genome_id = genome_id[:genome_id.rfind(".")]
            with open(path) as f:
                for line in f:
                    if line.startswith(">"):
                        line = ">" + genome_id + ":" + line[1:]            
                    fout.write(line)
        fout.close()

rule build_mmseqs_db:
    input:
        proteins = outdir + "/proteins/{genome_set}/db.faa"
    output:
        db = outdir + "/proteins/{genome_set}/db"
    shell:
        """
        export PATH=$PATH:/BioII/lulab_b/jinyunfan/miniforge3/envs/rna-analysis/bin
        mmseqs createdb --dbtype 1 {input.proteins} {output.db}
        """

rule homolog_search:
    input:
        proteins = outdir + "/proteins/"+genome_set+"/db.faa",
        db = outdir + "/proteins/"+genome_set+"/db",
        query = indir + "/proteins/{query_id}.faa" 
    output:
        hits = outdir + "/{query_id}/homologs/"+genome_set+"/best.hit.txt",
        ahits = outdir + "/{query_id}/homologs/"+genome_set+"/hits.tsv",
        dbt = outdir + "/{query_id}/homologs/"+genome_set+"/hits.dbtype" 
    log:
        log = outdir + "/{query_id}/homologs/"+genome_set+".log"
    params:
        od = outdir + "/{query_id}/homologs/" + genome_set,
        db = outdir + "/proteins/" + genome_set + "/db",
    threads: 5
    shell:
        """
        export PATH=$PATH:/BioII/lulab_b/jinyunfan/miniforge3/envs/rna-analysis/bin
        scripts/protein-homolog-search.py -q {input.query} -db {params.db} -p {input.proteins} -od {params.od} --rerun > {log.log} 2>&1
        """


rule collapse:
    input:
        hits = expand(outdir + "/{query_id}/homologs/"+genome_set+"/best.hit.txt",query_id=query_ids),
    output:
        table = outdir + "/proteins/"+genome_set+"/table.txt",
        collapse = outdir + "/proteins/"+genome_set+"/collapse.txt"
    params:
        wd = outdir,
        genome_set = genome_set,
        query_ids = ",".join(config["query-ids"])
    log:
        log = outdir + "/proteins/"+genome_set+"/collapse.log"
    shell:
        """
        scripts/merge-protein-homolog-hits.py -wd {params.wd} -qi {params.query_ids} -o {output.table} \
        -c {output.collapse} --genome-set {params.genome_set} > {log.log} 2>&1
        """

def get_score_requirements(wildcards):
    paths = []
    sRNA_id = wildcards.sRNA_id
    for g in genome_ids_by_srna[wildcards.sRNA_id]:
        paths.append(outdir + f"/{g}/combined.{hfqlabel}/{sRNA_id}.txt")
    return paths

rule collate_scores_by_homolog:
    input:
        srna = indir + "/" + hits + "/hits/{sRNA_id}.fa",
        score = get_score_requirements, 
        table = outdir + "/proteins/{genome_set}/" + f"table.txt",
        genome_ids = config["genome-ids"]
    output:
        scores = outdir + f"/scores-by-homolog-pair/{hfqlabel}--" + "{genome_set}/{sRNA_id}.txt"
    params:
        wd = outdir,
        genome_set = genome_set,
        hfqlabel = hfqlabel
    log:
        log = outdir + f"/scores-by-homolog-pair/{hfqlabel}--" + "{genome_set}/{sRNA_id}.log"
    shell:
        """
        scripts/collate-homologous-pairs.py -pg {input.table} --srna "{input.srna}" -hl {params.hfqlabel} \
        --input-directory {params.wd} --output "{output.scores}" -gi {input.genome_ids} > {log.log} 2>&1
        """ 

rule prepare_weights:
    input:
    output:
        weights = outdir + "/weights." + weights_string + ".json"
    params:
        weights = json.dumps(weights)
    run:
        with open(output.weights,"w") as f:
            f.write(params.weights)

ruleorder: extract_16S_rRNA_sequence > extract_marker_sequences
rule extract_16S_rRNA_sequence:
    input:
        fasta = indir +"/fasta/{genome_id}.fa"
    output:
        flag = indir + "/16S-rRNA/" + "{genome_id}.done"
    log:
        log = indir + "/16S-rRNA/" + "{genome_id}.log"
    params:
        fasta = indir + "/16S-rRNA/" + "{genome_id}.fa"
    shell:
        """
        scripts/extract-16S-rRNAs.py --fasta {input.fasta} --output {params.fasta} > {log.log} 2>&1 && touch {output.flag}
        """

rule extract_marker_sequences:
    input:
        genome = indir +"/fasta/{genome_id}.fa",
        CDS = indir + "/CDS/{genome_id}.bed",
        protein = indir + "/proteins/{genome_id}.faa",
        hmm = f"models/{marker}.hmm"
    output:
        flag = indir + f"/{marker}/" + "{genome_id}.done"
    params:
        fasta = indir + f"/{marker}/" + "{genome_id}.fa"
    log:
        log = indir + f"/{marker}/" + "{genome_id}.log"   
    shell:
        """
        scripts/extract-protein-marker-nuc-sequences.py -g {input.genome} -p {input.protein} -b {input.CDS} --hmm {input.hmm} -o {params.fasta} > {log.log} 2>&1 && touch {output.flag} 
        """

rule build_tree:
    input:
        fastas = expand(indir + f"/{marker}/" + "{genome_id}.done", genome_id=genome_ids),
        genome_ids = config["genome-ids"]
    output:
        tree = outdir + f"/phylogeny/{genome_set}/{marker}.nwk",
        msa = outdir + f"/phylogeny/{genome_set}/{marker}.afa"
    params:
        indir = indir + f"/{marker}",
        marker = marker,
        outdir = outdir + f"/phylogeny/{genome_set}"
    log:
        log =  outdir + f"/phylogeny/{genome_set}.log"
    shell:
        """
        scripts/build-tree.py -id {params.indir} -od {params.outdir} -gi {input.genome_ids} -m {params.marker} > {log.log} 2>&1
        """


rule extract_representative_genomes:
    input:
        msa = outdir + f"/phylogeny/{genome_set}/{marker}.afa"
    output:
        all_fasta = outdir + f"/phylogeny/{genome_set}/{marker}.fa",
        rep_fasta = outdir + f"/phylogeny/{genome_set}/{marker}.rep.fa",
        rep = outdir + f"/phylogeny/{genome_set}/{marker}-representative-genome-ids.txt"
    log:
        log = outdir + f"/phylogeny/{genome_set}/{marker}-representative-genome-ids.log"
    shell:
        """
        esl-reformat --informat afa fasta {input.msa} > {output.all_fasta} 2> {log.log}
        cd-hit-est -i {output.all_fasta} -o {output.rep_fasta} -c 0.9999999 -r 0 -d 1000 > {log.log} 2>&1
        cat {output.rep_fasta} | grep '>' | sed 's/>//g' > {output.rep}
        """

rule denosing:
    input:
        scores =  outdir + f"/scores-by-homolog-pair/{hfqlabel}" + "--{genome_set}/{sRNA_id}.txt",
        tree = outdir + f"/phylogeny/{genome_set}/"+f"{marker}.nwk",
        weights = outdir + "/weights." + weights_string + ".json"
    output:
        scores = outdir + f"/comparative-analysis-denoised-ml/{hfqlabel}" + "--{genome_set}/" + params_string + "/{sRNA_id}.txt"
    log:
        log = outdir + f"/comparative-analysis-denoised-ml/{hfqlabel}" + "--{genome_set}/" + params_string + "/{sRNA_id}.log"
    threads: 4
    params:
        srm = denoise.get("srm",10),
        nvm = denoise.get("nvm",10),
    shell:
        """
        scripts/comparative-scoring-ml.py --input {input.scores} --output {output.scores} --tree  {input.tree} -sr {params.srm} -nv {params.nvm}  --jobs 10 -w {input.weights} --reroot --normalize > {log.log} 2>&1
        """


rule extract_score:
    input:
        comparative_scores = outdir + f"/comparative-analysis-denoised-ml/{hfqlabel}" + "--{genome_set}/" + params_string + "/{sRNA_id}.txt", 
        genome_scores = outdir + f"/scores-by-homolog-pair/{hfqlabel}" + "--{genome_set}/{sRNA_id}.txt",        
        collapse = outdir + f"/proteins/{genome_set}/collapse.txt",
        weights = outdir + "/weights." + weights_string + ".json",
        rep = outdir + f"/phylogeny/{genome_set}/{marker}-representative-genome-ids.txt" 
    output:
        checkpoint = outdir + f"/checkpoints/{hfqlabel}" + ".{genome_set}/" + params_string + "/{sRNA_id}.txt"
    params:
        wd = outdir,
        hfqlabel = hfqlabel,
        name = "{sRNA_id}",
        tag = hfqlabel + "--" + genome_set + "--" + params_string,
        conservation_weight = config.get("conservation",0)
    log:
        log = outdir + f"/checkpoints/{hfqlabel}" + ".{genome_set}/" + params_string + "/{sRNA_id}.log"
    shell:
        """
        scripts/extract-final-score.py -hl {params.hfqlabel} --comparative-scores "{input.comparative_scores}" --genome-scores "{input.genome_scores}" -cw {params.conservation_weight} \
        --collapse {input.collapse} -w {input.weights} --working-directory {params.wd} \
        --srna-name "{params.name}" --tag {params.tag}".ml" -rgi {input.rep} > "{log.log}" 2>&1 && touch "{output.checkpoint}"
        """
