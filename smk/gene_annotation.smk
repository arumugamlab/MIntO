#!/usr/bin/env python

'''
Gene prediction and functional annotation step

Authors: Vithiagaran Gunalan, Carmen Saenz, Mani Arumugam
'''

# configuration yaml file
# import sys
import pathlib
from os import path

localrules: make_bed_for_genome, make_cd_transl_faa_for_genome, combine_annotation_outputs, make_merged_subset_bed

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, minto_dir, script_dir, metadata
include: 'include/cmdline_validator.smk'
include: 'include/config_parser.smk'

if config['map_reference'] in ("MAG", "reference_genome"):
    map_reference=config["map_reference"]
else:
    print('ERROR in ', config_path, ': map_reference variable is not correct. "map_reference" variable should be MAG or reference_genome.')

mag_omics = 'metaG'
if map_reference == 'MAG':
    if 'MAG_omics' in config and config['MAG_omics'] != None:
        mag_omics = config['MAG_omics']
    reference_dir="{wd}/{mag_omics}/8-1-binning/mags_generation_pipeline/unique_genomes".format(wd=working_dir, mag_omics=mag_omics)
    print('NOTE: MIntO is using "' + reference_dir + '" as PATH_reference variable')
elif map_reference == 'reference_genome':
    if config['PATH_reference'] is None:
        print('ERROR in ', config_path, ': PATH_reference variable is empty. Please, complete ', config_path)
    reference_dir=config["PATH_reference"]
    print('NOTE: MIntO is using "'+ reference_dir+'" as PATH_reference variable')

if path.exists(reference_dir) is False:
    print('ERROR in ', config_path, ': reference genome path ', reference_dir, ' does not exit. Please, complete ', config_path)

annot_list = list()
if 'ANNOTATION' in config:
    if config['ANNOTATION'] is None:
        print('ERROR in ', config_path, ': ANNOTATION list is empty. "ANNOTATION" variable should be dbCAN, kofam and/or eggNOG. Please, complete ', config_path)
    else:
        try:
            for annot in config["ANNOTATION"]:
                if annot in ['dbCAN', 'kofam', 'eggNOG']:
                    annot_list.append(annot)
                else:
                    raise TypeError
        except:
            print('ERROR in ', config_path, ': ANNOTATION variable is not correct. "ANNOTATION" variable should be dbCAN, kofam and/or eggNOG. Please, complete ', config_path)
else:
    print('ERROR in ', config_path, ': ANNOTATION list is empty. "ANNOTATION" variable should be dbCAN, kofam and/or eggNOG. Please, complete', config_path)

eggNOG_dbmem = True
if 'eggNOG' in annot_list and 'eggNOG_dbmem' in config:
    if not config['eggNOG_dbmem']:
        eggNOG_dbmem = False

# Define all the outputs needed by target 'all'

if map_reference == 'MAG':
    post_analysis_dir="9-MAG-genes-post-analysis"
    post_analysis_out="MAG-genes"
elif map_reference == 'reference_genome':
    post_analysis_dir="9-refgenome-genes-post-analysis"
    post_analysis_out="refgenome-genes"
elif map_reference == 'genes_db':
    post_analysis_dir="9-db-genes-post-analysis"
    post_analysis_out="db-genes"

if map_reference == 'genes_db':
    def merge_genes_output(): # CHECK THIS PART - do not output anything when map_reference == 'genes_db'
        result = expand()
        return(result)

def predicted_genes_collate_out():
    result = expand("{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET.annotations.tsv",
                    wd = working_dir,
                    post_analysis_dir = post_analysis_dir,
                    post_analysis_out = post_analysis_out)
    return(result)

rule all:
    input:
        predicted_genes_collate_out(),
        "{wd}/DB/{post_analysis_dir}/{post_analysis_out}_SUBSET.bed".format(wd = working_dir,
                    post_analysis_dir = post_analysis_dir,
                    post_analysis_out = post_analysis_out)


###############################################################################################
# Prepare predicted genes for annotation
## Combine faa and bed files from the different genomes
###############################################################################################

# Get a sorted list of genomes

def get_genomes_from_refdir(ref_dir):
    genomes = [ pathlib.Path(f).stem for f in os.scandir(ref_dir) if f.is_file() and f.name.endswith('.fna') ]
    return(sorted(genomes))

genomes = get_genomes_from_refdir(reference_dir)

# Generate unique locus ids that are not clashing with each other.
# Clashes are study-specific, so locus-ids will be reproducible for a study with same MAGs.
# But they might be different after resolution when compared between studies.

rule generate_locus_ids:
    input:
        expand("{ref_dir}/{genome}.fna", ref_dir = reference_dir, genome = genomes)
    output:
        "{wd}/DB/{post_analysis_dir}/genomes/locus_id_list.txt"
    localrule: True
    run:
        from hashlib import md5
        from re import sub
        import pathlib
        tr_tb = str.maketrans("0123456789", "GHIJKLMNOP")
        genome_locusids = {}
        for fna in input:
            g = pathlib.Path(fna).stem
            # Get a 10-char translation of md5 checksum
            lid = sub(r'(.)(\1+)', r'\1', md5(bytes(f"{g}\n", 'utf-8'), usedforsecurity=False).hexdigest().translate(tr_tb).upper())[:10]
            # Having FAIL in the locus_id will lead to barrnap failing due to a bug, so replace it
            lid = sub(r'FAIL', r'RAIL', lid)
            if lid not in genome_locusids:
                genome_locusids[lid] = g
            else:
                #clashing
                i = 9
                while lid in genome_locusids:
                    while lid[i] == "Z" and i > 1:
                        i -= 1
                    new_lid_char = chr(ord(lid[i]) + 1)
                    lid = lid[:i] + new_lid_char + lid[i+1:]
                genome_locusids[lid] = g
        with open(output[0], "w") as of:
            for lid, genome in genome_locusids.items():
                print(genome, lid, sep="\t", file=of)


########################
# Prokka on an fna file in reference_dir.
# To make it reproducible, we assign locus_tag based on the name of the MAG/genome.
# This makes testing and benchmarking easier.
########################

rule prokka_for_genome:
    input:
        fna=lambda wildcards: "{reference_dir}/{genome}.fna".format(reference_dir=reference_dir, genome=wildcards.genome),
        locus_ids=lambda wildcards: "{wd}/DB/{post_analysis_dir}/genomes/locus_id_list.txt".format(wd=wildcards.wd, post_analysis_dir=wildcards.post_analysis_dir)
    output:
        fna="{wd}/DB/{post_analysis_dir}/genomes/{genome}/{genome}.fna",
        faa="{wd}/DB/{post_analysis_dir}/genomes/{genome}/{genome}.faa",
        gff="{wd}/DB/{post_analysis_dir}/genomes/{genome}/{genome}.gff",
    log:
        "{wd}/logs/DB/{post_analysis_dir}/prokka/{genome}.log"
    resources:
        mem=10
    threads: 8
    conda:
        config["minto_dir"]+"/envs/gene_annotation.yml"
    shell:
        """
        rm -rf $(dirname {output})
        locus_tag=$(grep "{wildcards.genome}$(printf '\t')" {input.locus_ids} | cut -f 2 )
        prokka --outdir $(dirname {output.fna}) --prefix {wildcards.genome} --locustag $locus_tag --addgenes --cdsrnaolap --cpus {threads} --centre X --compliant {input} >& {log}
        """

# Dont combine FAA files ####

rule make_cd_transl_faa_for_genome:
    input: rules.prokka_for_genome.output.faa
    output: "{wd}/DB/{post_analysis_dir}/CD_transl/{genome}.faa"
    log:
        "{wd}/logs/DB/{post_analysis_dir}/{genome}.reformat_faa.log"
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #seqtk
    shell:
        """
        time (seqtk seq -l 0 {input} | sed -e "s/^>/>{wildcards.genome}|/g" > {output}) >& {log}
        """

# Convert GFF file into BED file ####

rule make_bed_for_genome:
    input: rules.prokka_for_genome.output.gff
    output:
        bed=                      "{wd}/DB/{post_analysis_dir}/GFF/{genome}.bed",
        bed_hdr_mod=              "{wd}/DB/{post_analysis_dir}/GFF/{genome}.header-modif.bed",
        bed_hdr_mod_coord_correct="{wd}/DB/{post_analysis_dir}/GFF/{genome}.header-modif.coord_correct.bed"
    log:
        "{wd}/logs/DB/{post_analysis_dir}/{genome}.make_bed_from_gff.log"
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #gff2bed
    shell:
        """
        time (\
            gff2bed < {input} > {output.bed}
            sed -e "s/^gnl|X|/{wildcards.genome}|/g" {output.bed} > {output.bed_hdr_mod}
            awk -v FS='\\t' -v OFS='\\t' '{{$2+=1; print $0}}' {output.bed_hdr_mod} > {output.bed_hdr_mod_coord_correct} \
        ) >& {log}
        """


rule gene_annot_subset:
    input:
        bed_file = "{wd}/DB/{post_analysis_dir}/GFF/{genome}.header-modif.coord_correct.bed",
        faa_cds = "{wd}/DB/{post_analysis_dir}/CD_transl/{genome}.faa"
    output:
        bed_subset="{wd}/DB/{post_analysis_dir}/subset/{genome}.{post_analysis_out}_SUBSET.bed",
        fasta_subset="{wd}/DB/{post_analysis_dir}/subset/{genome}.{post_analysis_out}_translated_cds_SUBSET.faa"
    log:
        "{wd}/logs/DB/{post_analysis_dir}/{genome}.{post_analysis_out}_subset_out.log"
    resources:
        mem=5
    threads: 1
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    shell:
        """
        remote_dir={wildcards.wd}/DB/{post_analysis_dir}/subset/
        file_prefix={wildcards.genome}.{post_analysis_out}
        time (Rscript {script_dir}/MAGs_genes_out_subset.R {input.bed_file} {input.faa_cds} ${{remote_dir}} ${{file_prefix}}) &> {log}
        """

def get_genome_bed(wildcards):
    #Collect the BED files for MAGs
    genomes = get_genomes_from_refdir(reference_dir)
    result = expand("{wd}/DB/{post_analysis_dir}/subset/{genome}.{post_analysis_out}_SUBSET.bed",
                    wd=wildcards.wd,
                    post_analysis_dir=wildcards.post_analysis_dir,
                    post_analysis_out=wildcards.post_analysis_out,
                    genome=genomes)
    return(result)

rule make_merged_subset_bed:
    input: get_genome_bed
    output:
        bed_file="{wd}/DB/{post_analysis_dir}/{post_analysis_out}_SUBSET.bed",
    log:
        "{wd}/logs/DB/{post_analysis_dir}/{post_analysis_out}.merge_bed.log"
    wildcard_constraints:
        post_analysis_out='MAG-genes|refgenome-genes'
    shell:
        """
        time (\
            cat {input} \
                    | awk -F'\\t' '{{gsub(/[;]/, "\\t", $10)}} 1' OFS='\\t' \
                    | cut -f 1-10 \
                    > {output.bed_file} \
        ) >& {log}
        """

################################################################################################
# Genes annotation using 3 different tools to retrieve functions from:
## kofam
## eggNOG
## dbCAN
################################################################################################

rule gene_annot_kofamscan:
    input:
        fasta_subset=rules.gene_annot_subset.output.fasta_subset,
        ko_list=lambda wildcards: "{minto_dir}/data/kofam_db/ko_list".format(minto_dir = minto_dir),
        prok_hal=lambda wildcards: "{minto_dir}/data/kofam_db/profiles/prokaryote.hal".format(minto_dir = minto_dir),
        module_map=lambda wildcards: "{minto_dir}/data/kofam_db/KEGG_Module2KO.tsv".format(minto_dir = minto_dir),
        pathway_map=lambda wildcards: "{minto_dir}/data/kofam_db/KEGG_Pathway2KO.tsv".format(minto_dir = minto_dir)
    output:
        "{wd}/DB/{post_analysis_dir}/annot/kofam/{genome}.{post_analysis_out}_translated_cds_SUBSET_kofam.tsv",
    shadow:
        "minimal"
    log:
        "{wd}/logs/DB/{post_analysis_dir}/{genome}.{post_analysis_out}_kofamscan.log"
    resources:
        mem=10
    threads: 16
    conda:
        config["minto_dir"]+"/envs/gene_annotation.yml"
    shell:
        """
        remote_dir=$(dirname {output})
        time (
            exec_annotation -k {input.ko_list} -p {input.prok_hal} --tmp-dir tmp --create-alignment -f mapper-one-line --cpu {threads} -o kofam_mapper.txt {input.fasta_subset}
            {script_dir}/kofam_hits.pl --pathway-map {input.pathway_map} --module-map {input.module_map} kofam_mapper.txt > {output}
        ) &> {log}
        """

rule gene_annot_dbcan:
    input:
        fasta_subset=rules.gene_annot_subset.output.fasta_subset
    output:
        "{wd}/DB/{post_analysis_dir}/annot/dbCAN/{genome}.{post_analysis_out}_translated_cds_SUBSET_dbCAN.tsv",
    shadow:
        "minimal"
    params:
        prefix="{post_analysis_out}_translated_cds_SUBSET",
        dbcan_db=lambda wildcards: "{minto_dir}/data/dbCAN_db/V12/".format(minto_dir = minto_dir)
    log:
        "{wd}/logs/DB/{post_analysis_dir}/{genome}.{post_analysis_out}_dbcan.log"
    resources:
        mem=10
    threads: 16
    conda:
        config["minto_dir"]+"/envs/gene_annotation.yml"
    shell:
        """
        if [[ ! -f {params.dbcan_db}/fam-substrate-mapping.tsv ]]; then
          echo "Error: dbCAN database needs to be re-built, run dependencies.smk" > {log};
          exit 1;
        fi

        time (run_dbcan {input.fasta_subset} protein --db_dir {params.dbcan_db} --dia_cpu {threads} --out_pre {params.prefix}_ --out_dir out
        {script_dir}/process_dbcan_overview.pl out/{params.prefix}_overview.txt > {output} )&> {log}
        """

rule gene_annot_eggnog:
    input:
        fasta_subset=rules.gene_annot_subset.output.fasta_subset
    output:
        "{wd}/DB/{post_analysis_dir}/annot/eggNOG/{genome}.{post_analysis_out}_translated_cds_SUBSET_eggNOG.tsv",
    shadow:
        "minimal"
    params:
        prefix="{post_analysis_out}_translated_cds_SUBSET",
        eggnog_inmem = lambda wildcards: "--dbmem" if eggNOG_dbmem else "",
        eggnog_db= lambda wildcards: "{minto_dir}/data/eggnog_data/data/".format(minto_dir = minto_dir) #config["EGGNOG_db"]
    log:
        "{wd}/logs/DB/{post_analysis_dir}/{genome}.{post_analysis_out}_eggnog.log"
    resources:
        mem=lambda wildcards: 50 if eggNOG_dbmem else 16
    threads: lambda wildcards: 24 if eggNOG_dbmem else 10
    conda:
        config["minto_dir"]+"/envs/gene_annotation.yml"
    shell:
        """
        time (
        mkdir out
        emapper.py --data_dir {params.eggnog_db} -o {params.prefix} \
                   --no_annot --no_file_comments --report_no_hits --override --output_dir out -m diamond -i {input.fasta_subset} --cpu {threads}
        emapper.py --annotate_hits_table out/{params.prefix}.emapper.seed_orthologs \
                   --data_dir {params.eggnog_db} -m no_search --no_file_comments --override -o {params.prefix} --output_dir out --cpu {threads} {params.eggnog_inmem}
        cut -f 1,5,12,13,14,21 out/{params.prefix}.emapper.annotations > out/{params.prefix}.eggNOG
        {script_dir}/process_eggNOG_OGs.pl out/{params.prefix}.eggNOG > out/{params.prefix}_eggNOG.tsv
        rm out/{params.prefix}.eggNOG
        sed -i 's/\#query/ID/' out/{params.prefix}_eggNOG.tsv
        sed 's/ko\://g' out/{params.prefix}_eggNOG.tsv > {output} )&> {log}
        """

################################################################################################
# Combine gene annotation results into one file
################################################################################################

def get_genome_annotation_tsvs(wildcards):
    #Collect the CDS faa files for MAGs
    genomes = get_genomes_from_refdir(reference_dir)
    inputs = expand("{wd}/DB/{post_analysis_dir}/annot/{annot}/{genome}.{post_analysis_out}_translated_cds_SUBSET_{annot}.tsv",
                    wd=wildcards.wd,
                    annot=wildcards.annot,
                    post_analysis_dir=wildcards.post_analysis_dir,
                    post_analysis_out=wildcards.post_analysis_out,
                    genome=genomes)
    return(inputs)

rule combine_annotation_outputs:
    input:
        get_genome_annotation_tsvs
    output:
        "{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET_{annot}.tsv"
    log:
        "{wd}/logs/DB/{post_analysis_dir}/{post_analysis_out}_combine_{annot}.log"
    wildcard_constraints:
        post_analysis_out='MAG-genes|refgenome-genes'
    shadow:
        "minimal"
    shell:
        """
        time (head -n 1 {input[0]} > {output};
        FOLDER=$(dirname {input[0]})
        tail -n +2 -q $FOLDER/*.tsv >> {output} ) &> {log}
        """

rule predicted_gene_annotation_collate:
    input:
        annot_out=expand("{{wd}}/DB/{{post_analysis_dir}}/annot/{{post_analysis_out}}_translated_cds_SUBSET_{annot}.tsv", annot=annot_list)
    output:
        annot_out="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_out}_translated_cds_SUBSET.annotations.tsv",
    params:
        prefix="{post_analysis_out}_translated_cds_SUBSET",
    log:
        "{wd}/logs/DB/{post_analysis_dir}/{post_analysis_out}_annot.log"
    resources:
        mem=10
    threads: 9
    conda:
        config["minto_dir"]+"/envs/mags.yml" # python with pandas
    shell:
        """
        time (python3 {script_dir}/collate.py {input} > {output}) >& {log}
        """
