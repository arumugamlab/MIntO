#!/usr/bin/env python

'''
Gene prediction and functional annotation step

Authors: Vithiagaran Gunalan, Carmen Saenz, Mani Arumugam
'''

# configuration yaml file
# import sys
import pathlib
from os import path

localrules: rename_prokka_sequences, combine_annotation_outputs, combine_individual_beds

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, minto_dir, script_dir, metadata
include: 'include/cmdline_validator.smk'
include: 'include/config_parser.smk'

module print_versions:
    snakefile:
        'include/versions.smk'
    config: config

use rule annotation_base from print_versions as version_*

snakefile_name = print_versions.get_smk_filename()

# MIntO mode and database-mapping

# Define the 2 modes ('catalog' mode is not allowed here)
valid_minto_modes = ['MAG', 'refgenome']

# Which database are we mapping reads to?
if 'MINTO_MODE' in config and config['MINTO_MODE'] != None:
    MINTO_MODE=config['MINTO_MODE']
else:
    raise Exception("ERROR in {}: 'MINTO_MODE' variable must be defined".format(config_path))

# Backward compatibility and common misnomers
if MINTO_MODE in ['reference_genome', 'reference-genome', 'reference', 'refgenomes']:
    MINTO_MODE = 'refgenome'
elif MINTO_MODE in ['MAGs', 'mag', 'mags']:
    MINTO_MODE = 'MAG'

if not MINTO_MODE in valid_minto_modes:
    raise Exception("ERROR in {}: 'MINTO_MODE' variable must be {}.".format(config_path, valid_minto_modes))

mag_omics = 'metaG'
if MINTO_MODE == 'MAG':
    if 'MAG_omics' in config and config['MAG_omics'] != None:
        mag_omics = config['MAG_omics']
    reference_dir="{wd}/{mag_omics}/8-1-binning/mags_generation_pipeline/unique_genomes".format(wd=working_dir, mag_omics=mag_omics)
    print('NOTE: MIntO is using "' + reference_dir + '" as PATH_reference variable')
elif MINTO_MODE == 'refgenome':
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

GENE_DB_TYPE = MINTO_MODE + '-genes'

rule all:
    input:
        "{wd}/DB/{subdir}/4-annotations/combined_annotations.tsv".format(
                    wd = working_dir,
                    subdir = MINTO_MODE),
        "{wd}/DB/{subdir}/genomes.marker_genes.table".format(
                    wd = working_dir,
                    subdir = MINTO_MODE),
        "{wd}/DB/{subdir}/{filename}.bed.mini".format(
                    wd = working_dir,
                    subdir = MINTO_MODE,
                    filename = GENE_DB_TYPE),
        print_versions.get_version_output(snakefile_name)
    default_target: True

######################
# 1. PREPARE GENOMES
######################

###############################################################################################
# Prepare predicted genes for annotation
## Combine faa and bed files from the different genomes
###############################################################################################

# Get a sorted list of genomes

def get_genomes_from_refdir(ref_dir):
    g= [ pathlib.Path(f).stem for f in os.scandir(ref_dir) if f.is_file() and f.name.endswith('.fna') ]
    return(sorted(g))

genomes = get_genomes_from_refdir(reference_dir)

# Generate unique locus ids that are not clashing with each other.
# Clashes are study-specific, so locus-ids will be reproducible for a study with same MAGs.
# But they might be different after resolution when compared between studies.

rule generate_locus_ids:
    input:
        expand("{ref_dir}/{genome}.fna", ref_dir = reference_dir, genome = genomes)
    output:
        "{wd}/DB/{minto_mode}/1-prokka/locus_id_list.txt"
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
            print('mag_id', 'locus_id', sep="\t", file=of)
            for lid, genome in genome_locusids.items():
                print(genome, lid, sep="\t", file=of)

######################
# 2. GENE PREDICTION
######################

########################
# Prokka on an fna file in reference_dir.
# To make it reproducible, we assign locus_tag based on the name of the MAG/genome.
# This makes testing and benchmarking easier.
########################

rule prokka_for_genome:
    input:
        fna=lambda wildcards: "{reference_dir}/{genome}.fna".format(reference_dir=reference_dir, genome=wildcards.genome),
        locus_ids="{wd}/DB/{minto_mode}/1-prokka/locus_id_list.txt"
    output:
        fna="{wd}/DB/{minto_mode}/1-prokka/{genome}/{genome}.fna",
        faa="{wd}/DB/{minto_mode}/1-prokka/{genome}/{genome}.faa",
        gff="{wd}/DB/{minto_mode}/1-prokka/{genome}/{genome}.gff",
        gbk="{wd}/DB/{minto_mode}/1-prokka/{genome}/{genome}.gbk",
        summary="{wd}/DB/{minto_mode}/1-prokka/{genome}/{genome}.summary",
    shadow:
        "minimal"
    log:
        "{wd}/logs/DB/{minto_mode}/1-prokka/{genome}.log"
    resources:
        mem=10
    threads: 8
    conda:
        config["minto_dir"]+"/envs/gene_annotation.yml"
    shell:
        """
        time (
            locus_tag=$(grep "{wildcards.genome}$(printf '\\t')" {input.locus_ids} | cut -f 2 )
            prokka --outdir out --prefix tmp --locustag $locus_tag --addgenes --cdsrnaolap --cpus {threads} --centre X --compliant {input.fna}
            rsync -a out/tmp.fna {output.fna}
            rsync -a out/tmp.faa {output.faa}
            rsync -a out/tmp.gff {output.gff}
            rsync -a out/tmp.gbk {output.gbk}
            rsync -a out/tmp.txt {output.summary}
        ) >& {log}
        """

# Dont combine FAA files ####

rule rename_prokka_sequences:
    input:
        gff=rules.prokka_for_genome.output.gff,
        fna=rules.prokka_for_genome.output.fna,
        faa=rules.prokka_for_genome.output.faa
    output:
        bed="{wd}/DB/{minto_mode}/2-postprocessed/{genome}.bed",
        fna="{wd}/DB/{minto_mode}/2-postprocessed/{genome}.fna",
        faa="{wd}/DB/{minto_mode}/2-postprocessed/{genome}.faa"
    shadow:
        "minimal"
    log:
        "{wd}/logs/DB/{minto_mode}/{genome}.rename_prokka.log"
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml" #gff2bed
    shell:
        """
        time (
            # Change sequence names in fna
            sed "s/^>gnl|X|/>/" {input.fna} > out.fna

            # Change sequence names in faa and keep only first word
            sed -e "s/\\s.*//" {input.faa} > out.faa

            # Change sequence names and make BED:
            # -----------------------------------
            # Ignore 'gene' lines
            # Call gff2bed
            # Modify header
            # Set ID based on ID=X
            # If no ID=X then set ID as <chr>_<start>_<stop>
            # When relevant, set feature_type to CRISPR (compliant with GFF3)

            cat {input.gff} \
                | awk -F '\\t' '$3 != "gene"' \
                | gff2bed \
                | perl -lan -F"\\t" \
                    -e '$F[0] =~ s/^gnl\\|X\\|//;' \
                    -e 'if ($F[9] =~ /ID=([^;]*);.*/) {{
                            $F[9] = "$1";
                        }} else {{
                            if ($F[9] =~ /note=CRISPR/) {{
                                $F[7] = "CRISPR";
                            }}
                            $F[9] = sprintf("%s_%s_%s", $F[0], $F[1]+1, $F[2]);
                        }}' \
                    -e 'print join("\\t", @F);' \
                > out.bed

            # Save output files
            rsync -a out.bed {output.bed}
            rsync -a out.fna {output.fna}
            rsync -a out.faa {output.faa}
        ) >& {log}
        """

# Get a list of genome bed files
def get_genome_bed(wildcards):
    #Collect the BED files for MAGs
    result = expand("{wd}/DB/{minto_mode}/2-postprocessed/{genome}.bed",
                    wd=wildcards.wd,
                    minto_mode=wildcards.minto_mode,
                    genome=genomes)
    return(result)

# Replace 'score' in 5th column with 'length' = stop - start. This is BED, so no need for 'stop-start+1'.
# Retain only first word in ID column
# Write full info into full BED file.
# Write only chr,start,stop,length,ID and skip CRISPR in mini BED file.

def process_genome_bed_list(input_list, output_full, output_mini, log_file):
    import pandas as pd
    import datetime

    def logme(stream, msg):
        print(datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), msg, file=stream)

    colnames = ['chr', 'start', 'stop', 'name', 'length', 'strand', 'source', 'feature', 'frame', 'ID']
    drop_cols = ['name', 'strand', 'source', 'feature', 'frame']
    with open(str(log_file), 'w') as f:
        with open(str(output_full), 'w') as full, open(str(output_mini), 'w') as mini:
            for i in range(0, len(input_list)):
                logme(f, "INFO: reading file {}".format(i))

                df = pd.read_csv(input_list[i], header=None, names=colnames, sep="\t", memory_map=True)
                df['length'] = df['stop'] - df['start']
                df = df.replace(to_replace={'ID': r'[;\s].*'}, value='', regex=True)
                df.to_csv(full, sep="\t", index=False, header=(i==0)) # Only write header the first time around

                df = df.query('feature != "CRISPR"').drop(drop_cols, axis='columns')
                df.to_csv(mini, sep="\t", index=False, header=False)

        logme(f, "INFO: done")

# Combine individual BEDs but only select chr,start,stop,strand,feature,ID
rule combine_individual_beds:
    input: get_genome_bed
    output:
        bed_full="{wd}/DB/{minto_mode}/{filename}.bed",
        bed_mini="{wd}/DB/{minto_mode}/{filename}.bed.mini",
    log:
        "{wd}/logs/DB/{minto_mode}/{filename}.merge_bed.log"
    shadow:
        "minimal"
    wildcard_constraints:
        minto_mode='MAG|refgenome'
    run:
        import shutil
        process_genome_bed_list(input, 'full.bed', 'mini.bed', log)
        shutil.copy2('full.bed', output.bed_full)
        shutil.copy2('mini.bed', output.bed_mini)

######################
# 3. MARKER GENES
######################

###############################################################################################
# Normalization of read counts to 10 marker genes (MG normalization)
## fetchMGs identifies 10 universal single-copy phylogenetic MGs
## (COG0012, COG0016, COG0018, COG0172, COG0215, COG0495, COG0525, COG0533, COG0541, and COG0552)
###############################################################################################

# Run fetchMGs
# fetchMG cannot handle '.' in gene names, so we replace '.' with '__MINTO_DOT__' in fasta headers; then change back in output.
rule fetchMG_genome_cds_faa:
    input:
        cds_faa      = "{wd}/DB/{minto_mode}/2-postprocessed/{genome}.faa",
        fetchMGs_dir = "{minto_dir}/data/fetchMGs-1.2".format(minto_dir = minto_dir)
    output: '{wd}/DB/{minto_mode}/fetchMGs/{genome}/{genome}.marker_genes.table'
    localrule: True
    shadow:
        "minimal"
    log: '{wd}/logs/DB/{minto_mode}/fetchMGs/{genome}.log'
    threads: 8
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    shell:
        """
        time (
            sed 's/\\s.*//;s/\\./__MINTO_DOT__/g' {input.cds_faa} > HEADER_FIXED.faa
            {input.fetchMGs_dir}/fetchMGs.pl -outdir out -protein_only -threads {threads} -x {input.fetchMGs_dir}/bin -m extraction HEADER_FIXED.faa
            sed 's/__MINTO_DOT__/./g' out/HEADER_FIXED.all.marker_genes_scores.table > {output}
        ) >& {log}
        """

def get_genome_MG_tables(wildcards):
    #Collect the CDS faa files for MAGs
    genomes = get_genomes_from_refdir(reference_dir)
    result = expand("{wd}/DB/{minto_mode}/fetchMGs/{genome}/{genome}.marker_genes.table",
                    wd=wildcards.wd,
                    minto_mode=wildcards.minto_mode,
                    genome=genomes)
    return(result)

rule merge_MG_tables:
    input: get_genome_MG_tables
    output: "{wd}/DB/{minto_mode}/genomes.marker_genes.table"
    localrule: True
    log: "{wd}/logs/DB/{minto_mode}/merge_marker_genes_scores.table.log"
    shell:
        """
        time (
                head -n 1 {input[0]} > {output}
                for file in {input}; do
                    awk 'FNR>1' ${{file}} >> {output}
                done
        ) >& {log}
        """

###########################
# 4. FUNCTIONAL ANNOTATION
###########################

################################################################################################
# Genes annotation using 3 different tools to retrieve functions from:
## kofam
## eggNOG
## dbCAN
################################################################################################

rule gene_annot_kofamscan:
    input:
        faa=rules.rename_prokka_sequences.output.faa,
        ko_list=lambda wildcards: "{minto_dir}/data/kofam_db/ko_list".format(minto_dir = minto_dir),
        prok_hal=lambda wildcards: "{minto_dir}/data/kofam_db/profiles/prokaryote.hal".format(minto_dir = minto_dir),
        module_map=lambda wildcards: "{minto_dir}/data/kofam_db/KEGG_Module2KO.tsv".format(minto_dir = minto_dir),
        pathway_map=lambda wildcards: "{minto_dir}/data/kofam_db/KEGG_Pathway2KO.tsv".format(minto_dir = minto_dir)
    output:
        "{wd}/DB/{minto_mode}/4-annotations/kofam/{genome}.kofam.tsv",
    shadow:
        "minimal"
    log:
        "{wd}/logs/DB/{minto_mode}/{genome}.kofamscan.log"
    resources:
        mem=10
    threads: 16
    conda:
        config["minto_dir"]+"/envs/gene_annotation.yml"
    shell:
        """
        remote_dir=$(dirname {output})
        time (
            exec_annotation -k {input.ko_list} -p {input.prok_hal} --tmp-dir tmp -f mapper-one-line --cpu {threads} -o kofam_mapper.txt {input.faa}
            echo -e "#Definitions downloaded\\t$(stat -c '%y' {input.ko_list} | cut -d' ' -f 1)\\t$(stat -c '%y' {input.module_map} | cut -d' ' -f 1)\\t$(stat -c '%y' {input.pathway_map} | cut -d' ' -f 1)" > kofam_processed.txt
            {script_dir}/kofam_hits.pl --pathway-map {input.pathway_map} --module-map {input.module_map} kofam_mapper.txt >>  kofam_processed.txt
            rsync -a  kofam_processed.txt {output}
        ) >& {log}
        """

rule gene_annot_dbcan:
    input:
        faa=rules.rename_prokka_sequences.output.faa,
        dbcan_db="{}/data/dbCAN_db/V12/fam-substrate-mapping.tsv".format(minto_dir)
    output:
        "{wd}/DB/{minto_mode}/4-annotations/dbCAN/{genome}.dbCAN.tsv",
    shadow:
        "minimal"
    log:
        "{wd}/logs/DB/{minto_mode}/{genome}.dbcan.log"
    resources:
        mem=10
    threads: 16
    conda:
        config["minto_dir"]+"/envs/gene_annotation.yml"
    shell:
        """
        time (
            run_dbcan {input.faa} protein --db_dir $(dirname {input.dbcan_db}) --dia_cpu {threads} --out_pre dbcan_ --out_dir out
            echo -e "#Database downloaded\\t$(conda list | sed -E 's|[[:space:]]+| |g' | cut -d' ' -f 1-2 | grep -P dbcan)\\t$(stat -c '%y' {input.dbcan_db} | cut -d' ' -f 1)" > dbcan_processed.txt
            {script_dir}/process_dbcan_overview.pl out/dbcan_overview.txt >> dbcan_processed.txt
            rsync -a dbcan_processed.txt {output}
        ) >& {log}
        """

rule gene_annot_eggnog:
    input:
        faa=rules.rename_prokka_sequences.output.faa,
        eggnog_db="{}/data/eggnog_data/data/eggnog.db".format(minto_dir)
    output:
        "{wd}/DB/{minto_mode}/4-annotations/eggNOG/{genome}.eggNOG.tsv",
    shadow:
        "minimal"
    params:
        eggnog_inmem = lambda wildcards: "--dbmem" if eggNOG_dbmem else ""
    log:
        "{wd}/logs/DB/{minto_mode}/{genome}.eggnog.log"
    resources:
        mem=lambda wildcards: 50 if eggNOG_dbmem else 16
    threads: lambda wildcards: 24 if eggNOG_dbmem else 10
    conda:
        config["minto_dir"]+"/envs/gene_annotation.yml"
    shell:
        """
        time (
            mkdir out
            emapper.py --data_dir $(dirname {input.eggnog_db}) -o tmp \
                       --no_annot --no_file_comments --report_no_hits --override --output_dir out -m diamond -i {input.faa} --cpu {threads}
            emapper.py --annotate_hits_table out/tmp.emapper.seed_orthologs \
                       --data_dir $(dirname {input.eggnog_db}) -m no_search --no_file_comments --override -o tmp --output_dir out --cpu {threads} {params.eggnog_inmem}
            cut -f 1,5,12,13,14,21 out/tmp.emapper.annotations > out/emapper.out
            echo -e "#Database version\\t$(emapper.py --data_dir $(dirname {input.eggnog_db}) -v | cut
 -d"/" -f3 | cut -d" " -f 6)" > out/eggNOG.tsv
            {script_dir}/process_eggNOG_OGs.pl out/emapper.out \
                    | sed 's/\#query/ID/; s/ko\://g' \
                    >> out/eggNOG.tsv
            rsync -a out/eggNOG.tsv {output}
        ) >& {log}
        """

################################################################################################
# Combine gene annotation results into one file
################################################################################################

def get_genome_annotation_tsvs(wildcards):
    #Collect the CDS faa files for MAGs
    inputs = expand("{wd}/DB/{minto_mode}/4-annotations/{annot}/{genome}.{annot}.tsv",
                    wd=wildcards.wd,
                    annot=wildcards.annot,
                    minto_mode=wildcards.minto_mode,
                    genome=genomes)
    return(inputs)

rule combine_annotation_outputs:
    input:
        get_genome_annotation_tsvs
    output:
        "{wd}/DB/{minto_mode}/4-annotations/{annot}.tsv"
    log:
        "{wd}/logs/DB/{minto_mode}/combine_{annot}.log"
    wildcard_constraints:
        annot="eggNOG|kofam|dbCAN"
    shadow:
        "minimal"
    shell:
        """
        time (
            head -n 2 {input[0]} > out.txt
            FOLDER=$(dirname {input[0]})
            tail -n +3 -q $FOLDER/*.tsv >> out.txt
            rsync -a out.txt {output}
        ) >& {log}
        """

rule predicted_gene_annotation_collate:
    input:
        annot_out=expand("{{wd}}/DB/{{minto_mode}}/4-annotations/{annot}.tsv", annot=annot_list)
    output:
        annot_out="{wd}/DB/{minto_mode}/4-annotations/combined_annotations.tsv",
    log:
        "{wd}/logs/DB/{minto_mode}/combine_annot.log"
    resources:
        mem=10
    threads: 1
    conda:
        config["minto_dir"]+"/envs/mags.yml" # python with pandas
    shell:
        """
        time (
            python3 {script_dir}/collate_gene_annotations.py {input} > {output}
        ) >& {log}
        """
