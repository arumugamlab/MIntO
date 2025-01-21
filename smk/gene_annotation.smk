#!/usr/bin/env python

'''
Gene prediction and functional annotation step

Authors: Vithiagaran Gunalan, Carmen Saenz, Judit Szarvas, Mani Arumugam
'''

import pathlib
import os.path

# Get common config variables
# These are:
#   config_path, project_id, omics, working_dir, minto_dir, script_dir, metadata

NEED_PROJECT_VARIABLES = True
include: 'include/cmdline_validator.smk'
include: 'include/config_parser.smk'

module print_versions:
    snakefile:
        'include/versions.smk'
    config: config

use rule annotation_base from print_versions as version_*

snakefile_name = print_versions.get_smk_filename()

# MIntO mode and database-mapping

# Which mode are we running?
MINTO_MODE = get_minto_mode(config)

# 'catalog' mode is not allowed here
valid_minto_modes = ['MAG', 'refgenome']
check_allowed_values('MINTO_MODE', MINTO_MODE, valid_minto_modes)

if MINTO_MODE == 'MAG':
    mag_omics = 'metaG'
    if (x := validate_optional_key(config, 'MAG_omics')):
        mag_omics = x
    reference_dir = f"{working_dir}/{mag_omics}/8-1-binning/mags_generation_pipeline/unique_genomes"
    print('NOTE: MIntO is using "' + reference_dir + '" as PATH_reference variable')
elif MINTO_MODE == 'refgenome':
    if (x := validate_required_key(config, 'PATH_reference')):
        reference_dir = x
        print('NOTE: MIntO is using "'+ reference_dir+'" as PATH_reference variable')

# Taxonomic annotation of genomes

taxonomies_versioned = list()
taxonomies = list()
run_taxonomy = False
if (x := validate_optional_key(config, 'RUN_TAXONOMY')):
    run_taxonomy = x

if run_taxonomy:
    TAXONOMY_CPUS   = validate_required_key(config, 'TAXONOMY_CPUS')
    TAXONOMY_memory = validate_required_key(config, 'TAXONOMY_memory')

    TAXONOMY = validate_required_key(config, 'TAXONOMY_NAME')
    allowed = ('phylophlan', 'gtdb')
    for x in TAXONOMY.split(","):
        check_allowed_values('TAXONOMY_NAME', x, allowed)

    taxonomies = TAXONOMY.split(",")
    for t in taxonomies:
        version="unknown"
        if t == "phylophlan":
            version = validate_required_key(config, "PHYLOPHLAN_TAXONOMY_VERSION")
        elif t == "gtdb":
            version = validate_required_key(config, "GTDB_TAXONOMY_VERSION")
        taxonomies_versioned.append(t+"."+version)
    print('NOTE: MIntO is running taxonomy labelling of the unique MAGs using [{}].'.format(", ".join(taxonomies_versioned)))

    if "phylophlan" in taxonomies:
        if "gtdb" in taxonomies:
            use rule annotation_base, annotation_phylophlan, annotation_gtdb from print_versions as version_*
        else:
            use rule annotation_base, annotation_phylophlan from print_versions as version_*
    else:
        use rule annotation_base, annotation_gtdb from print_versions as version_*

ANNOTATION = validate_required_key(config, 'ANNOTATION')
for x in ANNOTATION:
    check_allowed_values('ANNOTATION', x, ('dbCAN', 'kofam', 'eggNOG'))

eggNOG_dbmem = True
if 'eggNOG' in ANNOTATION:
    x = validate_optional_key(config, 'eggNOG_dbmem')
    if x is not None:
        eggNOG_dbmem = x

# Define all the outputs needed by target 'all'

GENE_DB_TYPE = MINTO_MODE + '-genes'

# Prepare module completeness information for KEGG-database annotations
# Currently only generated when 'kofam' annotation is available
def module_completeness():
    result = []
    if 'kofam' in ANNOTATION:
        result.append(
                    "{wd}/DB/{subdir}/4-annotations/kofam.per_mag.module_completeness.tsv".format(
                        wd = working_dir,
                        subdir = MINTO_MODE)
                    )
    return(result)

rule all:
    input:
        expand("{wd}/DB/{subdir}/3-taxonomy/taxonomy.{taxonomy}.tsv",
                    wd = working_dir,
                    subdir = MINTO_MODE,
                    taxonomy = taxonomies_versioned),
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
        module_completeness(),
        print_versions.get_version_output(snakefile_name),
        expand("{wd}/output/versions/annot_{taxonomy_method}.flag",
                    wd = working_dir,
                    taxonomy_method = taxonomies)
    default_target: True


######################
# 0. PREPARE GENOMES
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
        "{wd}/DB/{minto_mode}/{minto_mode}.locus_id_list.txt"
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
# 1. GENE PREDICTION
######################

########################
# Prokka on an fna file in reference_dir.
# To make it reproducible, we assign locus_tag based on the name of the MAG/genome.
# This makes testing and benchmarking easier.
########################

rule prokka_for_genome:
    input:
        fna=lambda wildcards: "{reference_dir}/{genome}.fna".format(reference_dir=reference_dir, genome=wildcards.genome),
        locus_ids="{wd}/DB/{minto_mode}/{minto_mode}.locus_id_list.txt"
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
        minto_dir + "/envs/gene_annotation.yml"
    shell:
        """
        time (
            locus_tag=$(grep "^{wildcards.genome}$(printf '\\t')" {input.locus_ids} | cut -f 2 )
            prokka --outdir out --prefix tmp --locustag $locus_tag --addgenes --cdsrnaolap --cpus {threads} --centre X --compliant {input.fna}
            rsync -a out/tmp.fna {output.fna}
            rsync -a out/tmp.faa {output.faa}
            rsync -a out/tmp.gff {output.gff}
            rsync -a out/tmp.gbk {output.gbk}
            rsync -a out/tmp.txt {output.summary}
        ) >& {log}
        """

######################
# 2a. RENAME SEQS
######################

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
    localrule: True
    conda:
        minto_dir + "/envs/MIntO_base.yml" #gff2bed from bedops
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

######################
# 2b. COMBINE BED FILES
######################

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
    localrule: True
    shadow:
        "minimal"
    wildcard_constraints:
        minto_mode = r'MAG|refgenome'
    run:
        import shutil
        process_genome_bed_list(input, 'full.bed', 'mini.bed', log)
        shutil.copy2('full.bed', output.bed_full)
        shutil.copy2('mini.bed', output.bed_mini)

######################
# 2c. MARKER GENES
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
        minto_dir + "/envs/r_pkgs.yml"
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

######################
# 3. TAXONOMY
######################

rule check_directory_for_taxonomic_annotation:
    input:
        rules.merge_MG_tables.output
    output:
        "{wd}/DB/{minto_mode}/2-postprocessed/all.done"
    localrule: True
    resources:
        mem=1
    threads: 1
    shell:
        """
        touch {output}
        """

########################
# PhyloPhlAn on fna files
# Back up the raw output
# Make a standard format output with:
# mag_id,kingdom,phylum,class,order,family,genus,species
# Use mash-distance cutoffs to assign taxonomy resolution:
#       d <  5% - species
#  5% < d < 10% - genus
# 10% < d < 20% - family
########################

rule phylophlan_taxonomy_for_genome_collection:
    input:
        genomes   = rules.check_directory_for_taxonomic_annotation.output,
        phylo_def = f"{minto_dir}/data/phylophlan/{{db_version}}.txt.bz2"
    output:
        orig  = "{wd}/DB/{minto_mode}/3-taxonomy/taxonomy.phylophlan.{db_version}.tsv.orig",
        fixed = "{wd}/DB/{minto_mode}/3-taxonomy/taxonomy.phylophlan.{db_version}.tsv"
    shadow:
        "minimal"
    log:
        "{wd}/logs/DB/{minto_mode}/3-taxonomy/taxonomy.phylophlan.{db_version}.log"
    params:
        db_folder=lambda wildcards: "{minto_dir}/data/phylophlan".format(minto_dir=minto_dir)
    resources:
        mem=TAXONOMY_memory
    threads:
        TAXONOMY_CPUS
    conda:
        minto_dir + "/envs/mags.yml"
    shell:
        """
        time (
            phylophlan_assign_sgbs --input $(dirname {input.genomes}) --input_extension fna --output_prefix taxonomy --nproc {threads} -d {wildcards.db_version} --database_folder {params.db_folder}
            echo -e "mag_id\\tkingdom\\tphylum\\tclass\\torder\\tfamily\\tgenus\\tspecies" > taxonomy.tsv.fixed
            cut -f1,2 taxonomy.tsv \
                    | tail -n +2 \
                    | sed "s/ /\\t/" \
                    | perl -lane '($id, undef, $tax, $dist) = split(/:/, $F[1]); @tax = map {{s/.__//; $_}} split(/\|/, $tax); pop(@tax); $, = "\\t"; print(defined($dist)?$dist:0.5, $F[0], @tax);' \
                    | perl -lane '$dist = shift(@F); $mag = shift(@F); pop(@F) if $dist>0.05; pop(@F) if $dist>0.1; pop(@F) if $dist>0.2; print join("\\t", $mag, @F)' \
                    >> taxonomy.tsv.fixed
        ) >& {log}
        rsync -a taxonomy.tsv {output.orig}
        rsync -a taxonomy.tsv.fixed {output.fixed}
        """

########################
# GTDB-tk on fna files
# Back up the raw output
# Make a standard format output with:
# mag_id,kingdom,phylum,class,order,family,genus,species
########################

rule gtdb_taxonomy_for_genome_collection:
    input:
        genomes  = rules.check_directory_for_taxonomic_annotation.output,
        gtdb_def = f"{minto_dir}/data/GTDB/{{db_version}}/taxonomy/gtdb_taxonomy.tsv"
    output:
        orig  = "{wd}/DB/{minto_mode}/3-taxonomy/taxonomy.gtdb.{db_version}.tsv.orig",
        fixed = "{wd}/DB/{minto_mode}/3-taxonomy/taxonomy.gtdb.{db_version}.tsv"
    shadow:
        "minimal"
    log:
        "{wd}/logs/DB/{minto_mode}/3-taxonomy/taxonomy.gtdb.{db_version}.log"
    params:
        db_folder=lambda wildcards: "{minto_dir}/data/GTDB/{db_version}".format(minto_dir=minto_dir, db_version=wildcards.db_version)
    resources:
        mem=70
    threads:
        TAXONOMY_CPUS
    conda:
        minto_dir + "/envs/gtdb.yml"
    shell:
        """
        time (
            export GTDBTK_DATA_PATH={params.db_folder}
            gtdbtk classify_wf --genome_dir $(dirname {input.genomes}) --extension fna --out_dir tmp --cpus {threads} --skip_ani_screen
            cat tmp/gtdbtk.*.summary.tsv | grep "user_genome" | head -1 > taxonomy.tsv
            cat tmp/gtdbtk.*.summary.tsv | grep -v "user_genome" >> taxonomy.tsv
            echo -e "mag_id\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies" > taxonomy.tsv.fixed
            cut -f1,2 taxonomy.tsv \
                    | tail -n +2 \
                    | sed "s/ /_/g" \
                    | perl -lane '@tax = map {{s/.__//; $_}} split(/;/, $F[1]); $, = "\t"; print($F[0], @tax);' \
                    >> taxonomy.tsv.fixed
        ) >& {log}
        rsync -a taxonomy.tsv {output.orig}
        rsync -a taxonomy.tsv.fixed {output.fixed}
        """

###########################
# 4. FUNCTIONAL ANNOTATION
###########################

###############################
# Prepare batches for annotation
###############################

# Get a list of genome faa files
def get_genome_faa(wildcards):
    #Collect the BED files for MAGs
    result = expand("{wd}/DB/{minto_mode}/2-postprocessed/{genome}.faa",
                    wd=wildcards.wd,
                    minto_mode=wildcards.minto_mode,
                    genome=genomes)
    return(result)

checkpoint prepare_genome_batches:
    input:
        faa = get_genome_faa
    output:
        batches = directory("{wd}/DB/{minto_mode}/4-annotations/batches")
    localrule: True
    threads: 1
    resources:
        mem = 2
    params:
        batch_size = 500
    run:
        import math

        # Make output dir if necessary
        if not os.path.exists(output.batches):
            os.makedirs(output.batches)

        # Get batch size
        n_genomes = len(input.faa)
        n_batches = math.ceil(n_genomes/params.batch_size)

        # Process batches one-by-one
        for batch in range(n_batches):

            # Open the batch file

            ofname = f"{output.batches}/{batch:04d}.faa"
            with open(ofname, 'w') as ofstream:

                # Get the genomes for this batch
                istart = batch*params.batch_size
                iend   = istart + params.batch_size

                # Go through each genome in this batch
                for ifname in input.faa[istart:iend]:
                    with open(ifname, "r") as ifstream:
                        for line in ifstream:
                            ofstream.write(line)


################################################################################################
# Genes annotation using 3 different tools to retrieve functions from:
## kofam
## eggNOG
## dbCAN
################################################################################################

rule gene_annot_kofamscan:
    input:
        faa="{wd}/DB/{minto_mode}/4-annotations/batches/{batch}.faa",
        ko_list=lambda wildcards: "{minto_dir}/data/kofam_db/ko_list".format(minto_dir = minto_dir),
        prok_hal=lambda wildcards: "{minto_dir}/data/kofam_db/profiles/prokaryote.hal".format(minto_dir = minto_dir),
        module_map=lambda wildcards: "{minto_dir}/data/kofam_db/KEGG_Module2KO.tsv".format(minto_dir = minto_dir),
        pathway_map=lambda wildcards: "{minto_dir}/data/kofam_db/KEGG_Pathway2KO.tsv".format(minto_dir = minto_dir)
    output:
        "{wd}/DB/{minto_mode}/4-annotations/kofam/{batch}.tsv",
    shadow:
        "minimal"
    log:
        "{wd}/logs/DB/{minto_mode}/{batch}.kofamscan.log"
    resources:
        mem=10
    threads: 16
    conda:
        minto_dir + "/envs/gene_annotation.yml"
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
        faa="{wd}/DB/{minto_mode}/4-annotations/batches/{batch}.faa",
        dbcan_db="{}/data/dbCAN_db/V12/fam-substrate-mapping.tsv".format(minto_dir)
    output:
        "{wd}/DB/{minto_mode}/4-annotations/dbCAN/{batch}.tsv",
    shadow:
        "minimal"
    log:
        "{wd}/logs/DB/{minto_mode}/{batch}.dbcan.log"
    resources:
        mem=10
    threads: 16
    conda:
        minto_dir + "/envs/gene_annotation.yml"
    shell:
        """
        export PATH="{script_dir:q}:$PATH"
        time (
            run_dbcan {input.faa} protein --db_dir $(dirname {input.dbcan_db}) --dia_cpu {threads} --hmm_cpu {threads} --out_pre dbcan_ --out_dir out
            echo -e "#Database downloaded\\t$(conda list | sed -E 's|[[:space:]]+| |g' | cut -d' ' -f 1-2 | grep -P dbcan)\\t$(stat -c '%y' {input.dbcan_db} | cut -d' ' -f 1)" > dbcan_processed.txt
            {script_dir}/process_dbcan_overview.pl out/dbcan_overview.txt >> dbcan_processed.txt
            rsync -a dbcan_processed.txt {output}
        ) >& {log}
        """

rule gene_annot_eggnog:
    input:
        faa="{wd}/DB/{minto_mode}/4-annotations/batches/{batch}.faa",
        eggnog_db="{}/data/eggnog_data/data/eggnog.db".format(minto_dir)
    output:
        "{wd}/DB/{minto_mode}/4-annotations/eggNOG/{batch}.tsv",
    shadow:
        "minimal"
    params:
        eggnog_inmem = lambda wildcards: "--dbmem" if eggNOG_dbmem else ""
    log:
        "{wd}/logs/DB/{minto_mode}/{batch}.eggnog.log"
    resources:
        mem=lambda wildcards: 50 if eggNOG_dbmem else 16
    threads: lambda wildcards: 24 if eggNOG_dbmem else 10
    conda:
        minto_dir + "/envs/gene_annotation.yml"
    shell:
        """
        time (
            mkdir out
            emapper.py --data_dir $(dirname {input.eggnog_db}) -o tmp \
                       --no_annot --no_file_comments --report_no_hits --override --output_dir out -m diamond -i {input.faa} --cpu {threads}
            emapper.py --annotate_hits_table out/tmp.emapper.seed_orthologs \
                       --data_dir $(dirname {input.eggnog_db}) -m no_search --no_file_comments --override -o tmp --output_dir out --cpu {threads} {params.eggnog_inmem}
            cut -f 1,5,12,13,14,21 out/tmp.emapper.annotations > out/emapper.out
            echo -e "#Database version\\t$(emapper.py --data_dir $(dirname {input.eggnog_db}) -v | cut -d'/' -f3 | cut -d' ' -f 6)" > out/eggNOG.tsv
            {script_dir}/process_eggNOG_OGs.pl out/emapper.out \
                    | sed 's/\#query/ID/; s/ko\://g' \
                    >> out/eggNOG.tsv
            rsync -a out/eggNOG.tsv {output}
        ) >& {log}
        """

###############################
# Get the batches
###############################

def get_annotation_for_genome_batches(wildcards):
    #Collect the genome bins from previous step
    checkpoint_output = checkpoints.prepare_genome_batches.get(**wildcards).output[0]
    result = expand("{wd}/DB/{minto_mode}/4-annotations/{annot}/{batch}.tsv",
                    wd=wildcards.wd,
                    minto_mode=wildcards.minto_mode,
                    annot=wildcards.annot,
                    batch=glob_wildcards(os.path.join(checkpoint_output, '{batch}.faa')).batch)
    return(result)

################################################################################################
# Combine gene annotation results into one file
################################################################################################

rule combine_annotation_batches:
    input:
        get_annotation_for_genome_batches
    output:
        "{wd}/DB/{minto_mode}/4-annotations/{annot}.tsv"
    log:
        "{wd}/logs/DB/{minto_mode}/combine_{annot}.log"
    localrule: True
    wildcard_constraints:
        annot = r'eggNOG|kofam|dbCAN'
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
        annot_out=expand("{{wd}}/DB/{{minto_mode}}/4-annotations/{annot}.tsv", annot=ANNOTATION)
    output:
        annot_out="{wd}/DB/{minto_mode}/4-annotations/combined_annotations.tsv",
    log:
        "{wd}/logs/DB/{minto_mode}/combine_annot.log"
    resources:
        mem=10
    threads: 1
    conda:
        minto_dir + "/envs/mags.yml" # python with pandas
    shell:
        """
        time (
            python3 {script_dir}/collate_gene_annotations.py {input} > {output}
        ) >& {log}
        """

###############################################################################################
# Estimate module completeness
###############################################################################################

# Prepare tsv file with format:
# <MAG>\t<KO1>\t<KO2>...
#
# Currently, only for 'kofam'

rule prepare_kos_per_mag:
    input:
        tsv = "{wd}/DB/{minto_mode}/4-annotations/{annot}.tsv"
    output:
        tsv = temp("{wd}/DB/{minto_mode}/4-annotations/{annot}/kos_per_mag.tsv"),
    log:
        "{wd}/logs/DB/{minto_mode}/combine_{annot}.log"
    wildcard_constraints:
        annot = r'kofam'
    localrule: True
    run:
        mag_kos = {}
        with open(input.tsv, "r") as fp:
            _ = fp.readline()
            if _.startswith("#"):
                _ = fp.readline()
            for line in fp:
                tmp = line.strip().split("\t")
                mag_id = tmp[0].split("|")[-1].split("_")[0]
                kos = tmp[1].split(",")
                if "-" not in kos:
                    if mag_id in mag_kos:
                        mag_kos[mag_id].extend(kos)
                    else:
                        mag_kos[mag_id] = kos
        with open(output.tsv, "w") as of:
            for m in mag_kos.keys():
                l = list(set(mag_kos.get(m)))
                print(m, "\t".join(l), file=of, sep="\t")

# Estimate completeness of each module per MAG
rule mag_completeness:
    input:
        tsv=rules.prepare_kos_per_mag.output.tsv,
        kpc_graph="{}/data/kofam_db/graphs.pkl".format(minto_dir),
        kpc_pathways="{}/data/kofam_db/all_pathways.txt".format(minto_dir),
        kpc_class="{}/data/kofam_db/all_pathways_class.txt".format(minto_dir),
        kpc_names="{}/data/kofam_db/all_pathways_names.txt".format(minto_dir),
    output:
        "{wd}/DB/{minto_mode}/4-annotations/{annot}.per_mag.module_completeness.tsv"
    log:
        "{wd}/logs/DB/{minto_mode}/{annot}.per_mag.module_completeness.log"
    shadow:
        "minimal"
    resources:
        mem = 4
    threads: 1
    conda:
        minto_dir + "/envs/gene_annotation.yml"
    shell:
        """
        time (
            give_completeness --input {input.tsv} \
                    --graphs {input.kpc_graph} \
                    --definitions {input.kpc_pathways} \
                    --names {input.kpc_names} \
                    --classes {input.kpc_class} \
                    --outprefix per_mag \
                    --add-per-contig
            sed -e 's|^contig\\b|ID|' per_mag_contigs.tsv > {output}
        ) >& {log}
        """
