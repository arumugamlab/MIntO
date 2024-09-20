#!/usr/bin/env python

'''
Print software and database versions for each module

Authors: Judit Szarvas
'''

import sys
from os import path

include: 'config_parser.smk'

def get_smk_filename():
    if '-s' in sys.argv:
        i = sys.argv.index('-s')
    elif '--snakefile' in sys.argv:
        i = sys.argv.index('--snakefile')
    smk_name = path.splitext(path.basename(sys.argv[i+1]))[0]
    return smk_name

def get_version_output(smkfilename):
    smkfiles = ["QC_0", "QC_1", "QC_2", "assembly", "binning_preparation", "mags_generation", "gene_annotation", "gene_abundance", "data_integration", "module_tester"]
    if smkfilename in smkfiles:
        result = "{wd}/output/versions/{smk}.flag".format(wd=working_dir, smk=smkfilename)
        return(result)
    else:
        print('WARNING:', smkfilename, 'not set up for printing software versions')
    return()

snakefile_name = get_smk_filename()

metaphlan_version = ""
if "metaphlan_version" in config and config['metaphlan_version'] is not None:
    metaphlan_version = config['metaphlan_version']
motus_version = ""
if "motus_version" in config and config['motus_version'] is not None:
    motus_version = config['motus_version']

spades_script = 'spades.py' # from conda environment
if 'METASPADES_custom_build' in config:
    spades_script = config['METASPADES_custom_build']

ppl_version = ""
gtdb_version = ""
if "PHYLOPHLAN_TAXONOMY_VERSION" in config and config["PHYLOPHLAN_TAXONOMY_VERSION"] is not None:
    ppl_version = config["PHYLOPHLAN_TAXONOMY_VERSION"]
if "GTDB_TAXONOMY_VERSION" in config and config["GTDB_TAXONOMY_VERSION"] is not None:
    gtdb_version = config["GTDB_TAXONOMY_VERSION"]

##################################
#  Rules for TESTING
##################################

rule test_base:
    output:
        temp("{wd}/output/versions/module_tester.flag")
    resources:
        mem=1
    threads:
        1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        echo "MIntO git commit $(cd {minto_dir} && git show --pretty=reference -q && cd - > /dev/null)" > {wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        touch {output}
        """

##################################
#  Rules for QC_0
##################################

rule QC_0_base:
    output:
        temp("{wd}/output/versions/QC_0_base.flag")
    resources:
        mem=1
    threads:
        1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        VOUT={wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        echo "MIntO git commit $(cd {minto_dir} && git show --pretty=reference -q && cd - > /dev/null)" > $VOUT
        fastqc --version >> $VOUT
        echo $(fastp --version 2>&1) >> $VOUT
        touch {output}
        """

rule QC_0_rpkg:
    input:
        "{wd}/output/versions/QC_0_base.flag"
    output:
        temp("{wd}/output/versions/QC_0.flag")
    resources:
        mem=1
    threads: 1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    shell:
        """
        Rscript --version >> {wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        touch {output}
        """

##################################
#  Rules for QC_1
##################################

rule QC_1_base:
    output:
        temp("{wd}/output/versions/QC_1_base.flag")
    resources:
        mem=1
    threads:
        1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        VOUT={wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        echo "MIntO git commit $(cd {minto_dir} && git show --pretty=reference -q && cd - > /dev/null)" > $VOUT
        echo "trimmomatic v$(trimmomatic -version)" >> $VOUT
        touch {output}
        """

rule QC_1_rpkg:
    input:
        "{wd}/output/versions/QC_1_base.flag"
    output:
        temp("{wd}/output/versions/QC_1.flag")
    resources:
        mem=1
    threads: 1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    shell:
        """
        Rscript --version >> {wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        touch {output}
        """

##################################
#  Rules for QC_2
##################################

rule QC_2_base:
    output:
        temp("{wd}/output/versions/QC_2_base.flag")
    resources:
        mem=1
    threads:
        1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        VOUT={wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        echo "MIntO git commit $(cd {minto_dir} && git show --pretty=reference -q && cd - > /dev/null)" > $VOUT
        seqkit version >> $VOUT
        echo "bwa-mem2 v$(bwa-mem2 version 2> /dev/null)" >> $VOUT
        echo "msamtools $(msamtools 2>&1 | grep "Version" | cut -d" " -f 2)" >> $VOUT
        echo "samtools $(samtools 2>&1 | grep "Version" | cut -d" " -f 2)" >> $VOUT
        echo "$(sortmerna 2>&1 | grep "Program" | cut -d" " -f 9-11)" >> $VOUT
        sourmash --version >> $VOUT
        touch {output}
        """

rule QC_2_mpl:
    input:
        "{wd}/output/versions/QC_2_base.flag"
    output:
        temp("{wd}/output/versions/QC_2_mpl.flag")
    resources:
        mem=1 
    threads: 1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/metaphlan.yml"
    shell:
        """
        VOUT={wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        metaphlan --version >> $VOUT
        echo "metaphlan database $(cat {minto_dir}/data/metaphlan/{metaphlan_version}/mpa_latest)" >> $VOUT
        touch {output}
        """

rule QC_2_motus:
    input:
        "{wd}/output/versions/QC_2_mpl.flag"
    output:
        temp("{wd}/output/versions/QC_2_motus.flag")
    resources:
        mem=1
    threads: 1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/motus_env.yml"
    shell:
        """
        (motus -db {minto_dir}/data/motus/{motus_version}/db_mOTU --version || motus --version) 2> /dev/null >> {wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        touch {output}
        """

rule QC_2_rpkg:
    input:
        "{wd}/output/versions/QC_2_motus.flag"
    output:
        temp("{wd}/output/versions/QC_2.flag")
    params:
        fn=snakefile_name
    resources:
        mem=1
    threads: 1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    shell:
        """
        Rscript --version >> {wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        touch {output}
        """

##################################
#  Rules for Assembly
##################################

rule assembly_base:
    output:
        temp("{wd}/output/versions/assembly.flag")
    resources:
        mem=1
    threads:
        1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        VOUT={wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        echo "MIntO git commit $(cd {minto_dir} && git show --pretty=reference -q && cd - > /dev/null)" > $VOUT
        $(which {spades_script}) --version >> $VOUT
        megahit --version >> $VOUT
        echo "flye v$(flye --version)" >> $VOUT
        touch {output}
        """

##################################
#  Rules for binning_preparation
##################################

rule binning_preparation_base:
    output:
        temp("{wd}/output/versions/binprep_base.flag")
    resources:
        mem=1
    threads:
        1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        VOUT={wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        echo "MIntO git commit $(cd {minto_dir} && git show --pretty=reference -q && cd - > /dev/null)" > $VOUT
        echo "bwa-mem2 v$(bwa-mem2 version 2> /dev/null)" >> $VOUT
        echo "msamtools $(msamtools 2>&1 | grep "Version" | cut -d" " -f 2)" >> $VOUT
        echo "samtools $(samtools 2>&1 | grep "Version" | cut -d" " -f 2)" >> $VOUT
        coverm --version >> $VOUT
        touch {output}
        """

rule binning_preparation_avamb:
    input:
        "{wd}/output/versions/binprep_base.flag"
    output:
        temp("{wd}/output/versions/binning_preparation.flag")
    resources:
        mem=1
    threads: 1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/avamb.yml"
    shell:
        """
        vamb --version >> {wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        touch {output}
        """

##################################
#  Rules for mags_generation
##################################

rule mags_base:
    output:
        temp("{wd}/output/versions/mags_base.flag")
    resources:
        mem=1
    threads:
        1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        VOUT={wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        echo "MIntO git commit $(cd {minto_dir} && git show --pretty=reference -q && cd - > /dev/null)" > $VOUT
        coverm --version >> $VOUT
        touch {output}
        """

rule mags_avamb:
    input:
        "{wd}/output/versions/mags_base.flag"
    output:
        temp("{wd}/output/versions/mags_avamb.flag")
    resources:
        mem=1
    threads: 1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/avamb.yml"
    shell:
        """
        vamb --version >> {wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        touch {output}
        """

rule mags_checkm2:
    input:
        "{wd}/output/versions/mags_base.flag"
    output:
        temp("{wd}/output/versions/mags_generation.flag")
    resources:
        mem=1
    threads: 1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/checkm2.yml"
    shell:
        """
        echo "checkm2 v$(checkm2 predict --version)" >> {wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        touch {output}
        """

##################################
#  Rules for gene_annotation
##################################

rule annotation_base:
    output:
        temp("{wd}/output/versions/gene_annotation.flag")
    params:
        p1="|".join(["dbcan"]),
        p2="|".join(["aragorn", "barrnap", "infernal"]),
        p3="|".join(["prodigal", "bedtools", "biopython", "blast", "hmmer", "parallel", "diamond", "mmseqs2"])
    resources:
        mem=1
    threads:
        1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/gene_annotation.yml"
    shell:
        """
        VOUT={wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        echo "MIntO git commit $(cd {minto_dir} && git show --pretty=reference -q && cd - > /dev/null)" > $VOUT
        prokka --version >> $VOUT 2>&1
        echo "kofamscan v$(exec_annotation --version | cut -d" " -f 2)" >> $VOUT
        conda list | sed -E 's|[[:space:]]+| |g' | cut -d" " -f 1-2 | grep -P "{params.p1}" >> $VOUT
        conda list | sed -E 's|[[:space:]]+| |g' | cut -d" " -f 1-2 | grep -P "{params.p2}" >> $VOUT
        conda list | sed -E 's|[[:space:]]+| |g' | cut -d" " -f 1-2 | grep -P "{params.p3}" | grep -v "perl" >> $VOUT

        echo "#----" >> $VOUT
        echo "dbCAN database files" >> $VOUT
        find $(dirname $(which dbcan_build))/.. -iname "dbcan_build.py" -type f -exec grep "dbCAN2/download/Databases" {{}} \; | grep -o -P "fam-substrate-mapping-.*.tsv|dbCAN-PUL_.*.xlsx|V12/CAZyDB.*.fa" >> $VOUT || echo "Error: dbcan_build.py not found"

        echo "#----" >> $VOUT
        emapper.py --data_dir {minto_dir}/data/eggnog_data/data -v | cut -d"/" -f 1,3 >> $VOUT

        echo "#----" >> $VOUT
        echo "KEGG/kofamscan files" >> $VOUT
        echo "KO profiles modification date $(stat -c '%y' {minto_dir}/data/kofam_db/ko_list | cut -d" " -f 1)" >> $VOUT
        echo "Module to KO file creation date $(stat -c '%y' {minto_dir}/data/kofam_db/KEGG_Module2KO.tsv | cut -d" " -f 1)" >> $VOUT
        echo "Pathway to KO file creation date $(stat -c '%y' {minto_dir}/data/kofam_db/KEGG_Pathway2KO.tsv | cut -d" " -f 1)" >> $VOUT

        touch {output}
        """

rule annotation_phylophlan:
    output:
        temp("{wd}/output/versions/annot_phylophlan.flag")
    params:
        db=ppl_version
    resources:
        mem=1
    threads: 1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/mags.yml"
    shell:
        """
        VOUT={wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        phylophlan --version >> $VOUT
        phylophlan_assign_sgbs --version >> $VOUT
        echo "phylophlan database {params.db}" >> $VOUT
        touch {output}
        """

rule annotation_gtdb:
    output:
        temp("{wd}/output/versions/annot_gtdb.flag")
    params:
        db=gtdb_version
    resources:
        mem=1
    threads: 1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/gtdb.yml"
    shell:
        """
        VOUT={wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        echo "gtdbtk v$(gtdbtk --version | cut -d" " -f 3 | head -n 1)" >> $VOUT
        echo "gtdbtk database {params.db}" >> $VOUT
        touch {output}
        """

##################################
#  Rules for gene_abundance
##################################

rule abundance_base:
    output:
        temp("{wd}/output/versions/abund_base.flag")
    resources:
        mem=1
    threads:
        1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """
        VOUT={wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        echo "MIntO git commit $(cd {minto_dir} && git show --pretty=reference -q && cd - > /dev/null)" > $VOUT
        echo "bwa-mem2 v$(bwa-mem2 version 2> /dev/null)" >> $VOUT
        echo "msamtools $(msamtools 2>&1 | grep "Version" | cut -d" " -f 2)" >> $VOUT
        echo "samtools $(samtools 2>&1 | grep "Version" | cut -d" " -f 2)" >> $VOUT
        bedtools --version >> $VOUT
        touch {output}
        """

rule abundance_rpkg:
    input:
        "{wd}/output/versions/abund_base.flag"
    output:
        temp("{wd}/output/versions/gene_abundance.flag")
    resources:
        mem=1
    threads: 1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    shell:
        """
        Rscript --version >> {wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        touch {output}
        """

##################################
#  Rules for data_integration
##################################

rule integration_rpkg:
    output:
        temp("{wd}/output/versions/data_integration.flag")
    resources:
        mem=1
    threads:
        1
    localrule: True
    conda:
        config["minto_dir"]+"/envs/r_pkgs.yml"
    shell:
        """
        VOUT={wildcards.wd}/output/versions/{snakefile_name}.$(date "+%Y-%m-%d").txt
        echo "MIntO git commit $(cd {minto_dir} && git show --pretty=reference -q && cd - > /dev/null)" > $VOUT
        Rscript --version >> $VOUT
        conda list | sed -E 's|[[:space:]]+| |g' | cut -d" " -f 1-2 | grep -P "^r-.*" >> $VOUT
        conda list | sed -E 's|[[:space:]]+| |g' | cut -d" " -f 1-2 | grep -P "^bioconductor-.*" >> $VOUT
        touch {output}
        """
