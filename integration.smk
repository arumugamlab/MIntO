#!/usr/bin/env python

'''
Gene and function profiling step

Authors: Carmen Saenz
'''

# snakemake --snakefile /emc/cbmr/users/rzv923/MIntO/integration.smk --cluster "qsub -pe smp {threads} -l h_vmem={resources.mem}G -N {name} -cwd" --latency-wait 60 --jobs 2 --configfile data_integration.yaml --use-conda --conda-prefix /emc/cbmr/users/rzv923/ibdmdb/tmp_porus/
# snakemake --snakefile /emc/cbmr/users/rzv923/MIntO/integration.smk --restart-times 0 --keep-going --latency-wait 30 --cluster "sbatch -J {name} --mem={resources.mem}G -c {threads} -e slurm-%x.e%A -o slurm-%x.o%A" --configfile integration.smk.yaml --use-conda --conda-prefix /emc/cbmr/users/rzv923/ibdmdb/tmp_ironman/ --jobs 8
# snakemake --snakefile /emc/cbmr/users/rzv923/MIntO/integration.smk --restart-times 0 --keep-going --latency-wait 30 --cluster "sbatch -J {name} --mem={resources.mem}G -c {threads} -e slurm-%x.e%A -o slurm-%x.o%A" --jobs 8 --configfile integration.smk.yaml --use-singularity

# Make list of illumina samples, if ILLUMINA in config
ilmn_samples = list()
if 'ILLUMINA' in config:
    print("Samples:")
    for ilmn in config["ILLUMINA"]:
            print(ilmn)
            ilmn_samples.append(ilmn)

local_dir = config['local_dir']
working_dir = config['working_dir']
omics = config['omics']

map_reference=config["map_reference"]
minto_dir=config["minto_dir"]
script_dir=config["minto_dir"]+"/scripts"
#script_dir=config["script_dir"]

normalization=config['abundance_normalization']
identity=config['alignment_identity']
metadata_file=config['METADATA']

if omics == 'metaG':
    omics_opt='metaG'
    omics_prof='A'
    #gene_abund="{wd}/{omics_opt}/6-mapping-profiles/BWA_reads-{post_analysis_out}/genes_abundances.p{identity}.{normalization}.csv".format(wd = working_dir, omics_opt = omics_opt, post_analysis_out = post_analysis_out, identity = identity, normalization = normalization),
    #gene_abund_file="metaG/6-mapping-profiles/BWA_reads-{map_reference}/genes_abundances.p{identity}.{normalization}.csv"
elif omics == 'metaT':
    omics_opt='metaT'
    omics_prof='T'
    #gene_abund_file="metaT/6-mapping-profiles/BWA_reads-{map_reference}/genes_abundances.p{identity}.{normalization}.csv"
elif omics == 'metaG_metaT':
    omics_opt=['metaG','metaT']
    omics_prof=['A','T','E']
    #gene_abund_file="metaG/6-mapping-profiles/BWA_reads-{map_reference}/genes_abundances.p{identity}.{normalization}.csv",
    #"metaT/6-mapping-profiles/BWA_reads-{map_reference}/genes_abundances.p{identity}.{normalization}.csv"

if map_reference == 'MAG':
    post_analysis_dir="9-MAGs-prokka-post-analysis"
    post_analysis_out="MAGs_genes"
    post_analysis_genome="MAGs_genes"
    annot_file="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_genome}_translated_cds_SUBSET.annotations.tsv".format(wd = working_dir,post_analysis_dir = post_analysis_dir, post_analysis_genome = post_analysis_genome)
    funct_opt=('eggNOG.OGs','KEGG_Pathway','KEGG_Module','merged_KO','PFAMs','dbCAN.mod','dbCAN.enzclass')
    funct_opt_list = ','.join(['"' + id + '"' for id in ('eggNOG OGs','KEGG_Pathway','KEGG_Module','merged_KO','PFAMs','dbCAN.mod','dbCAN.enzclass')])
elif map_reference == 'reference_genome':
    post_analysis_dir="9-reference-genes-post-analysis"
    post_analysis_out="reference_genes"
    post_analysis_genome="reference_genes"
    annot_file="{wd}/DB/{post_analysis_dir}/annot/{post_analysis_genome}_translated_cds_SUBSET.annotations.tsv".format(wd = working_dir,post_analysis_dir = post_analysis_dir, post_analysis_genome = post_analysis_genome)
    funct_opt=('eggNOG.OGs','KEGG_Pathway','KEGG_Module','merged_KO','PFAMs','dbCAN.mod','dbCAN.enzclass')
    funct_opt_list = ','.join(['"' + id + '"' for id in ('eggNOG OGs','KEGG_Pathway','KEGG_Module','merged_KO','PFAMs','dbCAN.mod','dbCAN.enzclass')])
elif map_reference == 'genes_db':
    post_analysis_dir="9-genes-db-post-analysis"
    post_analysis_out="db_genes"
    post_analysis_genome="None"
    annot_file=config['ANNOTATION_file']
    funct_opt=config['ANNOTATION_ids']
    funct_opt_list = ','.join(['"' + id + '"' for id in funct_opt])

if normalization == 'TPM' and (map_reference == 'MAG' or map_reference == 'reference_genome'):
    post_analysis_TPM=post_analysis_out
    post_analysis_other='None'
else:
    post_analysis_TPM="None"
    post_analysis_other=post_analysis_out


#for id in funct_opt:
#            print(id)
#            funct_opt_list.append(id)
print(funct_opt_list)
# if normalization == 'MG' and (map_reference == 'MAG' or map_reference == 'reference_genome'):
#     post_analysis_MG=post_analysis_out
# else:
#     post_analysis_MG="None"

def integration_merge_profiles():
    result = expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}.csv", 
            wd = working_dir,
            omics = omics,
            post_analysis_out = post_analysis_out,
            identity = identity,
            normalization = normalization)
    return(result)

def integration_gene_profiles():
    result = expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/G{omics_prof}.csv", 
            wd = working_dir,
            omics = omics,
            post_analysis_out = post_analysis_out,
            identity = identity,
            normalization = normalization,
            omics_prof = omics_prof),\
    expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/G{omics_prof}.rds", 
            wd = working_dir,
            omics = omics,
            post_analysis_out = post_analysis_out,
            identity = identity,
            normalization = normalization,
            omics_prof = omics_prof),\
    expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/plots/G{omics_prof}.PCoA-bray.pdf", 
            wd = working_dir,
            omics = omics,
            post_analysis_out = post_analysis_out,
            identity = identity,
            normalization = normalization,
            omics_prof = omics_prof)         
    return(result)

def integration_function_profiles():
    result = expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/F{omics_prof}.{funct_opt}.csv", 
            wd = working_dir,
            omics = omics,
            post_analysis_out = post_analysis_out,
            identity = identity,
            normalization = normalization,
            omics_prof = omics_prof,
            funct_opt = funct_opt),\
    expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/F{omics_prof}.{funct_opt}.rds", 
            wd = working_dir,
            omics = omics,
            post_analysis_out = post_analysis_out,
            identity = identity,
            normalization = normalization,
            omics_prof = omics_prof,
            funct_opt = funct_opt),\
    expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/plots/F{omics_prof}.{funct_opt}.PCoA-bray.pdf", 
            wd = working_dir,
            omics = omics,
            post_analysis_out = post_analysis_out,
            identity = identity,
            normalization = normalization,
            omics_prof = omics_prof,
            funct_opt = funct_opt)         
    return(result)

# Define all the outputs needed by target 'all'
rule all:
    input: 
        integration_merge_profiles(),
        integration_gene_profiles(),
        integration_function_profiles()


###############################################################################################
# Prepare gene profile
## Merge gene abundance and gene transcript profiles
###############################################################################################
rule integration_merge_profiles:
    input:
        gene_abund=expand("{wd}/{omics_opt}/6-mapping-profiles/BWA_reads-{post_analysis_out}/genes_abundances.p{identity}.{normalization}.csv", wd = working_dir, omics_opt = omics_opt, post_analysis_out = post_analysis_out, identity = identity, normalization = normalization),
    output: 
        gene_abund_merge="{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}.csv".format(wd = working_dir, omics = omics, post_analysis_out = post_analysis_out, identity = identity, normalization = normalization),
    log:
        "{wd}/logs/output/data_integration/{post_analysis_out}/integration_merge_profiles.log".format(wd = working_dir, post_analysis_out = post_analysis_out),
    resources:
        mem=config["MERGE_memory"]
    threads: config["MERGE_threads"]
    conda:
        config["minto_dir"]+"/envs/MIntO_base.yml"
    shell:
        """ 
        remote_dir=$(dirname {output.gene_abund_merge})
        time (if [ {omics} == 'metaG_metaT' ]
        then python3 {script_dir}/gene_abundances_merge_profiles.py {threads} {resources.mem} {working_dir} {omics} {post_analysis_out} {normalization} {identity}
        else 
        rsync {working_dir}/{omics}/6-mapping-profiles/BWA_reads-{post_analysis_out}/genes_abundances.p{identity}.{normalization}.csv ${{remote_dir}}/{omics}.genes_abundances.p{identity}.{normalization}.csv 
        fi ) &> {log}
        """
###############################################################################################
# Generate gene expression profile
###############################################################################################
# Change names to R scripts
# what do to with annot file
rule integration_gene_profiles:
    input:
        gene_abund_merge="{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}.csv".format(wd = working_dir, omics = omics, post_analysis_out = post_analysis_out, identity = identity, normalization = normalization),
        #gene_abund="{input_dir}/".format(input_dir=gene_abund_file),
        #gene_abund="{wd}/{omics}/6-mapping-profiles/BWA_reads-{map_reference}/genes_abundances.p{identity}.{normalization}.csv",
    output: 
        gene_abund_prof=expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/G{omics_prof}.csv", wd = working_dir, omics = omics, post_analysis_out = post_analysis_out, identity = identity, normalization = normalization, omics_prof = omics_prof),
        gene_abund_phyloseq=expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/G{omics_prof}.rds", wd = working_dir, omics = omics, post_analysis_out = post_analysis_out, identity = identity, normalization = normalization, omics_prof = omics_prof),
        gene_abund_plots=expand("{wd}/output/data_integration/{post_analysis_out}/{omics}.genes_abundances.p{identity}.{normalization}/plots/G{omics_prof}.PCoA-bray.pdf", wd = working_dir, omics = omics, post_analysis_out = post_analysis_out, identity = identity, normalization = normalization, omics_prof = omics_prof),
    params:
        annot_file={annot_file},
        funct_opt= funct_opt_list
    log:
        "{wd}/logs/output/data_integration/{post_analysis_out}/integration_gene_profiles.{omics}.{normalization}.log".format(wd = working_dir, post_analysis_out = post_analysis_out, normalization = normalization, omics = omics),
    resources:
        mem=config["MERGE_memory"]
    threads: config["MERGE_threads"]
    conda:
        config["minto_dir"]+"/envs/taxa_env.yml"
    shell:
        """ 
        time ( echo 'integration of gene profiles'
        if [[ {map_reference} == 'MAG' ]] || [[ {map_reference} == 'reference_genome' ]]
        then Rscript {script_dir}/gene_expression_profile_genome_based.R {threads} {resources.mem} {working_dir} {omics} {post_analysis_out} {normalization} {identity} {params.annot_file} {metadata_file} {input.gene_abund_merge} {params.funct_opt}
        elif [[ {map_reference} == 'genes_db' ]]
        then Rscript {script_dir}/gene_expression_profile_gene_based.R {threads} {resources.mem} {working_dir} {omics} {post_analysis_out} {normalization} {identity} {params.annot_file} {metadata_file} {input.gene_abund_merge} {params.funct_opt}
        fi ) &> {log}
        """

##################################################################################################
# Generate function expression profile for the different databases included in the annotation file
##################################################################################################
# Change names to R scripts
# what do to with annot file
# check that python3 {script_dir}/genes_abundances_merge_profiles.py works in taxa_env 
# Add options when individual or both omics
#if normalization == 'TPM' and (map_reference == 'MAG' or map_reference == 'reference_genome'):
    ## MAG or ref genome - TPM
    # Make exception to run genes_abundances_merge_profiles.py for raw counts
    # Reduce arguments in genes_abundances_merge_profiles.py
    # change name function_expression_profile_MAG.R 
    # modify def
rule integration_function_profiles_TPM:
    input:
        absolute_counts=expand("{wd}/{omics_opt}/6-mapping-profiles/BWA_reads-{post_analysis_TPM}/genes_abundances.p{identity}.bed", wd = working_dir, omics_opt = omics_opt, post_analysis_TPM = post_analysis_TPM, identity = identity),
        gene_abund_phyloseq=expand("{wd}/output/data_integration/{post_analysis_TPM}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/G{omics_prof}.rds", wd = working_dir, omics = omics, post_analysis_TPM = post_analysis_TPM, identity = identity, omics_prof = omics_prof, normalization = normalization),
    output: 
        absolute_counts_merge="{wd}/output/data_integration/{post_analysis_TPM}/{omics}.genes_abundances.p{identity}.bed".format(wd = working_dir, omics = omics, post_analysis_TPM = post_analysis_TPM, identity = identity),
        func_abund_prof=expand("{wd}/output/data_integration/{post_analysis_TPM}/{omics}.genes_abundances.p{identity}.TPM/F{omics_prof}.{funct_opt}.csv", wd = working_dir, omics = omics, post_analysis_TPM = post_analysis_TPM, identity = identity, omics_prof = omics_prof, funct_opt = funct_opt),
        func_abund_phyloseq=expand("{wd}/output/data_integration/{post_analysis_TPM}/{omics}.genes_abundances.p{identity}.TPM/phyloseq_obj/F{omics_prof}.{funct_opt}.rds", wd = working_dir, omics = omics, post_analysis_TPM = post_analysis_TPM, identity = identity, omics_prof = omics_prof, funct_opt = funct_opt),
        func_abund_plots=expand("{wd}/output/data_integration/{post_analysis_TPM}/{omics}.genes_abundances.p{identity}.TPM/plots/F{omics_prof}.{funct_opt}.PCoA-bray.pdf", wd = working_dir, omics = omics, post_analysis_TPM = post_analysis_TPM, identity = identity, omics_prof = omics_prof, funct_opt = funct_opt),
    params:
        annot_file={annot_file},
        funct_opt=funct_opt_list
    log:
        "{wd}/logs/output/data_integration/{post_analysis_TPM}/integration_function_profiles.{omics}.TPM.log".format(wd = working_dir, post_analysis_TPM = post_analysis_TPM, omics = omics),
    resources:
        mem=config["MERGE_memory"]
    threads: config["MERGE_threads"]
    conda:
        config["minto_dir"]+"/envs/taxa_env.yml" #R
    shell:
        """ 
        time (python3 {script_dir}/gene_abundances_merge_raw_profiles.py {threads} {resources.mem} {working_dir} {omics} {post_analysis_TPM} {normalization} {identity}
        Rscript {script_dir}/function_expression_profile_genome_TPM_based.R {threads} {resources.mem} {working_dir} {omics} {post_analysis_TPM} {normalization} {identity} {params.annot_file} {metadata_file} {minto_dir} {params.funct_opt}) &> {log}
        """

# what do to with annot file
# Add options when individual or both omics and which norm
#         ## MAG or ref genome - MAG
#         # Write script to generate function_expression_profile when MG norm
#         # modify def
#     ## gene_db - TPM
#     # Write script to generate function_expression_profile when TPM-db_genes norm with unique metaG or metaT
#     # modify def
# When no TPM + MAG/ref genomes
# include options to run in shell
rule integration_function_profiles_MG_TPM:
        input:
            gene_abund_phyloseq=expand("{wd}/output/data_integration/{post_analysis_other}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/G{omics_prof}.rds", wd = working_dir, omics = omics, post_analysis_other = post_analysis_other, identity = identity, omics_prof = omics_prof, normalization = normalization),
        output: 
            func_abund_prof=expand("{wd}/output/data_integration/{post_analysis_other}/{omics}.genes_abundances.p{identity}.{normalization}/F{omics_prof}.{funct_opt}.csv", wd = working_dir, omics = omics, post_analysis_other = post_analysis_other, identity = identity, omics_prof = omics_prof, normalization = normalization, funct_opt = funct_opt),
            func_abund_phyloseq=expand("{wd}/output/data_integration/{post_analysis_other}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/F{omics_prof}.{funct_opt}.rds", wd = working_dir, omics = omics, post_analysis_other = post_analysis_other, identity = identity, omics_prof = omics_prof, normalization = normalization, funct_opt = funct_opt),
            func_abund_plots=expand("{wd}/output/data_integration/{post_analysis_other}/{omics}.genes_abundances.p{identity}.{normalization}/plots/F{omics_prof}.{funct_opt}.PCoA-bray.pdf", wd = working_dir, omics = omics, post_analysis_other = post_analysis_other, identity = identity, omics_prof = omics_prof, normalization = normalization, funct_opt = funct_opt),
        params:
            annot_file={annot_file},
            funct_opt=funct_opt_list
        log:
            "{wd}/logs/output/data_integration/{post_analysis_other}/integration_funtion_profiles.{omics}.{normalization}.log".format(wd = working_dir, post_analysis_other = post_analysis_other, normalization = normalization, omics = omics),
        resources:
            mem=config["MERGE_memory"]
        threads: config["MERGE_threads"]
        conda:
            config["minto_dir"]+"/envs/taxa_env.yml" #R
        shell:
            """ 
            time ( echo 'integration of function profiles'
            if ( [[ {map_reference} == 'MAG' ]] || [[ {map_reference} == 'reference_genome' ]] ) && [[ {normalization} == 'MG' ]]
            then Rscript {script_dir}/function_expression_profile_genome_MG_based.R {threads} {resources.mem} {working_dir} {omics} {post_analysis_other} {normalization} {identity} {params.annot_file} {metadata_file} {minto_dir} {params.funct_opt}
            elif [[ {map_reference} == 'genes_db' ]] && [[ {normalization} == 'TPM' ]]
            then Rscript {script_dir}/function_expression_profile_gene_TPM_based.R {threads} {resources.mem} {working_dir} {omics} {post_analysis_other} {normalization} {identity} {params.annot_file} {metadata_file} {minto_dir} {params.funct_opt}
            fi) &> {log}
            """
# #if normalization == 'TPM' and (map_reference == 'MAG' or map_reference == 'reference_genome'):
#         ## MAG or ref genome - MAG
#         # Write script to generate function_expression_profile when MG norm
#         # modify def
# rule integration_function_profiles_MG_TPM:
#         input:
#             gene_abund_phyloseq=expand("{wd}/output/data_integration/{post_analysis_MG}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/G{omics_prof}.rds", wd = working_dir, omics = omics, post_analysis_MG = post_analysis_MG, identity = identity, omics_prof = omics_prof, normalization = normalization),
#         output: 
#             func_abund_prof=expand("{wd}/output/data_integration/{post_analysis_MG}/{omics}.genes_abundances.p{identity}.{normalization}/F{omics_prof}.{funct_opt}.csv", wd = working_dir, omics = omics, post_analysis_MG = post_analysis_MG, identity = identity, omics_prof = omics_prof, normalization = normalization, funct_opt = funct_opt),
#             func_abund_phyloseq=expand("{wd}/output/data_integration/{post_analysis_MG}/{omics}.genes_abundances.p{identity}.{normalization}/phyloseq_obj/F{omics_prof}.{funct_opt}.rds", wd = working_dir, omics = omics, post_analysis_MG = post_analysis_MG, identity = identity, omics_prof = omics_prof, normalization = normalization, funct_opt = funct_opt),
#             func_abund_plots=expand("{wd}/output/data_integration/{post_analysis_MG}/{omics}.genes_abundances.p{identity}.{normalization}/plots/F{omics_prof}.{funct_opt}.PCoA-bray.pdf", wd = working_dir, omics = omics, post_analysis_MG = post_analysis_MG, identity = identity, omics_prof = omics_prof, normalization = normalization, funct_opt = funct_opt),
#         params:
#             dir_annot=lambda wildcards: "{wd}/DB/{post_analysis_dir}/annot/reference_genes_translated_cds_SUBSET.annotations.tsv".format(wd = working_dir, post_analysis_dir=post_analysis_dir)
#         log:
#             "{wd}/logs/output/data_integration/{post_analysis_MG}/integration_gene_profiles.log".format(wd = working_dir, post_analysis_MG = post_analysis_MG),
#         resources:
#             mem=5
#         threads: 8
#         conda:
#             config["minto_dir"]+"/envs/taxa_env.yml" #R
#         shell:
#             """ 
#             time (Rscript {script_dir}/function_expression_profile_MAG.R {threads} {resources.mem} {working_dir} {omics} {post_analysis_MG} {normalization} {identity} {params.dir_annot} {metadata_file}) &> {log}
#             """
# rule integration_function_profiles_genes_db:
#     input:
#         gene_abund_phyloseq=expand("{{wd}}/output/data_integration/{{post_analysis_out}}/{{omics}}.genes_abundances.p{{identity}}.{{normalization}}/phyloseq_obj/G{omics_prof}.rds".format(omics_prof = omics_prof),
#     output: 
#         func_abund_prof=expand("{{wd}}/output/data_integration/{{post_analysis_out}}/{{omics}}.genes_abundances.p{{identity}}.{{normalization}}/F{omics_prof}.{funct_opt}.csv".format(omics_prof = omics_prof, funct_opt = funct_opt),
#         func_abund_phyloseq=expand("{{wd}}/output/data_integration/{{post_analysis_out}}/{{omics}}.genes_abundances.p{{identity}}.{{normalization}}/phyloseq_obj/F{omics_prof}.{funct_opt}.rds".format(omics_prof = omics_prof, funct_opt = funct_opt),
#         func_abund_plots=expand("{{wd}}/output/data_integration/{{post_analysis_out}}/{{omics}}.genes_abundances.p{{identity}}.{{normalization}}/plots/F{omics_prof}.{funct_opt}.PCoA-bray.pdf".format(omics_prof = omics_prof, funct_opt = funct_opt),
#     log:
#         "{wd}/logs/output/data_integration/{post_analysis_out}/integration_gene_profiles.log"#.format(wd = working_dir, post_analysis_out = post_analysis_out),
#     resources:
#         mem=5
#     threads: 8
#     conda:
#         config["minto_dir"]+"/envs/taxa_env_v2.yml"
#     shell:
#     ## gene_db - TPM
#     # Write script to generate function_expression_profile when TPM-db_genes norm with unique metaG or metaT
#     # modify def
#         """ 
#         time (Rscript {script_dir}/function_expression_profile_db_genes.R {threads} {resources.mem} {working_dir} {omics} {post_analysis_out} {normalization} {identity} {dir_annot} {metadata_file}) &> {log}
#         """