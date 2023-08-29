#!/usr/bin/env python

##############################################
# Clarify where the final QC_2 reads are located
##############################################
def get_qc2_output_location(omics):
    if omics == 'metaT':
        return '5-1-sortmerna'
    elif omics == 'metaG':
        return '4-hostfree'

def get_qc2_output_files_fwd_only(wildcards):
    files = expand("{wd}/{omics}/{location}/{sample}/{run}.{pair}.fq.gz",
                wd = wildcards.wd,
                omics = wildcards.omics,
                location = get_qc2_output_location(wildcards.omics),
                sample = wildcards.sample,
                run = get_runs_for_sample(wildcards),
                pair = '1')
    return(files)

def get_qc2_output_files_rev_only(wildcards):
    files = expand("{wd}/{omics}/{location}/{sample}/{run}.{pair}.fq.gz",
                wd = wildcards.wd,
                omics = wildcards.omics,
                location = get_qc2_output_location(wildcards.omics),
                sample = wildcards.sample,
                run = get_runs_for_sample(wildcards),
                pair = '2')
    return(files)
