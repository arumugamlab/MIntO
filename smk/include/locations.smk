#!/usr/bin/env python

##############################################
# Clarify where the final QC_2 reads are located
##############################################
def get_qc2_output_location(omics):
    if omics == 'metaT':
        return '5-1-sortmerna'
    elif omics == 'metaG':
        return '4-hostfree'

########################################
# Get a sorted list of runs for a sample
########################################
def get_runs_for_sample(wildcards):

    # Sequentially look for runs from the last step to first

    sample = None

    # In assembly, the wildcard value for sample is 'illumina'
    # In QC_2, it is 'sample'
    # Resolve this first
    if hasattr(wildcards, 'illumina'): # In assembly.smk
        sample = wildcards.illumina
    elif hasattr(wildcards, 'sample'): # In QC_2.smk
        sample = wildcards.sample
        # If it is composite sample, just return itself
        if sample in merged_illumina_samples:
            return(sample)
    else:
        raise Exception(f"Wildcard 'illumina' or 'sample' should be defined for me to resolve runs per sample")

    runs = []

    # If corrected runs are already present, take it from there
    # If not, look at QC2 outputs
    # This is to handle special cases where we delete qc2 output after error-correction to save space
    for loc in ['5-corrected-runs', '5-1-sortmerna', '4-hostfree', '3-minlength', '1-trimmed']:
        sample_dir = f"{wildcards.wd}/{wildcards.omics}/{loc}/{sample}"
        if path.exists(sample_dir):
            runs = [ re.sub("\.1\.fq\.gz$", "", path.basename(f)) for f in os.scandir(sample_dir) if f.is_file() and f.name.endswith(".1.fq.gz") ]
            if len(runs) > 0:
                break
            else:
                print(f"WARNING: Cannot find runs for sample={sample} in dir={sample_dir}")

    #print(runs)
    if len(runs) == 0:
        raise Exception(f"Cannot find runs for sample={sample}")
    return(sorted(runs))

########################################
# Get list of fwd reads for this sample
########################################
def get_qc2_output_files_fwd_only(wildcards):
    files = expand("{wd}/{omics}/{location}/{sample}/{run}.{pair}.fq.gz",
                wd = wildcards.wd,
                omics = wildcards.omics,
                location = get_qc2_output_location(wildcards.omics),
                sample = wildcards.sample,
                run = get_runs_for_sample(wildcards),
                pair = '1')
    return(files)

########################################
# Get list of rev reads for this sample
########################################
def get_qc2_output_files_rev_only(wildcards):
    files = expand("{wd}/{omics}/{location}/{sample}/{run}.{pair}.fq.gz",
                wd = wildcards.wd,
                omics = wildcards.omics,
                location = get_qc2_output_location(wildcards.omics),
                sample = wildcards.sample,
                run = get_runs_for_sample(wildcards),
                pair = '2')
    return(files)
