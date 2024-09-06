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
def get_runs_for_sample(wd, omics, sample):

    # If it is a merged_illumina_sample, then it will not exist until the last step.
    # Just return itself so that it is created after last step is done.
    if 'MERGE_ILLUMINA_SAMPLES' in config and config['MERGE_ILLUMINA_SAMPLES'] is not None:
        if sample in config['MERGE_ILLUMINA_SAMPLES']:
            return(sample)

    # Sequentially look for runs from the last step to first

    runs = []

    # If corrected runs are already present, take it from there
    # If not, look at QC2 outputs
    # This is to handle special cases where we delete qc2 output after error-correction to save space
    for loc in ['5-corrected-runs', '5-1-sortmerna', '4-hostfree', '3-minlength', '1-trimmed']:
        sample_dir = f"{wd}/{omics}/{loc}/{sample}"
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
# Get list of one-end reads for this sample
########################################
def get_qc2_output_files_one_end(wd, omics, sample, pair):
    files = expand("{wd}/{omics}/{location}/{sample}/{run}.{pair}.fq.gz",
                wd = wd,
                omics = omics,
                location = get_qc2_output_location(omics),
                sample = sample,
                run = get_runs_for_sample(wd, omics, sample),
                pair = pair)
    return(files)

########################################
# Get list of fwd reads for this sample
########################################
def get_qc2_output_files_fwd_only(wildcards):
    return(get_qc2_output_files_one_end(wildcards.wd, wildcards.omics, wildcards.sample, '1'))

########################################
# Get list of rev reads for this sample
########################################
def get_qc2_output_files_rev_only(wildcards):
    return(get_qc2_output_files_one_end(wildcards.wd, wildcards.omics, wildcards.sample, '2'))
