#!/usr/bin/env python

import os.path

##############################################
# Clarify where the final inputs are
# NOTE: Watch the difference between list.extend(list) and list.append(scalar)
##############################################
def get_ordered_fastq_dirs(omics, caller, reverse=False):
    # These are available for everyone
    locations = ['1-trimmed']

    # These are available for QC_2 onwards
    if caller not in ['QC_1']:
        locations.extend(['3-minlength', '4-hostfree'])
        if omics == 'metaT':
            locations.append('5-1-sortmerna')

    # error-corrected reads are available strictly after QC_2
    if caller not in ['QC_1', 'QC_2']:
        locations.extend(['5-corrected-runs', '6-corrected'])

    # Reverse if needed
    if reverse:
        locations.reverse()

    return(locations)

##############################################
# Clarify where the final QC_2 reads are located
##############################################
def get_qc2_output_location(omics):
    locations = get_ordered_fastq_dirs(omics, caller='QC_2', reverse=True)
    return(locations[0])

########################################
# Get a sorted list of runs for a sample
########################################
def get_runs_for_sample(wd, omics, sample, caller):

    # If it is a merged_illumina_sample, then it will not exist until the last step.
    # Just return itself so that it is created after last step is done.
    if 'MERGE_ILLUMINA_SAMPLES' in config and config['MERGE_ILLUMINA_SAMPLES'] is not None:
        if sample in config['MERGE_ILLUMINA_SAMPLES']:
            return(sample)

    # Sequentially look for runs from the last step to first

    runs = []

    # Look for individual runs in the output directories of each stage, from earliest-to-latest.
    # Sometimes, pipeline errors-out and stops in the middle, which means later stages may have
    # some runs missing. If we rerun the script, then sequentially looking from latest-to-earliest
    # will wrongly use the 'partial' list of runs. This will skip some runs without us realizing that things went wrong.
    #
    # If we delete earlier stages to save space, then this will still work, since run-names are
    # propagated as-is through all the directories listed in 'locations'.
    # If user wants to skip earlier steps and start with later step, it will also work.

    locations = get_ordered_fastq_dirs(omics, caller, reverse=False)

    # Now look for them
    for loc in locations:
        sample_dir = f"{wd}/{omics}/{loc}/{sample}"
        if os.path.exists(sample_dir):
            runs = list()
            for f in os.scandir(sample_dir):
                if f.is_file():
                    if f.name.endswith('.1.fq.gz') or f.name.endswith('_1.fq.gz'):
                        runs.append(re.sub("[\._]1\.fq\.gz", "", os.path.basename(f)))
            runs = sorted(runs)
            if len(runs) > 0:
                break
            else:
                print(f"WARNING: Cannot find runs for sample={sample} in dir={sample_dir}")

    #print(runs)
    if len(runs) == 0:
        raise Exception(f"Cannot find runs for sample={sample}")
    return(runs)

########################################
# Get list of one-end reads for this sample
########################################
def get_final_fastq_one_end(wd, omics, sample, stage, pair):

    files = list()

    locations = list()
    if stage == 'QC_1':
        # QC_1 has only one output
        locations.append('1-trimmed')
    elif stage == 'QC_2':
        # Check QC2 output but if that fails 6-corrected is admissible.
        locations.append(get_qc2_output_location(omics))
        locations.append('6-corrected')

    for rdir in locations:

        # Get runs
        for r in get_runs_for_sample(wd, omics, sample, stage):
            # Check for sample.1.fq.gz and sample_1.fq.gz, in that order
            for delim in ['.', '_']:
                f = f"{wd}/{omics}/{rdir}/{sample}/{r}{delim}{pair}.fq.gz"
                # Add the file to the list
                if os.path.exists(f):
                    files.append(f)

    if len(files) == 0:
        raise Exception(f"Final fastq files not found for sample={sample}")

    # Return the final list
    return(files)
