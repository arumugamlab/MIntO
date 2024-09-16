#!/usr/bin/env python

import os

################################################################################################
# Wrapper to create BWA index for any given fasta file, with extensions .fna or .fasta.
# If CLUSTER_NODES is defined, then index files should be distributed to all nodes to
# the local directory CLUSTER_LOCAL_DIR, using jobs submitted to CLUSTER_WORKLOAD_MANAGER.
#
# Steps involved:
# ---------------
# 1. Index generation:
#    input:  {somewhere}/{something}.fna
#    output: {somewhere}/BWA_index/{something}.fna.*
#    (All binary files from bwa-mem2 index are made.)
# 2. Sync'ing to clusters (staging):
#    input:  {somewhere}/BWA_index/{something}.fna.*
#    output: CLUSTER_LOCAL_DIR/{somewhere}/BWA_index/{something}.fna.*
# 3. Symlinking to local drive:
#    input:  CLUSTER_LOCAL_DIR/{somewhere}/BWA_index/{something}.fna.*
#    output: {somewhere}/BWA_index_local/{something}.fna.*
#
# Usage:
# ------
# When using this wrapper, a rule should figure out whether it will be executed locally or on a cluster.
# Then it should ask for bwaindex files as follows:
#   {somewhere}/BWA_index/{something}.fna.*:       if the files can be read directly without staging.
#   {somewhere}/BWA_index_local/{something}.fna.*: if the files need to be staged to all nodes in cluster.
#
# See rule 'genome_mapping_profiling' in gene_abundance.smk for usage example.
################################################################################################

if CLUSTER_WORKLOAD_MANAGER is not None:
    known_managers = ['SLURM', 'PBS', 'TORQUE', 'SGE', 'LSF']
    if CLUSTER_WORKLOAD_MANAGER.upper() not in known_managers:
        raise Exception(f"MIntO error: Unexpected value: CLUSTER_WORKLOAD_MANAGER='{CLUSTER_WORKLOAD_MANAGER}'. Must be one of {known_managers}")

if CLUSTER_NODES is not None:

    # Make sure other variables are defined
    if CLUSTER_WORKLOAD_MANAGER is None or CLUSTER_LOCAL_DIR is None:
        raise Exception(f"MIntO error: CLUSTER_WORKLOAD_MANAGER and CLUSTER_LOCAL_DIR must be defined, if CLUSTER_NODES is defined")

    # Baseline    : 5 GB
    # Size-based  : 22 byte per byte file size
    # New attempts: +40 GB each time
    rule BWA_index_contigs_nfs:
        input:
            fasta="{somewhere}/{something}.{fasta}"
        output:
            "{somewhere}/BWA_index/{something}.{fasta}.0123",
            "{somewhere}/BWA_index/{something}.{fasta}.amb",
            "{somewhere}/BWA_index/{something}.{fasta}.ann",
            "{somewhere}/BWA_index/{something}.{fasta}.bwt.2bit.64",
            "{somewhere}/BWA_index/{something}.{fasta}.pac",
        log:
            "{somewhere}/BWA_index/BWA_index.{something}.{fasta}.log"
        wildcard_constraints:
            fasta='fasta|fna'
        resources:
            mem = lambda wildcards, input, attempt: 5 + int(22*input.size_mb/1024) + 40*(attempt-1),
        threads: 4
        conda:
            config["minto_dir"]+"/envs/MIntO_base.yml" #bwa-mem2
        shell:
            """
            outdir=$(dirname {output[0]})
            umask 002 && mkdir -p {CLUSTER_LOCAL_DIR}/$outdir
            cd {CLUSTER_LOCAL_DIR}/$outdir
            prefix=$(basename {input})
            time (
                bwa-mem2 index {input.fasta} -p $prefix
                echo "NODE: $(hostname)"
                echo "SENDING:"
                echo "$(pwd)/$prefix --> $outdir"
                rsync -av --itemize-changes $prefix.* $outdir/
            ) >& {log}
            """

    checkpoint prepare_files_for_syncing:
        input:
            "{somewhere}/BWA_index/{something}.{fasta}.0123",
            "{somewhere}/BWA_index/{something}.{fasta}.amb",
            "{somewhere}/BWA_index/{something}.{fasta}.ann",
            "{somewhere}/BWA_index/{something}.{fasta}.bwt.2bit.64",
            "{somewhere}/BWA_index/{something}.{fasta}.pac",
        output:
            directory("{somewhere}/BWA_index/{something}.{fasta}.sync"),
        log:
            "{somewhere}/BWA_index/sync.{something}.{fasta}.log"
        wildcard_constraints:
            fasta='fasta|fna'
        localrule: True
        shell:
            """
            dir="{output}"
            mkdir -p $dir
            for node in {CLUSTER_NODES}; do
                touch $dir/$node
            done
            """

    def get_node_request_argument(node, cluster_workload_manager):
        if cluster_workload_manager.upper() == 'SLURM':
            return f"--nodelist {node}"
        elif cluster_workload_manager.upper() in['PBS', 'TORQUE']:
            return f"-l nodes={node}"
        elif cluster_workload_manager.upper() == 'SGE':
            return f"-l h={node}"
        elif cluster_workload_manager.upper() == 'LSF':
            return f"-m {node}"
        else:
            raise Exception(f"MIntO error: cannot handle cluster workload manager '{cluster_workload_manager}'")

    rule rsync_index_files_to_local_dir:
        input:
            "{somefasta}.sync/{node}"
        output:
            "{somefasta}.sync/{node}.done"
        log:
            "{somefasta}.sync/{node}.log"
        resources:
            mem = 2,
            qsub_args = lambda wildcards: get_node_request_argument(wildcards.node, CLUSTER_WORKLOAD_MANAGER)
        wildcard_constraints:
            node = '|'.join(CLUSTER_NODES)
        shell:
            """
            path=$(dirname {wildcards.somefasta})
            umask 002 && mkdir -p {CLUSTER_LOCAL_DIR}/$path
            time (
                echo "NODE: $(hostname)"
                echo "SYNCHING:"
                echo "{wildcards.somefasta} --> {CLUSTER_LOCAL_DIR}"
                rsync -av --itemize-changes {wildcards.somefasta}.{{0123,amb,ann,bwt.2bit.64,pac}} {CLUSTER_LOCAL_DIR}/$path/ && touch {output}
            ) >& {log}
            """

    def get_bwa_index_sync_flags(wildcards):
        checkpoint_output = checkpoints.prepare_files_for_syncing.get(**wildcards).output[0]
        result = expand("{somewhere}/BWA_index/{something}.{fasta}.sync/{node}.done",
                        somewhere=wildcards.somewhere,
                        something=wildcards.something,
                        fasta=wildcards.fasta,
                        node=CLUSTER_NODES)
        return(result)


    # Keeping the NFS bwa-index output in this rule so that temp() doesn't delete it before this step is finished, just in case.
    rule BWA_index_contigs:
        input:
            get_bwa_index_sync_flags
        output:
            "{somewhere}/BWA_index_local/{something}.{fasta}.0123",
            "{somewhere}/BWA_index_local/{something}.{fasta}.amb",
            "{somewhere}/BWA_index_local/{something}.{fasta}.ann",
            "{somewhere}/BWA_index_local/{something}.{fasta}.bwt.2bit.64",
            "{somewhere}/BWA_index_local/{something}.{fasta}.pac",
        log:
            "{somewhere}/BWA_index_local/symlink.{something}.{fasta}.log",
        wildcard_constraints:
            fasta='fasta|fna'
        localrule: True
        shell:
            """
            time (
                for out in {output}; do
                    target=$(echo $out | sed "s/BWA_index_local/BWA_index/")
                    ln -s --force {CLUSTER_LOCAL_DIR}/$target $out
                done
            ) >& {log}
            """

else:

    # Baseline    : 5 GB
    # Size-based  : 22 byte per byte file size
    # New attempts: +40 GB each time
    rule BWA_index_contigs:
        input:
            fasta="{somewhere}/{something}.{fasta}"
        output:
            "{somewhere}/BWA_index/{something}.{fasta}.0123",
            "{somewhere}/BWA_index/{something}.{fasta}.amb",
            "{somewhere}/BWA_index/{something}.{fasta}.ann",
            "{somewhere}/BWA_index/{something}.{fasta}.bwt.2bit.64",
            "{somewhere}/BWA_index/{something}.{fasta}.pac",
        log:
            "{somewhere}/BWA_index/BWA_index.{something}.{fasta}.log"
        wildcard_constraints:
            fasta='fasta|fna'
        shadow:
            "minimal"
        resources:
            mem = lambda wildcards, input, attempt: 5 + int(22*input.size_mb/1024) + 40*(attempt-1),
        threads: 4
        conda:
            config["minto_dir"]+"/envs/MIntO_base.yml" #bwa-mem2
        shell:
            """
            outdir=$(dirname {output[0]})
            prefix=$(basename {input})
            time (
                bwa-mem2 index {input.fasta} -p $prefix
                rsync -av --itemize-changes $prefix.* $outdir/
            ) >& {log}
            """
