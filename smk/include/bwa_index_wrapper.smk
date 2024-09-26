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
#
# Note on CLUSTER_NODES:
# ----------------------
#
# 1. What is happening behind the scenes?
#
# Values in the python array CLUSTER_NODES will be used as node names to distribute the files.
# This is easier to control when run with '--cluster' especially if we pass {resources.qsub_args}.
# Here is an example --cluster argument for slurm:
#   --cluster 'sbatch -J {name} --mem={resources.mem}G --gres=gpu:{resources.gpu} -c {threads} {resources.qsub_args}'
# To avoid snakemake errors due to missing values (e.g., not all rules define resources.mem), we should precede it with:
#   --default-resources gpu=0 mem=4 "qsub_args=''"
# With this setup, the rule corresponding the one node (via {wildcards.node}) will submit
# the job to the CLUSTER_WORKLOAD_MANAGER with a request to run only on that node.
#
# 2. What if I forget --cluster to snakemake?
#
# If, by mistake, it is run without '--cluster', this could cause trouble: now all jobs
# (multiple values of {wildcards.node}) are run where you invoked your snakemake command.
# And this node could wrongly create a 'status' file for {wildcards.node} corresponding to a
# different node that denotes sync has been done for the latter. We have a simple solution for this.
# The node names in CLUSTER_NODES should match the node name accessible by the command 'hostname --short'.
# The sync status files are named with the local node name. You ask for sync'ing in nodeX, and
# only nodeX can make that file. This avoids problems when this is NOT run using --cluster.
# There is anothe problem that multiple sync rules may run the same 'rsync source dest' command
# on a single server. This would likely cause a race condition. To avoid that, we create a lock
# using 'flock' so that only one rsync is running for a specific bwa index sync on one node.
#
# 3. What if I run on a single node (no cluster), but want to avoid NFS traffic for index files?
#
# Our implementation, as a bonus, also allows rsyncing an NFS bwa index onto a local disk
# when running without --cluster, if you define CLUSTER_NODES with just your exec node's name.
################################################################################################

#################
# First, make the index
#################

# Memory requirements:
# --------------------
# Baseline    : 5 GB
# Size-based  : 22 byte per byte file size
# New attempts: +40 GB each time
rule BWA_index_in_place:
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
        fasta     = 'fasta|fna',
        something = '[^/]+'
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
            echo "NODE: $(hostname)"
            echo "SENDING:"
            echo "$(pwd)/$prefix --> $outdir"
            rsync -av --itemize-changes $prefix.* $outdir/
        ) >& {log}
        """

####################################################################
# If 'BWA_index_local/*' is requested, then distribute to all nodes
####################################################################

if CLUSTER_WORKLOAD_MANAGER is not None:
    known_managers = ['SLURM', 'PBS', 'TORQUE', 'SGE', 'LSF']
    if CLUSTER_WORKLOAD_MANAGER.upper() not in known_managers:
        raise Exception(f"MIntO error: Unexpected value: CLUSTER_WORKLOAD_MANAGER='{CLUSTER_WORKLOAD_MANAGER}'. Must be one of {known_managers}")

if CLUSTER_NODES is not None:

    # Make sure other variables are defined
    if CLUSTER_WORKLOAD_MANAGER is None or CLUSTER_LOCAL_DIR is None:
        raise Exception(f"MIntO error: CLUSTER_WORKLOAD_MANAGER and CLUSTER_LOCAL_DIR must be defined, if CLUSTER_NODES is defined")

    # If the index files exist, create a clustersync file to register a node
    rule prepare_index_files_for_syncing:
        localrule: True
        input:
            "{somewhere}/BWA_index/{something}.{fasta}.0123",
            "{somewhere}/BWA_index/{something}.{fasta}.amb",
            "{somewhere}/BWA_index/{something}.{fasta}.ann",
            "{somewhere}/BWA_index/{something}.{fasta}.bwt.2bit.64",
            "{somewhere}/BWA_index/{something}.{fasta}.pac",
        output:
            touch("{somewhere}/BWA_index/{something}.{fasta}.clustersync/{node}")
        log:
            "{somewhere}/BWA_index/sync.{something}.{fasta}.{node}.log"
        wildcard_constraints:
            node      = '|'.join(CLUSTER_NODES),
            fasta     = 'fasta|fna',
            something = '[^/]+'

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

    # If a node is registered with clustersync file, then synch it to that node and update synch-done status
    rule rsync_index_files_to_local_dir:
        input:
            "{somefasta}.clustersync/{node}"
        output:
            temp("{somefasta}.clustersync/{node}.synched")
        log:
            "{somefasta}.clustersync/{node}.synch.log"
        resources:
            mem = 2,
            qsub_args = lambda wildcards: get_node_request_argument(wildcards.node, CLUSTER_WORKLOAD_MANAGER)
        wildcard_constraints:
            node = '|'.join(CLUSTER_NODES)
        shell:
            """
            path=$(dirname {wildcards.somefasta})
            actual_node=$(hostname --short)
            umask 002 && mkdir -p {CLUSTER_LOCAL_DIR}/$path
            time (
                echo "TARGET NODE: {wildcards.node}"
                echo "ACTUAL NODE: $actual_node"
                echo "SYNC'N FILE: {wildcards.somefasta}.* --> {CLUSTER_LOCAL_DIR}/$path/"
                flock --exclusive --nonblock {CLUSTER_LOCAL_DIR}/{wildcards.somefasta}.synch.lock rsync -av --itemize-changes {wildcards.somefasta}.{{0123,amb,ann,bwt.2bit.64,pac}} {CLUSTER_LOCAL_DIR}/$path/ && touch {wildcards.somefasta}.clustersync/${{actual_node}}.synched
                echo " DONE"
            ) >& {log}
            """

    # If synch-done status is available for all nodes, then symlink an NFS location to index files on local disk
    rule symlink_local_index_files_to_target:
        localrule: True
        input:
            lambda wildcards: expand("{somewhere}/BWA_index/{something}.{fasta}.clustersync/{node}.synched",
                                        somewhere=wildcards.somewhere,
                                        something=wildcards.something,
                                        fasta=wildcards.fasta,
                                        node=CLUSTER_NODES)
        output:
            "{somewhere}/BWA_index_local/{something}.{fasta}.0123",
            "{somewhere}/BWA_index_local/{something}.{fasta}.amb",
            "{somewhere}/BWA_index_local/{something}.{fasta}.ann",
            "{somewhere}/BWA_index_local/{something}.{fasta}.bwt.2bit.64",
            "{somewhere}/BWA_index_local/{something}.{fasta}.pac",
        log:
            "{somewhere}/BWA_index/symlink_local.{something}.{fasta}.log",
        wildcard_constraints:
            fasta     = 'fasta|fna',
            something = '[^/]+'
        shell:
            """
            time (
                echo "SYMLINKING:"
                for out in {output}; do
                    target=$(echo $out | sed "s/BWA_index_local/BWA_index/")
                    echo "{CLUSTER_LOCAL_DIR}/$target --> $out"
                    ln -s --force {CLUSTER_LOCAL_DIR}/$target $out
                done
            ) >& {log}
            """

    # If a node is registered with clustersync file, then clean it on that node and update clean-done status
    rule clean_index_files_in_local_dir:
        input:
            "{somefasta}.clustersync/{node}"
        output:
            temp("{somefasta}.clustersync/{node}.cleaned")
        log:
            "{somefasta}.clustersync/{node}.clean.log"
        resources:
            mem = 2,
            qsub_args = lambda wildcards: get_node_request_argument(wildcards.node, CLUSTER_WORKLOAD_MANAGER)
        wildcard_constraints:
            node = '|'.join(CLUSTER_NODES)
        shell:
            """
            time (
                flock --exclusive --nonblock {CLUSTER_LOCAL_DIR}/{wildcards.somefasta}.clean.lock sh -c '
                    actual_node=$(hostname --short)
                    echo "NODE - TARGET: {wildcards.node}"
                    echo "NODE - ACTUAL: $actual_node"
                    echo "CLEANING FILE: {CLUSTER_LOCAL_DIR}/{wildcards.somefasta}.*"
                    if [ -f "{CLUSTER_LOCAL_DIR}/{wildcards.somefasta}.0123" ]; then
                        echo -n "  Index files found; deleting: "
                        rm {CLUSTER_LOCAL_DIR}/{wildcards.somefasta}.{{0123,amb,ann,bwt.2bit.64,pac}} && echo -n "DONE"
                        echo ""
                    else
                        echo "  WARNING: Index files NOT found; doing nothing"
                    fi
                    echo "  Touching flag: {wildcards.somefasta}.clustersync/$actual_node.cleaned"
                    touch {wildcards.somefasta}.clustersync/$actual_node.cleaned
                '
                echo " DONE"
            ) >& {log}
            """

    # If clean-done status is available for all nodes, delete symlinks & make global clean-done status file.
    # Rules that want us to clean up should ask for:
    #  '{somewhere}/BWA_index/{something}.{fasta}.clustersync/cleaning.done'.
    # As a rule, this should be done in a separate run after mapping is finished, otherwise
    # cleanup will run at the same time as synching, and the index files won't be available for mapping.
    # See how QC_2.smk implements it via a config variable 'CLEAN_BWA_INDEX'
    # that can be set in a separate instance of snakemake after QC_2.smk finishes successfully.
    rule check_index_files_are_cleaned:
        localrule: True
        input:
            status = lambda wildcards: expand("{somewhere}/BWA_index/{something}.{fasta}.clustersync/{node}.cleaned",
                                            somewhere=wildcards.somewhere,
                                            something=wildcards.something,
                                            fasta=wildcards.fasta,
                                            node=CLUSTER_NODES)
        output:
            "{somewhere}/BWA_index/{something}.{fasta}.clustersync/cleaning.done",
        wildcard_constraints:
            fasta     = 'fasta|fna',
            something = '[^/]+'
        shell:
            """
            # Delete symlinks
            rm {wildcards.somewhere}/BWA_index_local/{wildcards.something}.{wildcards.fasta}.*

            # Delete directory, if empty
            if [ -z "$(ls -A '{wildcards.somewhere}/BWA_index_local')" ]; then
                rmdir {wildcards.somewhere}/BWA_index_local/
            fi

            # Touch output
            touch {output}
            """
