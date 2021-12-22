#!/usr/bin/env python

'''
MAGs recovery and annotation

1) Run the binning program  (vamb in different option)
2) Run Checkm on all the data
3) Copy the HQ genomes in a folder
4) Run Coverm on HQ (why coverm, becasue it is easier to add a new binner in the case)
5) Retrieving the score for the genomes
6) Retrieving the best and unique set of genomes (with old scored formula)
7) Run prokka on the genomes (prokka) [separate environment]
8) Run taxonomic label on the genomes (PhylopHlan Metagenomic) [separate environment] 

Authors: Eleonora Nigro, Mani Arumugam
'''

# snakemake --snakefile mags_generation.smk --cluster "qsub -pe smp {threads} -l h_vmem={resources.mem}G -N {name} -cwd" --latency-wait 60 --jobs 25 --configfile mags_generation.yaml --use-conda #--conda-prefix tmp
# snakemake  --cluster "sbatch -J try_vamb --gres=gpu:1  -c 8" --jobs 1 --use-conda
# snakemake --snakefile /emc/cbmr/users/rzv923/MIntO/mags_generation.smk --restart-times 0 --keep-going --latency-wait 30 --cluster "qsub -pe smp {threads} -l h_vmem={resources.mem}G -N {name} -cwd" --use-conda --conda-prefix /emc/cbmr/users/rzv923/ibdmdb_test/tmp_porus/ --configfile mags_generation.smk.yaml --jobs 10
# snakemake --snakefile /emc/cbmr/users/rzv923/MIntO/mags_generation_test.smk --restart-times 0 --keep-going --latency-wait 30 --cluster "sbatch -J {name} --mem={resources.mem}G --gres=gpu:1 -c {threads} -e slurm-%x.e%A -o slurm-%x.o%A"  --use-conda --conda-prefix /data/MIntO_snakemake_env/ --configfile mags_generation.smk.yaml --jobs 10

wdir = config['working_dir']
working_dir = config['working_dir'] #+ "/metaG/8-1-binning"
minto_dir=config["minto_dir"]
script_dir=config["minto_dir"]+"/scripts/"
#binner=config['BINNERS']
project = config['PROJECT']

print(config['BINNERS'])

## Configuration file
# This is given by Carmen (fasta files and contigs file)
#configfile: "/emc/cbmr/users/wdm117/github/mags_generation_pipeline/mags_generation_pipeline_snakemake/config_binning.yaml" # could be deleted 

##############################
def mags_recovery():
    result = expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/{binner}/clusters.tsv",
                    wd = working_dir,
                    binner = config['BINNERS']),\
    expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/{binner}/{binner}_discarded_genomes.txt", 
                    wd = working_dir,
                    binner = config['BINNERS']),\
	expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/all", 
                    wd = working_dir),\
	expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/copy_genomes_all_finished.txt", 
                    wd = working_dir),\
	expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/checkm_completed.txt", 
                    wd = working_dir),\
	expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/checkm/checkm-comprehensive.tsv", 
                    wd = working_dir),\
	expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/HQ_genomes_checkm.tsv", 
                    wd = working_dir),\
	expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/coverm_unique_cluster.tsv", 
                    wd = working_dir),\
	expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/HQ_genomes_checkm_scored.tsv", 
                    wd = working_dir),\
	expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/coverm_unique_cluster_scored.tsv", 
                    wd = working_dir),\
	expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/best_unique_genomes.txt", 
                    wd = working_dir),\
	expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/mags_generation_pipelined_ended!", 
                    wd = working_dir),\
	expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/prokka_info.txt", 
                    wd = working_dir),\
	expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/taxonomy_info.txt", 
                    wd = working_dir)
    return(result)

rule all:
    input: 
        mags_recovery()

##############################

### Run Vamb
rule run_vamb:
	input:
		contigs_file=lambda wildcards: "{wd}/metaG/8-1-binning/{project}_scaffolds.2500.fasta".format(wd = working_dir, project = project), 
		#contigs_file=config["CONTIGS_FILE"],
		depth_file=lambda wildcards: "{wd}/metaG/8-1-binning/{project}_scaffolds.2500.depth.txt".format(wd = working_dir, project = project),
		#depth_file=config["DEPTH_FILE"]		
	
	output:
		"{wd}/metaG/8-1-binning/mags_generation_pipeline/{binner}/clusters.tsv",#.format(wd = working_dir, binner = binner), 
		"{wd}/metaG/8-1-binning/mags_generation_pipeline/{binner}/log.txt",#.format(wd = working_dir, binner = binner), 
		"{wd}/metaG/8-1-binning/mags_generation_pipeline/{binner}/latent.npz",#.format(wd = working_dir, binner = binner), 
		"{wd}/metaG/8-1-binning/mags_generation_pipeline/{binner}/lengths.npz",#.format(wd = working_dir, binner = binner), 
		"{wd}/metaG/8-1-binning/mags_generation_pipeline/{binner}/mask.npz",#.format(wd = working_dir, binner = binner), 
		"{wd}/metaG/8-1-binning/mags_generation_pipeline/{binner}/model.pt",#.format(wd = working_dir, binner = binner), 
		"{wd}/metaG/8-1-binning/mags_generation_pipeline/{binner}/tnf.npz",#.format(wd = working_dir, binner = binner), 
	
	params:
		gpu=config["VAMB_GPU"],
		#run_vamb = "{script_dir}/run_vamb.sh"
	
	log:
		"{wd}/logs/metaG/mags_generation/run_vamb_{binner}.log"#.format(wdir = wdir, binner = binner)
	
	resources:
		mem=config['VAMB_memory']
	
	threads: 
		config["VAMB_THREADS"]
	
	conda:
		config["minto_dir"]+"/envs/vamb.yaml" #"/emc/cbmr/users/rzv923/vamb_orig.yaml" #f68d45ef
		#"/emc/cbmr/users/rzv923/vamb_orig.yaml" 90d82712
	
	shell:
		""" time (sh {script_dir}/run_vamb.sh {params.gpu} {wildcards.binner} {input.contigs_file} {input.depth_file} {threads} {wildcards.wd}/metaG/8-1-binning/mags_generation_pipeline/{wildcards.binner}
		rsync {wildcards.wd}/metaG/8-1-binning/mags_generation_pipeline/{wildcards.binner}/tmp/*  {wildcards.wd}/metaG/8-1-binning/mags_generation_pipeline/{wildcards.binner}
		rm -rf {wildcards.wd}/metaG/8-1-binning/mags_generation_pipeline/{wildcards.binner}/tmp) &> {log}"""

#####
		#"""gpu={params.gpu}; \
		#	binner={wildcards.binner};\
		#	contigs_file={input.contigs_file}; \
		#	depth_file={input.depth_file};\
		#	threads={threads}; \
		#	output={wildcards.initial_folder}/{wildcards.binner};\
		#	echo "gpu: ${{gpu}}" ;echo "binner: ${{binner}}"; echo "contigs file: ${{contigs_file}}"; echo "depth file:  ${{depth_file}}"; echo "threads: ${{threads}}"; echo "output: ${{output}}" ;\
		#	if [ $gpu == "yes" ]; then \
		#		if [ ${{binner}} == "vamb_256" ]; then \
		#			echo "Launching vamb 256" \
		#			vamb -l 16 -n 256 256  --fasta $contigs_file --jgi $depth_file -m 2500 -o _ -p $threads --cuda --outdir $output/tmp \
		#		elif [ ${{binner}} == "vamb_384" ]; then \
		#			vamb -l 24 -n 384 384  --fasta $contigs_file --jgi $depth_file -m 2500 -o _ -p $threads --cuda --outdir $output/tmp \
		#		elif [ ${{binner}} == "vamb_512" ]; then \
		#			vamb -l 32 -n 512 512 --fasta $contigs_file --jgi $depth_file -m 2500 -o _ -p $threads --cuda --outdir $output/tmp \
		#		elif [ $binner == "vamb_768" ]; then \
		#			vamb -l 40 -n 768 786 --fasta $contigs_file --jgi $depth_file -m 2500 -o _ -p $threads --cuda --outdir $output/tmp \
		#		else echo "Something went wrong" \
		#		fi \
		#	else \
		#		if [ ${{binner}} == "vamb_256" ]; then \
		#			vamb -l 16 -n 256 256  --fasta $contigs_file --jgi $depth_file -m 2500 -o _ -p $threads  --outdir $output/tmp \
		#		elif [ ${{binner}} == "vamb_384" ]; then \
		#			vamb -l 24 -n 384 384  --fasta $contigs_file --jgi $depth_file -m 2500 -o _ -p $threads  --outdir $output/tmp \
		#		elif [ ${{binner}} == "vamb_512" ]; then \
		#			vamb -l 32 -n 512 512  --fasta $contigs_file --jgi $depth_file -m 2500 -o _ -p $threads  --outdir $output/tmp \
		#		elif [ ${{binner}} == "vamb_768" ]; then \
		#			vamb -l 40 -n 768 768  --fasta $contigs_file --jgi $depth_file -m 2500 -o _ -p $threads  --outdir $output/tmp \
		#		else echo "Something went wrong" \
		#		fi \ 
		#	fi \
		#	rsync {wildcards.initial_folder}/{wildcards.binner}/tmp/* {wildcards.initial_folder}/{wildcards.binner} && rm -rf {wildcards.initial_folder}/{wildcards.binner}/tmp"""	

### Run take all genomes [put all the genomes in a folder "all" where CheckM will be launched] # this is on vamb, if there are other binners, depending on the output, the bins should be moved in all
rule take_all_genomes_for_each_run:
	input:
		# #name1 = configfile["NAME1"] # maybe they should be added 
		# #name2 = configfile["NAME2"] # maybe they should be added
		# vamb_cluster =  "{initial_folder}/{binner}/clusters.tsv", 
		# contigs_file = config["CONTIGS_FILE"]
		vamb_cluster =  "{wd}/metaG/8-1-binning/mags_generation_pipeline/{binner}/clusters.tsv",#.format(wd = working_dir, binner = binner),
		contigs_file = "{wd}/metaG/8-1-binning/{project}_scaffolds.2500.fasta".format(wd = working_dir, project = project),

	output:
		discarded_genomes = "{wd}/metaG/8-1-binning/mags_generation_pipeline/{binner}/{binner}_discarded_genomes.txt",#.format(wd = working_dir, binner = binner),
		tmp_folder = directory("{wd}/metaG/8-1-binning/mags_generation_pipeline/{binner}/tmp_folder"),
		
	params:
		min_fasta_length = config["MIN_FASTA_LENGTH"],
		#script_take_all_genomes = "{script_dir}/take_all_genomes.py",
		tmp_folder = directory("{wd}/metaG/8-1-binning/mags_generation_pipeline/{binner}/tmp_folder"),
	
	log:
		"{wd}/logs/metaG/mags_generation/{binner}.take_all_genomes_for_each_run.log"#.format(wdir = wdir, binner = binner)
	
	resources:
		mem=10
	
	threads: 
		8 # Decide number of threads
	
	shell:
		""" time (python {script_dir}/take_all_genomes.py --vamb_cluster_tsv {input.vamb_cluster} --contigs_file {input.contigs_file} --assembly_method_name {wildcards.binner} \
--min_fasta_length {params.min_fasta_length} --output_folder {params.tmp_folder} --discarded_genomes_info {output.discarded_genomes}) &> {log} """
		#--output_folder {wildcards.initial_folder}/{wildcards.binner}/tmp_genomes

### Run copy all the genomes and remove tmp folders
rule copy_genomes_in_all:
	input:
		#all_folder = expand('{initial_folder}/{binner}/tmp_genomes', initial_folder=config['INITIAL_FOLDER'], binner=config['BINNERS'])
		all_folder = expand("{wd}/metaG/8-1-binning/mags_generation_pipeline/{binner}/tmp_folder",wd = working_dir, binner = config['BINNERS'])

	output:
		all_genomes = directory("{wd}/metaG/8-1-binning/mags_generation_pipeline/all"), # remember to cancel it in the rule_all
		output = "{wd}/metaG/8-1-binning/mags_generation_pipeline/copy_genomes_all_finished.txt"
	
	log:
 		"{wd}/logs/metaG/mags_generation/copy_genomes_in_all.log"#.format(wdir = config['working_dir'])
	
	resources:
		mem=10
	
	threads: 
		8 # Decide number of threads

	shell:
		""" time (mkdir -p {output.all_genomes}
		cp -r {wildcards.wd}/metaG/8-1-binning/mags_generation_pipeline/vamb*/tmp_folder/*.fna {output.all_genomes}
		echo 'Finished to copy genomes in all' > {output.output}) &> {log} """
		#"find {input.all_folder} -name '*.fna' -type f -exec cp {} {output} \;"
		#"mkdir -p {wildcards.initial_folder}/all && cp -r {input.all_folder}/*.fna {output}"
		#"rsync -a {input.all_folder}/*.fna {output}" # this create {input.all_folder}/all/tmp_genomes
		#"cp -r {input.all_folder}/*.fna {output}"

## Run checkm on the genomes in all
rule run_checkm:
	input:
		all_genomes="{wd}/metaG/8-1-binning/mags_generation_pipeline/all"
	
	output:
		"{wd}/metaG/8-1-binning/mags_generation_pipeline/checkm_completed.txt"
	
	# params:
	# 	#checkm_threads=config["CHECKM_THREADS"], #Move to threads
	# 	#pplacer_threads=config["PPLACER_THREADS"],#Move to threads
	# 	making_batch="{script_dir}/making_batch_checkm.sh",#.format(config["SCRIPT_FOLDER"]),
	# 	run_checkm="{script_dir}/run_checkm.sh"#.format(config["SCRIPT_FOLDER"])
	
	log:
 		"{wd}/logs/metaG/mags_generation/run_checkm.log"#.format(wdir = config['working_dir'])

	resources:
		mem=config["CHECKM_memory"]
	
	threads: # Decide on of the two
 		config["CHECKM_THREADS"]
		#checkm_threads=config["CHECKM_THREADS"],
		#pplacer_threads=config["PPLACER_THREADS"],
	
	conda:
		config["minto_dir"]+"/envs/checkm.yaml"
	
	shell:
		# """ rsync {input.all_genomes}/*.fna {wildcards.initial_folder}/checkm
		# cd {wildcards.wd}/checkm
		# sh {params.making_batch}
		# sh {params.run_checkm} {params.checkm_threads} {params.pplacer_threads}
		# echo 'checkm finished!' >> {wildcards.wd}/checkm_completed.txt """
		""" time (rsync {input.all_genomes}/*.fna {wildcards.wd}/metaG/8-1-binning/mags_generation_pipeline/checkm
		cd {wildcards.wd}/metaG/8-1-binning/mags_generation_pipeline/checkm
		sh {script_dir}/making_batch_checkm.sh
		sh {script_dir}/run_checkm.sh {threads} {threads}
		echo 'checkm finished!' >> {wildcards.wd}/metaG/8-1-binning/mags_generation_pipeline/checkm_completed.txt) &> {log} """

## Create a comphrensive table with checkm 
rule make_comphrensive_table:
	input:
		"{wd}/metaG/8-1-binning/mags_generation_pipeline/checkm_completed.txt"
		
	output:
		checkm_total = "{wd}/metaG/8-1-binning/mags_generation_pipeline/checkm/checkm-comprehensive.tsv"
	
	params:
		checkm_tsv_tables = "{wd}/metaG/8-1-binning/mags_generation_pipeline/checkm",
		remove_intermediate_files_checkm = config["CLEAN_CHECKM"]
		#make_comprehensive="{}/make_comprehensive_checkm.py".format(config["SCRIPT_FOLDER"])
	
	log:
 		"{wd}/logs/metaG/mags_generation/make_comphrensive_table.log"#.format(wdir = config['working_dir'])

	resources:
		mem=10
	
	threads:
 		8

	#shell:
	#	"python {params.make_comprehensive --folder_with_checkm {input.checkm_tsv_tables} --output_file {output.checkm_total} --remove_intermediate_files {params.remove_intermediates_files_checkm}"

	run:	
		import glob
		import os
		import pandas as pd	
		import shutil
		all_checkm_files = glob.glob(params.checkm_tsv_tables + "/*.tsv")
		# concatenate all the .tsv file in the folder in order to create a comphresenive file
		li = []
		for filename in all_checkm_files:
			df = pd.read_csv(filename, index_col=None, header=0, sep = "\t")
			li.append(df)
		all_checkm_output = pd.concat(li, axis=0, ignore_index=True)
		# save the file with all the checkm in the same file
		all_checkm_output.to_csv("{}".format(output.checkm_total), sep = "\t", index = False)
		
		if params.remove_intermediate_files_checkm == "yes":
			folder=glob.glob(params.checkm_tsv_tables + "/*/") # (could be done os.remove(input.checkm_tsv_tables + "/*/) , but just to be sure not to take any .tsv)
			for f in folder:
				if not ".tsv" in f:
					print("[rule make_comphrensive_table]: removing intermediate files: {}".format(f))
					shutil.rmtree(f)

## Copy HQ genomes inside HQ_genomes folder
rule copy_HQ_genomes:
 	input:
 		checkm_total="{wd}/metaG/8-1-binning/mags_generation_pipeline/checkm/checkm-comprehensive.tsv",
		#all_genomes_folder = "{initial_folder}/all/tmp_genomes"
 	
	output:
 		HQ_table="{wd}/metaG/8-1-binning/mags_generation_pipeline/HQ_genomes_checkm.tsv",
		
	params:
		HQ_folder="{wd}/metaG/8-1-binning/mags_generation_pipeline/HQ_genomes",
		all_genomes_folder =  "{wd}/metaG/8-1-binning/mags_generation_pipeline/all/",
 		completeness = config["CHECKM_COMPLETENESS"],
 		contamination = config["CHECKM_CONTAMINATION"]
	
	log:
 		"{wd}/logs/metaG/mags_generation/copy_HQ_genomes.log"#.format(wdir = config['working_dir'])

	resources:
		mem=10
	
	threads:
 		8

	run:
		import subprocess
		import pandas as pd
		# open the checkm_comprhrensive table
		checkm_results=pd.read_csv(input.checkm_total, sep = "\t")
 		# take and save the HQ table
		HQ_checkm_results = checkm_results[(checkm_results["Completeness"] >= params.completeness) & (checkm_results["Contamination"] <= params.contamination)]
		HQ_checkm_results.to_csv(output.HQ_table, sep = "\t", index = False)
 		# create the path for copying the genomes
		try:
			os.mkdir(params.HQ_folder)
		except OSError:
			print("Creation of the directory {} failed!".format(params.HQ_folder))
 		# take the bins
		hq_bins = list(HQ_checkm_results["Bin Id"])
		for bin_id in hq_bins:
			source_file = params.all_genomes_folder +"/{}.fna".format(bin_id) 
			destination_file = params.HQ_folder  + "/{}.fna".format(bin_id)
			print("[rule copy_HQ_genomes] Copying {} to {}".format(source_file, destination_file))
			subprocess.call(["cp", source_file, destination_file] )

## Run coverm on HQ genomes to create the .tsv file
rule run_coverm:
	input:
		HQ_table="{wd}/metaG/8-1-binning/mags_generation_pipeline/HQ_genomes_checkm.tsv"
	
	output:
		coverm_output="{wd}/metaG/8-1-binning/mags_generation_pipeline/coverm_unique_cluster.tsv"  #unique-{}-cluster.tsv"
	
	params:
		#coverm_threads = config["COVERM_THREADS"], # Moved to threads
		HQ_folder="{wd}/metaG/8-1-binning/mags_generation_pipeline/HQ_genomes"
		
	log:
 		"{wd}/logs/metaG/mags_generation/run_coverm.log"

	resources:
		mem=config["COVERM_memory"]
	
	threads: 
		config["COVERM_THREADS"]

	conda:
		config["minto_dir"]+"/envs/coverm.yaml"

	shell:
		""" time (coverm cluster --genome-fasta-directory {params.HQ_folder} -x fna --ani 99  --output-cluster-definition {output.coverm_output} --threads {threads} --precluster-method finch) &> {log} """

#################### TESTED TILL HERE #############

## Run retrieving scored 
rule calculate_score_genomes:
	input:
		coverm_output = "{wd}/metaG/8-1-binning/mags_generation_pipeline/coverm_unique_cluster.tsv",
		HQ_table = "{wd}/metaG/8-1-binning/mags_generation_pipeline/HQ_genomes_checkm.tsv"

	output:
		scored_genomes = "{wd}/metaG/8-1-binning/mags_generation_pipeline/HQ_genomes_checkm_scored.tsv"
	
	params:
		HQ_folder="{wd}/metaG/8-1-binning/mags_generation_pipeline/HQ_genomes",
		#calculate_genomes_score="{script_dir}/calculate_genomes_score.py"
		score_method = config["SCORE_METHOD"]

	log:
 		"{wd}/logs/metaG/mags_generation/calculate_score_genomes.log"

	resources:
		mem=10
	
	threads: 
		8 # Decide number of threads

	#conda:
	#	config["minto_dir"]+"/envs/phylophlan3.yaml"

	shell:
		""" time (python {script_dir}calculate_genomes_score.py --checkm_output {input.HQ_table} --fasta_folder {params.HQ_folder} --output_file {output.scored_genomes} --score_method {params.score_method}) &> {log} """


## Run retrieved the best unique genomes
rule find_unique_and_best_genomes:
	input:
		scored_genomes = "{wd}/metaG/8-1-binning/mags_generation_pipeline/HQ_genomes_checkm_scored.tsv",
		coverm="{wd}/metaG/8-1-binning/mags_generation_pipeline/coverm_unique_cluster.tsv"
	
	output:
		scored = "{wd}/metaG/8-1-binning/mags_generation_pipeline/coverm_unique_cluster_scored.tsv",
		best_unique_genomes = "{wd}/metaG/8-1-binning/mags_generation_pipeline/best_unique_genomes.txt"

	log:
 		"{wd}/logs/metaG/mags_generation/find_unique_and_best_genomes.log"

	resources:
		mem=10
	
	threads: 
		8 # Decide number of threads

	run:
		import pandas as pd

		# read the table for the score
		score_table = pd.read_csv(input.scored_genomes, sep = "\t", index_col = "Bin_id", comment = "#") # we skip the first line with the --score_method

		# read coverm table
		coverm_table = pd.read_csv(input.coverm, sep = "\t", names = ["ref_cluster", "cluster_members"])


		# list of best genomes that should be written in the output
		best_genomes_list = []

		# create a dictionary of cluster
		d_cluster = {}

		for i in range(len(coverm_table)):

			ref_cluster = coverm_table["ref_cluster"][i]
			cluster_members = coverm_table["cluster_members"][i].split("/")[-1].replace(".fna", "") # it will append also the name of the genome withput the path

			if ref_cluster not in d_cluster:
				d_cluster[ref_cluster] = [cluster_members]

			else:
				d_cluster[ref_cluster].append(cluster_members)

		# now we take the best genome based on the score
		for cluster in d_cluster:

			genomes_in_the_cluster = d_cluster[cluster]
			dataframe_score = score_table.loc[genomes_in_the_cluster, ["Score"]].sort_values(by=["Score"], ascending=False) # we take the genomes from the score table
			best_genome = dataframe_score.index[0]
			#best_genomes = dataframe_score["Bin_id"][0] # we take the best genome
			best_genomes_list.append(best_genome)

		
		# Subset score table in order to have the best genomes only 
		best_genomes_scored = score_table.loc[best_genomes_list]

		# Save the file
		best_genomes_scored = best_genomes_scored.sort_values(by=["Score"], ascending=False)
		best_genomes_scored.to_csv(output.scored, sep = "\t", index = True)

		# Create the file
		best_unique_genomes_list= list(best_genomes_scored.index)

		with open(output.best_unique_genomes, "w") as fh:
			for genome in best_unique_genomes_list:
				fh.write("{}\n".format(genome))
		
## Run copy the best genomes
rule copy_best_genomes:
	input:
		best_unique_genomes = "{wd}/metaG/8-1-binning/mags_generation_pipeline/best_unique_genomes.txt"
	
	output:
		genomes_copied = "{wd}/metaG/8-1-binning/mags_generation_pipeline/mags_generation_pipelined_ended!"
	
	log:
 		"{wd}/logs/metaG/mags_generation/copy_best_genomes.log"

	resources:
		mem=10
	
	threads: 
		8 # Decide number of threads

	shell:
		""" time (mkdir {wildcards.wd}/metaG/8-1-binning/mags_generation_pipeline/unique_genomes
		while read line; do cp {wildcards.wd}/metaG/8-1-binning/mags_generation_pipeline/HQ_genomes/${{line}}.fna {wildcards.wd}/metaG/8-1-binning/mags_generation_pipeline/unique_genomes ; 
        done < {input.best_unique_genomes}
		touch {output.genomes_copied})&> {log}  """

## Run prokka
rule run_prokka:
	input:
		"{wd}/metaG/8-1-binning/mags_generation_pipeline/mags_generation_pipelined_ended!"
	
	output:
		prokka_ended = "{wd}/metaG/8-1-binning/mags_generation_pipeline/prokka_info.txt"
	
	params:
		unique_genomes_folder = "{wd}/metaG/8-1-binning/mags_generation_pipeline/unique_genomes",
		run_prokka = config["RUN_PROKKA"],
		#prokka_cpus = config["PROKKA_CPUS"], # Moved to threads
		#prokka_script = "{script_dir}/run_prokka.sh",
		#prokka_folder = config["INITIAL_FOLDER"] #Replaced by wd

	log:
 		"{wd}/logs/metaG/mags_generation/run_prokka.log"

	resources:
		mem=config["PROKKA_memory"]
	
	threads: 
		config["PROKKA_CPUS"]
	
	conda:
		config["minto_dir"]+"/envs/prokka.yaml"
	
	shell: 
		""" time (sh {script_dir}run_prokka.sh {params.run_prokka} {threads} {wildcards.wd}/metaG/8-1-binning/mags_generation_pipeline/ {params.unique_genomes_folder} {output.prokka_ended})&> {log} """


## Run taxonomy
rule run_taxonomy:
	input:
		"{wd}/metaG/8-1-binning/mags_generation_pipeline/mags_generation_pipelined_ended!"
	
	output:
		taxonomy_ended = "{wd}/metaG/8-1-binning/mags_generation_pipeline/taxonomy_info.txt"
	
	params:
		unique_genomes_folder = "{wd}/metaG/8-1-binning/mags_generation_pipeline/unique_genomes",
		output_phylophlan = "{wd}/metaG/8-1-binning/mags_generation_pipeline/taxonomy",
		run_taxonomy = config["RUN_TAXONOMY"],
		#taxonomy_cpus = config["TAXONOMY_CPUS"], # Moved to threads
		#taxonomy_script = "{script_dir}/run_taxonomy.sh",
		#taxonomy_folder = config["INITIAL_FOLDER"],
		taxonomy_database = config["TAXONOMY_DATABASE"],
		database_folder = config["DATABASE_FOLDER"]
	
	log:
 		"{wd}/logs/metaG/mags_generation/run_taxonomy.log"
	
	resources:
		mem=config["TAXONOMY_memory"]

	threads: 
		config["TAXONOMY_CPUS"]
	
	conda:
		config["minto_dir"]+"/envs/phylophlan3.yaml"

	shell: 
		"""time (sh {script_dir}run_taxonomy.sh {params.run_taxonomy} {threads} {params.unique_genomes_folder} {params.output_phylophlan} {params.database_folder} {output.taxonomy_ended} {params.taxonomy_database})&> {log} """

