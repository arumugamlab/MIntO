<h1 align="center"> MIntO </h1>
<h1 align="center"> A Modular and Scalable Pipeline for Microbiome Meta-omics Data Integration </h1>

MIntO (Microbiome Integrated Meta-Omics), a highly versatile pipeline that integrates metagenomic and metatranscriptomic data in a scalable way.
The distinctive feature of this pipeline is the computation of gene expression profile by taking into account the community turnover and gene expression variations as the underlying process that shapes the community transcript levels along the time and between conditions.
The integrated pipeline will be relevant to provide unique biochemical insights into the microbial ecology by linking functions to retrieved genomes and to examine gene expression variations. MIntO can be downloaded from https://github.com/arumugamlab/MIntO.


<p align="center"><img src="images/MIntO_logo.png" height="400" /></p>

MIntO implementation and automation are achieved by Snakemake (Mölder et al. 2021):
- User-friendly framework 
- Facilitates the scalability of the pipeline by optimizing the number of parallel processes from a single-core workstations to compute clusters. 
- Docker and Conda environments to ensure version control of the different libraries 
- Modular design to enable the user to choose between three available modes based on the input data and the experiment design, including de novo assembly.

## Dependencies

MIntO has two dependencies: Conda and FetchMGs

### FetchMGs installation 
FetchMGs can be downloaded and installed following the instruction at https://github.com/motu-tool/fetchMGs

### Anaconda installation
Anaconda can be downloaded and installed following the instructions at https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
Please note:
- We tested our pipeline on a 64-bit Linux system using Anaconda3. 
- In order to download and install Anaconda, there should be a minimum of 3 GB disk space.
- After installing conda, please close the console for the changes to be active. You might have to open a new terminal.

### Snakemake installation
Snakemake can be downloaded and installed using Conda, following the instructions at https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

### Download Databases
In addition, MIntO uses several databases that have to be downloaded by the user.
- Before running QC_2.smk for the first time, SortMeRNA database has to be downloaded and indexed.
  - SortMeRNA databases can be downloaded here:  https://github.com/biocore/sortmerna/tree/master/data/rRNA_databases
 
- Before running gene_annotation.smk for the first time, there are three databases to download/install: emapperdb, KOfam database and dbcan database
  - EggNOG-mapper tool - https://github.com/eggnogdb/eggnog-mapper/.
  
    Download emapperdb using download_eggnog_data.py available at https://github.com/eggnogdb/eggnog-mapper/
    
    The emapperdb should be downloaded to <PathToMIntO>/data/eggnog_data/data
    
    conda create env <PathToMIntO>/envs/py38_env.yml --name <CondaEnvName>
    conda activate <CondaEnvName>
    python3 download_eggnog_data.py --data_dir <PathToMIntO>/data/eggnog_data/data/ -P -M -f

  - KofamScan tool - https://github.com/takaram/kofam_scan.

    Download KOfam database from ftp://ftp.genome.jp/pub/db/kofam/ and decompress it. 
    
    The profile HMMs in profiles/ directory should be copied to MIntO/data/kofam_db/profiles and ko_list file should be copied to MIntO/data/kofam_db/

  - dbCAN2 tool - https://github.com/linnabrown/run_dbcan
  
    Follow instructions to install database in https://github.com/linnabrown/run_dbcan
  
    The generated files should be copied to MIntO/data/dbCAN_db/

## Three modes

The modular design of MIntO enables the user to run the pipeline using three available modes based on the input data and the experiment design: 
- Genome-based assembly-dependent mode recovers MAGs from metagenomics samples. The required inputs are metagenomic and/or metatranscriptomic raw reads and Snakemake configuration file for QC_1.smk
- Genome-based assembly-free mode uses publicly available genomes provided by the user. In addition to genome-based assembly-dependent requirements, the user should provide as input the genomes as fasta files, genome features as GFF files, and amino acid sequences of protein-coding genes as fasta files.
- Gene-catalog-based assembly-free requires a gene catalog provided by the user. In addition to genome-based assembly-dependent requirements, the user should provide as input the gene-based database.

## Structure
MIntO can be divided into seven major steps, which will be discussed in the next paragraphs using our analysis of example data (Figure 1A):
- Quality control and preprocessing (QC_1.smk and QC_2.smk) 
- Assembly-free taxonomy profiling (QC_2.smk)
- Recovery of MAGs and taxonomic annotation (assembly.smk, binning-preparation.smk and mags_generation.smk) 
  - only run in genome-based assembly-dependent mode
- Gene prediction and functional annotation (gene_annotation.smk)
  - only run in genome-based modes
- Alignment and normalization (gene_abundance.smk)
  - genome-based mode: recovered MAG or publicly available genomes
  - gene-based mode: gene catalog
- Integration: Gene and function profiling (integration.smk)
  - assembly-dependent or assembly-free modes
- Visualization and reporting

In the three modes, the pipeline workflow includes quality control and preprocessing; assembly-free taxonomy profiling of high quality metagenomic reads by identifying phylogenetic markers; alignment and normalization; integration: gene and function profiling; and visualization and reporting. The gene prediction and functional annotation step is run using the recovered MAGs or publicly available genomes. 

## Requirements

Before running MIntO, each pair of reads should be grouped per sample ID in a separate directory with the corresponding sample ID name. FASTQ input files should be gzipped compressed with the .1.fq.gz extension for forward reads and .2.fq.gz extension for reverse reads. The same sample ID in metaG and metaT is required to match the data and generate the gene and function expression profiles.

MIntO requires a configuration file as an input which includes several parameters used in each step by the .smk script. Each configuration file is generated when running the previous .smk script, but they should be filled in by the user.. A template of the required configuration file to run QC_1.smk can be found at MIntO/data/QC_1.yaml. 

## Configuration files parameters

### General settings used in most of the steps:

| Parameters     |   | Description                                                                                                                                                                               |
|-------------|---|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| PROJECT     | path  | Name of the project.                                                                                                                                                                               |
| working_dir |  metaG or metaT | Path to the directory where MIntO generates the output directories and files.                                                                                                                      |
| omics       | path | Type of reads used as input. ‘metaG’ for metagenomic reads. ‘metaT’ for metatranscriptomic reads.                                                                                                  |
| local_dir   | path  | Local directory where temporary files are generated while Snakemake rules are run. These temporary files are deleted when the rule is finished.                                                    |
| minto_dir   | path  | Location where MIntO has been downloaded.                                                                                                                                                          |
| METADATA    | path/file | Plain text metadata file with the first 3 columns as: sample ID, conditions and sample alias. This file is used to generate the output plots. If no metadata is provided, the sample IDs are used. |

### Parameters used in QC_1.smk should be included in QC_1.yaml
In the pre-processing of metaG and metaT data, the first step is the quality reads filtering. In order to keep as many sequences as possible, the raw reads are filtered first by sequencing adapters and low quality bases using Trimmomatic  v0.39. The distribution of the cumulative reads per base position on all the trimmed by quality paired-read fastq files can be checked in:

<working_dir>/output/1-trimmed/<omics>._cumulative_read_lenght_cutoff.pdf

| Parameters           |                            | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
|----------------------|----------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| raw_reads_dir        | path                       | Directory where the raw paired-reads are located. Each pair of reads should be grouped per sample ID in a separate directory with the corresponding sample ID name. FASTQ input files should be gzipped compressed with the .1.fq.gz extension for forward reads and .2.fq.gz extension for reverse reads.                                                                                                                                                                                                                                           |
| trimmomatic_adaptors | <adapters.fa, False or Skip> | adapters.fa, sequence adapters file path to be used by Trimmomatic. False, if there is no sequence adapters file, a custom script retrieves the adapters by selecting the most abundant index in the first 10,000 headers of the raw fastq files. The user should provide two files: a fasta file with the adapter sequences named as {omics}.adapters.fa.in and a file with the list of headers in {omics}.adapters.fa.in. This file should be called {omics}.palindrome.list. Skip, skip this step if adapter sequences have already been removed. |
| perc_remaining_reads | int [default 95]           | The MINLEN parameter in Trimmomatic (used in QC_2.smk) is used to remove the reads that are too short. This cutoff is estimated as the maximum length above which a predefined percentage of the reads from the previous step are retained. f the estimated read length cutoff is below 50bp, trimmomatic will use 50bp as the minimum sequence length.                                                                                                                                                                                              |
| TRIMMOMATIC_threads  | int                        | Number of threads used to run Trimmomatic.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| TRIMMOMATIC_memory   | int                        | Number of GBs used to run Trimmomatic per thread.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| ILLUMINA             |                            | List of sample IDs. The names should match with directories in the raw_reads_dir parameter.                                                                                                                                                                                                                                                                                                                                                                                                                                                          |

QC_1.smk script can be run like:

- snakemake --snakefile <PathToMIntO>/QC_1.smk --cluster "qsub -pe smp {threads} -l h_vmem={resources.mem}G -N {name} -cwd" --latency-wait 30 --jobs <int> --configfile <PathToConfigFile>QC1.yaml --use-conda --conda-prefix <PathToTmpCondaEnv>/tmp

### Parameters used in QC_2.smk should be included in QC_2.yaml
In the pre-processing of metaG and metaT data, the second step is the read length filtering. In order to keep as many sequences as possible, the raw reads are filtered by read length using Trimmomatic  v0.39 and the MINLEN parameter calculated in QC_1.smk. 
MetaG high-quality filtered reads can be profiled by MetaPhlAn3 or mOTUs2. To explore the similarities and dissimilarities of the data, the relative abundance of the species composition is used to generate two outputs: (1) the 15 most abundant genera across the samples, and (2) a principal-coordinate analysis (PCoA) using Bray Curtis distance generated in this directory:
<working_dir>/output/6-taxa_profile/

MIntO generates the required configuration file in <working_dir>/<omics>/QC_2.yaml

#### Pre-processing - length reads filtering. Trimmomatic filters sequences by read length

| Parameters          |     | Description                                                                                                                           |
|---------------------|-----|---------------------------------------------------------------------------------------------------------------------------------------|
| TRIMMOMATIC_threads | int | Number of threads used by Trimmomatic                                                                                                 |
| TRIMMOMATIC_memory  | int | Number of GBs used by Trimmomatic per thread                                                                                          |
| TRIMMOMATIC_minlen  | int | It is automatically filled. Notice that changing this parameter will affect the already established perc_remaining_reads in QC_1.yaml |

#### Pre-processing - host genome filtering. BWA maps reads to the host genome and the mapped reads are excluded from the FASTQ files by  Mseqtools.

| Parameters       |     | Description                           |
|------------------|-----|---------------------------------------|
| BWA_host_threads | int | Number of threads used by BWA.        |
| BWA_host_memory  | int | Number of GBs used by BWA per thread. |

#### Assembly-free taxonomy profiling. When <omics> = metaG, in <working_dir>/metaG/QC_2.yaml

| Parameters   |                                     | Description                                                                                                                                                                                                  |
|--------------|-------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| TAXA_threads | int                               | Number of threads used by the profiling tool                                                                                                                                                                 |
| TAXA_memory  | int                               | Number of GBs per thread used by the profiling tool                                                                                                                                                          |
| taxa_profile | <metaphlan, motus_rel or motus_raw> | High-quality filtered reads will be profiled by default with MetaPhlAn3 (metaphlan). An alternative tool is mOTUs2 in two different modes: relative abundance (motus_rel) or absolute abundance (motus_raw). |

#### Assembly-free taxonomy profiling. When <omics> = metaT, in <working_dir>/metaG/QC_2.yaml

| Parameters        |        | Description                                                                       |
|-------------------|--------|-----------------------------------------------------------------------------------|
| sortmeRNA_threads | int  | Number of threads used by sortMeRNA.                                              |
| sortmeRNA_memory  | int  | Number of GBs per thread used by sortMeRNA.                                       |
| sortmeRNA_db      | path | Location where the SortMeRNA databases will be downloaded          |
| sortmeRNA_db_idx  | path | Location where MIntO will output the indexed SortMeRNA databases and will be used |


QC_2.smk script can be run like:

- snakemake --snakefile <PathToMIntO>/QC_2.smk --cluster "qsub -pe smp {threads} -l h_vmem={resources.mem}G -N {name} -cwd" --latency-wait 30 --jobs <int> --configfile <working_dir>/<omics>/QC2.yaml --use-conda --conda-prefix <PathToTmpCondaEnv>/tmp

### Parameters used in assembly.smk should be included in assemblies.yaml

MIntO approach to reconstruct MAGs from high-quality host-free reads exploits metaG assembly of single samples as well as co-assembly of pre-defined sample groups followed by binning preparation and contig binning. 

The assembly can be as input short-reads, long-reads or a combination of both (hybrid).

MIntO generates the required configuration file in <working_dir>/<omics>/assemblies.yaml

#### MetaSPAdes settings - For the individual assembly of short-reads or hybrid reads.

| Parameters                |                            | Description                                                                                                             |
|---------------------------|----------------------------|-------------------------------------------------------------------------------------------------------------------------|
| METASPADES_qoffset        | <33 or 64> [default: auto] | PHRED quality offset in the input reads (33 or 64)                                                                      |
| METASPADES_threads        | int                      | Number of threads used by MetaSPADES                                                                                    |
| METASPADES_memory         | int                      | Number of GBs per thread used by MetaSPADES                                                                             |
| METASPADES_hybrid_max_k   | <33, 55, 77, 99 or 127>    | Maximum k-mer size (must be odd and less than 128). The minimum is set to 21. Option when running hybrid assembly.      |
| METASPADES_illumina_max_k | <33, 55, 77, 99 or 127>    | Maximum k-mer size (must be odd and less than 128). The minimum is set to 21. Option when running short-reads assembly. |

##### MEGAHIT settings - For the co-assembly

| Parameters      |                                 | Description                                                                                                                                                                                                                                                                                                       |
|-----------------|---------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| MEGAHIT_threads | int                          | Number of threads used by MEGAHIT                                                                                                                                                                                                                                                                                 |
| MEGAHIT_presets |  - meta-sensitive  - meta-large | List of MEGAHIT parameters to run per co-assembly. By default MEGAHIT is run with two preset parameters: meta-sensitive: '--min-count 1 --k-list 21,29,39,49,...,129,141'                                        meta-large: '--k-min 27 --k-max 127 --k-step 10'. For more details, please check MEGAHIT manual. |

Note:

Each MEGAHIT_presets parameter will be applied to each sample.			


#### MetaFlye settings - For the individual assembly of long-reads

| Parameters       |                                 | Description                                                                                                                                                                                                                                                                                                       |
|------------------|---------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| METAFLYE_presets | --meta --genome-size 3.0m       | List of Flye parameters to run per long-read assembly. By default Flye is run with one parameter: --meta --genome-size 3.0m For more details, please check Flye manual.                                                                                                                                           |

Notes:
1. Each METAFLYE_presets parameter will be applied to each sample. 
2. If nothing is listed in METAFLYE_presets, then MetaFlye won't be run.

#### Input data

| Parameters |                   | Description                                                                                                                                                                                                                                |
|------------|-------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| COASSEMBLY | <Subject1: I3+I4> | Illumina samples that will be co-assembled together using MEGAHIT. By default all the Illumina samples will be used. Please use the following definition: 'Subject1: I3+I4' will result in 1 co-assembly: 'Subject1' using I3 and I4 data. |
| NANOPORE   | < - N1 >      | List of nanopore samples that will be assembled individually using MetaFlye.                                                                                                                                                               |
| ILLUMINA   | < - I3 >      | List of Illumina samples that will be assembled individually using MetaSPADES.                                                                                                                                                             |

Note: 

Memory per coassembly is calculated to be 10G per sample in the coassembly. Please make sure that there is enough RAM on the server.

assembly.smk script can be run like:
- snakemake --snakefile <PathToMIntO>/assembly.smk --cluster "qsub -pe smp {threads} -l h_vmem={resources.mem}G -N {name} -cwd" --latency-wait 30 --jobs <int> --configfile <working_dir>/<omics>/assemblies.yam --use-conda --conda-prefix <PathToTmpCondaEnv>/tmp --use-singularity  --keep-going  --restart-times 1

### Parameters used in binning_preparation.smk should be included in assemblies.yaml

Contigs longer than 2500bp from all the combinations of assemblies in assembly.smk are combined together in preparation for binning. MetaG reads from individual short-read metagenomes are first mapped to this set of contigs using BWA aligner (Vasimuddin et al. 2019) 2.2.1 in paired-end mode. Sequencing depth of the contigs in each sample is estimated by jgi_summarize_bam_contig_depths program included in MetaBAT2 (Kang et al. 2019).

#### BWA settings - Used when mapping reads back to contigs

| Parameters  |       | Description                   |
|-------------|-------|-------------------------------|
| BWA_threads | int | Number of threads used by BWA |

#### samtools settings - Used when sorting bam files

| Parameters              |       | Description                               |
|-------------------------|-------|-------------------------------------------|
| SAMTOOLS_sort_threads   | int | Number of threads used by samtools        |
| SAMTOOLS_sort_memory_gb | int | Number of GBs per thread used by samtools |

binning_preparation.smk script can be run like:

- snakemake --snakefile <PathToMIntO>/binning_preparation.smk --cluster "qsub -pe smp {threads} -l h_vmem={resources.mem}G -N {name} -cwd" --latency-wait 30 --jobs <int> --configfile <working_dir>/<omics>/assemblies.yaml --use-conda --conda-prefix <PathToTmpCondaEnv>/tmp --use-singularity --keep-going  --restart-times 1

### Parameters used in mags_generation.smk should be included in mags_generation.yaml

Contig binning is then performed by executing VAMB (Nissen et al. 2021). When it is possible, GPU is highly recommended in order to speed up the binning process, especially if working with a large number of samples. By default, MIntO runs VAMB 4 times, each time with a different set of parameters. However the user can choose to perform just one run or a set of runs of their choice. Bins generated by VAMB are split into MAGs derived from individual metaG samples. Only the MAGs that passed quality control using CheckM (Parks et al. 2015) (completeness > 90% and contamination < 5% ) are kept. The MAGs are then subjected to cluster analysis performed with CoverM (v 0.6.0) (https://github.com/wwood/CoverM#usage, module cluster) in order to dereplicate them at 99% ANI (Jain et al. 2018). For each genome, a score is retrieved with the formula below.
assembly score = log10(longest contig length/ # contigs) + log10(N50/L50)
genome score = completeness - 2 * contamination
final score = 0.1 * genome score + assembly score 

Then for each cluster the genome with the highest score is chosen, generating a unique set of non-redundant MAGs which will be used in the next step.

Once the unique set of MAGs has been retrieved, taxonomy is assigned using the module phylophlan_metagenomic in PhyloPhlAn3 (Asnicar et al. 2020).

MIntO generates the required configuration file in <working_dir>/<omics>/mags_generation.yaml

| Parameters           |          | Description                                                                                                                                                                                          |
|----------------------|----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| BINNERS              | str    | The user can choose between: vamb_256, vamb_384, vamb_512 and vamb_768. default: run all of them                                                                                                     |
| VAMB_THREADS         | int    | The user can choose the number of threads that VAMB should use. If the GPU is not active, we recommend using a higher number of cores, to decrease the time of binning.  default: 8                  |
| VAMB_memory          | int    | number of GB to use for the memory default: 1                                                                                                                                                        |
| VAMB_GPU             | <yes/no> | if MIntO has to use the GPU or not. default: yes                                                                                                                                                     |
| MIN_FASTA_LENGTH     | str    | Minimum length of a bin to be kept in the next steps. default: “500,000”                                                                                                                             |
| CHECKM_THREADS       | int    | Threads to be used by checkm default: 8                                                                                                                                                              |
| PPLACER_THREADS      | int    | Threads to be used by pplacer  default: 8                                                                                                                                                            |
| CHECKM_memory        | int    | Gb used by checkm  default: 50                                                                                                                                                                       |
| CHECKM_COMPLETENESS  | int    | Only genomes with a completeness higher than this will be considered in the next steps default: 90                                                                                                   |
| CHECKM_CONTAMINATION | int    | Only genomes with a contamination lower than this will be considered default: 5                                                                                                                      |
| CLEAN_CHECKM         | yes/no> | Clean the checkm intermediates file.  default: yes                                                                                                                                                   |
| COVERM_THREADS       | int    | Cores to be used by coverm.  default: 10                                                                                                                                                             |
| COVERM_memory        | int    | memory to be used by Checkm default: 10                                                                                                                                                              |
| SCORE_METHOD         | <checkm> | This can only be checkm at the moment!                                                                                                                                                               |
| RUN_PROKKA           | <yes/no> | Run prokka on the unique genomes retrieved? default: yes Careful: if it is no, MIntO will not continue with the next rules, since it misses the input!                                               |
| PROKKA_CPUS          | int    | How many CPUs should be used by Prokka default: 10                                                                                                                                                   |
| PROKKA_memory        | int    | Gb Ram memory to be used by prokka default: 8                                                                                                                                                        |
| RUN_TAXONOMY         | <yes/no> | Run taxonomy labelling of the unique set of genomes using PhyloPhlAn3 default: yes                                                                                                                   |
| TAXONOMY_DATABASE    | str    | The database to download  or to be used default: SGB.Jan20                                                                                                                                           |
| TAXONOMY_CPUS        | int    | Cores to be used by PhyloPhlAn default: 10                                                                                                                                                           |
| TAXONOMY_memory      | int    | Gb Ram to be used default: 10                                                                                                                                                                        |
| DATABASE FOLDER      | str    | Path to the folder where to store the database, or where databases are already present.  We suggest having a folder with all the databases.  default: “.” ( same folder where MIntO has be launched) |

mags_generation.smk script can be run like:
- snakemake --snakefile <PathToMIntO>/mags_generation.smk --cluster "qsub -pe smp {threads} -l h_vmem={resources.mem}G -N {name} -cwd" --latency-wait 30 --jobs <int> --configfile <working_dir>/<omics>/mags_generation.yaml --use-conda --conda-prefix <PathToTmpCondaEnv>/tmp

### Parameters used in gene_annotation.smk should be included in mapping.yaml

This script should only be run with genome-based modes (map_reference = MAG or reference_genome).

It combines .faa and .bed files from the given/recovered genomes followed by function annotation. There are three tools used (EggNOG-mapper, KofamScan and dbCAN2) to annotate the genes with seven different functional databases: EggNOG (Yin et al. 2012; Huerta-Cepas et al. 2019), KEGG Pathways, Modules and KOs (Kanehisa and Goto 2000), dbCAN modules and enzyme classes (Yin et al. 2012), and Pfam (Mistry et al. 2021).

If map_reference = reference_genome (genome-based assembly-free mode) the user is using publicly available genomes. The required files are genome (.fna),  amino acid sequences of protein-coding genes (.faa) and genome features (.gff) files per genome. These files should be grouped per genome in a separate directory labeled as the corresponding genome name.

MIntO generates the required configuration file in <working_dir>/<omics>/mapping.yaml

| Parameters    |   | Description                                                                               |
|---------------|---|-------------------------------------------------------------------------------------------|
| reference_dir | path  | point to the directory where all the publicly available genome’s directories are located. |

gene_annotation.smk script can be run like:
- snakemake --snakefile <PathToMIntO>/gene_annotation.smk --cluster "qsub -pe smp {threads} -l h_vmem={resources.mem}G -N {name} -cwd" --latency-wait 30 --jobs <int> --configfile <working_dir>/<omics>/mapping.yaml --use-conda --conda-prefix <PathToTmpCondaEnv>/tmp

### Parameters used in gene_abundance.smk should be included in mapping.yaml

The high-quality filtered (host-free for metaG and host- and rRNA-free for metaT) reads are used to generate the function profiles following several steps: metaG and metaT read alignments and read count normalization.

#### BWA-MEM2 Alignment

| Parameters              |     | Description                                                                                                                  |
|-------------------------|-----|------------------------------------------------------------------------------------------------------------------------------|
| msamtools_filter_length | int | min. length of alignment [default: 50]                                                                                       |
| alignment_identity      | int | min. sequence identity of alignment, in percentage, integer between 0 and 100; requires NM field to be present [default: 95] |

#### Normalization approach

| Parameters              |             | Description                                                     |
|-------------------------|-------------|-----------------------------------------------------------------|
| abundance_normalization | <TPM or MG> | Type of abundance normalization (gene abundance or transcripts) |

#### Map reads to reference
| Parameters       |                                      | Description                                                                               |
|------------------|--------------------------------------|-------------------------------------------------------------------------------------------|
| map_reference    | <MAG, reference_genome  or genes_db> | Select the mode that MIntO should be run depending on the input given by the user:        |
| reference_dir    | path                               | point to the directory where all the publicly available genome’s directories are located. |
| BWAindex_threads | int                                | Number of threads used by BWA to generate the index.                                      |
| BWAindex_memory: | int                                | Number of GBs per thread used by BWA to generate the index.                               |
| BWA_threads      | int                                | Number of threads used by BWA.                                                            |
| BWA_memory       | int                                | Number of GBs per thread used by BWA.                                                     |
| ILLUMINA         |                                      | List of Illumina samples that will be aligned individually to the reference..             |

gene_abundance.smk script can be run like:
- snakemake --snakefile <PathToMIntO>/gene_abundance.smk --cluster "qsub -pe smp {threads} -l h_vmem={resources.mem}G -N {name} -cwd" --latency-wait 30 --jobs <int> --configfile <working_dir>/<omics>/mapping.yaml --use-conda --conda-prefix <PathToTmpCondaEnv>/tmp

### Parameters used in integration.smk should be included in data_integration.yaml

| Parameters              |                                      | Description                                                                |
|-------------------------|--------------------------------------|----------------------------------------------------------------------------|
| alignment_identity      | int                                |                                                                            |
| abundance_normalization | <TPM or MG>                          |                                                                            |
| map_reference           | <MAG, reference_genome  or genes_db> |                                                                            |
| MERGE_memory            | int                                | Number of threads used to generate the gene and function profiles.         |
| MERGE_threads           | int                                | Number of GBs per thread used  to generate the gene and function profiles. |
| ANNOTATION_file         | path                               |                                                                            |
| ANNOTATION_ids          |                                      | list                                                                       |

integration.smk script can be run like:
- snakemake --snakefile <PathToMIntO>/integration.smk --cluster "qsub -pe smp {threads} -l h_vmem={resources.mem}G -N {name} -cwd" --latency-wait 30 --jobs <int> --configfile <working_dir>/data_integration.yaml --use-conda --conda-prefix <PathToTmpCondaEnv>/tmp

Further analyses can be done using the output files. MIntO generates three different types of table: 
- (1) assembly-free and assembly-based taxonomic profiles
- (2) gene profiles, including the gene IDs (generated by Prokka (Beghini et al. 2021; Seemann 2014) when selecting assembly-dependent mode or sequence IDs when choosing assembly-free mode) and normalized gene abundance, transcript or expression
- (3) function profiles per database, including the function IDs, function description and function abundance, transcript or expression normalized counts. For an easier downstream analysis of these data, phyloseq objects are generated for the taxonomic, gene and function profiles.

MIntO also outputs several plots as preliminary results to help the user in the downstream analysis: 
- (1) relative abundance of the top 15 genera from the assembly-free taxonomic profile (Figure 2A)
- (2) PCoA plots using the taxonomy (Figure 2B), gene profiles (Figure 4A) or function profiles (Figure 4B) from the different databases
- (3) summary of genes and features per function database that are expressed in the data (Figure 3).


## Contributors

Carmen Saenz (1), Eleonora Nigro (1), Vithiagaran Gunalan (1), Manimozhiyan Arumugam (1)*

1 Novo Nordisk Foundation Center for Basic Metabolic Research, Faculty of Health and Medical Sciences, University of Copenhagen, Copenhagen, Denmark
*Correspondence: arumugam@sund.ku.dk ; Tel.: +45 35 33 75 81
