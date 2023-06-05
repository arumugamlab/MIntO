#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(data.table)
library(tidyr) #
library(Biostrings)
library(stringr) #
library(rlang)

threads_n <- args[1]
dir_out <- args[2]
file_name <- args[3]

##########################  ** Load directories **  ########################## 
# threads_n <- 8
# bed_file <- "/emc/cbmr/users/rzv923/ibdmdb/DB/9-MAGs-prokka-post-analysis/GFF/MAGs_genes.bed"
# cd_transl_file <- "/emc/cbmr/users/rzv923/ibdmdb/DB/9-MAGs-prokka-post-analysis/CD_transl/MAGs_genes_translated_cds.faa"
# dir_out <- "/emc/cbmr/users/rzv923/ibdmdb/DB/9-MAGs-prokka-post-analysis/"
# file_name <- "MAGs_genes"
# path_out <- "/emc/cbmr/users/rzv923/ibdmdb/DB/9-MAGs-prokka-post-analysis/"

setDTthreads(threads = threads_n)
##########################  ** Read file **  ########################## 
# Before running fetchMG (to identify marker genes in genomes)
## From CDS tranls fasta file with ALL-SUBSET sequences, separate by genome ####
inFasta <- readAAStringSet(paste0(dir_out, file_name, "_translated_cds_SUBSET.faa"))

name_list <- sort(list.files(paste0(dir_out, 'fasta/')))
name_list <- name_list[name_list %like% ".fna"]

for (cds_file in (1:length(name_list))){
    cds_file_name <- name_list[cds_file]
    cds_file_genome <- gsub('.fna', '', cds_file_name)
    ## subset fasta file by matching seq names from fasta and genome ID
    inFasta_sub <- inFasta[names(inFasta) %like% cds_file_genome]

    ## Modify headers - replace '.' by '-'
    fa_given_names <- names(inFasta_sub)
    ## prepare data frame
    inFasta_sub_df <- data.frame(names(inFasta_sub))
    names(inFasta_sub_df) <- c('header')
    inFasta_sub_df$header_repl <- gsub('\\.', '\\-', inFasta_sub_df$header)
    inFasta_sub_df <- subset(inFasta_sub_df, select=c('header', 'header_repl'))

    ## assign new seq names by mapping fasta seq name to data frame names
    names(inFasta_sub) <- inFasta_sub_df[match(fa_given_names , inFasta_sub_df$header) , "header_repl"]

    writeXStringSet(inFasta_sub , paste0(dir_out, "CD_transl/",cds_file_genome, '_translated_cds_SUBSET.faa'))
}
cat(paste0("Ready to run FetchMGs"),file=paste0(dir_out, "CD_transl/", file_name, "_fetchMGs_translated_cds_SUBSET.txt"),sep="\n")
