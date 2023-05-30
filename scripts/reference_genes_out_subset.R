#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(data.table)
library(tidyr) #
library(Biostrings)
library(stringr) #
library(rlang)

threads_n <- args[1]
bed_file <- args[2]
bed_file_short <- args[3]
cd_transl_file <- args[4]
dir_out <- args[5]
file_name <- args[6]

##########################  ** Load directories **  ########################## 
# threads_n <- 8
# #bed_file <- "/emc/cbmr/data/microbiome/processed/SHIME/Cdiff_sterile1/DB/9-reference-genes-post-analysis/GFF/reference_genes.bed"
# bed_file <- "/emc/cbmr/data/microbiome/processed/SHIME/Cdiff_sterile1/DB/9-reference-genes-post-analysis/GFF/reference_genes.header-modif.coord_correct.bed"
# bed_file_short <- "/emc/cbmr/data/microbiome/processed/SHIME/Cdiff_sterile1/DB/9-reference-genes-post-analysis/GFF/reference_genes.bed"
# cd_transl_file <- "/emc/cbmr/data/microbiome/processed/SHIME/Cdiff_sterile1/DB/9-reference-genes-post-analysis/CD_transl/reference_genes_translated_cds.faa"
# dir_out <- "/emc/cbmr/data/microbiome/processed/SHIME/Cdiff_sterile1/DB/9-reference-genes-post-analysis/"
# file_name <- "reference_genes"
# normalization <- "MG"
# path_out <- "/emc/cbmr/data/microbiome/processed/SHIME/Cdiff_sterile1/DB/9-reference-genes-post-analysis/"

setDTthreads(threads = threads_n)
##########################  ** Read file **  ########################## 
# Genomes genes - bed file ####

# Load bed file - with full information column ####
file_bed <- read.csv(bed_file, sep ="\t", header = F, stringsAsFactors = F)
# columns to not paste together
cols <- c( 'V1' , 'V2' , 'V3' ,'V4' ,'V5' ,'V6' , 'V7' ,'V8' ,'V9')
# create a new column `x` with the columns collapsed together that don't match cols
file_bed$x <- apply( file_bed[,!names(file_bed) %in% cols] , 1 , paste , collapse = " " )
# remove the unnecessary columns
file_bed <- file_bed[,names(file_bed) %in% c(cols, 'x')]

# Select only locus_tag ID
bed_info <- c("chr","start","stop","name","score","strand","source","feature","frame","info2")
names(file_bed) <- bed_info
file_bed$locus_tag <- unlist(lapply(strsplit(file_bed$info, "locus\\_tag\\="), `[`, 2))
file_bed$locus_tag <- unlist(lapply(strsplit(file_bed$locus_tag, "\\;"), `[`, 1))
# Returns string without trailing white space
file_bed$locus_tag <- gsub("\\s+$", "", file_bed$locus_tag)

# Load bed file - with only ID in information column ####
file <- read.csv(bed_file_short, sep ="\t", header = F, stringsAsFactors = F)
bed_info <- c("chr","start","stop","name","score","strand","source","feature","frame","info")
names(file) <- bed_info
# merge file_bed and file to replace ID by locus_tag ID when locus_tag is not NA 
# (In this way, it will be easier to match locus_tag IDs to aa seq headers in fasta file later for annotation)
file <- merge(file, file_bed, by=c("chr","start","stop","name","score","strand","source","feature","frame"))
file$name[!is.na(file$locus_tag)] <- file$locus_tag[!is.na(file$locus_tag)]

file$coord <- paste0(file$chr, '_',file$start, '_', file$stop)
feature_list <- unique(file$feature)
level.order <- c('gene','CDS','region','exon','rRNA','tRNA',
                 'binding_site','SRP_RNA','signal_peptide_region_of_CDS',
                 'sequence_feature','ribozyme','pseudogene','ncRNA',
                 'RNase_P_RNA','tmRNA' ,'Homology' ,'direct_repeat',
                 'riboswitch')
level.feature <-c(intersect(level.order, feature_list), setdiff(feature_list,level.order))

file_order <- file[order(match(file$feature, level.feature)),]
file_dplyr <- as.data.frame(file_order %>%
                              dplyr::group_by(coord)%>%
                              dplyr::summarise(feature = paste(feature, collapse=';'),
                                               name = paste(name, collapse=',')))

file_dplyr_gene <- file_dplyr[file_dplyr$feature %like% '^\\gene',]
file_dplyr_gene$name <- unlist(lapply(strsplit(file_dplyr_gene$name, "\\,"), `[`, 1))
file_dplyr_gene$coord_name <- paste(file_dplyr_gene$coord, file_dplyr_gene$name, sep='|')
file$coord_name <- paste(file$coord, file$name, sep='|')
genomes_df <- file[file$coord_name %in% file_dplyr_gene$coord_name,]
genomes_df$coord <- NULL
genomes_df$coord_name <- NULL
genomes_df$name <- gsub('gene-', '', genomes_df$name )
genomes_df$info <- gsub('gene-', '', genomes_df$info )
genomes_df$info2 <- NULL
genomes_df$locus_tag <- NULL

genomes_df$coord <- paste0(genomes_df$chr, '_',genomes_df$start, '_', genomes_df$stop)
genomes_df$coord <- gsub('-', '_', genomes_df$coord)
genomes_df$genome <- unlist(lapply(strsplit(as.character(genomes_df$chr),'\\|'), `[`, 1))
genomes_df$contig <- unlist(lapply(strsplit(as.character(genomes_df$chr),'\\|'), `[`, 2))
genomes_df$name2 <- genomes_df$name

genomes_df$ID_gene <- paste0(genomes_df$genome, '|',genomes_df$name2)

genomes_df$ID_gene[genomes_df$name == '.'] <- paste0(genomes_df$chr[genomes_df$name == '.'], '|',genomes_df$start[genomes_df$name == '.'], '_',genomes_df$stop[genomes_df$name == '.'])

genomes_df$name <- genomes_df$name2
genomes_df$info <- genomes_df$ID_gene
bed_info <- c("chr","start","stop","name","score","strand","source","feature","frame","info")

genomes_sub_df <- subset(genomes_df, select=bed_info)

write.table(genomes_sub_df,paste0(dir_out,"/GFF/", file_name,"_names_modif.bed"),sep="\t",row.names=FALSE, col.names = F, quote = F)
rm(genomes_sub_df)

# Read bed file (modified names)
genomes_df <- as.data.frame(fread(paste0(dir_out, "/GFF/", file_name, "_names_modif.bed"), header=F), stringsAsFactors = F)
bed_info <- c("chr","start","stop","name","score","strand","source","feature","frame","info")
names(genomes_df) <- bed_info

# Group duplicated features
genomes_df$coord <- paste0(genomes_df$chr, '_',genomes_df$start, '_', genomes_df$stop)

genomes_sub <- subset(genomes_df, select=c('coord', 'name', 'source', 'feature'))

genomes_sub <- genomes_sub[order(genomes_sub$source),]
genomes_sub$source <- as.character(genomes_sub$source)
genomes_sub_dplyr <- as.data.frame(genomes_sub %>%
                                     dplyr::group_by(coord)%>%
                                     dplyr::summarise(name = paste(unique(name), collapse=';'),
                                                      source = paste(unique(source), collapse=','),
                                                      feature = paste(unique(feature), collapse='|')))

genomes_sub2 <- subset(genomes_df, select=c('coord','chr', 'start', 'stop'))
genomes_sub2_u <- genomes_sub2[!duplicated(genomes_sub2[,'coord']),]

genomes_df_u <- merge(genomes_sub2_u, genomes_sub_dplyr, by='coord')
#unique(genomes_df_u$source)

#length(unique(genomes_df_u$coord))==nrow(genomes_df_u)

genomes_sub_df <- genomes_df_u
#unique(genomes_sub_df$source)
#names(genomes_sub_df)

genomes_sub_df$score <- '.'
genomes_sub_df$frame <- '.'

genomes_sub_df$source2 <- unlist(lapply(strsplit(as.character(genomes_sub_df$source),'\\,'), `[`, 1))
genomes_sub_df$genome <- unlist(lapply(strsplit(as.character(genomes_sub_df$chr),'\\|'), `[`, 1))
genomes_sub_df$contig <- unlist(lapply(strsplit(as.character(genomes_sub_df$chr),'\\|'), `[`, 2))
genomes_sub_df$name[genomes_sub_df$name == '.'] <- paste0(genomes_sub_df$contig[genomes_sub_df$name == '.'], '_',genomes_sub_df$start[genomes_sub_df$name == '.'], '_',genomes_sub_df$stop[genomes_sub_df$name == '.'])
genomes_sub_df$name2 <- unlist(lapply(strsplit(as.character(genomes_sub_df$name),'\\;'), `[`, 1))

#genomes_sub_df$info <- paste0(genomes_sub_df$genome, '|',genomes_sub_df$name2) #******** genomes_sub_df$info <- paste0(genomes_sub_df$chr, '|',genomes_sub_df$name2)
genomes_sub_df$info <- paste0(genomes_sub_df$chr, '|',genomes_sub_df$name2)

genomes_sub_df$strand <- '.'

genomes_sub_df <- as.data.frame(genomes_sub_df %>%
                                  dplyr::select(chr,start,stop,name,score,strand,source,feature,frame,info))

write.table(genomes_sub_df,paste0(dir_out, file_name, "_SUBSET.bed"),sep="\t",row.names=FALSE, col.names = F, quote = F)

rm(genomes_df, genomes_df_u, genomes_sub_df, genomes_sub, genomes_sub_dplyr, genomes_sub2, genomes_sub2_u,
   file, file_bed, file_dplyr, file_dplyr_gene, file_order)

#### #### ####  CD_transl file #### #### ####
## Load bed file SUBSET
genomes_sub_df <- as.data.frame(fread(paste0(dir_out, file_name, "_SUBSET.bed"), header=F), stringsAsFactors = F)
bed_info <- c("chr","start","stop","name","score","strand","source","feature","frame","info")
names(genomes_sub_df) <- bed_info

# genomes_sub_df$coord <- paste0(genomes_sub_df$chr, '_',genomes_sub_df$start, '_', genomes_sub_df$stop)
# genomes_sub_df$coord <- gsub('-', '_', genomes_sub_df$coord)
# genomes_sub_df$source2 <- unlist(lapply(strsplit(as.character(genomes_sub_df$source),'\\,'), `[`, 1))

#genomes_sub_df$coord_name <- paste(genomes_sub_df$chr, genomes_sub_df$name, sep='|')

## Load CDs transl sequences
fastafile_REFgenes<- readAAStringSet(file = cd_transl_file)

## get sequence names from fasta
fa_given_names <- names(fastafile_REFgenes)

## prepare data frame
## Get locus_tag ID 
## (In this way, it will be easier to match locus_tag IDs from bed file to aa seq headers)
fastafile_REFgenes_df <- data.frame(names(fastafile_REFgenes))
names(fastafile_REFgenes_df) <- c('header')
fastafile_REFgenes_df$ID <- unlist(lapply(strsplit(fastafile_REFgenes_df$header, "locus_tag\\="), `[`, 1))
fastafile_REFgenes_df$info <- unlist(lapply(strsplit(fastafile_REFgenes_df$header, "locus_tag\\="), `[`, 2))
fastafile_REFgenes_df$ID <- unlist(lapply(strsplit(as.character(fastafile_REFgenes_df$ID),' '), `[`, 1))
fastafile_REFgenes_df$ID <- unlist(lapply(strsplit(as.character(fastafile_REFgenes_df$ID),'\\_prot\\_'), `[`, 1))
fastafile_REFgenes_df$info <- unlist(lapply(strsplit(as.character(fastafile_REFgenes_df$info),'\\]'), `[`, 1))

fastafile_REFgenes_df$header_repl <- paste(fastafile_REFgenes_df$ID, fastafile_REFgenes_df$info, sep='|')
fastafile_REFgenes_df <- subset(fastafile_REFgenes_df, select=c('header', 'header_repl'))

## assign new seq names  by mapping fasta sequence name to data frame names
names(fastafile_REFgenes) <- fastafile_REFgenes_df[match(fa_given_names , fastafile_REFgenes_df$header) , "header_repl"]
## write data to fasta file with updated names
writeXStringSet(fastafile_REFgenes , paste0(dir_out,'CD_transl/', file_name, '_translated_cds.header-modif.faa'))

fastafile_REFgenes2<- readAAStringSet(file = paste0(dir_out,'CD_transl/', file_name, '_translated_cds.header-modif.faa'))

## get seq names from fasta
fa_given_names2 <- names(fastafile_REFgenes2)

## prepare data frame
fastafile_REFgenes2_df <- data.frame(names(fastafile_REFgenes2))
names(fastafile_REFgenes2_df) <- c('header')

## Select fasta sequence headers in bed file
fastafile_REFgenes2_sub <- fastafile_REFgenes2[names(fastafile_REFgenes2) %in% unique(genomes_sub_df$info)]

fastafile_genes2 <- c(fastafile_REFgenes2_sub)

writeXStringSet(fastafile_genes2 , paste0(dir_out, file_name, "_translated_cds_SUBSET.faa"))

# # remove duplicated headers from fasta file
# # bash
# cd /emc/cbmr/data/microbiome/processed/SHIME/Diet1/DB/9-MAGs-post-analysis-prokka/CD_transl
# awk '/^>/ { f = !($0 in a); a[$0] } f==1 { print $0 }' genes-PBGC-SBGC-smORFs_SUBSET_translated_cds.faa > genes-PBGC-SBGC-smORFs_SUBSET_translated_cds_no-dupl.faa

# if (normalization == 'MG'){
#   # Before running fetchMG (to identify marker genes in genomes)
#   ## From CDS tranls fasta file with ALL-SUBSET sequences, separate by genome ####
#   inFasta <- readAAStringSet(paste0(dir_out, file_name, "_translated_cds_SUBSET.faa"))
  
#   name_list <- sort(list.files(paste0(dir_out, 'fasta/')))
#   name_list <- name_list[name_list %like% ".fna"]
  
#   for (cds_file in (1:length(name_list))){
#     cds_file_name <- name_list[cds_file]
#     cds_file_genome <- gsub('.fna', '', cds_file_name)
#     ## subset fasta file by matching seq names from fasta and genome ID
#     inFasta_sub <- inFasta[names(inFasta) %like% cds_file_genome]
    
#     ## Modify headers - replace '.' by '-'
#     fa_given_names <- names(inFasta_sub)
#     ## prepare data frame
#     inFasta_sub_df <- data.frame(names(inFasta_sub))
#     names(inFasta_sub_df) <- c('header')
#     inFasta_sub_df$header_repl <- gsub('\\.', '\\-', inFasta_sub_df$header)
#     inFasta_sub_df <- subset(inFasta_sub_df, select=c('header', 'header_repl'))
    
#     ## assign new seq names by mapping fasta seq name to data frame names
#     names(inFasta_sub) <- inFasta_sub_df[match(fa_given_names , inFasta_sub_df$header) , "header_repl"]
    
#     writeXStringSet(inFasta_sub , paste0(dir_out, "CD_transl/",cds_file_genome, '_translated_cds_SUBSET.faa'))
#   }
#   cat(paste0("Ready to run FetchMGs"),file=paste0(out_dir, "CD_transl/fetchMGs_translated_cds_SUBSET.txt"),sep="\n")
# }
