#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(data.table)
library(tidyr)
library(Biostrings)
library(stringr)
library(rlang)

threads_n <- args[1]
bed_file <- args[2]
cd_transl_file <- args[3]
dir_out <- args[4]
file_name <- args[5]

##########################  ** Load directories **  ########################## 

setDTthreads(threads = threads_n)

##########################  ** Read file **  ########################## 
# Prokka genes - bed file ####

bed_info <- c("chr","start","stop","name","score","strand","source","feature","frame","info")
genomes_df = read.csv(bed_file, sep ="\t", header = F, stringsAsFactors = F)
names(genomes_df) <- bed_info
genomes_df$coord <- paste0(genomes_df$chr, '_',genomes_df$start, '_', genomes_df$stop)
genomes_df$coord <- gsub('-', '_', genomes_df$coord)
genomes_df$genome <- unlist(lapply(strsplit(as.character(genomes_df$chr),'\\|'), `[`, 1))
genomes_df$contig <- unlist(lapply(strsplit(as.character(genomes_df$chr),'\\|'), `[`, 2))
genomes_df$name2 <- genomes_df$name
#genomes_df$name2[genomes_df$source == 'macrel']  <-  paste0(genomes_df$contig[genomes_df$source  == 'macrel'], '_', unlist(lapply(strsplit(as.character(genomes_df$name[genomes_df$source == 'macrel']),'\\_'), `[`, 2)))

genomes_df$ID_gene <- paste0(genomes_df$genome, '|',genomes_df$name2)
# genomes_df$ID_gene[genomes_df$source %in% c('antismash_SBGC')] <- paste0(genomes_df$chr[genomes_df$source %in% c('antismash_SBGC')], '|',genomes_df$name2[genomes_df$source %in% c('antismash_SBGC')])
# genomes_df$ID_gene[genomes_df$source %in% c('gutsmash_PBGC')] <- paste0(genomes_df$chr[genomes_df$source %in% c('gutsmash_PBGC')], '|',genomes_df$name2[genomes_df$source %in% c('gutsmash_PBGC')])
# #genomes_df$ID_gene[genomes_df$source2 %in% c('macrel')] <- paste0(genomes_df$genome[genomes_df$source2 %in% c('macrel')], '|',genomes_df$name2[genomes_df$source2 %in% c('macrel')])
genomes_df$ID_gene[genomes_df$name == '.'] <- paste0(genomes_df$chr[genomes_df$name == '.'], '|',genomes_df$start[genomes_df$name == '.'], '_',genomes_df$stop[genomes_df$name == '.'])

genomes_df$name <- genomes_df$name2
genomes_df$info <- genomes_df$ID_gene

genomes_sub_df <- subset(genomes_df, select=bed_info)

write.table(genomes_sub_df,paste0(dir_out,"/", file_name,"_names_modif.bed"),sep="\t",row.names=FALSE, col.names = F, quote = F)

## Keep all antiSMASH and macrel genes
genomes_df <- as.data.frame(fread(paste0(dir_out, "/", file_name, "_names_modif.bed"), header=F), stringsAsFactors = F)
names(genomes_df) <- bed_info

# Group duplicated features
genomes_df$coord <- paste0(genomes_df$chr, '_',genomes_df$start, '_', genomes_df$stop)

genomes_sub <- subset(genomes_df, select=c('coord', 'name', 'source', 'feature'))

#unique(genomes_sub$source)

source_list <- unique(genomes_sub$source)
level.order <- c('Prodigal:002006','Aragorn:001002','barrnap:0.9',
                 'minced:0.4.2', 'prokka')
level.source <-c(intersect(level.order, source_list), setdiff(source_list,level.order))

genomes_sub$source <- factor(genomes_sub$source, levels = level.source)

genomes_sub <- genomes_sub[order(match(genomes_sub$source, level.source)),]
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

genomes_sub_df$info <- paste0(genomes_sub_df$genome, '|',genomes_sub_df$name2) #******** genomes_sub_df$info <- paste0(genomes_sub_df$chr, '|',genomes_sub_df$name2)

genomes_sub_df$strand <- '.'

genomes_sub_df <- as.data.frame(genomes_sub_df %>%
                                  dplyr::select(chr,start,stop,name,score,strand,source,feature,frame,info)) %>%
                                  dplyr::arrange(chr, start, stop)

write.table(genomes_sub_df,paste0(dir_out, file_name, "_SUBSET.bed"),sep="\t",row.names=FALSE, col.names = F, quote = F)

rm(genomes_df, genomes_df_u, genomes_sub_df, genomes_sub, genomes_sub_dplyr, genomes_sub2, genomes_sub2_u)

#### #### ####  CD_transl file #### #### ####
## Load bed file SUBSET
genomes_sub_df <- as.data.frame(fread(paste0(dir_out, file_name, "_SUBSET.bed"), header=F), stringsAsFactors = F)
names(genomes_sub_df) <- bed_info

#genomes_sub_df$coord <- paste0(genomes_sub_df$chr, '_',genomes_sub_df$start, '_', genomes_sub_df$stop)
#genomes_sub_df$coord <- gsub('-', '_', genomes_sub_df$coord)
#genomes_sub_df$source2 <- unlist(lapply(strsplit(as.character(genomes_sub_df$source),'\\,'), `[`, 1))

## Load CDs transl sequences
fastafile_MAGgenes<- readAAStringSet(file = cd_transl_file)

## get seq names from fasta
fa_given_names <- names(fastafile_MAGgenes)

## prepare data frame
fastafile_MAGgenes_df <- data.frame(names(fastafile_MAGgenes))
names(fastafile_MAGgenes_df) <- c('header')
fastafile_MAGgenes_df$ID <- unlist(lapply(strsplit(fastafile_MAGgenes_df$header, " "), `[`, 1))
# fastafile_MAGgenes_df$ID <- gsub('-', '_', fastafile_MAGgenes_df$ID) #*****

## assign new seq names by mapping fasta seq name to data frame names
names(fastafile_MAGgenes) <- fastafile_MAGgenes_df[match(fa_given_names , fastafile_MAGgenes_df$header) , "ID"]

fastafile_genes <- c(fastafile_MAGgenes)

fastafile_genes_sub <- fastafile_genes[names(fastafile_genes) %in% genomes_sub_df$info]
fastafile_sub <- c(fastafile_genes_sub)

writeXStringSet(fastafile_sub , paste0(dir_out, file_name, "_translated_cds_SUBSET.faa"))
