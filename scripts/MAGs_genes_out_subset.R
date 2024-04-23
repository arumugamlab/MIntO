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

# Group duplicated features
genomes_df$coord <- paste0(genomes_df$chr, '_',genomes_df$start, '_', genomes_df$stop)

genomes_sub <- subset(genomes_df, select=c('coord', 'name', 'source', 'feature'))

# rank feature sources
source_list <- unique(genomes_sub$source)
level.order <- c('Prodigal:002006','Aragorn:001002','barrnap:0.9',
                 'minced:0.4.2', 'prokka')
# include all sources
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
genomes_sub2 <- genomes_sub2[!duplicated(genomes_sub2[,'coord']),]

# join the coords details with the features
genomes_sub_df <- merge(genomes_sub2, genomes_sub_dplyr, by='coord')

genomes_sub_df$score <- '.'
genomes_sub_df$frame <- '.'
genomes_sub_df$strand <- '.'

#genomes_sub_df$source2 <- unlist(lapply(strsplit(as.character(genomes_sub_df$source),'\\,'), `[`, 1))
genomes_sub_df$genome <- unlist(lapply(strsplit(as.character(genomes_sub_df$chr),'\\|'), `[`, 1))
genomes_sub_df$contig <- unlist(lapply(strsplit(as.character(genomes_sub_df$chr),'\\|'), `[`, 2))
genomes_sub_df$name[genomes_sub_df$name == '.'] <- paste0(genomes_sub_df$contig[genomes_sub_df$name == '.'], '_',genomes_sub_df$start[genomes_sub_df$name == '.'], '_',genomes_sub_df$stop[genomes_sub_df$name == '.'])
genomes_sub_df$name2 <- unlist(lapply(strsplit(as.character(genomes_sub_df$name),'\\;'), `[`, 1))

genomes_sub_df$info <- paste0(genomes_sub_df$genome, '|',genomes_sub_df$name2)

genomes_sub_df <- as.data.frame(genomes_sub_df %>%
                                  dplyr::select(chr,start,stop,name,score,strand,source,feature,frame,info)) %>%
                                  dplyr::arrange(chr, start, stop)

write.table(genomes_sub_df,paste0(dir_out, file_name, "_SUBSET.bed"),sep="\t",row.names=FALSE, col.names = F, quote = F)

# delete from use
rm(genomes_df, genomes_sub, genomes_sub_dplyr, genomes_sub2)


#### #### ####  CD_transl file #### #### ####
## Load CDs transl sequences
fastafile_MAGgenes<- readAAStringSet(file = cd_transl_file)

## prepare data frame
fastafile_MAGgenes_df <- data.frame(names(fastafile_MAGgenes))
names(fastafile_MAGgenes_df) <- c('header')
fastafile_MAGgenes_df$ID <- unlist(lapply(strsplit(fastafile_MAGgenes_df$header, " "), `[`, 1))

# order can have changed so just plop ID back in
names(fastafile_MAGgenes) <- fastafile_MAGgenes_df$ID

# subselect features in BED
fastafile_genes_sub <- fastafile_MAGgenes[names(fastafile_MAGgenes) %in% genomes_sub_df$info]

writeXStringSet(fastafile_genes_sub , paste0(dir_out, file_name, "_translated_cds_SUBSET.faa"))
