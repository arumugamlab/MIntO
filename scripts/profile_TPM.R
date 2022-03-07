#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if("dplyr" %in% rownames(installed.packages()) == FALSE) {install.packages("dplyr")}
library(dplyr)
if("data.table" %in% rownames(installed.packages()) == FALSE) {install.packages("data.table")}
library(data.table)

gene_abund_bed <- args[1]
gene_tpm_csv <- args[2]
omics <- args[3]
read_n <- args[4]
read_n <- as.numeric(read_n)


## Bed file colnames
gene_info <- c("chr","start","stop","name","score","strand","source","feature","frame","info")

# GENE ABUNDANCES per SAMPLE
gene_abund_bed_df <- as.data.frame(fread(gene_abund_bed, header=T), stringsAsFactors = F)
## Filter number of mapped reads bellow the threashold
gene_abund_bed_coord_df <- gene_abund_bed_df
gene_abund_bed_coord_df[,!colnames(gene_abund_bed_coord_df) 
                        %in% gene_info][gene_abund_bed_coord_df[,!colnames(gene_abund_bed_coord_df)
                                                                %in% gene_info] <=read_n] <- 0
gene_abund_bed_coord_df$coord <- paste0(gene_abund_bed_coord_df$chr, '_',gene_abund_bed_coord_df$start, '_', gene_abund_bed_coord_df$stop)
#gene_abund_bed_coord_df$coord[!gene_abund_bed_coord_df$name == '.'] <- gene_abund_bed_coord_df$name[!gene_abund_bed_coord_df$name == '.']
gene_abund_bed_coord_df$coord <- gsub('-', '_', gene_abund_bed_coord_df$coord)

gene_abund_bed_coord_df$gene_lenght <- abs(gene_abund_bed_coord_df$stop- gene_abund_bed_coord_df$start) +1
## Subset df by colname
### Gene info
gene_info_bed_df <- gene_abund_bed_coord_df[colnames(gene_abund_bed_coord_df) %in% c(gene_info, "coord", "gene_lenght")]
colnames(gene_info_bed_df) <- c("chr","start","stop","name","score","strand","source","feature","frame","info", "coord", "gene_lenght")
### Gene abundances
gene_abund_samples_df <- gene_abund_bed_coord_df[!colnames(gene_abund_bed_coord_df) %in% gene_info]
#### Filter gene abundances df by gene coordinates
gene_abund_samples_u_df <- gene_abund_samples_df[!duplicated(gene_abund_samples_df[,'coord']),]
rownames(gene_abund_samples_u_df) <- gene_abund_samples_u_df$coord
gene_abund_samples_u_df$coord <- NULL

## Calculate RPK from gene abundances
gene_rpk_samples_u_df <- gene_abund_samples_u_df[1:ncol(gene_abund_samples_u_df)-1]/gene_abund_samples_u_df$gene_lenght
### Traslocate df to merge it later with total_reads_df
gene_rpk_samples_u_t_df <- as.data.frame(t(gene_rpk_samples_u_df))

## RPK sum per sample
rpk_sum <- as.data.frame(colSums(gene_rpk_samples_u_df))
names(rpk_sum) <- 'total_sum'

# # TOTAL READS per SAMPLE
# total_reads_df <- as.data.frame(fread(total_reads, header=F), stringsAsFactors = F)
# colnames(total_reads_df) <- c('sample', 'total_reads')
# # colnames(total_reads_df) <- c('sample', 'reads_info')
# # total_reads_df$total_reads <- as.numeric(gsub(" +.*", "", total_reads_df$reads_info))
# ## Calculate RPK and SCALING FACTOR from TOTAL READS
# total_reads_df$scaling_factor <- total_reads_df$total_reads/1e6
# rownames(total_reads_df) <- total_reads_df$sample
# total_reads_df$sample <- NULL
# total_reads_df$reads_info <- NULL
# total_reads_df$total_reads <- NULL

# Merge total_reads_df and gene_rpk_samples_u_t_df by rownames
gene_rpk_samples_total_reads_t_df <- cbind(gene_rpk_samples_u_t_df, rpk_sum)
# ## Calculate TPM from RPK and SCALING FACTOR
# gene_tpm_samples_u_t_df <- gene_rpk_samples_total_reads_t_df[1:ncol(gene_rpk_samples_total_reads_t_df)-1]/gene_rpk_samples_total_reads_t_df$scaling_factor

## Calculate TPM from RPK and RPK sum per sample
gene_tpm_samples_u_t_df <- gene_rpk_samples_total_reads_t_df[1:ncol(gene_rpk_samples_total_reads_t_df)-1]/gene_rpk_samples_total_reads_t_df$total_sum
gene_tpm_samples_u_t_df <- gene_tpm_samples_u_t_df*1e6

gene_tpm_samples_u_noNA_df <- as.data.frame(t(gene_tpm_samples_u_t_df))
gene_tpm_samples_u_noNA_df[is.na(gene_tpm_samples_u_noNA_df)] <- 0

colSums(gene_tpm_samples_u_noNA_df)
gene_tpm_samples_u_noNA_df$coord <- rownames(gene_tpm_samples_u_noNA_df)
rownames(gene_tpm_samples_u_noNA_df) <- NULL

gene_tpm_samples_u_noNA_df <- as.data.frame(gene_tpm_samples_u_noNA_df %>% 
                                              dplyr::rename_at(vars(!matches("coord")), ~ paste0(omics,".", .)))

# Merge with gene_info
gene_tpm_bed_df <- merge(gene_info_bed_df, gene_tpm_samples_u_noNA_df, by='coord')
print(dim(gene_info_bed_df))
print(dim(gene_tpm_samples_u_noNA_df))
print(dim(gene_tpm_bed_df))

write.table(gene_tpm_bed_df, gene_tpm_csv, row.names = FALSE, col.names = T, sep = "\t", quote = F)
