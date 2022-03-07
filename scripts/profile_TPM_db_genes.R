#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if("dplyr" %in% rownames(installed.packages()) == FALSE) {install.packages("dplyr")}
library(dplyr)
if("data.table" %in% rownames(installed.packages()) == FALSE) {install.packages("data.table")}
library(data.table)

gene_abund_txt <- args[1]
gene_tpm_csv <- args[2]
omics <- args[3]
read_n <- args[4]
read_length_out <- args[5]

## Bed file colnames
gene_info <- c('Genes','gene_lenght')

# GENE ABUNDANCES per SAMPLE
gene_abund_txt_df <- as.data.frame(fread(gene_abund_txt, header=T), stringsAsFactors = F)
## Filter number of mapped reads bellow the threashold
gene_abund_txt_coord_df <- gene_abund_txt_df
gene_abund_txt_coord_df[,2:ncol(gene_abund_txt_coord_df)][gene_abund_txt_coord_df[,2:ncol(gene_abund_txt_coord_df)] <=read_n] <- 0

#gene_abund_txt<- "/emc/cbmr/users/rzv923/ibdmdb_test/metaG/6-mapping-profiles/BWA_reads-db_genes/genes_abundances.p95.TPM.csv"
#read_length_out <- "/emc/cbmr/users/rzv923/Databases/IGC/output_count_fasta_v2.txt"
# Read file with length of fasta sequences 
read_length_df <- as.data.frame(fread(read_length_out, header=F, fill=T, stringsAsFactors = F, check.names = F, strip.white = T, sep = ''), stringsAsFactors = F)
read_length_df <- data.frame(do.call('rbind', str_split(string = read_length_df$V1, pattern = '__')))
names(read_length_df) <- c('Genes','gene_lenght')
read_length_df$gene_lenght <- as.numeric(read_length_df$gene_lenght)

read_length_df$Genes <- unlist(lapply(strsplit(as.character(read_length_df$Genes),' '), `[`, 1))

## Subset df by colname
### Gene info
gene_info_bed_df <- read_length_df

### Gene abundances
gene_abund_samples_df <- merge(gene_abund_txt_coord_df, read_length_df, by ='Genes')
#### Filter gene abundances df by gene coordinates
gene_abund_samples_u_df <- gene_abund_samples_df[!duplicated(gene_abund_samples_df[,'Genes']),]
rownames(gene_abund_samples_u_df) <- gene_abund_samples_u_df$Genes
gene_abund_samples_u_df$Genes <- NULL

## Calculate RPK from gene abundances
gene_rpk_samples_u_df <- gene_abund_samples_u_df[1:ncol(gene_abund_samples_u_df)-1]/gene_abund_samples_u_df$gene_lenght
### Traslocate df to merge it later with total_reads_df
gene_rpk_samples_u_t_df <- as.data.frame(t(gene_rpk_samples_u_df))

## RPK sum per sample
rpk_sum <- as.data.frame(colSums(gene_rpk_samples_u_df))
names(rpk_sum) <- 'total_sum'

# Merge total_reads_df and gene_rpk_samples_u_t_df by rownames
gene_rpk_samples_total_reads_t_df <- cbind(gene_rpk_samples_u_t_df, rpk_sum)

## Calculate TPM from RPK and RPK sum per sample
gene_tpm_samples_u_t_df <- gene_rpk_samples_total_reads_t_df[1:ncol(gene_rpk_samples_total_reads_t_df)-1]/gene_rpk_samples_total_reads_t_df$total_sum
gene_tpm_samples_u_t_df <- gene_tpm_samples_u_t_df*1e6

gene_tpm_samples_u_noNA_df <- as.data.frame(t(gene_tpm_samples_u_t_df))
gene_tpm_samples_u_noNA_df[is.na(gene_tpm_samples_u_noNA_df)] <- 0

colSums(gene_tpm_samples_u_noNA_df)
gene_tpm_samples_u_noNA_df$coord <- rownames(gene_tpm_samples_u_noNA_df)
rownames(gene_tpm_samples_u_noNA_df) <- NULL

gene_tpm_samples_u_noNA_df <- as.data.frame(gene_tpm_samples_u_noNA_df %>% 
                                              dplyr::rename_at(vars(!matches("Genes")), ~ paste0(omics,".", .)))

# Merge with gene_info
gene_tpm_bed_df <- merge(gene_info_bed_df, gene_tpm_samples_u_noNA_df, by='Genes')
print(dim(gene_info_bed_df))
print(dim(gene_tpm_samples_u_noNA_df))
print(dim(gene_tpm_bed_df))

write.table(gene_tpm_bed_df, gene_tpm_csv, row.names = FALSE, col.names = T, sep = "\t", quote = F)
