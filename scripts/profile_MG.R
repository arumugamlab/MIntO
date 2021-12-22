#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


threads_n <- args[1]
memory_lim <- args[2]
gene_abund_bed <- args[3]
fetchMG_table <- args[4]
gene_MG_csv <- args[5]
omics <- args[6]

##########################  ** Load libraries **  ########################## 
#if("dplyr" %in% rownames(installed.packages()) == FALSE) {install.packages("dplyr")}
library(dplyr)
#if("data.table" %in% rownames(installed.packages()) == FALSE) {install.packages("data.table")}
library(data.table)

setDTthreads(threads = threads_n)

# #################
# #** Load directories **
# dir_path <- "/emc/cbmr/users/rzv923/ibdmdb/" #"/emc/cbmr/data/microbiome/processed/SHIME/Diet1/"
# dir_fetchMG <- paste0(dir_path,"DB/9-MAGs-post-analysis-prokka/fetchMG/")
# omics <- 'metaT'
# gene_abund_bed <- paste0(dir_path,omics,"/6-mapping-profiles/BWA_reads-MAGs_genes/genes_abundances.p95.bed")
# fetchMG_table <- paste0(dir_path,"/DB/9-MAGs-prokka-post-analysis//all.marker_genes_scores.table")
# gene_MG_csv <- paste0(dir_path,omics,"/6-mapping-profiles/BWA_reads-MAGs_genes/genes_abundances.p95.MG.csv")
# 
# # Load data - fetchMG output
# fetchMG_table <- paste0(dir_fetchMG,"all.marker_genes_scores.table")
# #################

fetchMG_df <- as.data.frame(fread(fetchMG_table, header=T), stringsAsFactors = F)
names(fetchMG_df) <- c('header', 'HMM_score', 'COG', 'taxid.projectid')

# replace '-' by '.'
fetchMG_df$header <- gsub('\\-','\\_', fetchMG_df$header)

fetchMG_sub_df <- subset(fetchMG_df, select=c('header','COG')) #c('ID', 'ID_MAG', 'ID_gene', 'COG'))
# Subset COGs
fetchMG_sub_df <- fetchMG_sub_df[fetchMG_sub_df$COG %in% c('COG0012', 'COG0016', 'COG0018', 'COG0172', 'COG0215', 
                                                           'COG0495','COG0525', 'COG0533', 'COG0541', 'COG0552'),]

# Load data - raw counts
## GENE ABUNDANCES per SAMPLE
gene_abund_bed_df <- as.data.frame(fread(gene_abund_bed, header=T), stringsAsFactors = F)
gene_abund_bed_coord_df <- gene_abund_bed_df
gene_abund_bed_coord_df$coord <- paste0(gene_abund_bed_coord_df$chr, '_',gene_abund_bed_coord_df$start, '_', gene_abund_bed_coord_df$stop)
gene_abund_bed_coord_df$coord <- gsub('\\-', '\\_', gene_abund_bed_coord_df$coord)

### Calculate gene length
gene_abund_bed_coord_df$gene_lenght <- abs(gene_abund_bed_coord_df$stop- gene_abund_bed_coord_df$start) +1
### Subset df by colname
gene_info <- c("chr","start","stop","name","score","strand","source","feature","frame","info")

### Gene abundances
gene_abund_samples_df <- gene_abund_bed_coord_df[!colnames(gene_abund_bed_coord_df) %in% gene_info]
#### Filter gene abundances df by gene coordinates
gene_abund_samples_u_df <- gene_abund_samples_df[!duplicated(gene_abund_samples_df[,'coord']),]
rownames(gene_abund_samples_u_df) <- gene_abund_samples_u_df$coord
gene_abund_samples_u_df$coord <- NULL

## Calculate RPK from gene abundances (gene abundances normalized by gene length)
gene_rpk_samples_u_df <- gene_abund_samples_u_df[1:ncol(gene_abund_samples_u_df)-1]/gene_abund_samples_u_df$gene_lenght
gene_rpk_samples_u_df$coord <- rownames(gene_rpk_samples_u_df)
### Gene info
gene_info_bed_df <- gene_abund_bed_coord_df[colnames(gene_abund_bed_coord_df) %in% c(gene_info, "coord", "gene_lenght")]
colnames(gene_info_bed_df) <- c("chr","start","stop","name","score","strand","source","feature","frame","info", "coord", "gene_lenght")
gene_rpk_bed_df <- merge(gene_info_bed_df, gene_rpk_samples_u_df, by='coord')
#dim(gene_rpk_bed_df)
gene_rpk_bed_df$header <- gsub('\\-', '\\_', gene_rpk_bed_df$info)

gene_rpk_bed_df$header[gene_rpk_bed_df$name=='.'] <- gene_rpk_bed_df$coord[gene_rpk_bed_df$name=='.']

######## Merge fetchMG_sub_df and gene_rpk_bed_df
gene_rpk_bed_df$ID_MAG <- unlist(lapply(stringr::str_split(pattern = '\\|', string = gene_rpk_bed_df$header), `[`, 1))
gene_rpk_bed_df$ID_MAG <- gsub('\\.', '\\_', gene_rpk_bed_df$ID_MAG)
gene_rpk_bed_df$header <- gsub('\\.', '\\_', gene_rpk_bed_df$header)
fetchMG_sub_df$ID_MAG <- unlist(lapply(stringr::str_split(pattern = '\\|', string = fetchMG_sub_df$header), `[`, 1))
gene_rpk_fetchMG <- merge(fetchMG_sub_df, gene_rpk_bed_df, by=c('ID_MAG','header'))
#setdiff(fetchMG_sub_df$ID_gene,gene_rpk_fetchMG$ID_gene)


gene_rpk_fetchMG_sub <- dplyr::select(gene_rpk_fetchMG, -c("header","COG","coord","chr", "start","stop","name","score","strand","source","feature","frame","info", "gene_lenght"))
# Calculate median abundance of the 10 single-copy marker genes
gene_rpk_fetchMG_sub_median <- as.data.frame(gene_rpk_fetchMG_sub %>%
                                               dplyr::group_by(ID_MAG) %>%
                                               dplyr::summarise(across(everything(),median)))

gene_rpk_fetchMG_sub_median$header <- 'MG'
gene_rpk_fetchMG_sub_median <- as.data.frame(gene_rpk_fetchMG_sub_median %>%
                                               dplyr::select(ID_MAG, header, everything()))

gene_rpk_bed_sub_df <-  dplyr::select(gene_rpk_bed_df, -c("coord","chr", "start","stop","name","score","strand","source","feature","frame","info", "gene_lenght"))
gene_rpk_bed_sub_df <- as.data.frame(gene_rpk_bed_sub_df %>%
                                       dplyr::select(ID_MAG, header, everything()))

# Merge median abundance of the 10 single-copy marker genes and transcript abundance
gene_rpk_bed_fetchMG <- rbind(gene_rpk_fetchMG_sub_median,gene_rpk_bed_sub_df)

names_file <-  names(gene_rpk_bed_fetchMG)[!names(gene_rpk_bed_fetchMG) %in% c('ID_MAG', 'header')]

gene_rpk_bed_fetchMG_melt <- reshape2::melt(gene_rpk_bed_fetchMG, id.vars = c("ID_MAG", "header"))
gene_rpk_bed_fetchMG_melt_percell <- as.data.frame(gene_rpk_bed_fetchMG_melt %>% 
                                                     group_by(ID_MAG) %>%
                                                     mutate(value = value/value[header == 'MG']))
# Remove MG rows
gene_rpk_bed_fetchMG_melt_percell <- gene_rpk_bed_fetchMG_melt_percell[gene_rpk_bed_fetchMG_melt_percell$header != 'MG',] 
# Replace NaN and infinite values by 0
gene_rpk_bed_fetchMG_melt_percell$value[is.infinite(gene_rpk_bed_fetchMG_melt_percell$value)] <- 0
gene_rpk_bed_fetchMG_melt_percell$value[is.nan(gene_rpk_bed_fetchMG_melt_percell$value)] <- 0

gene_rpk_bed_fetchMG_percell <- reshape2::dcast(gene_rpk_bed_fetchMG_melt_percell, ID_MAG + header ~ variable,  value.var="value") #fun.aggregate = list,

gene_rpk_bed_fetchMG_percell_norm <- gene_rpk_bed_fetchMG_percell
gene_rpk_bed_fetchMG_percell_norm$ID_MAG <- NULL
#rownames(gene_rpk_bed_fetchMG_percell_norm) <- gene_rpk_bed_fetchMG_percell_norm$header
#gene_rpk_bed_fetchMG_percell_norm$header <- NULL
#gene_rpk_bed_fetchMG_percell_norm[] <-  lapply(gene_rpk_bed_fetchMG_percell_norm, function(X) (X/max(X))*1e6) # (X/sum(X))*1e6 or (X/max(X))*1e6
#gene_rpk_bed_fetchMG_percell_norm$header <- rownames(gene_rpk_bed_fetchMG_percell_norm)

gene_rpk_bed_fetchMG_percell_norm <- as.data.frame(gene_rpk_bed_fetchMG_percell_norm %>% 
                                                dplyr::rename_at(vars(!matches(c("header"))), ~ paste0(omics,".", .)))

# Merge with gene_info
gene_info_bed_df <- gene_rpk_bed_df[colnames(gene_rpk_bed_df) %in% c('header','coord', gene_info, 'gene_lenght')]
gene_rpk_bed_df <- merge(gene_info_bed_df, gene_rpk_bed_fetchMG_percell_norm, by='header')
gene_rpk_bed_df$header <- NULL
# print(dim(gene_info_bed_df))
# print(dim(gene_rpk_bed_fetchMG_percell_norm))
# print(dim(gene_rpk_bed_df))

write.table(gene_rpk_bed_df, gene_MG_csv, row.names = F, col.names = T, sep = "\t", quote = F)

