#!/usr/bin/env Rscript

# '''
# Generate gene expression profile from genome-based mode gene profiles 
#
# Authors: Carmen Saenz, Mani Arumugam
#
# '''

args = commandArgs(trailingOnly=TRUE)

##########################  ** Load arguments **  ########################## 
threads_n <- args[1]
integration_dir <- args[2]
omics <- args[3] #'metaG_metaT'
annot_file <- args[4]
metadata_file <- args[5]
input_file <- args[6]
annot_arg  <- args[7]
annot_names <- unlist(strsplit(annot_arg, split = "\\,")[[1]])
print(annot_names)

##########################  ** Load libraries **  ########################## 
library(dplyr)
library(data.table)
library(phyloseq)
library(KEGGREST)
library(ggplot2)
library(ggrepel)

#library(unix)

setDTthreads(threads = threads_n)
set.seed(1234)
##########################  ** Load functions **  ########################## 
# plot_PCoA <- function(distance_lab, data_phyloseq, color, label){ #out_name,title_name,
#   label2 <- noquote(label)
#   #### PCoA
#   ord<- ordinate(data_phyloseq, method = "PCoA", distance = distance_lab)
#   PcoA_Sample_site_abundance <- plot_ordination(data_phyloseq, ord, color = color, shape=NULL, label=label) #type = type,
#   PcoA_Sample_site_abundance <- PcoA_Sample_site_abundance +
#     ggtitle(title_name, distance_lab)  +
#     geom_point(size = 2.5) +
#     theme_bw() +
#     theme(legend.position = "top")
#   PcoA_Sample_site_abundance$layers[[1]] <- NULL
#   PcoA_Sample_site_abundance$layers[[2]] <- NULL
#   PcoA_Sample_site_abundance$layers[[3]] <- NULL
#   PcoA_Sample_site_abundance_2<- PcoA_Sample_site_abundance +
#     geom_text_repel(aes(label = sample_alias), size = 3.0, segment.alpha = 0.5, max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
#     geom_point(size = 2, shape = 16)
#   PcoA_Sample_site_abundance_2$layers[[1]] <- NULL
#   return(PcoA_Sample_site_abundance_2)
# }

plot_PCA <- function(data_phyloseq, color, label){ 
  ##### metadata:
  sample_data_df <- data.frame(sample_data(data_phyloseq), stringsAsFactors = F)
  ### Counts
  otu_table_df <-  as.data.frame(unclass(otu_table(data_phyloseq)), stringsAsFactors = F)
  ### Counts - Log transformed
  otu_table_df_rowsum_0 = c(which(rowSums(otu_table_df)==0))
  otu_table_df_rowsum_0_list <- rownames(otu_table_df[otu_table_df_rowsum_0,])
  otu_table_df_rowsum_not_0 <- otu_table_df[! rownames(otu_table_df) %in% otu_table_df_rowsum_0_list,]
  otu_table_df_log <- log(otu_table_df_rowsum_not_0)
  otu_table_df_log$ID <- rownames(otu_table_df_log)
  
  ##### Replace infinite values by log2(min value * 1e-2) for the PCA
  replace_val <- log((min(otu_table_df[otu_table_df >0 & !is.na(otu_table_df)]))*1e-2)
  otu_table_df_log_noinf2 <- otu_table_df_log %>% dplyr::mutate_all(function(x) ifelse(is.infinite(x), replace_val, x))
  rownames(otu_table_df_log_noinf2) <- otu_table_df_log_noinf2$ID
  otu_table_df_log_noinf2$ID <- NULL
  
  min(otu_table_df_log_noinf2[!is.na(otu_table_df_log_noinf2)])
  
  #Plotting scores of PC1 and PC2 with log transformation
  gene_expression_pca <- prcomp(t(otu_table_df_log_noinf2), center = T, sca=T)
  
  # ***************** ggplot way - w coloring *****************
  #library(ggplot2)
  dtp <- data.frame('sample_alias' = sample_data_df$sample_alias, 
                    'condition' = sample_data_df$condition, 
                    gene_expression_pca$x[,1:2]) # the first two componets are selected (NB: you can also select 3 for 3D plottings or 3+)
  df_out <- as.data.frame(gene_expression_pca$x)
  percentage <- round(gene_expression_pca$sdev / sum(gene_expression_pca$sdev) * 100, 2)
  percentage <- paste0( colnames(df_out), " (", paste0(as.character(percentage), "%", ")"))
  
  PCA_Sample_site_abundance <- (ggplot(data = dtp, aes(x = PC1, y = PC2, color = condition)) +
                                  geom_text_repel(aes(label = sample_alias),nudge_x = 0.04, size = 3.5, segment.alpha = 0.5) +
                                  geom_point(size = 2, shape = 16)+
                                  xlab(percentage[1]) + ylab(percentage[2]) +
                                  #labs(title = title) +
                                  theme_bw()+
                                  theme(plot.title = element_text(size=10), legend.position="bottom"))#+ stat_ellipse(type = "norm", linetype = 2))
  return(PCA_Sample_site_abundance)
}
# # *******************************************************************************************
# # **adonis/adonis2, Permutational Multivariate Analysis of Variance Using Distance Matrix**:
# library(vegan)
# #library(ztable)
# metadata_adonis <- as(sample_data(metaG_tpm_profile_rowsum_not_0_phyloseq), "data.frame")
# set.seed(2)
# adonis_Status<-vegan::adonis(phyloseq::distance(metaG_tpm_profile_rowsum_not_0_phyloseq, method= distance_lab) ~ Condition+time, data = metadata_adonis)
# #ztable(data.frame(adonis_Status$aov.tab),digits = 3)
# adonis_Status_df <- data.frame(adonis_Status$aov.tab)
# adonis_Status_df$Coefficients <- rownames(adonis_Status_df)
# write.table(adonis_Status_df, paste0(dir_diff_analysis_plot, db , '.FE.',filename3,'_',compartment,"_no0in_metaG.PCoA-",distance_lab,"_PERMANOVA.csv"), row.names = F, col.names = T, sep = "\t", quote = F)
# # *******************************************************************************************

#####
# Gene abundance profile
profiles_tpm <- input_file

print('#################################### Paths ####################################')
# Generate output directories ####
dir.create(file.path(integration_dir), showWarnings = FALSE)
visual_dir=paste0(integration_dir,'/plots/')
dir.create(file.path(visual_dir), showWarnings = FALSE)
phyloseq_dir=paste0(integration_dir,'/phyloseq_obj/')
dir.create(file.path(phyloseq_dir), showWarnings = FALSE)

if (omics == 'metaG_metaT'){
  print('#################################### DB profiles ####################################')
  tpm_profile_df = as.data.frame(fread(profiles_tpm, header=T), stringsAsFactors = F)
  tpm_profile_df$info[tpm_profile_df$name=='.'] <- tpm_profile_df$coord[tpm_profile_df$name=='.']
  bed_colnames <- c("chr","start","stop","name","score","strand","source","feature","frame","info")
  tpm_profile_bed_df <- tpm_profile_df[colnames(tpm_profile_df) %in% c('coord',bed_colnames)]
  tpm_profile_sub_df <- tpm_profile_df[!colnames(tpm_profile_df) %in% c(bed_colnames, 'gene_length')]
  tpm_profile_sub2 <- tpm_profile_sub_df[!duplicated(tpm_profile_sub_df[,'coord']),]
  rm(tpm_profile_sub_df)
  rm(tpm_profile_df)
  rownames(tpm_profile_sub2) <- NULL
  tpm_profile <- tpm_profile_sub2 %>%
    mutate(across(-1,  ~ as.numeric(replace(., . == '', 0))))
  tpm_profile[is.na(tpm_profile)] <- 0
  rm(tpm_profile_sub2)
  
  ## Bed file ####
  ## Subset df by colname
  gene_info <- c("coord","chr","start","stop", "name", "source","feature")
  ## Gene abundances
  tpm_profile_bed_df$coord <- paste0(tpm_profile_bed_df$chr, '_',tpm_profile_bed_df$start, '_', tpm_profile_bed_df$stop)
  #tpm_profile_bed_df$coord <- gsub('-', '_', tpm_profile_bed_df$coord)
  tpm_profile_bed_df$source2 <- unlist(lapply(strsplit(as.character(tpm_profile_bed_df$source),'\\,'), `[`, 1))
  tpm_profile_bed_df$name2 <- unlist(lapply(strsplit(as.character(tpm_profile_bed_df$name),'\\;'), `[`, 1))
  tpm_profile_bed_df$genome <- unlist(lapply(strsplit(as.character(tpm_profile_bed_df$chr),'\\|'), `[`, 1))
  tpm_profile_bed_df$contig <- unlist(lapply(strsplit(as.character(tpm_profile_bed_df$chr),'\\|'), `[`, 2))
  
  tpm_profile_bed_df$ID_gene <- tpm_profile_bed_df$info
  tpm_profile_bed_df$name3 <- unlist(lapply(stringr::str_split(pattern = '\\|', string = tpm_profile_bed_df$ID_gene), `[`, 3))
  tpm_profile_bed_df$name3[is.na(tpm_profile_bed_df$name3)] <- unlist(lapply(stringr::str_split(pattern = '\\|', string = tpm_profile_bed_df$ID_gene[is.na(tpm_profile_bed_df$name3)]), `[`, 2))
  tpm_profile_bed_df$name <- tpm_profile_bed_df$name3
  
  # Annotation ####
  profiles_annot_df <- as.data.frame(fread(annot_file,header=T), stringsAsFactors = F, row.names = T)
  
  ## Replace , by ; in annotation IDs
  profiles_annot_df[,2:ncol(profiles_annot_df)] <- data.frame(lapply(profiles_annot_df[,2:ncol(profiles_annot_df)], function(x) {
    gsub("\\,", "\\;", x)
  }))
  
  profiles_annot_df$ID_gene <- profiles_annot_df$ID
  profiles_annot_df$ID <- NULL
  profiles_annot_df <- subset(profiles_annot_df, select=unique(c("ID_gene", c(annot_names))))
  
  #### Bed file + Annotations ####
  ## Subset df by colname
  gene_info <- c("coord","chr","start","stop", "name", "source","feature", "ID_gene")
  ## Gene abundances
  tpm_profile_bed_sub_df <- tpm_profile_bed_df[colnames(tpm_profile_bed_df) %in% gene_info]
  profiles_annot_tpm_df <- merge(profiles_annot_df, tpm_profile_bed_sub_df, by='ID_gene', all.y=T)
  rownames(profiles_annot_tpm_df) <- profiles_annot_tpm_df$coord
  profiles_annot_tpm_df$coord <- NULL
  
  names_list <- names(profiles_annot_df)
  profiles_annot_tpm_df_sub <- subset(profiles_annot_tpm_df, select= c('name',names_list))#"KEGG_ko","kofam_KO"
  mappings_ncol <- ncol(profiles_annot_tpm_df_sub)
  # Delete data
  rm(profiles_annot_df, tpm_profile_bed_sub_df, tpm_profile_bed_df)
  
  # Split df into metaG and metaT ####
  rownames(tpm_profile) <- tpm_profile$coord
  tpm_profile$coord <- NULL
  
  # Subset metaG and metaT - keep only common samples
  ## Subset profile file by metaG and metaT
  metaG_tpm_profile = tpm_profile[grepl("metaG", colnames(tpm_profile))]
  colnames(metaG_tpm_profile) <- gsub(x = colnames(metaG_tpm_profile), pattern = "metaG\\.", replacement = "")
  metaT_tpm_profile = tpm_profile[grepl("metaT", colnames(tpm_profile))]
  colnames(metaT_tpm_profile) <- gsub(x = colnames(metaT_tpm_profile), pattern = "metaT\\.", replacement = "")
  ## Keep only samples in common meta G and metaT
  metaG_samples <- colnames(metaG_tpm_profile)
  metaT_samples <- colnames(metaT_tpm_profile)
  samples_intersect <- intersect(metaG_samples, metaT_samples)
  ## Subset metaG and metaT - keep only common samples
  ### Replace column names witht the same extention
  metaG_tpm_profile_sub <- metaG_tpm_profile[, colnames(metaG_tpm_profile) %in% samples_intersect]
  metaT_tpm_profile_sub <- metaT_tpm_profile[, colnames(metaT_tpm_profile) %in% samples_intersect]
  # Delete data
  rm(metaG_tpm_profile, metaT_tpm_profile, tpm_profile)
  
  ## Rows that sum 0 in metaG and metaT ####
  # metaG
  metaG_tpm_profile_rowsum_0 = c(which(rowSums(metaG_tpm_profile_sub)==0))
  metaG_tpm_profile_rowsum_0_list <- rownames(metaG_tpm_profile_sub[metaG_tpm_profile_rowsum_0,])
  # metaT
  metaT_tpm_profile_rowsum_0 = c(which(rowSums(metaT_tpm_profile_sub)==0))
  metaT_tpm_profile_rowsum_0_list <- rownames(metaT_tpm_profile_sub[metaT_tpm_profile_rowsum_0,])
  
  if (length(intersect(metaG_tpm_profile_rowsum_0_list,metaT_tpm_profile_rowsum_0_list)) != 0) {
    tpm_profile_rowsum_0 <- c(intersect(metaG_tpm_profile_rowsum_0_list,metaT_tpm_profile_rowsum_0_list))
    # metaG
    metaG_tpm_profile_rowsum_not_0 <- metaG_tpm_profile_sub[! rownames(metaG_tpm_profile_sub) %in% tpm_profile_rowsum_0,]
    # metaT
    metaT_tpm_profile_rowsum_not_0 <- metaT_tpm_profile_sub[! rownames(metaT_tpm_profile_sub) %in% tpm_profile_rowsum_0,]
  }else{
    metaG_tpm_profile_rowsum_not_0 <- metaG_tpm_profile_sub
    metaT_tpm_profile_rowsum_not_0 <- metaT_tpm_profile_sub
  }
  # Save GA and GT data ####
  # metaG
  metaG_tpm_profile_rowsum_not_0$info <- rownames(metaG_tpm_profile_rowsum_not_0)
  metaG_tpm_profile_rowsum_not_0 <- as.data.frame(metaG_tpm_profile_rowsum_not_0 %>%
                                                    dplyr::select(info, everything()))
  write.csv(metaG_tpm_profile_rowsum_not_0,paste0(integration_dir, '/GA.csv'), row.names = F, quote = F)
  metaG_tpm_profile_rowsum_not_0$info <- NULL
  
  # metaT
  metaT_tpm_profile_rowsum_not_0$info <- rownames(metaT_tpm_profile_rowsum_not_0)
  metaT_tpm_profile_rowsum_not_0 <- as.data.frame(metaT_tpm_profile_rowsum_not_0 %>%
                                                    dplyr::select(info, everything()))
  write.csv(metaT_tpm_profile_rowsum_not_0,paste0(integration_dir, '/GT.csv'), row.names = F, quote = F)
  metaT_tpm_profile_rowsum_not_0$info <- NULL
  
  ################################## GENE EXPRESSION  ################################## 
  print('#################################### GENE EXPRESSION  ####################################')
  profiles_annot_sub <- profiles_annot_tpm_df_sub[rownames(profiles_annot_tpm_df_sub) %in% rownames(metaG_tpm_profile_rowsum_not_0),]
  # Abundances
  all(rownames(metaG_tpm_profile_rowsum_not_0) %in% rownames(profiles_annot_sub))
  ## [1] TRUE
  all(rownames(metaG_tpm_profile_rowsum_not_0) == rownames(profiles_annot_sub))
  ## [1] FALSE
  profiles_annot_sub <- profiles_annot_sub[rownames(metaG_tpm_profile_rowsum_not_0),]
  all(rownames(metaG_tpm_profile_rowsum_not_0) == rownames(profiles_annot_sub))
  ## [1] TRUE
  
  write.csv(profiles_annot_sub,paste0(integration_dir, '/Annotations.csv'), row.names = F, quote = F)
  
  # Delete data
  rm(metaG_tpm_profile_sub, metaT_tpm_profile_sub, profiles_annot_tpm_df, profiles_annot_tpm_df_sub)
  
  #### phyloseq object ####
  # Metadata
  if (metadata_file=='None'){
    metadata_df <- data.frame(sample=names(metaG_tpm_profile_rowsum_not_0), condition='condition', sample_alias=names(metaG_tpm_profile_rowsum_not_0))
  } else{
    metadata_df <- as.data.frame(fread(metadata_file,  header = T), stringsAsFactors = F)
    names(metadata_df)[[1]] <- "sample"
    names(metadata_df)[[2]] <- "condition"
    names(metadata_df)[[3]] <- "sample_alias"
  }
  samp <- sample_data(metadata_df)
  rownames(samp) <- samp$sample
  
  ## metaG
  metaG_tpm_profile_rowsum_not_0_phyloseq <- phyloseq(otu_table(as.matrix(metaG_tpm_profile_rowsum_not_0), taxa_are_rows = T), 
                                                      tax_table(as.matrix(profiles_annot_sub)), samp)
  saveRDS(metaG_tpm_profile_rowsum_not_0_phyloseq, file = paste0(phyloseq_dir, '/GA.rds'))
  
  ## metaT
  metaT_tpm_profile_rowsum_not_0_phyloseq <- phyloseq(otu_table(as.matrix(metaT_tpm_profile_rowsum_not_0), taxa_are_rows = T), 
                                                      tax_table(as.matrix(profiles_annot_sub)), samp)
  saveRDS(metaT_tpm_profile_rowsum_not_0_phyloseq, file = paste0(phyloseq_dir, '/GT.rds'))
  
  # # Normalization
  # When metaG is 0, just get metaT value (replace metaG 0 values by 1) ####
  metaG_tpm_profile_rowsum_not_0_1 <- metaG_tpm_profile_rowsum_not_0
  metaG_tpm_profile_rowsum_not_0_1[metaG_tpm_profile_rowsum_not_0_1==0] <- 1
  
  metaG_tpm_profile_pseud <- metaG_tpm_profile_rowsum_not_0_1
  metaT_tpm_profile_pseud <- metaT_tpm_profile_rowsum_not_0
  
  # # metaG and metaT matrix same function order
  gene_ids = rownames(metaG_tpm_profile_pseud)
  relat_num_transc_cell_tpm = metaT_tpm_profile_pseud[match(rownames(metaG_tpm_profile_pseud), rownames(metaT_tpm_profile_pseud)),]
  relat_num_functions_cell_tpm = metaG_tpm_profile_pseud
  
  # Gene Expression
  gene_expression_tpm=(relat_num_transc_cell_tpm)/(relat_num_functions_cell_tpm)
  
  gene_expression_tpm_noinf <- as.data.frame(gene_expression_tpm %>%
                                               tibble::rownames_to_column('rownames') %>%
                                               dplyr::mutate_all(function(x) ifelse(is.infinite(x), NA, x)) %>%
                                               tibble::column_to_rownames('rownames'))
  
  # Save GE data ####
  gene_expression_tpm_noinf$info <- rownames(gene_expression_tpm_noinf)
  gene_expression_tpm_noinf <- as.data.frame(gene_expression_tpm_noinf %>%
                                               dplyr::select(info, everything()))
  write.csv(gene_expression_tpm_noinf,paste0(integration_dir, '/GE.csv'), row.names = F, quote = F)
  gene_expression_tpm_noinf$info <- NULL
  
  #### phyloseq object
  gene_expression_tpm_nonbinary_phyloseq <- phyloseq(otu_table(as.matrix(gene_expression_tpm_noinf), taxa_are_rows = T), 
                                                     tax_table(as.matrix(profiles_annot_sub)), samp)
  head(taxa_names(gene_expression_tpm_nonbinary_phyloseq))
  saveRDS(gene_expression_tpm_nonbinary_phyloseq, file = paste0(phyloseq_dir, '/GE.rds'))
  
  # Delete data
  rm(gene_expression_tpm_nonbinary_phyloseq, relat_num_transc_cell_tpm, relat_num_functions_cell_tpm, 
     metaG_tpm_profile_pseud, metaT_tpm_profile_pseud, metaG_tpm_profile_rowsum_not_0_phyloseq, metaT_tpm_profile_rowsum_not_0_phyloseq,
     gene_expression_tpm, metadata_df, metaG_tpm_profile_rowsum_not_0_1)
  # Free memory
  gc()
  print('#################################### PCA  ####################################')#####
  manual_plot_colors =c('#9D0208', '#264653','#e9c46a','#D8DCDE','#B6D0E0',
                        '#FFC87E','#F4A261','#E34F33','#E9C46A',
                        '#A786C9','#D4C0E2','#975773','#6699FF','#000066',
                        '#7AAFCA','#006699','#A9D181','#2F8475','#264445')  
  #distance_lab='bray'
  # PCA - metaG ####
  ## Prepare data - replace NA values by 0
  metaG_tpm_profile_rowsum_noinf <- as.data.frame(metaG_tpm_profile_rowsum_not_0 %>%
                                                    tibble::rownames_to_column('rownames') %>%
                                                    dplyr::mutate_all(function(x) ifelse(is.infinite(x), 0, x)) %>%
                                                    tibble::column_to_rownames('rownames'))
  metaG_tpm_profile_rowsum_noinf <- as.data.frame(metaG_tpm_profile_rowsum_noinf %>%
                                                    tibble::rownames_to_column('rownames') %>%
                                                    dplyr::mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
                                                    tibble::column_to_rownames('rownames'))
  
  metaG_tpm_profile_rowsum_noinf_plot <- phyloseq(otu_table(as.matrix(metaG_tpm_profile_rowsum_noinf), taxa_are_rows = T), tax_table(as.matrix(profiles_annot_sub)), samp)
  
  # title_name <- 'PCoA - gene abundance'
  # out_name <- paste0(visual_dir, 'GA.PCoA-',distance_lab,'.pdf')
  # plot_PCoA_out <- plot_PCoA(distance_lab, metaG_tpm_profile_rowsum_noinf_plot, "condition", "sample_alias")
  # 
  # pdf(out_name,width=8,height=8,paper="special" )
  # print(plot_PCoA_out  + 
  #         scale_color_manual(values=manual_plot_colors, name='Conditions') +
  #         #coord_fixed() +
  #         theme(legend.position="bottom")+
  #         ggtitle(title_name, 'Bray Curtis'))
  # print(plot_PCoA_out  + 
  #         facet_wrap(.~condition)+
  #         scale_color_manual(values=manual_plot_colors, name='Conditions') +
  #         #coord_fixed() +
  #         theme(legend.position="bottom")+
  #         ggtitle(title_name, 'Bray Curtis'))
  # dev.off()
  
  title_name <- 'PCA - gene abundance'
  out_name <- paste0(visual_dir, 'GA.PCA.pdf')
  plot_PCA_out <- plot_PCA(metaG_tpm_profile_rowsum_noinf_plot, "condition", "sample_alias")
  pdf(out_name,width=8,height=8,paper="special" )
  print(plot_PCA_out  + scale_color_manual(values=manual_plot_colors, name='Condition') +
          #coord_fixed() +
          theme(legend.position="bottom")+
          ggtitle(title_name))
  print(plot_PCA_out  + 
          facet_wrap(.~condition)+
          scale_color_manual(values=manual_plot_colors, name='Condition') +
          #coord_fixed() +
          theme(legend.position="bottom")+
          ggtitle(title_name))
  dev.off()
  
  # Delete data
  rm(metaG_tpm_profile_rowsum_noinf_plot, metaG_tpm_profile_rowsum_noinf)
  
  # PCA - metaT ####
  ## Prepare data - replace NA values by 0
  metaT_tpm_profile_rowsum_noinf <- as.data.frame(metaT_tpm_profile_rowsum_not_0 %>%
                                                    tibble::rownames_to_column('rownames') %>%
                                                    dplyr::mutate_all(function(x) ifelse(is.infinite(x), 0, x)) %>%
                                                    tibble::column_to_rownames('rownames'))
  metaT_tpm_profile_rowsum_noinf <- as.data.frame(metaT_tpm_profile_rowsum_noinf %>%
                                                    tibble::rownames_to_column('rownames') %>%
                                                    dplyr::mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
                                                    tibble::column_to_rownames('rownames'))
  
  metaT_tpm_profile_rowsum_noinf_plot <- phyloseq(otu_table(as.matrix(metaT_tpm_profile_rowsum_noinf), taxa_are_rows = T), tax_table(as.matrix(profiles_annot_sub)), samp)
  
  # title_name <- 'PCoA - gene transcript'
  # out_name <- paste0(visual_dir, 'GT.PCoA-',distance_lab,'.pdf')
  # plot_PCoA_out <- plot_PCoA(distance_lab, metaT_tpm_profile_rowsum_noinf_plot, "condition", "sample_alias")
  # 
  # pdf(out_name,width=8,height=8,paper="special" )
  # print(plot_PCoA_out  + 
  #         scale_color_manual(values=manual_plot_colors, name='Conditions') +
  #         #coord_fixed() +
  #         theme(legend.position="bottom")+
  #         ggtitle(title_name, 'Bray Curtis'))
  # print(plot_PCoA_out  + 
  #         facet_wrap(.~condition)+
  #         scale_color_manual(values=manual_plot_colors, name='Conditions') +
  #         #coord_fixed() +
  #         theme(legend.position="bottom")+
  #         ggtitle(title_name, 'Bray Curtis'))
  # dev.off()
  
  title_name <- 'PCA - gene transcript'
  out_name <- paste0(visual_dir, 'GT.PCA.pdf')
  plot_PCA_out <- plot_PCA(metaT_tpm_profile_rowsum_noinf_plot, "condition", "sample_alias")
  pdf(out_name,width=8,height=8,paper="special" )
  print(plot_PCA_out  + scale_color_manual(values=manual_plot_colors, name='Condition') +
          #coord_fixed() +
          theme(legend.position="bottom")+
          ggtitle(title_name))
  print(plot_PCA_out  + 
          facet_wrap(.~condition)+
          scale_color_manual(values=manual_plot_colors, name='Condition') +
          #coord_fixed() +
          theme(legend.position="bottom")+
          ggtitle(title_name))
  dev.off()
  
  # Delete data
  rm(metaT_tpm_profile_rowsum_noinf_plot, metaT_tpm_profile_rowsum_noinf)
  
  # PCA - GE ####
  ## Prepare data - replace NA values by 0
  gene_expression_tpm_profile_rowsum_noinf <- as.data.frame(gene_expression_tpm_noinf %>%
                                                              tibble::rownames_to_column('rownames') %>%
                                                              dplyr::mutate_all(function(x) ifelse(is.infinite(x), 0, x)) %>%
                                                              tibble::column_to_rownames('rownames'))
  gene_expression_tpm_profile_rowsum_noinf <- as.data.frame(gene_expression_tpm_profile_rowsum_noinf %>%
                                                              tibble::rownames_to_column('rownames') %>%
                                                              dplyr::mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
                                                              tibble::column_to_rownames('rownames'))
  
  gene_expression_tpm_profile_rowsum_noinf_plot <- phyloseq(otu_table(as.matrix(gene_expression_tpm_profile_rowsum_noinf), taxa_are_rows = T), tax_table(as.matrix(profiles_annot_sub)), samp)
  
  # title_name <- 'PCoA - gene expression'
  # out_name <- paste0(visual_dir, 'GE.PCoA-',distance_lab,'.pdf')
  # plot_PCoA_out <- plot_PCoA(distance_lab, gene_expression_tpm_profile_rowsum_noinf_plot, "condition", "sample_alias")
  # 
  # pdf(out_name,width=8,height=8,paper="special" )
  # print(plot_PCoA_out  + 
  #         scale_color_manual(values=manual_plot_colors, name='Conditions') +
  #         #coord_fixed() +
  #         theme(legend.position="bottom")+
  #         ggtitle(title_name, 'Bray Curtis'))
  # print(plot_PCoA_out  + 
  #         facet_wrap(.~condition)+
  #         scale_color_manual(values=manual_plot_colors, name='Conditions') +
  #         #coord_fixed() +
  #         theme(legend.position="bottom")+
  #         ggtitle(title_name, 'Bray Curtis'))
  # dev.off()
  
  title_name <- 'PCA - gene expression'
  out_name <- paste0(visual_dir, 'GE.PCA.pdf')
  plot_PCA_out <- plot_PCA(gene_expression_tpm_profile_rowsum_noinf_plot, "condition", "sample_alias")
  pdf(out_name,width=8,height=8,paper="special" )
  print(plot_PCA_out  + scale_color_manual(values=manual_plot_colors, name='Condition') +
          #coord_fixed() +
          theme(legend.position="bottom")+
          ggtitle(title_name))
  print(plot_PCA_out  + 
          facet_wrap(.~condition)+
          scale_color_manual(values=manual_plot_colors, name='Condition') +
          #coord_fixed() +
          theme(legend.position="bottom")+
          ggtitle(title_name))
  dev.off()
  
  # Delete data
  rm(gene_expression_tpm_profile_rowsum_noinf_plot, gene_expression_tpm_profile_rowsum_noinf)
  
  print('Done!')
  
} else{
  print('#################################### DB profiles ####################################')
  tpm_profile_df = as.data.frame(fread(profiles_tpm, header=T), stringsAsFactors = F)
  tpm_profile_df$info[tpm_profile_df$name=='.'] <- tpm_profile_df$coord[tpm_profile_df$name=='.']
  bed_colnames <- c("chr","start","stop","name","score","strand","source","feature","frame","info")
  tpm_profile_bed_df <- tpm_profile_df[colnames(tpm_profile_df) %in% c('coord',bed_colnames)]
  tpm_profile_sub_df <- tpm_profile_df[!colnames(tpm_profile_df) %in% c(bed_colnames, 'gene_length')]
  tpm_profile_sub2 <- tpm_profile_sub_df[!duplicated(tpm_profile_sub_df[,'coord']),]
  rm(tpm_profile_sub_df)
  rm(tpm_profile_df)
  rownames(tpm_profile_sub2) <- NULL
  tpm_profile <- tpm_profile_sub2 %>%
    mutate(across(-1,  ~ as.numeric(replace(., . == '', 0))))
  tpm_profile[is.na(tpm_profile)] <- 0
  rm(tpm_profile_sub2)
  
  ## Bed file ####
  ## Subset df by colname
  gene_info <- c("coord","chr","start","stop", "name", "source","feature")
  ## Gene abundances
  tpm_profile_bed_df$coord <- paste0(tpm_profile_bed_df$chr, '_',tpm_profile_bed_df$start, '_', tpm_profile_bed_df$stop)
  #tpm_profile_bed_df$coord <- gsub('-', '_', tpm_profile_bed_df$coord)
  tpm_profile_bed_df$source2 <- unlist(lapply(strsplit(as.character(tpm_profile_bed_df$source),'\\,'), `[`, 1))
  tpm_profile_bed_df$name2 <- unlist(lapply(strsplit(as.character(tpm_profile_bed_df$name),'\\;'), `[`, 1))
  tpm_profile_bed_df$genome <- unlist(lapply(strsplit(as.character(tpm_profile_bed_df$chr),'\\|'), `[`, 1))
  tpm_profile_bed_df$contig <- unlist(lapply(strsplit(as.character(tpm_profile_bed_df$chr),'\\|'), `[`, 2))
  
  tpm_profile_bed_df$ID_gene <- tpm_profile_bed_df$info
  tpm_profile_bed_df$name3 <- unlist(lapply(stringr::str_split(pattern = '\\|', string = tpm_profile_bed_df$ID_gene), `[`, 3))
  tpm_profile_bed_df$name3[is.na(tpm_profile_bed_df$name3)] <- unlist(lapply(stringr::str_split(pattern = '\\|', string = tpm_profile_bed_df$ID_gene[is.na(tpm_profile_bed_df$name3)]), `[`, 2))
  tpm_profile_bed_df$name <- tpm_profile_bed_df$name3
  
  # Annotation ####
  profiles_annot_df <- as.data.frame(fread(annot_file,header=T), stringsAsFactors = F, row.names = T)
  
  ## Replace , by ; in annotation IDs
  profiles_annot_df[,2:ncol(profiles_annot_df)] <- data.frame(lapply(profiles_annot_df[,2:ncol(profiles_annot_df)], function(x) {
    gsub("\\,", "\\;", x)
  }))
  
  profiles_annot_df$ID_gene <- profiles_annot_df$ID
  profiles_annot_df$ID <- NULL
  profiles_annot_df <- subset(profiles_annot_df, select=unique(c("ID_gene", c(annot_names))))
  
  #### Bed file + Annotations ####
  ## Subset df by colname
  gene_info <- c("coord","chr","start","stop", "name", "source","feature", "ID_gene")
  ## Gene abundances
  tpm_profile_bed_sub_df <- tpm_profile_bed_df[colnames(tpm_profile_bed_df) %in% gene_info]
  profiles_annot_tpm_df <- merge(profiles_annot_df, tpm_profile_bed_sub_df, by='ID_gene', all.y=T)
  rownames(profiles_annot_tpm_df) <- profiles_annot_tpm_df$coord
  profiles_annot_tpm_df$coord <- NULL
  
  names_list <- names(profiles_annot_df)
  profiles_annot_tpm_df_sub <- subset(profiles_annot_tpm_df, select= c('name',names_list))#"KEGG_ko","kofam_KO"
  mappings_ncol <- ncol(profiles_annot_tpm_df_sub)
  # Delete data
  rm(profiles_annot_df, tpm_profile_bed_sub_df, tpm_profile_bed_df)
  
  # Subset data - exclude gene abundances = 0 for all samples####
  rownames(tpm_profile) <- tpm_profile$coord
  tpm_profile$coord <- NULL
  omics_tpm_profile = tpm_profile[grepl(omics, colnames(tpm_profile))]
  colnames(omics_tpm_profile) <- gsub(x = colnames(omics_tpm_profile), pattern = paste0(omics,"\\."), replacement = "")
  # Delete data
  rm(tpm_profile)
  
  ## Delete rows that sum 0 ####
  omics_tpm_profile_rowsum_0 = c(which(rowSums(omics_tpm_profile)==0))
  omics_tpm_profile_rowsum_0_list <- rownames(omics_tpm_profile[omics_tpm_profile_rowsum_0,])
  
  if (length(omics_tpm_profile_rowsum_0_list) != 0) {
    tpm_profile_rowsum_0 <- c(omics_tpm_profile_rowsum_0_list)
    omics_tpm_profile_rowsum_not_0 <- omics_tpm_profile[! rownames(omics_tpm_profile) %in% tpm_profile_rowsum_0,]
  }else{
    omics_tpm_profile_rowsum_not_0 <- omics_tpm_profile
  }
  # Save gene abundance data ####
  omics_tpm_profile_rowsum_not_0$info <- rownames(omics_tpm_profile_rowsum_not_0)
  omics_tpm_profile_rowsum_not_0 <- as.data.frame(omics_tpm_profile_rowsum_not_0 %>%
                                                    dplyr::select(info, everything()))
  if (omics == 'metaG'){
    omics_label <-'GA'
    omics_label_plot <- 'gene abundance'
  }else if (omics == 'metaT'){
    omics_label <-'GT'
    omics_label_plot <- 'gene transcript'
  }
  write.csv(omics_tpm_profile_rowsum_not_0,paste0(integration_dir, '/', omics_label, '.csv'), row.names = F, quote = F)
  omics_tpm_profile_rowsum_not_0$info <- NULL
  
  #### Annotations ####
  profiles_annot_sub <- profiles_annot_tpm_df_sub[rownames(profiles_annot_tpm_df_sub) %in% rownames(omics_tpm_profile_rowsum_not_0),]
  # Abundances
  all(rownames(omics_tpm_profile_rowsum_not_0) %in% rownames(profiles_annot_sub))
  ## [1] TRUE
  all(rownames(omics_tpm_profile_rowsum_not_0) == rownames(profiles_annot_sub))
  ## [1] TRUE
  profiles_annot_sub <- profiles_annot_sub[rownames(omics_tpm_profile_rowsum_not_0),]
  all(rownames(omics_tpm_profile_rowsum_not_0) == rownames(profiles_annot_sub))
  ## [1] TRUE
  
  write.csv(profiles_annot_sub,paste0(integration_dir, '/Annotations.csv'), row.names = F, quote = F)
  
  # Metadata
  if (metadata_file=='None'){
    metadata_df <- data.frame(sample=names(omics_tpm_profile_rowsum_not_0), condition='condition', sample_alias=names(omics_tpm_profile_rowsum_not_0))
  } else{
    metadata_df <- as.data.frame(fread(metadata_file,  header = T), stringsAsFactors = F)
    names(metadata_df)[[1]] <- "sample"
    names(metadata_df)[[2]] <- "condition"
    names(metadata_df)[[3]] <- "sample_alias"
  }
  samp <- sample_data(metadata_df)
  rownames(samp) <- samp$sample
  
  #### phyloseq object ####
  omics_tpm_profile_rowsum_not_0_phyloseq <- phyloseq(otu_table(as.matrix(omics_tpm_profile_rowsum_not_0), taxa_are_rows = T), 
                                                      tax_table(as.matrix(profiles_annot_sub)), samp)
  saveRDS(omics_tpm_profile_rowsum_not_0_phyloseq, file = paste0(phyloseq_dir, '/', omics_label, '.rds'))
  
  print('#################################### PCA  ####################################')#####
  manual_plot_colors =c('#9D0208', '#264653','#e9c46a','#D8DCDE','#B6D0E0',
                        '#FFC87E','#F4A261','#E34F33','#E9C46A',
                        '#A786C9','#D4C0E2','#975773','#6699FF','#000066',
                        '#7AAFCA','#006699','#A9D181','#2F8475','#264445')
  #distance_lab='bray'
  # PCA ####
  ## Prepare data - replace NA values by 0
  omics_tpm_profile_rowsum_noinf <- as.data.frame(omics_tpm_profile_rowsum_not_0 %>%
                                                    tibble::rownames_to_column('rownames') %>%
                                                    dplyr::mutate_all(function(x) ifelse(is.infinite(x), 0, x)) %>%
                                                    tibble::column_to_rownames('rownames'))
  omics_tpm_profile_rowsum_noinf <- as.data.frame(omics_tpm_profile_rowsum_noinf %>%
                                                    tibble::rownames_to_column('rownames') %>%
                                                    dplyr::mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
                                                    tibble::column_to_rownames('rownames'))
  
  omics_tpm_profile_rowsum_noinf_plot <- phyloseq(otu_table(as.matrix(omics_tpm_profile_rowsum_noinf), taxa_are_rows = T), tax_table(as.matrix(profiles_annot_sub)), samp)
  
  # title_name <- paste0('PCoA - ', omics_label_plot)
  # out_name <- paste0(visual_dir, omics_label, '.PCoA-',distance_lab,'.pdf')
  # plot_PCoA_out <- plot_PCoA(distance_lab, omics_tpm_profile_rowsum_noinf_plot, "condition", "sample_alias")
  # 
  # pdf(out_name,width=8,height=8,paper="special" )
  # print(plot_PCoA_out  + 
  #         scale_color_manual(values=manual_plot_colors, name='Conditions') +
  #         #coord_fixed() +
  #         theme(legend.position="bottom")+
  #         ggtitle(title_name, 'Bray Curtis'))
  # print(plot_PCoA_out  + 
  #         facet_wrap(.~condition)+
  #         scale_color_manual(values=manual_plot_colors, name='Conditions') +
  #         #coord_fixed() +
  #         theme(legend.position="bottom")+
  #         ggtitle(title_name, 'Bray Curtis'))
  # dev.off()
  
  title_name <- paste0('PCA -  ', omics_label_plot)
  out_name <- paste0(visual_dir, omics_label, '.PCA.pdf')
  plot_PCA_out <- plot_PCA(omics_tpm_profile_rowsum_noinf_plot, "condition", "sample_alias")
  pdf(out_name,width=8,height=8,paper="special" )
  print(plot_PCA_out  + scale_color_manual(values=manual_plot_colors, name='Condition') +
          #coord_fixed() +
          theme(legend.position="bottom")+
          ggtitle(title_name))
  print(plot_PCA_out  + 
          facet_wrap(.~condition)+
          scale_color_manual(values=manual_plot_colors, name='Condition') +
          #coord_fixed() +
          theme(legend.position="bottom")+
          ggtitle(title_name))
  dev.off()
  
  print('Done!')
}
