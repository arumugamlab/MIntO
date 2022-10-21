#!/usr/bin/env Rscript

# ---
# title: "Taxonomic profile"
# author: "Carmen Saenz"
# date: "09/10/2021"
# Generate Visualization plots from assembly-free taxonomy profile (from MetaPhlan or mOTUs2):
#   - PCoA
#   - Relative abundances of 15 most abundant genera in each sample
# ---

# Load libraries
#if("data.table" %in% rownames(installed.packages()) == FALSE) {install.packages("data.table")}
library(data.table)
#if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse")}
library(tidyverse)
#if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2")}
library(ggplot2)
#if("ggrepel" %in% rownames(installed.packages()) == FALSE) {install.packages("ggrepel")}
library(ggrepel)
#if("phyloseq" %in% rownames(installed.packages()) == FALSE) {install.packages("phyloseq")}
library(phyloseq)
#if("reshape2" %in% rownames(installed.packages()) == FALSE) {install.packages("reshape2")}
library(reshape2)
library(tidyr)
library(vegan)

set.seed(1234)

# Arguments
# args = commandArgs(trailingOnly=TRUE)
# profile_dir = args[1]
# profile_param = args[2]
# out_dir = args[3]
# metadata_file = args[4]

# # Arguments - For 16S rRNA
path_dir <- "/emc/cbmr/data/microbiome/processed/SHIME/Cdiff_sterile1/16SrRNA/"  
profile_dir <- 
profile_param <- 'raw' #'RA'
out_dir <- "/emc/cbmr/data/microbiome/processed/SHIME/Cdiff_sterile1/diff_analysis/genes-PBGCs-SBGC-smORFs-Cdiff_Cscind.p95.TPM.locus-tag.csv_nofilter/DC/"
metadata_file <- "/emc/cbmr/data/microbiome/processed/SHIME/Cdiff_sterile1/16SrRNA/RNA_samples_SUBSET_short.csv"
out_phyloseq = paste0(out_dir, '/16SrRNA.ASVs.Silva_SUBSET.Prev-Relative_copy_number.rds') # '16SrRNA.ASVs.Silva_SUBSET.Prev-Relative_copy_number_RA.rds'
out_name <- '16SrRNA.ASVs.Silva_SUBSET.Prev-Relative_copy_number_RA'

# # Arguments - For GE: metaT (RPK) relative to 16S rRNA (normalized by copy number)
path_dir <- "/emc/cbmr/data/microbiome/processed/SHIME/Cdiff_sterile1/"
profile_param <- 'RA' 
out_dir <- "/emc/cbmr/data/microbiome/processed/SHIME/Cdiff_sterile1/data_integration/paper_genes_abundances.p95/TPM/"
metadata_file <- "/emc/cbmr/data/microbiome/processed/SHIME/Cdiff_sterile1/16SrRNA/RNA_samples_SUBSET_short.csv"
out_phyloseq = paste0(out_dir, 'phyloseq_obj/GT-RPK.16SrRNA.ASV.Silva.Relative_copy_number.rds')
out_name <- 'GT-RPK.16SrRNA.ASV.Silva.Relative_copy_number'

##########################  ** Function - PCoA **  ########################## 

plot_PCoA <- function(distance_lab, data_phyloseq, color, label){ #out_name,title_name,
  label2 <- noquote(label)
  #### PCoA
  ord<- ordinate(data_phyloseq, method = "PCoA", distance = distance_lab)
  PcoA_Sample_site_abundance <- plot_ordination(data_phyloseq, ord, color = color, shape=NULL, label=label) #type = type,
  PcoA_Sample_site_abundance <- PcoA_Sample_site_abundance +
    ggtitle(title_name, distance_lab)  +
    geom_point(size = 2.5) +
    theme_bw() +
    theme(legend.position = "top")
  PcoA_Sample_site_abundance$layers[[1]] <- NULL
  PcoA_Sample_site_abundance$layers[[2]] <- NULL
  PcoA_Sample_site_abundance$layers[[3]] <- NULL
  PcoA_Sample_site_abundance_2<- PcoA_Sample_site_abundance +
    geom_text_repel(aes(label = time), size = 3.0, segment.alpha = 0.5, max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
    #geom_point(aes(fill=participant_ID2), colour="black", size = 2.5, shape = 21, alpha=0.7)
    geom_point(size = 2.1, shape = 16, alpha=1)
  PcoA_Sample_site_abundance_2$layers[[1]] <- NULL
  return(PcoA_Sample_site_abundance_2)
}


# Phyloseq object
if (profile_param %like% 'raw') {
  profile_phyloseq <- readRDS(out_phyloseq)
  #### Transform counts to Relative Abundance
  profile_phyloseq_ra<-transform_sample_counts(profile_phyloseq, function(x){x/sum(x)})
} else{
  profile_phyloseq_ra <- readRDS(out_phyloseq)
} 


# PCoA - Bray-Curtis ##### ALL samples ####
distance_lab = 'bray'
title_name <- c(paste0("PCoA - 16S rRNA"), "Bray-Curtis")
out_name_plot <- paste0(out_dir, out_name, ".PCoA.Bray_Curtis.pdf")
title_name_pval <- paste0("Bray-Curtis")
metadata_df <- as(sample_data(profile_phyloseq_ra), "data.frame")
# 16S
if(length(unique(metadata_df$Condition))>1){
  #**adonis/adonis2, Permutational Multivariate Analysis of Variance Using Distance Matrix**: ####
  library(vegan)
  #adonis_list$bray
  otu_table_sort <- as.data.frame(unclass(t(otu_table(profile_phyloseq_ra))), stringsAsFactors = F)
  metadata_df <- data.frame(sample_data(profile_phyloseq_ra), stringsAsFactors = F)
  rownames(metadata_df) <- metadata_df$Sample
  metadata <-metadata_df[match(names(otu_table_sort), rownames(metadata_df)),]
  #adonis_bray_Status<-adonis(as.dist(vegdist(t(otu_table_sort), method="bray", na.rm = T)) ~ condition, data = metadata)
  #adonis_bray_Status$aov.tab
  adonis_bray_Status2<-adonis2(as.dist(vegdist(t(otu_table_sort), method="bray", na.rm = T)) ~ Condition, data = metadata)
  r2_value <- format(round(adonis_bray_Status2$R2[1],3), nsmall = 3)
  p_value <- adonis_bray_Status2$`Pr(>F)`[1]
  
  title_name_pval <- paste0("Bray-Curtis; R2=",r2_value,"; pval=",p_value)
}
distance_lab = 'jaccard'
title_name <- c(paste0("PCoA - gene expression"), "Bray-Curtis")
# metaT
if(length(unique(metadata_df$Condition))>1){
  #**adonis/adonis2, Permutational Multivariate Analysis of Variance Using Distance Matrix**: ####
  library(vegan)
  #adonis_list$bray
  otu_table_sort <- as.data.frame(unclass(t(otu_table(profile_phyloseq_ra))), stringsAsFactors = F)
  metadata_df <- data.frame(sample_data(profile_phyloseq_ra), stringsAsFactors = F)
  rownames(metadata_df) <- metadata_df$Sample
  metadata <-metadata_df[match(rownames(otu_table_sort), rownames(metadata_df)),]
  #adonis_bray_Status<-adonis(as.dist(vegdist(t(otu_table_sort), method="bray", na.rm = T)) ~ condition, data = metadata)
  #adonis_bray_Status$aov.tab
  adonis_bray_Status2<-adonis2(as.dist(vegdist(otu_table_sort, method="bray", na.rm = T)) ~ Condition, data = metadata)
  r2_value <- format(round(adonis_bray_Status2$R2[1],3), nsmall = 3)
  p_value <- adonis_bray_Status2$`Pr(>F)`[1]
  
  title_name_pval <- paste0("Bray-Curtis; R2=",r2_value,"; pval=",p_value)
}

#plot_PCoA_out <- plot_PCoA(distance_lab, profile_phyloseq_ra, title_name, out_name)

plot_PCoA_out <- plot_PCoA(distance_lab, profile_phyloseq_ra, "Condition", "time")

# manual_plot_colors =c("nonIBD" = "#ADB6BA", "CD" = "#264653", "UC" = "#9D0208")
# manual_plot_colors =c("nIBD1" = "#798083", "nIBD2" = "#CAD0D2",
#                       "CD1" = "#1C4151", "CD2" = "#4B6C7A",
#                       "UC1" = "#610508", "UC2" = "#A24B4E")

#manual_plot_colors =c("Cdiff" = "#D46C4E", "Cdiff_Csind" = "#264D59") #, "Cdiff_4str" = "#43978D"
manual_plot_colors =c("Cdiff" = "#92392F", "Cdiff_Csind" = "#11698E") 

# manual_plot_colors =c("nIBD1" = "#798083", "nIBD2" = "#CAD0D2",
#                       "CD1" = "#0C5B56", "CD2" = "#5CA5A0",
#                       "UC1" = "#DEAB34", "UC2" = "#DA5B41")


# Save ####
pdf(out_name_plot,width=3.2,height=4,paper="special" )
print(plot_PCoA_out  + 
        scale_color_manual(values=manual_plot_colors, name='Condition') +
        #coord_fixed() +
        theme(legend.position="bottom")+
        ggtitle(title_name, title_name_pval))
print(plot_PCoA_out  + 
        facet_wrap(.~Condition)+
        scale_color_manual(values=manual_plot_colors, name='Condition') +
        #coord_fixed() +
        theme(legend.position="bottom")+
        ggtitle(title_name, title_name_pval))
dev.off()

#**adonis/adonis2, Permutational Multivariate Analysis of Variance Using Distance Matrix**: ####
library(vegan)
adonis_list$bray
metadata <- as(sample_data(profile_phyloseq_ra), "data.frame")
adonis_bray_Status<-adonis2(as.dist(vegdist(t(otu_table_df), method="bray", na.rm = T)) ~ condition+week_n, data = metadata)
adonis_bray_Status$aov.tab
adonis_bray_Status2<-adonis2(as.dist(vegdist(t(otu_table_df), method="bray", na.rm = T)) ~ condition+week_n, data = metadata)
adonis_bray_Status2

#adonis_bray_Status$`Pr(>F)`
#ztable(as.data.frame(adonis_bray_Status$aov.tab),digits = 3)

# Barchart - Top 15 genera####
##### metadata:
sample_data_df <- data.frame(sample_data(profile_phyloseq_ra), stringsAsFactors = F)

### OTUs - db
taxa_table_df <- as.data.frame(unclass(tax_table(profile_phyloseq_ra)), stringsAsFactors = F)
taxa_table_df$taxa_ID <- rownames(taxa_table_df)
rownames(taxa_table_df) <- NULL

### Counts
otu_table_df <-  as.data.frame(unclass(t(otu_table(profile_phyloseq_ra))), stringsAsFactors = F)
otu_table_df$taxa_ID <- rownames(otu_table_df)
rownames(otu_table_df) <- NULL
#otu_table_df$taxa_ID <- gsub('value\\.', '', otu_table_df$taxa_ID )
#otu_table_df$taxa_ID <- gsub('\\.', ' ', otu_table_df$taxa_ID )
#otu_table_df$taxa_ID[otu_table_df$taxa_ID %like% 'Ralstonia'] <- gsub(' ', '/', otu_table_df$taxa_ID[otu_table_df$taxa_ID %like% 'Ralstonia'])
#otu_table_df$taxa_ID[otu_table_df$taxa_ID %like% 'Ralstonia'] <- sub('/', ' ', otu_table_df$taxa_ID[otu_table_df$taxa_ID %like% 'Ralstonia'])

otu_taxa_merge <- merge(otu_table_df, taxa_table_df, by = 'taxa_ID' , all.x = T)

# Filter by mOTUs present in the data
otu_taxa_filt_df = otu_taxa_merge[rowSums(otu_taxa_merge[, 2:ncol(otu_table_df)])>0,]
otu_taxa_filt_df = otu_taxa_filt_df %>% dplyr::select(1:ncol(otu_table_df))

# # Prevalence 10% 
# otu_taxa_filt_df = otu_taxa_merge[rowSums(otu_taxa_merge[, 2:ncol(otu_table_df)])>0.01,]
# otu_taxa_filt_df = otu_taxa_filt_df %>% select(1:ncol(otu_table_df)) #"kingdom", "phylum","class","order","family","genus","mOTU","short_name"

# # After filtering by prevalence 10%, samples don't sum to 1 anymore, make it sum to 1
# colSums(otu_taxa_filt_df[2:ncol(otu_table_df)])
# myNumCols <- which(unlist(lapply(otu_taxa_filt_df, is.numeric)))
# otu_taxa_filt_df[(nrow(otu_taxa_filt_df) + 1), myNumCols] <- 1- colSums(otu_taxa_filt_df[, myNumCols], na.rm=TRUE)
# otu_taxa_filt_df[is.na(otu_taxa_filt_df)] <- '-1'

otu_taxa_melt <- melt(otu_taxa_filt_df, id.vars=c("taxa_ID"))
otu_taxa_melt <- merge(taxa_table_df, otu_taxa_melt, by="taxa_ID")
#colnames(otu_taxa_melt) <- c('taxa_ID',"kingdom", "phylum","class","order","family","genus","mOTU","short_name","sample", "RA"   )

otu_taxa_metadata =  merge(sample_data_df,otu_taxa_melt,  by.x = 'Sample', by.y = 'variable')

otu_taxa_metadata$value <- as.numeric(otu_taxa_metadata$value)

#### STATS ######
# Cdiff and Cscind -summary 
summary(otu_taxa_metadata$value[which(otu_taxa_metadata$Species == 'difficile' & otu_taxa_metadata$group == 'Cdiff')])
summary(otu_taxa_metadata$value[which(otu_taxa_metadata$Species == 'difficile' & otu_taxa_metadata$group == 'Cdiff_Csind')])
summary(otu_taxa_metadata$value[which(otu_taxa_metadata$Species == 'scindens' & otu_taxa_metadata$group == 'Cdiff_Csind')])
# Cdiff and Cscind
otu_taxa_metadata$value[which(otu_taxa_metadata$Species == 'difficile' & otu_taxa_metadata$group == 'Cdiff' & otu_taxa_metadata$time == 4)]
otu_taxa_metadata$value[which(otu_taxa_metadata$Species == 'difficile' & otu_taxa_metadata$group == 'Cdiff' & otu_taxa_metadata$time == 49)]
otu_taxa_metadata$value[which(otu_taxa_metadata$Species == 'difficile' & otu_taxa_metadata$group == 'Cdiff_Csind' & otu_taxa_metadata$time == 4)]
otu_taxa_metadata$value[which(otu_taxa_metadata$Species == 'difficile' & otu_taxa_metadata$group == 'Cdiff_Csind' & otu_taxa_metadata$time == 49)]
otu_taxa_metadata$value[which(otu_taxa_metadata$Species == 'scindens' & otu_taxa_metadata$group == 'Cdiff_Csind' & otu_taxa_metadata$time == 4)]
otu_taxa_metadata$value[which(otu_taxa_metadata$Species == 'scindens' & otu_taxa_metadata$group == 'Cdiff_Csind' & otu_taxa_metadata$time == 49)]

# contaminants
summary(otu_taxa_metadata$value[which(otu_taxa_metadata$Genus == 'Ralstonia' & otu_taxa_metadata$group == 'Cdiff')])
summary(otu_taxa_metadata$value[which(otu_taxa_metadata$Genus == 'Ralstonia' & otu_taxa_metadata$group == 'Cdiff_Csind')])
summary(otu_taxa_metadata$value[which(otu_taxa_metadata$Genus == 'Pseudolabrys' & otu_taxa_metadata$group == 'Cdiff')])
summary(otu_taxa_metadata$value[which(otu_taxa_metadata$Genus == 'Pseudolabrys' & otu_taxa_metadata$group == 'Cdiff_Csind')])
summary(otu_taxa_metadata$value[which(is.na(otu_taxa_metadata$Genus) & otu_taxa_metadata$group == 'Cdiff')])
summary(otu_taxa_metadata$value[which(is.na(otu_taxa_metadata$Genus) & otu_taxa_metadata$group == 'Cdiff_Csind')])

otu_taxa_metadata_contamin_df <- otu_taxa_metadata[which(!otu_taxa_metadata$Species %in% c('difficile', 'scindens')),]
otu_taxa_metadata_contamin_sum_df <- data.frame(otu_taxa_metadata_contamin_df %>% 
                                           dplyr::group_by(Sample, group) %>% 
                                           dplyr::summarise(RA_count = sum(value)) %>% 
                                           dplyr::arrange(desc(RA_count))) 
summary(otu_taxa_metadata_contamin_sum_df$RA_count[which(otu_taxa_metadata_contamin_sum_df$group == 'Cdiff')])
summary(otu_taxa_metadata_contamin_sum_df$RA_count[which(otu_taxa_metadata_contamin_sum_df$group == 'Cdiff_Csind')])

#################
otu_taxa_metadata_top15_df <- data.frame(otu_taxa_metadata %>% 
                                           dplyr::group_by(Genus) %>% 
                                           dplyr::summarise(RA_count = sum(value)) %>% 
                                           dplyr::arrange(desc(RA_count))) 

otu_taxa_metadata_top15_list = otu_taxa_metadata_top15_df$Genus
otu_taxa_metadata_top15_list[is.na(otu_taxa_metadata_top15_list)] <- 'Unknown'
#otu_taxa_metadata_top15_list <- otu_taxa_metadata_top15_list[otu_taxa_metadata_top15_list != "Unknown"]
#otu_taxa_metadata_top15_list <- gsub('(.*)\\[(.*)\\]', "\\2", otu_taxa_metadata_top15_list)

otu_taxa_metadata_top15 <- otu_taxa_metadata
otu_taxa_metadata_top15$Genus[is.na(otu_taxa_metadata_top15$Genus)] <- 'Unknown'
otu_taxa_metadata_top15$Genus[!otu_taxa_metadata_top15$Genus %in% c(otu_taxa_metadata_top15_list, "Unknown")] <- 'Other'
#otu_taxa_metadata_top15$Genus[is.na(otu_taxa_metadata_top15$Genus)] <- 'Other'

# Plot A
otu_taxa_metadata_top15_sum <- data.frame(otu_taxa_metadata_top15 %>% 
                                            dplyr::group_by(Sample, Condition, time, Genus) %>% 
                                            dplyr::summarise(RA_count = sum(value)))

otu_taxa_metadata_top15_sum$Genus<- factor(otu_taxa_metadata_top15_sum$Genus, levels = rev(c(otu_taxa_metadata_top15_list)))
otu_taxa_metadata_top15_sum$Condition <- factor(otu_taxa_metadata_top15_sum$Condition, levels = c('Cdiff', 'Cdiff_Csind'))
##otu_taxa_metadata_top15_sum$sample_alias <- factor(otu_taxa_metadata_top15_sum$sample_alias, levels = c('pA', 'pB', 'pC','pD', 'pE', 'pF'))
library(gtools)
#otu_taxa_metadata_top15_sum$sample_alias <- factor(otu_taxa_metadata_top15_sum$sample_alias, levels=mixedsort(as.character(unique(otu_taxa_metadata_top15_sum$sample_alias))))
otu_taxa_metadata_top15_sum$time <- as.character(otu_taxa_metadata_top15_sum$time)
day_order <- unique(otu_taxa_metadata_top15_sum$time)[order(as.numeric(unique(otu_taxa_metadata_top15_sum$time)))]
otu_taxa_metadata_top15_sum$time <- factor(otu_taxa_metadata_top15_sum$time, levels = day_order, ordered = TRUE)

#variables = unique(otu_taxa_metadata$species[order(-otu_taxa_metadata$value)])
#otu_taxa_metadata$species<- factor(otu_taxa_metadata$species, levels = variables)
colors_kit<-c("Clostridioides"="#1C897E", 
              "Lachnoclostridium" = "#CA5C2D", 
              "Ralstonia" = "#EEDBC4", 
              "Pseudolabrys" = "#DB9B3C", 
              "Unknown" = "#CAE0E6")

out_name_plot <- paste0(out_dir, out_name, '.Top15genera.pdf')
pdf(out_name_plot,width=5,height=3.5,paper="special" )
print(ggplot(data=otu_taxa_metadata_top15_sum, aes(x=time, group = Genus)) +
        geom_bar(aes(y=RA_count, fill = Genus),stat="identity", alpha=.9) +
        theme_minimal() + theme(axis.text = element_text(size = 8),panel.grid.minor = element_blank()) + labs(x = "Time (h)", y = "Relative abundance") +
        theme(title = element_text(size = 10),
              axis.text.x = element_text(color = "grey20", size = 10,angle = 0, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
              axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = 0.5, vjust = 1, face = "plain"),
              axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"), 
              panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())+
        guides(fill=guide_legend(ncol= 1))+ 
        facet_grid(~Condition, scales = "free_x")+
        scale_fill_manual(values = colors_kit, name="Genera"))
dev.off()

# plotB
otu_taxa_metadata_top15 <- otu_taxa_metadata
otu_taxa_metadata_top15_sub <- otu_taxa_metadata_top15[otu_taxa_metadata_top15$Species %in% c('difficile', 'scindens'),]
unique(otu_taxa_metadata_top15_sub$Species)

otu_taxa_metadata_top15_sub$species <- paste(otu_taxa_metadata_top15_sub$Genus, otu_taxa_metadata_top15_sub$Species)
otu_taxa_metadata_top15_sub$species<- factor(otu_taxa_metadata_top15_sub$species, levels = c("Clostridioides difficile", "Lachnoclostridium scindens"))
otu_taxa_metadata_top15_sub$Condition <- factor(otu_taxa_metadata_top15_sub$Condition, levels = c('Cdiff', 'Cdiff_Csind'))
##otu_taxa_metadata_top15_sum$sample_alias <- factor(otu_taxa_metadata_top15_sum$sample_alias, levels = c('pA', 'pB', 'pC','pD', 'pE', 'pF'))
library(gtools)
#otu_taxa_metadata_top15_sum$sample_alias <- factor(otu_taxa_metadata_top15_sum$sample_alias, levels=mixedsort(as.character(unique(otu_taxa_metadata_top15_sum$sample_alias))))
otu_taxa_metadata_top15_sub$time <- as.character(otu_taxa_metadata_top15_sub$time)
day_order <- unique(otu_taxa_metadata_top15_sub$time)[order(as.numeric(unique(otu_taxa_metadata_top15_sub$time)))]
otu_taxa_metadata_top15_sub$time <- factor(otu_taxa_metadata_top15_sub$time, levels = day_order, ordered = TRUE)

#variables = unique(otu_taxa_metadata$species[order(-otu_taxa_metadata$value)])
#otu_taxa_metadata$species<- factor(otu_taxa_metadata$species, levels = variables)


otu_taxa_metadata_top15_sub_plot <- otu_taxa_metadata_top15_sub
otu_taxa_metadata_top15_sub_plot$value[otu_taxa_metadata_top15_sub_plot$value == 0] <- NA
otu_taxa_metadata_top15_sub_plot$time <- as.numeric(as.character(otu_taxa_metadata_top15_sub_plot$time))
out_name_plot <- paste0(out_dir, out_name, '.Cdiff_Cscind_LOESS.pdf')
pdf(out_name_plot,width=2.8,height=4,paper="special" )
colors_kit<-c("Clostridioides difficile"="#1C897E", 
              "Lachnoclostridium scindens" = "#CA5C2D")
print(ggplot(data=otu_taxa_metadata_top15_sub_plot, aes(x=time, value, group=species, color=species)) +
        geom_line() +
        #geom_point(size = 1.5, alpha=.9)  +
        theme_bw()+
        theme(title = element_text(size = 10),
              axis.text.x = element_text(color = "grey20", size = 10,angle = 0, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
              axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = 0.5, vjust = 1, face = "plain"),
              axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"), 
              panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), legend.position="bottom") + 
        guides(color=guide_legend(nrow= 1))+ 
        xlab(NULL)+ylab("Relative abundance")+facet_wrap(.~Condition, ncol=2) + geom_line(size=1.2, alpha=.9) +
        scale_color_manual(values = colors_kit, name="Species"))

test_plot <- ggplot(data=otu_taxa_metadata_top15_sub_plot, aes(x=time, y=value, group = species, color = species)) +
  geom_point(alpha=0.5) + 
  #geom_bar(aes(y=RA_sum, fill = taxa_ID),stat="identity", alpha=.7) +
  theme_bw() + 
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(color = "grey20", size = 10,angle = 0, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = 0.5, vjust = 1, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"), 
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), legend.position="bottom") + 
  guides(color=guide_legend(nrow= 1))+ 
  xlab(NULL)+ylab("Relative abundance")+
  facet_grid(.~Condition)+
  scale_color_manual(values = colors_kit, name="Conditions")
                     
# print(test_plot +geom_smooth(method="lm", se=FALSE, fill=NA,
#                             formula=y ~ poly(x, 3, raw=TRUE)))
print(test_plot +geom_smooth(method="gam"))

colors_kit =c("Cdiff" = "#92392F", "Cdiff_Csind" = "#11698E") 
print(ggplot(data=otu_taxa_metadata_top15_sub_plot, aes(x=time, value, group=Condition, color=Condition)) +
        #geom_point(size = 1.5, alpha=.9)  +
        theme_bw()+
        theme(title = element_text(size = 10),
              axis.text.x = element_text(color = "grey20", size = 10,angle = 0, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
              axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = 0.5, vjust = 1, face = "plain"),
              axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"), 
              panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), legend.position="bottom") + 
        guides(color=guide_legend(nrow= 1))+ 
        xlab(NULL)+ylab("Relative abundance")+facet_wrap(.~species, ncol=2) + geom_line(size=1.2, alpha=.9) +
        scale_color_manual(values = colors_kit, name="Condition"))

test_plot <- ggplot(data=otu_taxa_metadata_top15_sub_plot, aes(x=time, y=value, group = Condition, color = Condition)) +
  geom_point(alpha=0.5) + 
  #geom_bar(aes(y=RA_sum, fill = taxa_ID),stat="identity", alpha=.7) +
  theme_bw() + 
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(color = "grey20", size = 10,angle = 60, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = 0.5, vjust = 1, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"), 
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), legend.position="bottom") + 
  guides(color=guide_legend(nrow= 1))+ 
  xlab(NULL)+ylab("Relative abundance")+
  facet_grid(.~species)+
  scale_color_manual(values = colors_kit, name="Conditions")

# print(test_plot +geom_smooth(method="lm", se=FALSE, fill=NA,
#                        formula=y ~ poly(x, 3, raw=TRUE)))
print(test_plot +geom_smooth(method="loess"))

dev.off()

# Calculate fold change
# TcdA/ tcdB per cell
## Cdiff
otu_taxa_metadata_top15_sub$value[(otu_taxa_metadata_top15_sub$Condition=='Cdiff' & otu_taxa_metadata_top15_sub$time=='4'& otu_taxa_metadata_top15_sub$Species=='difficile')]/otu_taxa_metadata_top15_sub$value[(otu_taxa_metadata_top15_sub$Condition=='Cdiff' & otu_taxa_metadata_top15_sub$time=='49'& otu_taxa_metadata_top15_sub$Species=='difficile')] #  0.9987745
## Cdiff_Csind
otu_taxa_metadata_top15_sub$value[(otu_taxa_metadata_top15_sub$Condition=='Cdiff_Csind' & otu_taxa_metadata_top15_sub$time=='4'& otu_taxa_metadata_top15_sub$Species=='difficile')]/otu_taxa_metadata_top15_sub$value[(otu_taxa_metadata_top15_sub$Condition=='Cdiff_Csind' & otu_taxa_metadata_top15_sub$time=='49'& otu_taxa_metadata_top15_sub$Species=='difficile')] # 3.524017



