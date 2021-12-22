#!/usr/bin/env Rscript

# ---
# title: "Taxonomic profile"
# author: "Carmen Saenz"
# date: "09/10/2021"
# Generate Visualization plots from assembly-free taxonomy profile (from MetaPhlan or mOTUs2):
#   - PCoA
#   - Relative abundances of 15 most abundant genera across the samples
# ---

# Load libraries
#if("data.table" %in% rownames(installed.packages()) == FALSE) {install.packages("data.table")}
library(data.table)
#if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse")}
library(tidyverse)
library(tidyr)
#if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2")}
library(ggplot2)
#if("ggrepel" %in% rownames(installed.packages()) == FALSE) {install.packages("ggrepel")}
library(ggrepel)
#if("phyloseq" %in% rownames(installed.packages()) == FALSE) {install.packages("phyloseq")}
library(phyloseq)
#if("reshape2" %in% rownames(installed.packages()) == FALSE) {install.packages("reshape2")}
library(reshape2)
library(vegan)

set.seed(1234)

# Arguments
args = commandArgs(trailingOnly=TRUE)
profile_dir = args[1]
profile_param = args[2]
out_dir = args[3]
metadata_file = args[4]
out_phyloseq = paste0(out_dir, profile_param,'_phyloseq.rds')

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
    geom_text_repel(aes(label = sample_alias), size = 3.0, segment.alpha = 0.5, max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
    geom_point(size = 2, shape = 16)
  PcoA_Sample_site_abundance_2$layers[[1]] <- NULL
  return(PcoA_Sample_site_abundance_2)
}

# Generate phyloseq object - mOUTs2 output ####
if (profile_param %like% 'motus') {
  # **********************************           mOTUs output          ********************************
  files <- list.files(profile_dir, recursive = T)
  files <- files[grepl(paste0("\\.", profile_param,'$'), files)]
  
  countslist <- list()
  for(i in files){
    #readLines(paste0(profile_dir, i)) %>% head(10)
    infile <- as.data.frame(fread(paste0(profile_dir, i),  header = F, skip = 3, sep = '\t'), stringsAsFactors = F)
    names(infile)[1] <- c("mOTUs_id")
    names(infile)[2] <- c("consensus_taxonomy")
    names(infile)[3] <- unlist(lapply(strsplit(gsub(i, pattern=paste0("\\.", profile_param), replacement=""),'\\/'), `[`, 2))
    countslist[[i]] <- infile
  }
  
  countsdf <- data.frame(countslist %>% purrr::reduce(full_join, by = c("mOTUs_id", "consensus_taxonomy")))
  #countsdf[1:10,1:4]
  
  #### Taxa table
  taxa_df <- subset(countsdf, select=c("mOTUs_id","consensus_taxonomy"))
  
  taxa_df <- as.data.frame(taxa_df %>%
                             tidyr::separate(consensus_taxonomy, c("levels1", "levels2"), "k__"))
  taxa_df <- as.data.frame(taxa_df %>%
                             tidyr::separate(levels2, c("kingdom", "levels3"), "\\|p__"))
  taxa_df <- as.data.frame(taxa_df %>%
                             tidyr::separate(levels3, c("phylum", "levels4"), "\\|c__"))
  taxa_df <- as.data.frame(taxa_df %>%
                             tidyr::separate(levels4, c("class", "levels5"), "\\|o__"))
  taxa_df <- as.data.frame(taxa_df %>%
                             tidyr::separate(levels5, c("order", "levels6"), "\\|f__"))
  taxa_df <- as.data.frame(taxa_df %>%
                             tidyr::separate(levels6, c("family", "levels7"), "\\|g__"))
  taxa_df <- as.data.frame(taxa_df %>%
                             tidyr::separate(levels7, c("genus", "species_long"), "\\|s__"))
  taxa_df$levels1 <- NULL
  taxa_df$species <- unlist(lapply(strsplit(taxa_df$species_long,'\\['), `[`, 1))
  taxa_df[taxa_df$mOTUs_id=='-1',] <- c('-1','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown')
  taxa_df$species[taxa_df$species==""] <- taxa_df$species_long[taxa_df$species==""]
  taxa_df$species <- gsub(pattern='species incertae sedis', replacement='', taxa_df$species)
  #taxa_df$species[which(taxa_df$species == "")] <- taxa_df$species[which(taxa_df$species == "")]
  taxa_df$species <- gsub(pattern="\\[", replacement="", taxa_df$species)
  taxa_df$species <- gsub(pattern="\\]", replacement="", taxa_df$species)
  taxa_df$species <- paste(unlist(lapply(strsplit(taxa_df$species,'\\ '), `[`, 1)), unlist(lapply(strsplit(taxa_df$species,'\\ '), `[`, 2)))
  taxa_df$genus <- gsub(pattern="\\ gen\\. incertae sedis", replacement="", taxa_df$genus)
  taxa_df$genus <- unlist(lapply(strsplit(taxa_df$genus,' '), `[`, 1))
  #taxa_df[taxa_df$mOTUs_id=='-1',] <- c('-1','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown')

  #### OTU table
  otu_table <- subset(countsdf, select= -c(consensus_taxonomy))
  rownames(otu_table) <- otu_table$mOTUs_id
  otu_table$mOTUs_id <- NULL
  #### Taxa table
  rownames(taxa_df) <- taxa_df$mOTUs_id
}

# Generate phyloseq object - metaphlan output ####
if (profile_param == 'metaphlan') {
  # **********************************           metaphlan output          ********************************
  # files <- list.files(profile_dir, recursive = T)
  # files <- files[grepl(paste0("\\.", profile_param), files)]
  
  all_taxa <- as.data.frame(fread(paste0(profile_dir, '/merged_abundance_table.txt'),  header = T), stringsAsFactors = F)
  all_taxa_sub <- subset(all_taxa, select = c('clade_name', 'clade_taxid'))
  all_taxa_sub$clade_name <- paste0('|', all_taxa_sub$clade_name)
  all_taxa_sub2 <- as.data.frame(all_taxa_sub %>% tidyr::separate(clade_name, c("delete","kingdom", "phylum", "class", "order", "family", "genus", "species"), "\\|[kpcofgs]__"))
  all_taxa_sub2 <- all_taxa_sub2[!is.na(all_taxa_sub2$species),]
  all_taxa_sub2$delete <- NULL
  rownames(all_taxa_sub2) <- all_taxa_sub2$species
  
  #### OTU table
  sp_taxa <- as.data.frame(fread(paste0(profile_dir, '/merged_abundance_table_species.txt'),  header = T), stringsAsFactors = F)
  names(sp_taxa) <- gsub('\\.metaphlan','', names(sp_taxa))
  otu_table <- sp_taxa
  rownames(otu_table) <- otu_table$body_site
  otu_table$body_site <- NULL
  otu_table <- otu_table/100
  
  #### Taxa table
  taxa_df <- all_taxa_sub2
  rownames(taxa_df) <- taxa_df$species
  taxa_df <- as.data.frame(taxa_df)
}

# **********************************                                             ********************************
# **********************************          Generate phyloseq object           ********************************
# **********************************                                             ********************************

# OTUs
otu_table <- as.data.frame(otu_table)
# TAX table
## Sort TAX table row names as OTUs table row names
otu_table_sort <- otu_table[order(rownames(otu_table)),]
# Sort original data frame by the order of the new data frame
taxa_df_sort <- taxa_df[match(rownames(otu_table_sort), rownames(taxa_df)),]
# # Metadata
if (metadata_file=='None'){
  metadata_df <- data.frame(sample=names(otu_table_sort), condition='condition', sample_alias=names(otu_table_sort))
} else{
  metadata_df <- as.data.frame(fread(metadata_file,  header = T), stringsAsFactors = F)
  names(metadata_df)[[1]] <- "sample"
  names(metadata_df)[[2]] <- "condition"
  names(metadata_df)[[3]] <- "sample_alias"
}
samp <- sample_data(metadata_df)
rownames(samp) <- samp$sample

# phyloseq object
profile_phyloseq <- phyloseq(otu_table(as.matrix(otu_table_sort), taxa_are_rows = T), tax_table(as.matrix(taxa_df_sort)), samp)

head(tax_table(profile_phyloseq))
head(taxa_names(profile_phyloseq))
# Save taxonomic profile as phyloseq object
saveRDS(profile_phyloseq, file = out_phyloseq)

# Phyloseq object
if (profile_param %like% 'motus_raw') {
  profile_phyloseq <- readRDS(out_phyloseq)
  #### Transform counts to Relative Abundance
  profile_phyloseq_ra<-transform_sample_counts(profile_phyloseq, function(x){x/sum(x)})
} else{
  profile_phyloseq_ra <- readRDS(out_phyloseq)
} 

# PCoA - Bray-Curtis ##### ALL samples ####
distance_lab = 'bray'
title_name <- c(paste0("PCoA - Taxonomic profile - ", profile_param), "Bray-Curtis")
out_name <- paste0(out_dir, profile_param, ".PCoA.Bray_Curtis.pdf")
title_name_pval <- paste0("Bray-Curtis")
if(length(unique(metadata_df$condition))>1){
  #**adonis/adonis2, Permutational Multivariate Analysis of Variance Using Distance Matrix**: ####
  library(vegan)
  #adonis_list$bray
  otu_table_sort <- as.data.frame(unclass(otu_table(profile_phyloseq_ra)), stringsAsFactors = F)
  metadata_df <- data.frame(sample_data(profile_phyloseq_ra), stringsAsFactors = F)
  rownames(metadata_df) <- metadata_df$sample
  metadata <-metadata_df[match(names(otu_table_sort), rownames(metadata_df)),]
  #adonis_bray_Status<-adonis(as.dist(vegdist(t(otu_table_sort), method="bray", na.rm = T)) ~ condition, data = metadata)
  #adonis_bray_Status$aov.tab
  adonis_bray_Status2<-adonis2(as.dist(vegdist(t(otu_table_sort), method="bray", na.rm = T)) ~ condition, data = metadata)
  r2_value <- format(round(adonis_bray_Status2$R2[1],3), nsmall = 3)
  p_value <- adonis_bray_Status2$`Pr(>F)`[1]
  
  title_name_pval <- paste0("Bray-Curtis; R2=",r2_value,"; pval=",p_value)
}

#plot_PCoA_out <- plot_PCoA(distance_lab, profile_phyloseq_ra, title_name, out_name)

plot_PCoA_out <- plot_PCoA(distance_lab, profile_phyloseq_ra, "condition", "sample_alias")

manual_plot_colors =c('#9D0208', '#264653','#e9c46a','#D8DCDE','#B6D0E0',
                      '#FFC87E','#F4A261','#E34F33','#E9C46A',
                      '#A786C9','#D4C0E2','#975773','#6699FF','#000066',
                      '#7AAFCA','#006699','#A9D181','#2F8475','#264445') 

# Save ####
pdf(out_name,width=8,height=8,paper="special" )
print(plot_PCoA_out  + 
        scale_color_manual(values=manual_plot_colors, name='Condition') +
        coord_fixed() +
        theme(legend.position="bottom")+
        ggtitle(title_name, title_name_pval))
print(plot_PCoA_out  + 
        facet_wrap(.~condition)+
        scale_color_manual(values=manual_plot_colors, name='Condition') +
        coord_fixed() +
        theme(legend.position="bottom")+
        ggtitle(title_name, title_name_pval))
dev.off()


# Barchart - Top 15 genera####
##### metadata:
sample_data_df <- data.frame(sample_data(profile_phyloseq_ra), stringsAsFactors = F)

### OTUs - db
taxa_table_df <- as.data.frame(unclass(tax_table(profile_phyloseq_ra)), stringsAsFactors = F)
taxa_table_df$taxa_ID <- rownames(taxa_table_df)
rownames(taxa_table_df) <- NULL

### Counts
otu_table_df <- as.data.frame(unclass(otu_table(profile_phyloseq_ra)), stringsAsFactors = F)
otu_table_df$taxa_ID <- rownames(otu_table_df)
rownames(otu_table_df) <- NULL
otu_taxa_merge <- merge(otu_table_df, taxa_table_df, by = 'taxa_ID' , all.x = T)
write.csv(otu_taxa_merge,paste0(out_dir, profile_param, ".csv"), row.names = F, quote = F)

# Filter by mOTUs present in the data
otu_taxa_filt_df = otu_taxa_merge[rowSums(otu_taxa_merge[, 2:ncol(otu_table_df)])>0,]
otu_taxa_filt_df = otu_taxa_filt_df %>% select(1:ncol(otu_table_df))

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

otu_taxa_metadata =  merge(sample_data_df,otu_taxa_melt,  by.x = 'sample', by.y = 'variable')

otu_taxa_metadata$value <- as.numeric(otu_taxa_metadata$value)
otu_taxa_metadata_top15_df <- data.frame(otu_taxa_metadata %>% 
                                     dplyr::group_by(genus) %>% 
                                     dplyr::summarise(RA_count = sum(value)) %>% 
                                     dplyr::arrange(desc(RA_count))) 

otu_taxa_metadata_top15_list = otu_taxa_metadata_top15_df$genus[1:16]
otu_taxa_metadata_top15_list <- otu_taxa_metadata_top15_list[otu_taxa_metadata_top15_list != "Unknown"][1:15]
#otu_taxa_metadata_top15_list <- gsub('(.*)\\[(.*)\\]', "\\2", otu_taxa_metadata_top15_list)

otu_taxa_metadata_top15 <- otu_taxa_metadata
otu_taxa_metadata_top15$genus[!otu_taxa_metadata_top15$genus %in% c(otu_taxa_metadata_top15_list, "Unknown")] <- 'Other'
otu_taxa_metadata_top15$genus[is.na(otu_taxa_metadata_top15$genus)] <- 'Other'

# Plot
otu_taxa_metadata_top15_sum <- data.frame(otu_taxa_metadata_top15 %>% 
                                            dplyr::group_by(sample, condition, sample_alias, genus) %>% 
                                            dplyr::summarise(RA_count = sum(value)))

otu_taxa_metadata_top15_sum$genus<- factor(otu_taxa_metadata_top15_sum$genus, levels = rev(c(otu_taxa_metadata_top15_list, 'Other', 'Unknown')))
#variables = unique(otu_taxa_metadata$species[order(-otu_taxa_metadata$value)])
#otu_taxa_metadata$species<- factor(otu_taxa_metadata$species, levels = variables)

colors_kit<-c('#D8DCDE', '#B6D0E0',
              '#9D0208','#FFC87E','#F4A261','#E34F33','#E9C46A',
              '#A786C9','#D4C0E2','#975773','#6699FF','#000066',
              '#7AAFCA','#006699','#A9D181','#2F8475','#264445')

out_name <- paste0(out_dir, profile_param, '.Top15genera.pdf')
pdf(out_name,width=15,height=6,paper="special" )
print(ggplot(data=otu_taxa_metadata_top15_sum, aes(x=sample_alias, group = genus)) +
        geom_bar(aes(y=RA_count, fill = genus),stat="identity", alpha=.7) +
        theme_minimal() + theme(axis.text = element_text(size = 8),panel.grid.minor = element_blank()) + labs(x = "Samples", y = "Relative abundance") +
        theme(title = element_text(size = 10),
              axis.text.x = element_text(color = "grey20", size = 10,angle = 60, vjust =1,face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
              axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = 0.5, vjust = 1, face = "plain"),
              axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"), 
              panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())+
        guides(fill=guide_legend(ncol= 1))+ 
        facet_wrap(.~condition, scales = "free_x")+
        scale_fill_manual(values = colors_kit, name="Top 15 genera"))
dev.off()
