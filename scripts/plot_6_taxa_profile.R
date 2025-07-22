#!/usr/bin/env Rscript

# ---
# title: "Taxonomic profile"
# author: "Carmen Saenz, Mani Arumugam"
# date: "09/10/2021"
# Generate Visualization plots from assembly-free taxonomy profile (from MetaPhlan or mOTUs3):
#   - PCoA
#   - Relative abundances of 15 most abundant genera across the samples
#   - Alpha diversity plot
# ---

# Parse command line arguments
library(optparse)
option_list = list(
                make_option(c("--table"),    type="character", default=NULL, help="taxonomic profile table", metavar="character"),
                make_option(c("--profiler"), type="character", default=NULL, help="name of the taxonomic profiler", metavar="character"),
                make_option(c("--metadata"), type="character", default=NULL, help="metadata file", metavar="character"),
                make_option(c("--factor"),   type="character", default=NULL, help="name of key factor from metadata file", metavar="character"),
                make_option(c("--factor2"),  type="character", default=NULL, help="name of 2nd factor from metadata file", metavar="character"),
                make_option(c("--time"),     type="character", default=NULL, help="name of time variable from metadata file", metavar="character"),
                make_option(c("--outdir"),   type="character", default=NULL, help="output directory", metavar="character")
                )
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

profile_file = opt$table
profile_param = opt$profiler
out_dir = opt$outdir
metadata_file = opt$metadata
out_phyloseq = paste0(out_dir, '/', profile_param,'.phyloseq.rds')

if (any(is.null(c(opt$table, opt$profiler, opt$metadata, opt$factor, opt$outdir)))) {
  print_help(opt_parser)
  stop("Missing required arguments\n", call.=FALSE)
}

# Load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(phyloseq)
library(reshape2)
library(vegan)

set.seed(1234)

##########################  ** Function - PCoA **  ########################## 
plot_PCoA <- function(distance_lab, data_phyloseq, color, label, shape=NULL){ #out_name,title_name,
  label2 <- noquote(label)
  #### PCoA
  my_physeq = subset_taxa(data_phyloseq, kingdom != 'Unknown')
  ord<- ordinate(my_physeq, method = "PCoA", distance = distance_lab)
  PcoA_Sample_site_abundance <- plot_ordination(my_physeq, ord, color = color, shape = as.factor(shape), label = label) #type = type,
  PcoA_Sample_site_abundance <- PcoA_Sample_site_abundance +
    ggtitle(title_name, distance_lab)  +
    geom_point(size = 2.5) +
    theme_bw() +
    theme(legend.position = "top")
  PcoA_Sample_site_abundance$layers[[1]] <- NULL
  PcoA_Sample_site_abundance$layers[[2]] <- NULL
  PcoA_Sample_site_abundance$layers[[3]] <- NULL
  PcoA_Sample_site_abundance_2<- PcoA_Sample_site_abundance +
    geom_text_repel(aes(label = get(label)), 
                          size = 3.0, 
                          segment.alpha = 0.5, 
                          max.overlaps = getOption("ggrepel.max.overlaps", default = 20)
                    )
  if (!is.null(shape)) {
    PcoA_Sample_site_abundance_2 <- PcoA_Sample_site_abundance_2 + geom_point(aes(shape = as.factor(get(shape))), size = 2)
  } else {
    PcoA_Sample_site_abundance_2 <- PcoA_Sample_site_abundance_2 + geom_point(shape = 16, size = 2)
  }
  PcoA_Sample_site_abundance_2$layers[[1]] <- NULL
  return(PcoA_Sample_site_abundance_2)
}

# Generate phyloseq object - metaphlan/mOTUs3 output ####
# Common workflow, since we reformat mOTUs3 output like MetaPhlAn output
  
  species_table <- as.data.frame(fread(profile_file, header = T), stringsAsFactors = F) %>%
    filter(grepl('s__|^Unknown', clade_name)) %>%
    filter(!grepl('t__', clade_name)) %>%
    mutate(across('clade_name', \(x) str_replace_all(x, '[kpcofgs]__', ''))) %>%
    tidyr::separate(clade_name, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), "[\\|]") %>%
    mutate_all(~replace_na(., "Unknown"))

  if (profile_param %like% 'motus') {
    rnames <- str_extract(species_table$species, "\\[(\\S+)\\]$", group=1)
    rnames[is.na(rnames)] <- 'Unknown'
    rownames(species_table) <- rnames
  } else {
    rownames(species_table) <- paste0(species_table$species)
  }

  #### OTU table
  otu_table <- species_table %>%
    select(-c("kingdom", "phylum", "class", "order", "family", "genus", "species"))
  if (profile_param %like% 'metaphlan') {
    otu_table <- otu_table/100
  }
  otu_table <- as.data.frame(otu_table)

  #### Taxa table
  taxa_df <- species_table %>%
    select(c("kingdom", "phylum", "class", "order", "family", "genus", "species"))
  taxa_df <- as.data.frame(taxa_df)

  # Filter by OTUs present in the data
  taxa_df   =   taxa_df[rowSums(otu_table)>0 | rownames(otu_table) == 'Unknown',]
  otu_table = otu_table[rowSums(otu_table)>0 | rownames(otu_table) == 'Unknown',]

  # Catch the edge-case in motus_raw where sample has 0 readcount for all taxa including unknown
  # By updating its 'Unknown' to 1, we set Unknown=100% in relative abundance world
  to_update = colSums(otu_table) == 0
  otu_table[c("Unknown"), to_update] = 1

# **********************************                                             ********************************
# **********************************          Generate phyloseq object           ********************************
# **********************************                                             ********************************

# # Metadata
if (is.null(metadata_file)) {
  metadata_df <- data.frame(sample=names(otu_table), factor='default', sample_alias=names(otu_table))
  names(metadata_df)[2] <- opt$factor
} else{
  metadata_df <- as.data.frame(fread(metadata_file,  header = T), stringsAsFactors = F)
}
samp <- sample_data(metadata_df)
rownames(samp) <- samp$sample


# estimating number of facets/grids for pdf sizing
n_factor <- length(unique(metadata_df[[opt$factor]]))
if (!is.null(opt$factor2)) {
  n_factor <- length(unique(metadata_df[[opt$factor2]]))
}
n_facet_row <- ceiling(n_factor / 10)

# phyloseq object
profile_phyloseq <- phyloseq(otu_table(as.matrix(otu_table), taxa_are_rows = T), tax_table(as.matrix(taxa_df)), samp)

#head(tax_table(profile_phyloseq))
#head(taxa_names(profile_phyloseq))

# Save taxonomic profile as phyloseq object
saveRDS(profile_phyloseq, file = out_phyloseq)


# Read the phyloseq object
profile_phyloseq <- readRDS(out_phyloseq)

#### metadata:
sample_data_df <- data.frame(sample_data(profile_phyloseq), stringsAsFactors = F)

### OTUs - db
taxa_table_df <- as.data.frame(unclass(tax_table(profile_phyloseq)), stringsAsFactors = F)
taxa_table_df$taxa_ID <- rownames(taxa_table_df)
rownames(taxa_table_df) <- NULL

### Counts
otu_table_df <- as.data.frame(unclass(otu_table(profile_phyloseq)), stringsAsFactors = F)

otu_table_df$taxa_ID <- rownames(otu_table_df)
rownames(otu_table_df) <- NULL

##########################################
# Output 1: Merged profile file with abundance and taxonomic annotation
##########################################

# Merge taxonomy and abundance
# The logical in arrange() puts 'Unknown' as the first row
otu_taxa_merge <- merge(otu_table_df, taxa_table_df, by = 'taxa_ID' , all.x = T) %>%
                    dplyr::arrange(taxa_ID != 'Unknown')

# Write output
fwrite(otu_taxa_merge, file = paste0(out_dir, '/', profile_param, ".tsv"), sep = '\t', row.names = F, quote = F)

# Convert to relative abundance if necessary
# From now on, phyloseq object is only in RA mode
if (profile_param %like% 'motus_raw') {
    profile_phyloseq <-transform_sample_counts(profile_phyloseq, function(x){x/sum(x)})
}

# Remove 'Unknown'
profile_phyloseq = subset_taxa(profile_phyloseq, kingdom != 'Unknown')

# Remove samples without taxa
profile_phyloseq <- prune_samples(sample_sums(profile_phyloseq)>0.1, profile_phyloseq)

otu_table_df <- as.data.frame(unclass(otu_table(profile_phyloseq)), stringsAsFactors = F)
otu_table_df$taxa_ID <- rownames(otu_table_df)
rownames(otu_table_df) <- NULL

##########################################
# Output 2: Beta diversity - PCoA
##########################################

distance_lab = 'bray'
title_name <- c(paste0("PCoA - Taxonomic profile - ", profile_param), "Bray-Curtis")
out_name <- paste0(out_dir, '/', profile_param, ".PCoA.Bray_Curtis.pdf")
title_name_pval <- paste0("Metric: Bray-Curtis")
if(length(unique(metadata_df[[opt$factor]]))>1){
  #**adonis/adonis2, Permutational Multivariate Analysis of Variance Using Distance Matrix**: ####
  library(vegan)
  #adonis_list$bray
  otu_table <- as.data.frame(unclass(otu_table(profile_phyloseq )), stringsAsFactors = F)
  metadata_df <- data.frame(sample_data(profile_phyloseq ), stringsAsFactors = F)
  rownames(metadata_df) <- metadata_df$sample
  metadata <-metadata_df[match(names(otu_table), rownames(metadata_df)),]
  dist <- as.dist(vegdist(t(otu_table), method="bray", na.rm = T))
  adonis_bray_Status<-adonis2(as.formula(paste("dist", "~", opt$factor)), data = metadata)
  r2_value <- format(round(adonis_bray_Status$R2[1],3), nsmall = 3)
  p_value <- adonis_bray_Status$`Pr(>F)`[1]
  
  title_name_pval <- paste0("Metric: Bray-Curtis; PERMANOVA on ", opt$factor, ": R2=", r2_value, ", pval=", p_value)
}

#plot_PCoA_out <- plot_PCoA(distance_lab, profile_phyloseq , title_name, out_name)

plot_PCoA_out <- plot_PCoA(distance_lab, profile_phyloseq , color=opt$factor2, label = if (!is.null(opt$time)) opt$time else "sample", shape=opt$factor)

manual_plot_colors =c('#9D0208', '#264653','#e9c46a','#D8DCDE','#B6D0E0',
                      '#FFC87E','#F4A261','#E34F33','#E9C46A',
                      '#A786C9','#D4C0E2','#975773','#6699FF','#000066',
                      '#7AAFCA','#006699','#A9D181','#2F8475','#264445') 

# Save ####
pdf_size <- max(round(n_facet_row * 1.8, 0), 8)
pdf(out_name,width=pdf_size,height=pdf_size,paper="special" )
print(plot_PCoA_out  + 
        #scale_color_manual(values=manual_plot_colors, name='Condition') +
        coord_fixed() +
        theme(legend.position="bottom")+
        ggtitle(title_name, title_name_pval))
print(plot_PCoA_out  + 
        facet_wrap(as.formula(paste(".", "~", opt$factor)))+
        #scale_color_manual(values=manual_plot_colors, name='Condition') +
        coord_fixed() +
        theme(legend.position="bottom")+
        ggtitle(title_name, title_name_pval))
dev.off()

##########################################
# Output 3: Barchart - Top 15 genera
##########################################

# # Prevalence 10% 
# otu_taxa_filt_df = otu_taxa_merge[rowSums(otu_taxa_merge[, 2:ncol(otu_table_df)])>0.01,]
# otu_taxa_filt_df = otu_taxa_filt_df %>% select(1:ncol(otu_table_df)) #"kingdom", "phylum","class","order","family","genus","mOTU","short_name"

# # After filtering by prevalence 10%, samples don't sum to 1 anymore, make it sum to 1
# colSums(otu_taxa_filt_df[2:ncol(otu_table_df)])
# myNumCols <- which(unlist(lapply(otu_taxa_filt_df, is.numeric)))
# otu_taxa_filt_df[(nrow(otu_taxa_filt_df) + 1), myNumCols] <- 1- colSums(otu_taxa_filt_df[, myNumCols], na.rm=TRUE)
# otu_taxa_filt_df[is.na(otu_taxa_filt_df)] <- '-1'

otu_taxa_melt <- melt(otu_table_df, id.vars=c("taxa_ID"))
otu_taxa_melt <- merge(taxa_table_df, otu_taxa_melt, by="taxa_ID")
#colnames(otu_taxa_melt) <- c('taxa_ID', "kingdom", "phylum", "class", "order", "family", "genus", "mOTU", "short_name", "sample", "RA")

otu_taxa_metadata =  merge(sample_data_df, otu_taxa_melt, by.x = 'sample', by.y = 'variable')

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
sample_var <- if (!is.null(opt$time)) opt$time else "sample"
group_by_vars <- c(opt$factor, sample_var) 
if (!is.null(opt$factor2)) {
    group_by_vars <- c(opt$factor, opt$factor2, sample_var) 
}
otu_taxa_metadata_top15_sum <- data.frame(otu_taxa_metadata_top15 %>% 
                                            dplyr::group_by(across(all_of(group_by_vars)), genus, sample) %>% 
                                            dplyr::summarise(RA_count = sum(value), .groups="drop_last") %>%
                                            ungroup())
if (!is.null(opt$time)){
  otu_taxa_metadata_top15_sum <- otu_taxa_metadata_top15_sum %>% 
    dplyr::group_by(across(all_of(group_by_vars)), genus, sample) %>%
    dplyr::summarise(RA_count = mean(RA_count), .groups="drop_last")
}


otu_taxa_metadata_top15_sum$genus<- factor(otu_taxa_metadata_top15_sum$genus, levels = rev(c(otu_taxa_metadata_top15_list, 'Other', 'Unknown')))
#variables = unique(otu_taxa_metadata$species[order(-otu_taxa_metadata$value)])
#otu_taxa_metadata$species<- factor(otu_taxa_metadata$species, levels = variables)

colors_kit<-c('#D8DCDE', '#B6D0E0',
              '#9D0208','#FFC87E','#F4A261','#E34F33','#E9C46A',
              '#A786C9','#D4C0E2','#975773','#6699FF','#000066',
              '#7AAFCA','#006699','#A9D181','#2F8475','#264445')

out_name <- paste0(out_dir, '/', profile_param, '.Top15genera.pdf')
pdf_size <- max(round(n_factor * 0.55, 0), 15)
print(paste(pdf_size, "15top"))
pdf(out_name,width=pdf_size * 0.80,height=pdf_size,paper="special" )
plot_genera_out <- ggplot(data=otu_taxa_metadata_top15_sum, aes(x = as.factor(.data[[sample_var]]), group = genus)) +
        geom_bar(aes(y=RA_count, fill = genus), stat="identity", alpha=.7) +
        theme_minimal() + 
        theme(axis.text = element_text(size = 8), panel.grid.minor = element_blank()) + 
        labs(x = "Samples", y = "Relative abundance") +
        theme(title = element_text(size = 10),
              axis.text.x = element_text(color = "grey20", size = 10, angle = 60, hjust = 1.00, vjust = 1.00, face = "plain"),
              axis.text.y = element_text(color = "grey20", size = 10, angle = 00, hjust = 1.00, vjust = 0.00, face = "plain"),
              axis.title.x = element_text(color = "grey20", size = 12, angle = 00, hjust = 0.5, vjust = 1.0, face = "plain"),
              axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5, face = "plain"), 
              panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
              ) +
        guides(fill=guide_legend(ncol= 1)) +
#        theme(panel.margin.y = unit(0, "lines")) +
        scale_fill_manual(values = colors_kit, name="Top 15 genera")

if (!is.null(opt$factor2)) {
    plot_genera_out <- plot_genera_out + facet_grid(as.formula(paste(opt$factor, "~", opt$factor2)), scales = "free")
} else {
    plot_genera_out <- plot_genera_out + facet_wrap(as.formula(paste( opt$factor, "~", ".")), ncol = 5, scales = "free")
}
print(plot_genera_out)
dev.off()

##########################################
# Output 4: Alpha diversity - richness
##########################################

# Remove Unknown/Unassigned
richness_df <- filter(otu_table_df, taxa_ID != "Unknown")

# Make df
richness_df <- as.data.frame(colSums(richness_df>0, na.rm = TRUE))
names(richness_df) <- c("richness")
richness_df["sample"] <- row.names(richness_df)
richness_df <- inner_join(richness_df, sample_data_df, by="sample")


# Plot
out_name <- paste0(out_dir, '/', profile_param, '.richness.pdf')
pdf_size <- max(round(n_facet_row * 1.8, 0), 10)
pdf(out_name, width=pdf_size, height=pdf_size, paper="special" )

if (!is.null(opt$time)) {
    group_var = if (!is.null(opt$factor2)) opt$factor2 else opt$factor
    richness_plot <-ggplot(data=richness_df, aes(x=.data[[opt$time]], y=richness, group=.data[[group_var]])) +
                        geom_line(aes(color=.data[[group_var]])) +
                        geom_point() +
                        ylim(0, NA) +
                        theme(legend.position = "top") +
                        theme(axis.text = element_text(size = 8), panel.grid.minor = element_blank()) +
                        labs(x = opt$time, y = "Richness") +
                        theme(title = element_text(size = 10),
                              panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                            ) +
                        facet_wrap(as.formula(paste(".", "~", opt$factor)), scales = "free_x")
} else if (!is.null(opt$factor2)) {
    richness_plot <-ggplot(data=richness_df, aes(x=.data[[opt$factor]], y=richness, group=.data[[opt$factor2]])) +
                        geom_line(aes(color=.data[[opt$factor2]])) +
                        geom_point() +
                        ylim(0, NA) +
                        theme(legend.position = "top") +
                        theme(axis.text = element_text(size = 8), panel.grid.minor = element_blank()) +
                        labs(x = opt$factor, y = "Richness") +
                        theme(title = element_text(size = 10),
                              panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                            )
} else {
    richness_plot <-ggplot(data=richness_df, aes(x=.data[[opt$factor]], y=richness)) +
                        geom_boxplot() +
                        geom_point() +
                        ylim(0, NA) +
                        theme(legend.position = "top") +
                        theme(axis.text = element_text(size = 8), panel.grid.minor = element_blank()) +
                        labs(x = opt$factor, y = "Richness") +
                        theme(title = element_text(size = 10),
                              panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                            ) 
}
print(richness_plot)
dev.off()
