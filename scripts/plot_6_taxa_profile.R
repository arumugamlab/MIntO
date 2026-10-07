#!/usr/bin/env Rscript

# ---
# title: "Taxonomic profile"
# author: "Carmen Saenz, Mani Arumugam"
# Generate Visualization plots from assembly-free taxonomy profile (from MetaPhlan or mOTUs3):
#   - PCoA
#   - Stacked barplot of relative abundances of N most abundant taxa (at requested rank) across the samples
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
  make_option(c("--taxrank"),  type="character", default="genus", help="taxonomic rank for plotting top taxa", metavar="character"),
  make_option(c("--toptaxcount"), type="integer", default=15, help="number of top taxa at this rank to plot", metavar="character"),
  make_option(c("--outdir"),   type="character", default=NULL, help="output directory", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Catch missing required args

required_args <- list(
  table = opt$table,
  profiler = opt$profiler,
  metadata = opt$metadata,
  factor = opt$factor,
  outdir = opt$outdir
)

missing_args <- names(required_args)[vapply(required_args, is.null, logical(1))]

if (length(missing_args) > 0) {
  print_help(opt_parser)
  stop(
    "Missing required arguments: ",
    paste(missing_args, collapse = ", "),
    "\n",
    call. = FALSE
  )
}

profile_file = opt$table
profile_param = opt$profiler
out_dir = opt$outdir
metadata_file = opt$metadata
factor = opt$factor
factor2 = opt$factor2
req_rank = opt$taxrank
num_taxa = opt$toptaxcount
out_phyloseq = paste0(out_dir, '/', profile_param,'.phyloseq.rds')

# Validate rank request
valid_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
if (!req_rank %in% valid_ranks) {
  stop("Invalid --taxrank: ", req_rank,
       ". Must be one of: ", paste(valid_ranks, collapse = ", "))
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
  
  # Create the ordination object
  my_physeq = subset_taxa(data_phyloseq, kingdom != 'Unknown')
  ord <- ordinate(my_physeq, method = "PCoA", distance = distance_lab)
  p <- plot_ordination(my_physeq, ord, color = color, shape = shape, label = label)
  
  # Keep the ordination object/axis labels, but remove phyloseq's default layers
  p$layers <- list()
  
  # Set shape/color aesthetics
  point_aes <- ggplot2::aes()
  
  if (!is.null(color)) {
    p$data$.plot_colour <- p$data[[color]]
    point_aes$colour <- rlang::sym(color)
  }
  
  if (!is.null(shape)) {
    p$data$.plot_shape <- as.factor(p$data[[shape]])
    point_aes$shape <- rlang::sym(shape)
  }
  
  # Plot
  p <- p +
    ggplot2::geom_point(
      mapping = point_aes,
      size = 2
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = .data[[label]]),
      size = 3.0,
      segment.alpha = 0.5
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top")
  
  # Add shapes
  if (!is.null(shape)) {
    p <- p + ggplot2::scale_shape_manual(values = c(16:18, 15, 21:25, 0:14))
  }
  
  return(p)
}

# Generate phyloseq object - metaphlan/mOTUs3 output ####
# Common workflow, since we reformat mOTUs3 output like MetaPhlAn output

species_table <- as.data.frame(fread(profile_file, header = T), stringsAsFactors = F) %>%
  filter(grepl('s__|^Unknown', clade_name)) %>%
  filter(!grepl('t__', clade_name)) %>%
  mutate(across('clade_name', \(x) str_replace(x, '^Unknown$', paste(rep('Unknown', 7), collapse='|')))) %>%
  mutate(across('clade_name', \(x) str_replace_all(x, '[kpcofgs]__', ''))) %>%
  tidyr::separate(clade_name,
                  into = valid_ranks,
                  sep = "[\\|]",
                  fill = "right",
                  extra = "drop") %>%
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
  select(-any_of(valid_ranks))
if (profile_param %like% 'metaphlan') {
  otu_table <- otu_table/100
}
otu_table <- as.data.frame(otu_table)

#### Taxa table
taxa_df <- species_table %>%
  select(all_of(valid_ranks))
taxa_df <- as.data.frame(taxa_df)

# Filter by OTUs present in the data
taxa_df   =   taxa_df[rowSums(otu_table)>0 | rownames(otu_table) == 'Unknown',]
otu_table = otu_table[rowSums(otu_table)>0 | rownames(otu_table) == 'Unknown',]

# Catch the edge-case in motus_raw where sample has 0 readcount for all taxa including unknown
# By updating its 'Unknown' to 1, we set Unknown=100% in relative abundance world
empty_sample <- colSums(otu_table) == 0
if (any(empty_sample)) {
  if (!"Unknown" %in% rownames(otu_table)) {
    otu_table["Unknown", ] <- 0
    taxa_df["Unknown", ] <- "Unknown"
  }
  otu_table["Unknown", empty_sample] <- 1
}

# **********************************                                             ********************************
# **********************************          Generate phyloseq object           ********************************
# **********************************                                             ********************************

# # Metadata
metadata_df <- as.data.frame(fread(metadata_file,  header = T), stringsAsFactors = F)

# Validate metadata
if (!"sample" %in% names(metadata_df)) {
  stop("metadata must contain a 'sample' column")
}
if (anyDuplicated(metadata_df$sample)) {
  stop("Duplicate sample IDs in metadata: ",
       paste(unique(metadata_df$sample[duplicated(metadata_df$sample)]), collapse = ", "))
}
if (!factor %in% names(metadata_df)) {
  stop("--factor column not found in metadata: ", factor)
}
if (!is.null(factor2) && !factor2 %in% names(metadata_df)) {
  stop("--factor2 column not found in metadata: ", factor2)
}
if (!is.null(opt$time) && !opt$time %in% names(metadata_df)) {
  stop("--time column not found in metadata: ", opt$time)
}

abund_samples <- colnames(otu_table)
missing_meta <- setdiff(abund_samples, metadata_df$sample)
if (length(missing_meta) > 0) {
  stop("Samples missing from metadata: ", paste(missing_meta, collapse = ", "))
}
metadata_df <- metadata_df[match(abund_samples, metadata_df$sample), , drop = FALSE]
rownames(metadata_df) <- metadata_df$sample

# estimating number of facets/grids for pdf sizing
n_factor <- length(unique(metadata_df[[factor]]))
if (!is.null(factor2)) {
  n_factor <- length(unique(metadata_df[[factor2]]))
}
n_facet_row <- ceiling(n_factor / 10)

# phyloseq object
profile_phyloseq <- phyloseq(otu_table(as.matrix(otu_table), taxa_are_rows = T), tax_table(as.matrix(taxa_df)), sample_data(metadata_df))

# Save taxonomic profile as phyloseq object
saveRDS(profile_phyloseq, file = out_phyloseq)

# Read the phyloseq object
profile_phyloseq <- readRDS(out_phyloseq)

#### metadata:
sample_data_df <- data.frame(sample_data(profile_phyloseq), stringsAsFactors = F)

### OTUs - db
taxa_table_df <- as.data.frame(unclass(tax_table(profile_phyloseq)), stringsAsFactors = F) %>%
  tibble::rownames_to_column("taxa_ID")

### Counts
otu_table_df <- as.data.frame(unclass(otu_table(profile_phyloseq)), stringsAsFactors = F) %>%
  tibble::rownames_to_column("taxa_ID")

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
profile_phyloseq <- prune_samples(sample_sums(profile_phyloseq)>0, profile_phyloseq)

otu_table_df <- as.data.frame(unclass(otu_table(profile_phyloseq)), stringsAsFactors = F) %>%
  tibble::rownames_to_column("taxa_ID")

##########################################
# Output 2: Beta diversity - PCoA (always at species level)
##########################################

distance_lab = 'bray'
title_name <- paste0("PCoA - Taxonomic profile - ", profile_param)
out_name <- paste0(out_dir, '/', profile_param, ".species.PCoA.Bray_Curtis.pdf")
title_name_pval <- paste0("Metric: Bray-Curtis")
if (length(unique(metadata_df[[factor]]))>1) {
  #**adonis/adonis2, Permutational Multivariate Analysis of Variance Using Distance Matrix**: ####
  library(vegan)
  #adonis_list$bray
  otu_table <- as.data.frame(unclass(otu_table(profile_phyloseq )), stringsAsFactors = F)
  metadata_df <- data.frame(sample_data(profile_phyloseq ), stringsAsFactors = F)
  rownames(metadata_df) <- metadata_df$sample
  metadata <-metadata_df[match(names(otu_table), rownames(metadata_df)),]
  dist <- as.dist(vegdist(t(otu_table), method="bray", na.rm = T))
  
  # Create the formula safely and run permanova
  adonis_formula <- stats::as.formula(
    paste0("dist ~ `", factor, "`")
  )
  adonis_bray_Status <- vegan::adonis2(adonis_formula, data = metadata)
  
  r2_value <- format(round(adonis_bray_Status$R2[1],3), nsmall = 3)
  p_value <- adonis_bray_Status$`Pr(>F)`[1]
  
  title_name_pval <- paste0("Metric: Bray-Curtis; PERMANOVA on ", factor, ": R2=", r2_value, ", pval=", p_value)
}

plot_PCoA_out <- plot_PCoA(
  distance_lab,
  profile_phyloseq,
  color = if (!is.null(factor2)) factor2 else factor,
  label = if (!is.null(opt$time)) opt$time else "sample",
  shape = if (!is.null(factor2)) factor else NULL
)

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
        facet_wrap(as.formula(paste(".", "~", factor)))+
        #scale_color_manual(values=manual_plot_colors, name='Condition') +
        coord_fixed() +
        theme(legend.position="bottom")+
        ggtitle(title_name, title_name_pval))
dev.off()

##########################################
# Some common elements for top taxa and richness
##########################################

# Plot variables
sample_var <- if (!is.null(opt$time)) opt$time else "sample"
group_by_vars <- c(factor, sample_var) 
if (!is.null(factor2)) {
  group_by_vars <- c(factor, factor2, sample_var) 
}
group_by_vars <- unique(group_by_vars)

# # Prevalence 10% 
# otu_taxa_filt_df = otu_taxa_merge[rowSums(otu_taxa_merge[, 2:ncol(otu_table_df)])>0.01,]
# otu_taxa_filt_df = otu_taxa_filt_df %>% select(1:ncol(otu_table_df)) #"kingdom", "phylum","class","order","family","genus","mOTU","short_name"

# # After filtering by prevalence 10%, samples don't sum to 1 anymore, make it sum to 1
# colSums(otu_taxa_filt_df[2:ncol(otu_table_df)])
# myNumCols <- which(unlist(lapply(otu_taxa_filt_df, is.numeric)))
# otu_taxa_filt_df[(nrow(otu_taxa_filt_df) + 1), myNumCols] <- 1- colSums(otu_taxa_filt_df[, myNumCols], na.rm=TRUE)
# otu_taxa_filt_df[is.na(otu_taxa_filt_df)] <- '-1'

otu_taxa_melt <- reshape2::melt(otu_table_df, id.vars=c("taxa_ID"))
taxa_table_df <- taxa_table_df %>%
  dplyr::select(taxa_ID, dplyr::all_of(req_rank))
otu_taxa_melt <- merge(taxa_table_df, otu_taxa_melt, by="taxa_ID")

otu_taxa_metadata <-  merge(sample_data_df, otu_taxa_melt, by.x = 'sample', by.y = 'variable') %>%
  dplyr::rename(all_of(c(taxon = req_rank))) %>%
  dplyr::mutate(value = as.numeric(value))

##########################################
# Output 3: Barchart - Top N taxa
##########################################

# Get a list of top N taxa
top_taxa_list <- otu_taxa_metadata %>% 
  dplyr::filter(taxon != "Unknown") %>% 
  dplyr::group_by(taxon) %>% 
  dplyr::summarise(RA_count = sum(value)) %>% 
  dplyr::slice_max(order_by = RA_count, n = num_taxa, with_ties = FALSE) %>%
  dplyr::select(taxon) %>%
  pull()

# Consolidate all taxa not in top N into "Other" category (except "Unknown")
top_taxa_df <- otu_taxa_metadata %>%
  dplyr::mutate(taxon = if_else(taxon %in% c(top_taxa_list, "Unknown"), taxon, "Other"))

top_taxa_df_grouped_summed <- top_taxa_df %>% 
                                           dplyr::group_by(across(all_of(group_by_vars)), taxon, sample) %>% 
                                           dplyr::summarise(RA_count = sum(value), .groups="drop")

# Stack order: most abundant taxa at bottom, then less abundant, then Other/Unknown on top
top_taxa_df_grouped_summed$taxon<- factor(top_taxa_df_grouped_summed$taxon, levels = rev(c(top_taxa_list, 'Other', 'Unknown')))

colors_kit <- rev(c('#9D0208','#FFC87E','#F4A261','#E34F33','#264445',
                    '#E9C46A','#A786C9','#D4C0E2','#975773','#6699FF',
                    '#000066','#7AAFCA','#006699','#A9D181','#2F8475'))
if (num_taxa > 15) {
  colors_kit <- rep(colors_kit,
                    length.out = num_taxa)
}

# Colour order: rank 1 gets colour 1, rank 2 gets colour 2, etc.
taxon_colours <- stats::setNames(
  colors_kit[seq_along(top_taxa_list)],
  top_taxa_list
)
reserved_colours <- c(
  Other = "#BDBDBD",
  Unknown = "#FFFFFF"
)
fill_colours <- c(taxon_colours, reserved_colours)

out_name <- paste0(out_dir, '/', profile_param, '.', req_rank, '.', 'top', num_taxa, '.pdf')
pdf_size <- max(round(n_factor * 0.55, 0), 15)
print(paste(pdf_size, "15top"))
pdf(out_name,width=pdf_size * 0.80,height=pdf_size,paper="special" )
plot_genera_out <- ggplot(
  data=top_taxa_df_grouped_summed,
  aes(x = .data[[sample_var]], y=RA_count, fill = taxon)) +
  geom_bar(stat="identity", alpha=.7) +
  ylim(0, 1) +
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
  scale_fill_manual(values = fill_colours,
                    breaks = rev(c(top_taxa_list, "Other")), # Add "Unknown" in the end if you want it to be explicit
                    name=paste0("Top ", num_taxa, " taxa (", req_rank, ")"),
                    drop = FALSE
                    )

if (!is.null(factor2)) {
  plot_genera_out <- plot_genera_out + facet_grid(as.formula(paste(factor, "~", factor2)), scales = "free")
} else {
  plot_genera_out <- plot_genera_out + facet_wrap(as.formula(paste(factor, "~", ".")), ncol = 5, scales = "free")
}
print(plot_genera_out)
dev.off()

##########################################
# Output 4: Alpha diversity - richness
##########################################

# Remove Unknown/Unassigned
# Count the number of distinct non-Unknown taxa at the requested rank
# Merge with metadata
richness_df <- otu_taxa_metadata %>%
  dplyr::filter(taxon != "Unknown", value > 0) %>%
  dplyr::group_by(across(all_of(group_by_vars)), sample, taxon) %>%
  dplyr::summarise(cumulative = sum(value), .groups="keep") %>%
  dplyr::group_by(across(all_of(group_by_vars)), sample) %>%
  dplyr::summarize(richness = n(), .groups="drop")

# Plot
out_name <- paste0(out_dir, '/', profile_param, '.', req_rank, '.richness.pdf')
pdf_size <- max(round(n_facet_row * 1.8, 0), 10)
pdf(out_name, width=pdf_size, height=pdf_size, paper="special" )

if (!is.null(opt$time)) {
  group_var = if (!is.null(factor2)) factor2 else factor
  richness_plot <-ggplot(data=richness_df, aes(x=.data[[opt$time]], y=richness, group=.data[[group_var]])) +
    geom_line(aes(color=.data[[group_var]])) +
    geom_point() +
    ylim(0, NA) +
    theme(legend.position = "top") +
    theme(axis.text = element_text(size = 8), panel.grid.minor = element_blank()) +
    labs(x = opt$time, y = paste(req_rank, " richness")) +
    theme(title = element_text(size = 10),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    ) +
    facet_wrap(as.formula(paste(".", "~", factor)), scales = "free_x")
} else if (!is.null(factor2)) {
  richness_plot <-ggplot(data=richness_df, aes(x=.data[[factor]], y=richness, group=.data[[factor2]])) +
    geom_line(aes(color=.data[[factor2]])) +
    geom_point() +
    ylim(0, NA) +
    theme(legend.position = "top") +
    theme(axis.text = element_text(size = 8), panel.grid.minor = element_blank()) +
    labs(x = factor, y = paste(req_rank, " richness")) +
    theme(title = element_text(size = 10),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    )
} else {
  richness_plot <-ggplot(data=richness_df, aes(x=.data[[factor]], y=richness)) +
    geom_boxplot() +
    geom_point() +
    ylim(0, NA) +
    theme(legend.position = "top") +
    theme(axis.text = element_text(size = 8), panel.grid.minor = element_blank()) +
    labs(x = factor, y = paste(req_rank, " richness")) +
    theme(title = element_text(size = 10),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    ) 
}
print(richness_plot)
dev.off()
