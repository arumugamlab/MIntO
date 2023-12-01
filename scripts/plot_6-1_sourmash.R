#!/usr/bin/env Rscript

# ---
# title: "Sourmash clustering"
# author: "Judit Szarvas, Carmen Saenz, Mani Arumugam"
# date: "23/11/2023"
# Generate Visualization plots from k-mer/minhash based similarity matrix (sourmash):
#   - UMAP
#   - Average linkage hierarchical clustering dendogram
#   - Relative abundances of 15 most abundant genera across the samples in clusters
# ---

# Parse command line arguments
library(optparse)
option_list = list(
                make_option(c("--csv"),    type="character", default=NULL, help="sourmash compare csv", metavar="character"),
                make_option(c("--cutoff"), type="double", default=0.6, help="cutoff for hierarchical clustering", metavar="double"),
                make_option(c("--table"),    type="character", default=NULL, help="merged taxonomy profile", metavar="character"),
                make_option(c("--metadata"), type="character", default=NULL, help="metadata file", metavar="character"),
                make_option(c("--factor"),   type="character", default=NULL, help="name of key factor from metadata file", metavar="character"),
                make_option(c("--factor2"),  type="character", default=NULL, help="name of 2nd factor from metadata file", metavar="character"),
                make_option(c("--time"),     type="character", default=NULL, help="name of time variable from metadata file", metavar="character"),
                make_option(c("--outdir"),   type="character", default=NULL, help="output directory", metavar="character")
                )
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (any(is.null(c(opt$table, opt$metadata, opt$factor, opt$outdir)))) {
  print_help(opt_parser)
  stop("Missing required arguments\n", call.=FALSE)
}

# Load libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(umap)
library(dendextend)

set.seed(1234)

##########################  ** Functions **  ########################## 
plot_umap_layout <- function(dist_matrix, label_data, sampleid_v, color_var, shape_var){
  set.seed(1234)
  umap_est <- umap(dist_matrix, input="dist", preserve.seed = T)
  t.label <-  left_join(data.frame(sample=sampleid_v), label_data) %>% select(sample, {{ shape_var }}, {{ color_var }})
  umap_layout_df <- data.frame(umap_est$layout)
  umap_layout_df$sample <- rownames(umap_layout_df)
  umap_layout_df <- left_join(umap_layout_df, t.label, by = join_by(sample))

  p <- ggplot(umap_layout_df, aes(x=X1, y=X2))
  if (! is.null(shape_var)){
    p <- p + 
      geom_point(aes(color = .data[[color_var]], shape = as.factor(.data[[shape_var]])))
  } else {
    p <- p + 
      geom_point(aes(color = .data[[color_var]]))
  }

  p <- p +
      xlab("Dim 1") +
      ylab("Dim 2") +
      ggtitle("UMAP") +
      theme(plot.margin = margin(1, 4, 1, 1, "cm"))
  return(p)
}

plot_branch_colored_dnd <- function(dist_matrix, cut_height = 0.60, label_data, sampleid_v, linkage_method = "average", factor, factor2){
  t.label <- left_join(data.frame(sample=sampleid_v), label_data)

  hc <- hclust(as.dist(dist_matrix), method = linkage_method)
  dend <- as.dendrogram(hc)
  
  # create subclusters and color braches of those with more than one tip
  t_subclusters <- cutree(hc, h = cut_height, order_clusters_as_data = F)
  v <- data.frame(x=t_subclusters) %>% group_by(x) %>% count() %>% filter(n==1)
  t_subclusters[which(t_subclusters %in% v$x)] <- 0
  dend <- dendextend::color_branches(dend, clusters = t_subclusters) %>% set("branches_lwd", 2) %>% set("labels_cex", 0.1)
  plot(dend, horiz = F)
  
  # color bar below dendogram
  bcolor <- as.numeric(as.factor(t.label[[factor]]))
  rL <- c(factor)
  if (! is.null(factor2)){
    bcolor_1 <- as.numeric(as.factor(t.label[[factor2]]))
    bcolor <- cbind(bcolor_1, bcolor)
    rL <- c(factor2, factor)
  }
  colored_bars(colors = bcolor, dend = dend, rowLabels = rL, y_shift = 0)

  abline(h=cut_height,col="red")
}

get_hcl_subclucters <- function(dist_matrix, cut_height = 0.60, linkage_method = "average"){
  hc <- hclust(as.dist(dist_matrix), method = linkage_method)
  dend <- as.dendrogram(hc)
  
  t_subclusters <- cutree(hc, h = cut_height, order_clusters_as_data = F)
  v <- data.frame(x=t_subclusters) %>% group_by(x) %>% count() %>% filter(n==1)
  t_subclusters <- t_subclusters[which(! t_subclusters %in% v$x)]
  t_subclusters_df <- data.frame(clustering=t_subclusters)
  t_subclusters_df$sample <- rownames(t_subclusters_df)
  return(t_subclusters_df)
}

get_minsized_clusters <- function(sc_df, min_size = 3){
  cl_v <- {sc_df %>% group_by(clustering) %>% count() %>% filter(n > min_size-1)}$clustering
  return(cl_v)
}

plot_cluster_top15_genera <- function(top15_sum_df = otu_taxa_metadata_top15_sum, sample_var = "sample", row_factor_var, col_factor_var = NULL, cluster_mmbs, sample_metadata_df = NULL){
  if (! is.null(sample_metadata_df)){
    top15_sum_df <- left_join(top15_sum_df, sample_metadata_df, by = join_by(sample)) %>% 
      select({{ sample_var }}, genus, RA_count, {{ row_factor_var }}, {{ col_factor_var}})
  }
  top15_sum_df <- top15_sum_df %>% filter(sample %in% cluster_mmbs)
  plot_genera_out <- ggplot(data=top15_sum_df, aes(x = .data[[sample_var]], group = genus)) +
        geom_col(aes(y=RA_count, fill = genus), alpha=.7) +
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
        scale_fill_manual(values = colors_kit, name="Top 15 genera")

  if (!is.null(col_factor_var) & !is.null(row_factor_var)) {
    plot_genera_out <- plot_genera_out + facet_grid(as.formula(paste(row_factor_var, "~", col_factor_var)), scales = "free_x", space='free') }
  else if (!is.null(row_factor_var)) {
    plot_genera_out <- plot_genera_out + facet_wrap(vars({{ row_factor_var }}), ncol = 1, strip.position = "right")
  } else {
    plot_genera_out <- plot_genera_out
  }
    return(plot_genera_out)
}

read_smash_sim_csv <- function(csv_file){
  similarity_df <-as.data.frame(data.table::fread(csv_file, header = T), stringsAsFactors = F)
  rownames(similarity_df) <- colnames(similarity_df)
  similarity_matrix <- as.matrix(similarity_df)
  return(similarity_matrix)
}

# **********************************                                             ********************************
# **********************************          Generate top15 taxa                ********************************
# **********************************                                             ********************************


# Generate phyloseq object - metaphlan/mOTUs3 output ####
# Common workflow, since we reformat mOTUs3 output like MetaPhlAn output

output_prefix <- str_remove(basename(opt$table), ".merged_abundance_table.txt")
profile_param <- str_split_1(basename(opt$table), "\\.")[2]

species_table <- as.data.frame(data.table::fread(opt$table, header = T), stringsAsFactors = F) %>%
  filter(grepl('s__', clade_name)) %>%
  filter(!grepl('t__', clade_name)) %>%
  mutate(across('clade_name', \(x) str_replace_all(x, '[kpcofgs]__', ''))) %>%
   tidyr::separate(clade_name, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), "[\\|]")
    
if (str_starts(profile_param, "motus")) {
    rownames(species_table) <- str_extract(species_table$species, "\\[(\\S+)\\]$", group=1)
  } else {
    rownames(species_table) <- paste0(species_table$species)
  }
  
#### OTU table
otu_table <- species_table %>%
    select(-c("kingdom", "phylum", "class", "order", "family", "genus", "species"))
if (str_starts(profile_param, "metaphlan")) {
    otu_table <- otu_table/100
  }
  
#### Taxa table
taxa_df <- species_table %>%
    select(c("kingdom", "phylum", "class", "order", "family", "genus", "species"))

# # Metadata
if (is.null(opt$metadata)) {
  metadata_df <- data.frame(sample=names(otu_table), factor='default', sample_alias=names(otu_table))
  names(metadata_df)[2] <- opt$factor
} else{
  metadata_df <- as.data.frame(data.table::fread(opt$metadata,  header = T), stringsAsFactors = F)
  #metadata_df <- metadata_df %>% rename(sample="gid_wgs")
}

# transform raw to rel. abundance
if (str_starts(profile_param, "motus_raw")) {
  #### Transform counts to Relative Abundance
  rav <- function(x)(x/sum(x))
  otu_table_ra<- otu_table %>% mutate_all(rav)
} else{
  otu_table_ra <- otu_table
}

#### metadata
sample_data_df <- metadata_df

### OTUs - db
taxa_table_df <- taxa_df
taxa_table_df$taxa_ID <- rownames(taxa_table_df)
rownames(taxa_table_df) <- NULL

### Counts
otu_table_df <- otu_table_ra
otu_table_df$taxa_ID <- rownames(otu_table_df)
rownames(otu_table_df) <- NULL

otu_taxa_merge <- merge(otu_table_df, taxa_table_df, by = 'taxa_ID' , all.x = T)

# Filter by mOTUs present in the data
otu_taxa_filt_df = otu_taxa_merge[rowSums(otu_taxa_merge[, 2:ncol(otu_table_df)])>0,]
otu_taxa_filt_df = otu_taxa_filt_df %>% select(1:ncol(otu_table_df))

otu_taxa_melt <- melt(otu_taxa_filt_df, id.vars=c("taxa_ID"))
otu_taxa_melt <- merge(taxa_table_df, otu_taxa_melt, by="taxa_ID")

otu_taxa_metadata <-  merge(sample_data_df, otu_taxa_melt, by.x = 'sample', by.y = 'variable')

otu_taxa_metadata$value <- as.numeric(otu_taxa_metadata$value)
otu_taxa_metadata_top15_df <- data.frame(otu_taxa_metadata %>% 
                                     dplyr::group_by(genus) %>% 
                                     dplyr::summarise(RA_count = sum(value)) %>% 
                                     dplyr::arrange(desc(RA_count))) 

otu_taxa_metadata_top15_list = otu_taxa_metadata_top15_df$genus[1:16]
otu_taxa_metadata_top15_list <- otu_taxa_metadata_top15_list[otu_taxa_metadata_top15_list != "Unknown"][1:15]

otu_taxa_metadata_top15 <- otu_taxa_metadata
otu_taxa_metadata_top15$genus[!otu_taxa_metadata_top15$genus %in% c(otu_taxa_metadata_top15_list, "Unknown")] <- 'Other'
otu_taxa_metadata_top15$genus[is.na(otu_taxa_metadata_top15$genus)] <- 'Other'


otu_taxa_metadata_top15_sum <- data.frame(otu_taxa_metadata_top15 %>% 
                                            dplyr::group_by(sample, genus) %>% 
                                            dplyr::summarise(RA_count = sum(value), .groups="drop_last"))

otu_taxa_metadata_top15_sum$genus<- factor(otu_taxa_metadata_top15_sum$genus, levels = rev(c(otu_taxa_metadata_top15_list, 'Other', 'Unknown')))

colors_kit<-c('#D8DCDE', '#B6D0E0',
              '#9D0208','#FFC87E','#F4A261','#E34F33','#E9C46A',
              '#A786C9','#D4C0E2','#975773','#6699FF','#000066',
              '#7AAFCA','#006699','#A9D181','#2F8475','#264445')


rm(taxa_table_df, sample_data_df, otu_table_df, otu_table_ra, otu_table)
rm(otu_taxa_metadata_top15, otu_taxa_metadata, otu_taxa_metadata_top15_df, otu_taxa_filt_df, otu_taxa_merge, otu_taxa_melt)

##########################################
# Clustering
##########################################

smash_cosine_matrix <- read_smash_sim_csv(opt$csv)
smash_cosine_matrix.ids <- rownames(smash_cosine_matrix)
smash_dist_mat <- 1-smash_cosine_matrix

smash_subcl <- get_hcl_subclucters(smash_dist_mat, cut_height = opt$cutoff, linkage_method = "average")

out_name <- paste0(opt$outdir, '/', output_prefix, '.sourmash_clusters.tsv')
write.table(smash_subcl, file = out_name, sep = "\t", row.names = FALSE, col.names = TRUE)

##########################################
# Dendogram
##########################################

out_name <- paste0(opt$outdir, '/', str_split_1(output_prefix, "\\.")[1], '.sourmash_plots.pdf')
pdf(out_name, width=12, height=6, paper="special" )

print(plot_umap_layout(smash_dist_mat, metadata_df, smash_cosine_matrix.ids, opt$factor, opt$factor2))

print(plot_branch_colored_dnd(smash_dist_mat, cut_height = opt$cutoff, metadata_df, smash_cosine_matrix.ids, linkage_method = "average", opt$factor, opt$factor2))
dev.off()

##########################################
# Barplots
##########################################

out_name <- paste0(opt$outdir, '/', output_prefix, '.clusters.pdf')
pdf(out_name, width=5, height=4, paper="special" )

t_clv <- get_minsized_clusters(smash_subcl, min_size = 2)
for (i in 1:length(t_clv)){
  clustermmb <- {smash_subcl %>% filter(clustering == t_clv[i])}$sample
  if (length(clustermmb) > 0){
    print(plot_cluster_top15_genera(otu_taxa_metadata_top15_sum, "sample", opt$time, opt$factor, clustermmb, metadata_df))
  } else {
    print(ggplot(smash_subcl) + geom_col(aes(clustering)))
  }
}

dev.off()

