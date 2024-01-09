#!/usr/bin/env Rscript

# '''
# Generate function expression profile from genome-based mode gene profiles (MG normalized)
#
# Authors: Carmen Saenz, Mani Arumugam
#
# '''


##########################  ** Load libraries **  ##########################
library(data.table)
library(phyloseq)
library(ggplot2)
library(ggrepel)
library(PFAM.db)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(optparse)

# Parse command line arguments
opt_list <- list(
                make_option("--threads", type="integer", default=4, help="number of threads [default: %default]"),
                make_option("--outdir", type="character", default=NULL, help="output directory to write normalized counts", metavar="directory"),
                make_option("--omics", type="character", default=NULL, help="which omics to summarize: metaG, metaT or metaG_metaT", metavar="string"),
                make_option("--normalization", type="character", default=NULL, help="normalization type: MG or TPM"),
                make_option("--funcat-name", type="character", default=NULL, help="name of functional category to summarize: e.g., eggNOG_OGs, dbCAN.EC", metavar="string"),
                make_option("--funcat-desc", type="character", default=NULL, help="file containing descriptions of functions", metavar="file"),
                make_option("--genome-weights-metaG", type="character", default=NULL, help="file with metaG profiles", metavar="file"),
                make_option("--genome-weights-metaT", type="character", default=NULL, help="file with metaT profiles", metavar="file"),
                make_option("--main-factor", type="character", default=NULL, help="main factor variable to use when plotting data", metavar="string")
                )
opts <- parse_args(OptionParser(option_list=opt_list))


##########################  ** Load arguments **  ##########################
threads_n <- opts$threads
output_dir <- opts$outdir
omics <- opts$omics #'metaG_metaT'
normalization <- opts$normalization
funcat_name  <- opts[['funcat-name']]
funcat_desc_file <- opts[['funcat-desc']]
main_factor <- opts[['main-factor']]
metaG_profile_file <- opts[['genome-weights-metaG']]
metaT_profile_file <- opts[['genome-weights-metaT']]

# genome profiles to use as weights

metaG_genome_weights <- NULL
metaT_genome_weights <- NULL

if (normalization == 'MG') {
    if (omics %in% c('metaG', 'metaG_metaT')) {
        if (is.null(metaG_profile_file)) {
            stop("--genome-weights-metaG must be provided with --omics=metaG")
        } else {
            metaG_genome_weights <- fread(metaG_profile_file, header=T) %>%
                                as.data.frame(stringsAsFactors = F) %>%
                                dplyr::rename(MAG = ID) %>%
                                mutate(Funct = "REL")
            colnames(metaG_genome_weights) <- gsub(x = colnames(metaG_genome_weights), pattern = "metaG\\.", replacement = "")
        }
    }

    if (omics %in% c('metaT', 'metaG_metaT')) {
        if (is.null(metaT_profile_file)) {
            stop("--genome-weights-metaT must be provided with --omics=metaT")
        } else {
            metaT_genome_weights <- fread(metaT_profile_file, header=T) %>%
                                as.data.frame(stringsAsFactors = F) %>%
                                dplyr::rename(MAG = ID) %>%
                                mutate(Funct = "REL")
            colnames(metaT_genome_weights) <- gsub(x = colnames(metaT_genome_weights), pattern = "metaT\\.", replacement = "")
        }
    }
}


setDTthreads(threads = as.numeric(threads_n))
set.seed(1234)

##########################  ** Load functions **  ##########################

plot_PCA <- function(data_phyloseq, color, label) {

  library(phyloseq)
  library(ggplot2)

  ##### metadata:
  sample_data_df <- data.frame(sample_data(data_phyloseq), stringsAsFactors = F)

  ### Counts - non-zero/inf
  ## Prepare data - replace NA values by 0
  otu_table_df <-  unclass(otu_table(data_phyloseq)) %>%
                      as.data.frame(stringAsFactors = F) %>%
                      tibble::rownames_to_column('rownames') %>%
                      dplyr::mutate_all(function(x) ifelse(is.infinite(x) | is.na(x), 0, x)) %>%
                      tibble::column_to_rownames('rownames')

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
  dtp <- data.frame('sample' = sample_data_df[[label]],
                    'group' = sample_data_df[[color]],
                    gene_expression_pca$x[,1:2]) # the first two componets are selected (NB: you can also select 3 for 3D plottings or 3+)
  df_out <- as.data.frame(gene_expression_pca$x)
  percentage <- round(gene_expression_pca$sdev / sum(gene_expression_pca$sdev) * 100, 2)
  percentage <- paste0( colnames(df_out), " (", paste0(as.character(percentage), "%", ")"))

  PCA_Sample_site_abundance <- (ggplot(data = dtp, aes(x = PC1, y = PC2, color = group)) +
                                  geom_text_repel(aes(label = sample),nudge_x = 0.04, size = 3.5, segment.alpha = 0.5) +
                                  geom_point(size = 2, shape = 16)+
                                  xlab(percentage[1]) + ylab(percentage[2]) +
                                  #labs(title = title) +
                                  theme_bw()+
                                  theme(plot.title = element_text(size=10), legend.position="bottom"))#+ stat_ellipse(type = "norm", linetype = 2))
  return(PCA_Sample_site_abundance)
}

prepare_PCA <- function(physeq, type, database, color) {
    # Colors
    manual_plot_colors =c('#9D0208', '#264653','#e9c46a','#D8DCDE','#B6D0E0',
                          '#FFC87E','#F4A261','#E34F33','#E9C46A',
                          '#A786C9','#D4C0E2','#975773','#6699FF','#000066',
                          '#7AAFCA','#006699','#A9D181','#2F8475','#264445')

    # Title
    title_name <- paste0('PCA - function ', type, ' - ', database)

    # File names
    # abundance  -> FA
    # transcript -> FT
    # expression -> FE
    out_name <- paste0(visual_dir, 'F', str_sub(toupper(type), 1, 1), '.', database ,'.PCA.pdf')

    plot_PCA_out <- plot_PCA(physeq, color=color, label="sample_alias")
    pdf(out_name,width=8,height=8,paper="special" )
    print(plot_PCA_out  + scale_color_manual(values=manual_plot_colors, name=color) +
          #coord_fixed() +
          theme(legend.position="bottom")+ggtitle(title_name))
    print(plot_PCA_out  +
          facet_wrap(. ~ group)+
          scale_color_manual(values=manual_plot_colors, name=color) +
          #coord_fixed() +
          theme(legend.position="bottom")+
          ggtitle(title_name))
    dev.off()
}


get_annotation_descriptions <- function(db_name) {
    ### Annotation descriptions ####
    # library(purrr)
    # library(magrittr)
    if  (db_name %in% c('eggNOG_OGs', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_KO', 'kofam_Pathway', 'kofam_Module', 'kofam_KO', 'merged_KO', 'dbCAN.EC')){
        file_name <- db_name %>%
                            gsub('kofam_', 'KEGG_', .) %>%
                            gsub('merged_', 'KEGG_', .)
        def_list <- as.data.frame(fread(funcat_desc_file, header=T), stringsAsFactors = F)
        def_list_sub <- subset(def_list, select=c('Funct', 'Description'))
        def_list_sub$Description <- iconv(def_list_sub$Description, from = "ISO-8859-1", to = "UTF-8")
        keyMap_funct_desc <- merge(keyMap_funct, def_list_sub, by='Funct', all=T)
        keyMap_funct_desc2 <- as.data.frame(keyMap_funct_desc %>%
                                              dplyr::group_by(Funct)%>%
                                              dplyr::summarise(Description = paste(sort(unique(Description)), collapse=';')))
        annot_df <- keyMap_funct_desc2[!duplicated(keyMap_funct_desc2),]
        rm(keyMap_funct_desc,keyMap_funct_desc2)
    }
    else if  (db_name %in% c('PFAMs')){
        x <- PFAMID
        mapped_keys <- mappedkeys(x)
        xx <- as.list(x[mapped_keys])
        test <- as.data.frame(do.call(rbind, xx))
        names(test) <- 'Pfam_family'
        test$Pfam_id <- rownames(test)
        x <- PFAMDE
        mapped_keys <- mappedkeys(x)
        xx <- as.list(x[mapped_keys])
        test2 <- as.data.frame(do.call(rbind, xx))
        names(test2) <- 'Description'
        test2$Pfam_id <- rownames(test2)
        test3 <- merge(test, test2, by='Pfam_id')
        keyMap_funct_desc <- merge(keyMap_funct, test3, by.x='Funct' , by.y='Pfam_family', all.x=T)
        keyMap_funct_desc$Pfam_id <- NULL
        keyMap_funct_desc2 <- as.data.frame(keyMap_funct_desc %>%
                                              dplyr::group_by(Funct)%>%
                                              dplyr::summarise(Description = paste(sort(unique(Description)), collapse=';')))
        annot_df <- keyMap_funct_desc2[!duplicated(keyMap_funct_desc2),]
        rm(test, test2, test3, x, xx, keyMap_funct_desc,keyMap_funct_desc2)
    }
    else if  (db_name %in% c('dbCAN.module', 'dbCAN.enzclass', 'CAZy')){
        x <- PFAMCAZY
        mapped_keys <- mappedkeys(x)
        xx <- as.list(x[mapped_keys])
        test <- as.data.frame(do.call(rbind, xx))
        #data.frame(matrix((xx), nrow=length(xx), byrow=F))
        #names(test) <- 'CAZY_id'
        test$Pfam_id <- rownames(test)
        test_cazy <- as.data.frame(test %>% group_by(Pfam_id) %>%
          mutate(CAZY_id = paste(unique(c(V1, V2, V3, V4)), collapse = ',')))
        test_cazy$V1 <- test_cazy$V2 <- test_cazy$V3 <- test_cazy$V4 <- NULL
        test_cazy.singleKeys <- strsplit(test_cazy$CAZY_id, split = "\\  |\\,")
        if (length(test_cazy.singleKeys) > 0){
          test_cazy.keyMap <- data.frame(
            Pfam_id = rep(test_cazy$Pfam_id, sapply(test_cazy.singleKeys, length)),
            Funct = unlist(test_cazy.singleKeys),
            stringsAsFactors = FALSE
          )
        }
        x <- PFAMDE
        mapped_keys <- mappedkeys(x)
        xx <- as.list(x[mapped_keys])
        test2 <- as.data.frame(do.call(rbind, xx))
        names(test2) <- 'Description'
        test2$Pfam_id <- rownames(test2)
        test3 <- merge(test_cazy.keyMap, test2, by='Pfam_id')
        keyMap_funct_desc <- merge(keyMap_funct, test3, by='Funct', all.x=T)
        keyMap_funct_desc$Pfam_id <- NULL
        keyMap_funct_desc2 <- as.data.frame(keyMap_funct_desc %>%
                                              dplyr::group_by(Funct)%>%
                                              dplyr::summarise(Description = paste(sort(unique(Description)), collapse=';')))
        annot_df <- keyMap_funct_desc2[!duplicated(keyMap_funct_desc2),]
        rm(test, test2, test3, test_cazy, test_cazy.keyMap, test_cazy.singleKeys, x, xx, keyMap_funct_desc,keyMap_funct_desc2)
    }
    else {
        keyMap_funct_desc <- keyMap_funct
        keyMap_funct_desc$Description <- '-'
        keyMap_funct_desc2 <- as.data.frame(keyMap_funct_desc %>%
                                              dplyr::group_by(Funct)%>%
                                              dplyr::summarise(Description = paste(sort(unique(Description)), collapse=';')))
        annot_df <- keyMap_funct_desc2[!duplicated(keyMap_funct_desc2),]
        rm(keyMap_funct_desc,keyMap_funct_desc2)
    }
    #annot_df$Description <- gsub(",", ";", annot_df$Description)
    return(annot_df)
}

make_profile_files_old <- function(keys, profile, file_label, database, annotations, metadata) {
      counts_by_gene <- merge(keys, profile, by='coord', all=T)
      counts_by_gene$coord <- NULL
      counts_by_gene$ID_gene <- NULL
      counts_by_gene$name <- NULL

      # From count-per-gene, make count-per-function

      counts_by_func <- counts_by_gene %>%
                             replace_na(list(Funct = 'Unknown')) %>%
                             mutate(Funct = ifelse(Funct == '-' | Funct == '', 'Unknown', Funct)) %>%
                             dplyr::group_by(Funct) %>%
                             dplyr::summarise(across(everything(),sum)) %>%
                             as.data.frame(stringsAsFactors = F) %>%
                             dplyr::select(Funct, everything())
      names(counts_by_func) <- gsub('X', '', names(counts_by_func))

      counts_by_func_annotated <- merge(annotations, counts_by_func, by='Funct', all.y=T)
      fwrite(counts_by_func_annotated, file=paste0(output_dir, '/', file_label, '.', database ,'.tsv'), sep='\t', row.names = F, quote = F)

      #### phyloseq object ####

      rownames(counts_by_func) <- counts_by_func$Funct
      counts_by_func$Funct <- NULL

      ## metadata
      samp <- sample_data(metadata_df)
      rownames(samp) <- samp[["sample_alias"]]

      # taxa
      func_desc <- subset(counts_by_func_annotated, select=c("Funct", "Description"))
      rownames(func_desc) <- func_desc$Funct

      physeq <- phyloseq(otu_table(as.matrix(counts_by_func), taxa_are_rows = T),
                                                     tax_table(as.matrix(func_desc)), samp)
      saveRDS(physeq, file = paste0(phyloseq_dir, '/', file_label, '.', database ,'.rds'))

      return(physeq)
}

make_profile_files <- function(keys, profile, file_label, database, annotations, metadata, weights) {

      # Get abundance information by gene and select only relevant columns.
      counts_by_gene <- inner_join(keys, profile, by='coord') %>%
                            mutate(MAG = stringr::str_split(coord, '\\|') %>% map_chr(.,1)) %>%
                            dplyr::select(-coord)

      # From count-per-gene, make count-per-function
      # If necessary, weight the functions by REL of genomes
      if (!is.null(weights)) {

          # Summarize by MAG,Funct
          counts_by_func <- counts_by_gene %>%
                                dplyr::select(MAG, Funct, everything()) %>%
                                bind_rows(., weights) %>%
                                dplyr::summarise(across(everything(), sum),
                                                 .by = c(MAG, Funct))

          # melt to make one gene rpk per row for operational convenience
          cf_melt <- reshape2::melt(counts_by_func, id.vars = c("MAG", "Funct"))

          # Normalize by weighting by 'weights'
          # Remove REL rows
          # Replace NaN and infinite values by 0
          cf_melt_normalized <- cf_melt %>%
            mutate(value = value * value[Funct == 'REL'],
                   .by = c(MAG, variable)) %>%
            filter(Funct != 'REL') %>%
            dplyr::select(-MAG) %>%
            dplyr::summarise(value = sum(value),
                             .by = c(Funct, variable))


          # unmelt to bring back table
          counts_by_func <- reshape2::dcast(cf_melt_normalized, Funct ~ variable,  value.var="value")
          rm(cf_melt, cf_melt_normalized)
          gc()

      } else {

          # Summarize by Funct
          counts_by_func <- counts_by_gene %>%
                                dplyr::select(-any_of(c("MAG"))) %>%         # remove MAG if exists
                                dplyr::select(Funct, everything()) %>%
                                dplyr::summarise(across(everything(), sum),
                                                 .by = c(Funct))
      }

      # Clean up memory
      rm(counts_by_gene)
      gc()

      # Add annotations
      counts_by_func_annotated <- right_join(annotations, counts_by_func, by='Funct')

      # Write annotated functional profile as tsv
      fwrite(counts_by_func_annotated %>%
                filter(rowSums(across(-c(Funct, Description))) > 0),
             file=paste0(output_dir, '/', file_label, '.', database ,'.tsv'),
             sep='\t',
             row.names = F,
             quote = F)

      #### Create phyloseq object ####

      rownames(counts_by_func) <- counts_by_func$Funct
      counts_by_func$Funct <- NULL

      ## metadata
      samp <- sample_data(metadata_df)
      rownames(samp) <- samp[["sample_alias"]]

      # taxa
      func_desc <- subset(counts_by_func_annotated, select=c("Funct", "Description"))
      rownames(func_desc) <- func_desc$Funct

      physeq <- phyloseq(otu_table(as.matrix(counts_by_func), taxa_are_rows = T),
                                                     tax_table(as.matrix(func_desc)), samp)

      # Write phyloseq object
      saveRDS(physeq, file = paste0(phyloseq_dir, '/', file_label, '.', database ,'.rds'))

      return(physeq)
}

process_phyloseq_obj <- function(omics_label) {
    physeq <- readRDS(paste0(phyloseq_dir, '/', omics_label, '.rds'))
    counts <- unclass(otu_table(physeq)) %>%
                                   as.data.frame(stringsAsFactors = F) %>%
                                   tibble::rownames_to_column('rownames') %>%
                                   dplyr::mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
                                   tibble::column_to_rownames('rownames') %>%
                                   as.data.frame(stringsAsFactors = F)

    taxa <- as.data.frame(unclass(tax_table(physeq)), stringsAsFactors = F)
    metadata <- as.data.frame(unclass(sample_data(physeq)), stringsAsFactors = F)
    rownames(metadata) <- metadata[["sample_alias"]]

    return(list(counts=counts, taxa=taxa, metadata=metadata))
}

# *************************** Use GA and GT raw counts to generate FE profile *************************** ####
# The genes will be clustered by function IDs.
# TPM normalization - sum all TPMs for the member genes
#                     (there will be redundancy and total might be over 1 million)
# MG normalization  - weight individual member genes by the rel. abundance of the genome and sum.
# FE profile will be generated: (FT_norm) / (FA_norm)
#Raw gene abundances #####

print('#################################### Paths ####################################')

# Generate output directories ####
dir.create(file.path(output_dir), showWarnings = FALSE)
visual_dir=paste0(output_dir,'/plots/')
dir.create(file.path(visual_dir), showWarnings = FALSE)
phyloseq_dir=paste0(output_dir,'/phyloseq_obj/')
dir.create(file.path(phyloseq_dir), showWarnings = FALSE)


print('#################################### DB profiles ####################################')

# Initialize some variables
ge_norm_taxa_genes <- NULL
metadata_df <- NULL
profiles_annot_norm_df_sub <- NULL
ga_fa_df <- NULL
gt_ft_df <- NULL
ge_fe_df <- NULL

# GE
if (omics == 'metaG_metaT'){
    res <- process_phyloseq_obj("GE")
    ge_norm_count <- res$counts
    ge_norm_taxa <- res$taxa
    ge_norm_taxa_genes <- rownames(ge_norm_taxa)
    profiles_annot_norm_df_sub <- ge_norm_taxa

    # metadata
    metadata_df <- res$metadata

    ge_fe_df <- data.frame(DB='genes', feature_n=length(ge_norm_taxa_genes), feature = 'Genes')

    rm(ge_norm_count, ge_norm_taxa)
}

# metaG
if (omics %like% "metaG") {
    res <- process_phyloseq_obj("GA")
    ga_norm_count <- res$counts
    ga_norm_taxa <- res$taxa
    ga_norm_taxa_genes <- rownames(ga_norm_taxa)
    ga_fa_df <- data.frame(DB='genes', feature_n=length(ga_norm_taxa_genes), feature = 'Genes')

    if (omics == 'metaG') {
        metaG_norm_profile <- ga_norm_count
        profiles_annot_norm_df_sub <- ga_norm_taxa
        metadata_df <- res$metadata
    } else {
        metaG_norm_profile <- ga_norm_count[rownames(ga_norm_count) %in% ge_norm_taxa_genes,]
    }
    metaG_norm_profile <- metaG_norm_profile %>%
                                            tibble::rownames_to_column('coord') %>%
                                            dplyr::select(coord, everything())
    rm(ga_norm_count, ga_norm_taxa, ga_norm_taxa_genes)
}

# metaT
if (omics %like% "metaT") {
    res <- process_phyloseq_obj("GT")
    gt_norm_count <- res$counts
    gt_norm_taxa <- res$taxa
    gt_norm_taxa_genes <- rownames(gt_norm_taxa)
    gt_ft_df <- data.frame(DB='genes', feature_n=length(gt_norm_taxa_genes), feature = 'Genes')

    if (omics == 'metaT') {
        metaT_norm_profile <- gt_norm_count
        profiles_annot_norm_df_sub <- gt_norm_taxa
        metadata_df <- res$metadata
    } else {
        metaT_norm_profile <- gt_norm_count[rownames(gt_norm_count) %in% ge_norm_taxa_genes,]
    }
    metaT_norm_profile <- metaT_norm_profile %>%
                                            tibble::rownames_to_column('coord') %>%
                                            dplyr::select(coord, everything())
    rm(gt_norm_count, gt_norm_taxa, gt_norm_taxa_genes)
}

# Delete data
rm(ge_norm_taxa_genes)

# By now, we will have:
#   metadata_df
#   profiles_annot_norm_df_sub
#   ge_fe_df


#################################### FUNCTION EXPRESSION  ####################################
print('#################################### FUNCTION EXPRESSION  ####################################')

# Filter profile based on non-zero features if GE. Otherwise, use themselves

if (omics == 'metaT') {
    filter_profile <- metaT_norm_profile
} else {
    filter_profile <- metaG_norm_profile
}
profiles_annot_norm_df_sub <- profiles_annot_norm_df_sub %>%
                                  mutate(coord = rownames(.)) %>%
                                  dplyr::select(coord, everything()) %>%
                                  filter(coord %in% filter_profile$coord)
rm(filter_profile)


# Get entries for this funcat_name
dbMap <- profiles_annot_norm_df_sub %>%
            dplyr::select(c(coord, all_of(funcat_name))) %>%
            filter(!is.na(!!as.symbol(funcat_name)) & !!as.symbol(funcat_name) != "" & !!as.symbol(funcat_name) != "-")

if (dim(dbMap)[1] > 0) {

    message(format(Sys.time(), digits=0), " - ", funcat_name, " started")

    ################################################
    # Prepare annotations and their descriptions
    ################################################

    # Replicate the gene annotation table by having N rows for gene with N annotations.

    singleKeys <- strsplit(dbMap[,2], split = "\\; |\\;|\\, |\\,")
    keyMap <- data.frame(coord = rep(dbMap$coord, sapply(singleKeys, length)),
                         Funct = unlist(singleKeys),
                         stringsAsFactors = FALSE) %>%
                 distinct()

    # Since KEGG annotations can have map12345 and ko12345 for the same pathway,
    # remove such redundancy.
    if (funcat_name == "KEGG_Pathway"){
        keyMap$Funct <- gsub('ko', 'map', keyMap$Funct)
        keyMap <- distinct(keyMap)
    }

    keyMap_funct <- data.frame(Funct = c(unique(keyMap$Funct)))

    # Get annotation descriptions from external sources

    annot_df <- get_annotation_descriptions(funcat_name)


    #####################################################
    # Make profiles and write them to 'output' directory
    #####################################################

    # metaG or metaG_metaT
    if (omics %like% "metaG") {
        metaG_physeq <- make_profile_files(keys=keyMap,
                                           profile=metaG_norm_profile,
                                           file_label="FA",
                                           database=funcat_name,
                                           annotations=annot_df,
                                           metadata=metadata_df,
                                           weights=metaG_genome_weights)
    }

    # metaT or metaG_metaT
    if (omics %like% "metaT") {
        metaT_physeq <- make_profile_files(keys=keyMap,
                                           profile=metaT_norm_profile,
                                           file_label="FT",
                                           database=funcat_name,
                                           annotations=annot_df,
                                           metadata=metadata_df,
                                           weights=metaT_genome_weights)
    }

    # metaG_metaT
    if (omics == 'metaG_metaT') {
        FE_physeq <- make_profile_files(keys=keyMap,
                                           profile=metaT_norm_profile,
                                           file_label="FT",
                                           database=funcat_name,
                                           annotations=annot_df,
                                           metadata=metadata_df,
                                           weights=metaT_genome_weights)
        # Expression -
        function_ids <- taxa_names(metaG_physeq)
        metaG_values <- as.data.frame(otu_table(metaG_physeq))
        metaT_values <- as.data.frame(otu_table(metaT_physeq))
        metaT_values <- metaT_values[match(rownames(metaT_values), rownames(metaT_values)),]

        # Function Expression
        function_expression_norm=(metaT_values)/(metaG_values)

        # Get rid of rows with Inf or 0-sum
        function_expression_norm <- as.data.frame(function_expression_norm %>%
                                                  tibble::rownames_to_column('rownames') %>%
                                                  dplyr::mutate_all(function(x) ifelse(is.infinite(x), NA, x)) %>%
                                                  drop_na() %>%
                                                  filter(rowSums(across(where(is.numeric), ~.x != 0))>0) %>%
                                                  mutate(Funct = rownames) %>%
                                                  tibble::column_to_rownames('rownames')) %>%
                                                  dplyr::select(Funct, everything())

        function_expression_norm_desc <- right_join(annot_df, function_expression_norm, by='Funct')
        fwrite(function_expression_norm_desc, file=paste0(output_dir, '/FE.', funcat_name ,'.tsv'), sep='\t', row.names = F, quote = F)

        # Delete data
        rm(dbMap, keyMap_funct, singleKeys,
         metaG_values, metaT_values)

        #### phyloseq object - metaT relative to metaG -####

        rownames(function_expression_norm) <- function_expression_norm$Funct
        function_expression_norm$Funct <- NULL
        #min(function_expression_norm[function_expression_norm >0 & !is.na(function_expression_norm)])

        ## metadata
        samp <- sample_data(metadata_df)
        rownames(samp) <- samp[["sample_alias"]]

        # taxa
        FE_annot_desc <- subset(function_expression_norm_desc, select=c("Funct","Description"))
        rownames(FE_annot_desc) <- FE_annot_desc$Funct

        FE_physeq <- phyloseq(otu_table(as.matrix(function_expression_norm), taxa_are_rows = T),
                                               tax_table(as.matrix(FE_annot_desc)),
                                               samp)
        saveRDS(FE_physeq, file = paste0(phyloseq_dir, '/FE.', funcat_name ,'.rds'))
    }


    ###########
    # Plot PCA
    ###########

    # PCA - metaG ####
    if (omics %like% 'metaG') {
        prepare_PCA(metaG_physeq, type="abundance", database=funcat_name, color=main_factor)
    }

    # PCA - metaT ####
    if (omics %like% 'metaT') {
        prepare_PCA(metaT_physeq, type="transcript", database=funcat_name, color=main_factor)
    }

    # PCA - GE ####
    if (omics == 'metaG_metaT') {
        prepare_PCA(FE_physeq, type="expression", database=funcat_name, color=main_factor)
    }

    # Create df to plot number of genes and functions ####
    if (!is.null(ga_fa_df)) {
        x <- metaG_physeq
        n_row <- unclass(tax_table(x)) %>%
                    as.data.frame(stringsAsFactors = F) %>%
                    dplyr::filter(Funct != 'Unknown') %>%
                    nrow()
        ga_fa_df <- rbind(ga_fa_df, c(funcat_name, n_row, 'Functions'))
    }
    if (!is.null(gt_ft_df)) {
        x <- metaT_physeq
        n_row <- unclass(tax_table(x)) %>%
                    as.data.frame(stringsAsFactors = F) %>%
                    dplyr::filter(Funct != 'Unknown') %>%
                    nrow()
        gt_ft_df <- rbind(gt_ft_df, c(funcat_name, n_row, 'Functions'))
    }
    if (!is.null(ge_fe_df)) {
        x <- FE_physeq
        n_row <- unclass(tax_table(x)) %>%
                    as.data.frame(stringsAsFactors = F) %>%
                    dplyr::filter(Funct != 'Unknown') %>%
                    nrow()
        ge_fe_df <- rbind(ge_fe_df, c(funcat_name, n_row, 'Functions'))
    }

    message(format(Sys.time(), digits=0), " - ", funcat_name, " finished!")

} else {
    # create .tsv, .rds, .pdf
    file.create(paste0(output_dir, '/FA.', funcat_name ,'.tsv'))
    file.create(paste0(output_dir, '/FT.', funcat_name ,'.tsv'))
    file.create(paste0(output_dir, '/FE.', funcat_name ,'.tsv'))
    file.create(paste0(phyloseq_dir, 'FA.', funcat_name ,'.rds'))
    file.create(paste0(phyloseq_dir, 'FT.', funcat_name ,'.rds'))
    file.create(paste0(phyloseq_dir, 'FE.', funcat_name ,'.rds'))
    file.create(paste0(visual_dir, 'FA.', funcat_name ,'.PCA.pdf'))
    file.create(paste0(visual_dir, 'FT.', funcat_name ,'.PCA.pdf'))
    file.create(paste0(visual_dir, 'FE.', funcat_name ,'.PCA.pdf'))
    ge_fe_df <- rbind(ge_fe_df, c(funcat_name, 0, 'Functions'))
}

# Free memory
gc()

#####################################
# Write out feature count statistics
#####################################

print('#################################### FEATURE COUNTS  ####################################')#####

# Write counts
if (!is.null(ga_fa_df)) {
    fwrite(ga_fa_df, file=paste0(output_dir, '/GA_FA_features.', funcat_name, '.tsv'), sep='\t', row.names = F, quote = F)
}

if (!is.null(gt_ft_df)) {
    fwrite(gt_ft_df, file=paste0(output_dir, '/GT_FT_features.', funcat_name, '.tsv'), sep='\t', row.names = F, quote = F)
}

if (!is.null(ge_fe_df)) {
    fwrite(ge_fe_df, file=paste0(output_dir, '/GE_FE_features.', funcat_name, '.tsv'), sep='\t', row.names = F, quote = F)
}

print('done!')
