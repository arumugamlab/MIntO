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
library(tidyr)
library(stringr)
library(purrr)
library(optparse)
library(qs)

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

setDTthreads(threads = as.numeric(threads_n))
set.seed(1234)

# genome profiles to use as weights

metaG_genome_weights <- NULL
metaT_genome_weights <- NULL

if (normalization == 'MG') {
    if (omics %in% c('metaG', 'metaG_metaT')) {
        if (is.null(metaG_profile_file)) {
            stop("--genome-weights-metaG must be provided with --omics=metaG")
        } else {
            metaG_genome_weights <- (
                                     fread(metaG_profile_file, header=T)
                                     [, Funct := 'REL']
                                    )
            setnames(metaG_genome_weights, "ID", "MAG")
            setnames(metaG_genome_weights, colnames(metaG_genome_weights),
                     gsub(x = colnames(metaG_genome_weights), pattern = "metaG\\.", replacement = ""))
        }
    }

    if (omics %in% c('metaT', 'metaG_metaT')) {
        if (is.null(metaT_profile_file)) {
            stop("--genome-weights-metaT must be provided with --omics=metaT")
        } else {
            metaT_genome_weights <- (
                                     fread(metaT_profile_file, header=T)
                                     [, Funct := 'REL']
                                    )
            setnames(metaT_genome_weights, "ID", "MAG")
            setnames(metaT_genome_weights, colnames(metaT_genome_weights),
                     gsub(x = colnames(metaT_genome_weights), pattern = "metaT\\.", replacement = ""))
        }
    }
}



##########################  ** Load functions **  ##########################

plot_PCA <- function(data_phyloseq, color, label) {

  library(ggplot2)
  library(ggrepel)
  library(dplyr)

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
        annotations <- (fread(funcat_desc_file, header=T)
                        [, .(Funct, Description)]
                        [, Description := iconv(Description, from = "ISO-8859-1", to = "UTF-8")]
                       )
        annotations <- unique(merge(keyMap_funct, annotations, by='Funct', all=T)
                              [, Description := paste(sort(unique(Description)), collapse=';'), by=Funct]
                             )
    }
    else if  (db_name %in% c('PFAMs')){

        library(PFAM.db)

        # Get PFAM IDs
        pfam_names <- as.list(PFAMID[mappedkeys(PFAMID)])
        setDT(pfam_names, keep.rownames = TRUE)
        pfam_names <- data.table::transpose(pfam_names, keep.names = "Pfam_id")
        setnames(pfam_names, "V1", "Funct")

        # Get PFAM descriptions
        pfam_desc <- as.list(PFAMDE[mappedkeys(PFAMDE)])
        setDT(pfam_desc, keep.rownames = TRUE)
        pfam_desc <- data.table::transpose(pfam_desc, keep.names = "Pfam_id")
        setnames(pfam_desc, "V1", "Description")

        # Merge based on PFAM_ID and delete that join key
        pfam_table <- merge(pfam_names, pfam_desc, by='Pfam_id')[, Pfam_id := NULL]

        # Merge descriptions for the relevant PFAMs
        annotations <- unique(merge(keyMap_funct, pfam_table, by='Funct', all.x=T)
                                    [,Description := paste(sort(unique(Description)), collapse=';'), by='Funct']
                             )
        rm(pfam_names, pfam_desc, pfam_table, x, xx, keyMap_funct_desc,keyMap_funct_desc2)
    }
    else if  (db_name %in% c('dbCAN.module', 'dbCAN.enzclass', 'CAZy')){

        library(PFAM.db)
        library(dplyr)

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
          test_cazy.keyMap <- data.table(
            Pfam_id = rep(test_cazy$Pfam_id, sapply(test_cazy.singleKeys, length)),
            Funct = unlist(test_cazy.singleKeys)
          )
        }

        # Get PFAM descriptions
        pfam_desc <- as.list(PFAMDE[mappedkeys(PFAMDE)])
        setDT(pfam_desc, keep.rownames = TRUE)
        pfam_desc <- data.table::transpose(pfam_desc, keep.names = "Pfam_id")
        setnames(pfam_desc, "V1", "Description")

        # Merge based on PFAM_ID and delete that join key
        cazy_table <- merge(test_cazy.keyMap, pfam_desc, by='Pfam_id')[, Pfam_id := NULL]

        # Merge descriptions for the relevant PFAMs
        annotations <- unique(merge(keyMap_funct, cazy_table, by='Funct', all.x=T)
                              [,Description := paste(sort(unique(Description)), collapse=';'), by='Funct']
                             )
    }
    else {
        annotations <- unique(keyMap_funct[, Description := '-'])
    }
    return(annotations)
}

make_profile_files <- function(keys, profile, file_label, database, annotations, metadata, weights) {

      # Make a list of sample-names
      sample_cols <- setdiff(colnames(profile), c('coord'))

      # Get abundance information by gene and select only relevant columns.
      message(format(Sys.time(), digits=0), " - total genes ", dim(profile)[1])
      counts <- (merge(keys, profile, by='coord', all=FALSE)
                 [, MAG := stringr::str_split(coord, "\\|") %>% map_chr(., 1)] # Get the first field
                 [, coord := NULL]
                 [Reduce(`+`, mget(sample_cols)) > 0]
                )
      message(format(Sys.time(), digits=0), " - clean genes ", dim(counts)[1])

      # From count-per-gene, make count-per-function
      # If necessary, weight the functions by REL of genomes
      if (!is.null(weights)) {

          # Summarize by MAG,Funct
          # And prepend with weights
          counts <- rbind(weights,
                          counts[, lapply(.SD, sum, na.rm=TRUE), by = c('MAG', 'Funct')]
                         )
          message(format(Sys.time(), digits=0), " - MAGxfunctions ", dim(counts)[1])

          # weigh by the given weights per MAG
          weight_MAG_by_proportions <- function(x){
            if (is.character(x[1])){
              return(x)
            }else{
              x = x*x[1]
              return(x)
            }
          }

          # Normalize by weighting by 'weights'
          # Remove REL rows
          # Summarize by Funct
          counts <- (
                      counts
                      [, lapply(.SD, weight_MAG_by_proportions), by = MAG]
                      [! Funct == 'REL']
                      [, MAG := NULL]
                      [, lapply(.SD, sum, na.rm=TRUE), .SDcols = sample_cols, by = Funct]
                    )
      } else {

          # Summarize by Funct
          counts <- (
                      counts
                      [, MAG := NULL]
                      [, lapply(.SD, sum, na.rm=TRUE), .SDcols = sample_cols, by = Funct]
                    )
      }

      message(format(Sys.time(), digits=0), " - functions ", dim(counts)[1])

      # Add annotations
      counts <- merge(annotations, counts, by='Funct', all.x=FALSE, all.y=TRUE)

      # Write annotated functional profile as tsv
      fwrite(counts,
             file=paste0(output_dir, '/', file_label, '.', database ,'.tsv'),
             sep='\t',
             row.names = F,
             quote = F)

      #### Create phyloseq object ####

      my_otu_table <- as.matrix(counts[, -c('Funct', 'Description')])
      rownames(my_otu_table) <- counts$Funct

      ## metadata
      samp <- sample_data(metadata_df)
      rownames(samp) <- samp[["sample_alias"]]

      # taxa
      my_tax_table <- as.matrix(counts[, .(Funct, Description)])
      rownames(my_tax_table) <- counts$Funct

      physeq <- phyloseq(otu_table(my_otu_table, taxa_are_rows = T),
                         tax_table(my_tax_table),
                         samp)

      # Write phyloseq object
      qsave(physeq,
            preset = "balanced",
            file = paste0(phyloseq_dir, '/', file_label, '.', database ,'.qs'))

      return(physeq)
}

process_phyloseq_obj <- function(physeq_file) {
    physeq <- qread(physeq_file, nthreads=threads_n)
    counts <- (
                data.table(unclass(otu_table(physeq)), keep.rownames=TRUE)
                [, lapply(.SD, function(x) ifelse(is.na(x), 0, x))]
              )
    setnames(counts, "rn", "coord")

    taxa <- data.table(unclass(tax_table(physeq)), keep.rownames=TRUE)
    setnames(taxa, "rn", "coord")

    metadata <- (
                  data.table(as.data.frame(unclass(sample_data(physeq))), keep.rownames=TRUE)
                  [, rn := sample_alias]
                )

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
ge_gene_names <- NULL
metadata_df <- NULL
gene_annotation <- NULL
ga_fa_df <- NULL
gt_ft_df <- NULL
ge_fe_df <- NULL

# GE
if (omics == 'metaG_metaT'){
    message(format(Sys.time(), digits=0), " - began reading phyloseq")
    res <- process_phyloseq_obj(paste0(phyloseq_dir, '/GE.qs'))
    message(format(Sys.time(), digits=0), " - done  reading phyloseq")

    ge_fe_df <- data.frame(DB='genes', feature_n=dim(res$taxa)[1], feature = 'Genes')

    ge_gene_names <- res$counts$coord
    gene_annotation <- res$taxa

    # metadata
    metadata_df <- res$metadata

    rm(res)
}

# metaG
if (omics %like% "metaG") {
    message(format(Sys.time(), digits=0), " - began reading phyloseq")
    res <- process_phyloseq_obj(paste0(phyloseq_dir, '/GA.qs'))
    message(format(Sys.time(), digits=0), " - done  reading phyloseq")

    ga_fa_df <- data.frame(DB='genes', feature_n=dim(res$taxa)[1], feature = 'Genes')

    if (omics == 'metaG') {
        metaG_norm_profile <- res$counts
        gene_annotation <- res$taxa[coord %in% metaG_norm_profile$coord, ]
        metadata_df <- res$metadata
    } else {
        metaG_norm_profile <- res$counts[coord %in% ge_gene_names, ]
    }

    rm(res)
}

# metaT
if (omics %like% "metaT") {
    message(format(Sys.time(), digits=0), " - began reading phyloseq")
    res <- process_phyloseq_obj(paste0(phyloseq_dir, '/GT.qs'))
    message(format(Sys.time(), digits=0), " - done  reading phyloseq")

    gt_ft_df <- data.frame(DB='genes', feature_n=dim(res$taxa)[1], feature = 'Genes')

    if (omics == 'metaT') {
        metaT_norm_profile <- res$counts
        gene_annotation <- res$taxa[coord %in% metaT_norm_profile$coord, ]
        metadata_df <- res$metadata
    } else {
        metaT_norm_profile <- res$counts[coord %in% ge_gene_names, ]
    }

    rm(res)
}

# Delete data
rm(ge_gene_names)

# By now, we will have:
#   metadata_df
#   gene_annotation
#   ge_fe_df


#################################### FUNCTION EXPRESSION  ####################################
print('#################################### FUNCTION EXPRESSION  ####################################')

# Filter tax.table based on features in metaG, if metaG_metaT. metaG/metaT have been filtered earlier.

if (omics == 'metaG_metaT') {
    gene_annotation <- gene_annotation[coord %in% metaG_norm_profile$coord, ]
}


# Get entries for this funcat_name
gene_annotation <- (gene_annotation
          [, c('coord', funcat_name), with = FALSE]
          [!get(funcat_name) %in% c(NA, "", "-")]
         )

if (dim(gene_annotation)[1] > 0) {

    message(format(Sys.time(), digits=0), " - ", funcat_name, " started")

    ################################################
    # Prepare annotations and their descriptions
    ################################################

    # Replicate the gene annotation table by having N rows for gene with N annotations.

    singleKeys <- strsplit(gene_annotation[[funcat_name]], split = "\\; |\\;|\\, |\\,")
    keyMap <- unique(data.table(coord = rep(gene_annotation$coord, sapply(singleKeys, length)),
                                Funct = unlist(singleKeys),
                                key = "coord"))


    # Since KEGG annotations can have map12345 and ko12345 for the same pathway,
    # remove such redundancy.
    if (funcat_name %in% c("KEGG_Pathway", "kofam_Pathway")) {
        keyMap <- unique(keyMap[, Funct := gsub('ko', 'map', Funct)])
    }

    keyMap_funct <- data.table(Funct = c(unique(keyMap$Funct)))

    # Get annotation descriptions from external sources

    annot_dt <- get_annotation_descriptions(funcat_name)


    #####################################################
    # Make profiles and write them to 'output' directory
    #####################################################

    # metaG or metaG_metaT
    if (omics %like% "metaG") {
        metaG_physeq <- make_profile_files(keys=keyMap,
                                           profile=metaG_norm_profile,
                                           file_label="FA",
                                           database=funcat_name,
                                           annotations=annot_dt,
                                           metadata=metadata_df,
                                           weights=metaG_genome_weights)
    }

    # metaT or metaG_metaT
    if (omics %like% "metaT") {
        metaT_physeq <- make_profile_files(keys=keyMap,
                                           profile=metaT_norm_profile,
                                           file_label="FT",
                                           database=funcat_name,
                                           annotations=annot_dt,
                                           metadata=metadata_df,
                                           weights=metaT_genome_weights)
    }

    # metaG_metaT
    if (omics == 'metaG_metaT') {

        library(dplyr)

        # Expression -
        function_ids <- taxa_names(metaG_physeq)
        metaG_values <- as.data.frame(otu_table(metaG_physeq))
        metaT_values <- as.data.frame(otu_table(metaT_physeq))
        feature_ids <- intersect(rownames(metaG_values), rownames(metaT_values))
        metaG_values <- metaG_values[feature_ids,]
        metaT_values <- metaT_values[feature_ids,]

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

        function_expression_norm_desc <- right_join(as.data.frame(annot_dt), function_expression_norm, by='Funct')
        fwrite(function_expression_norm_desc, file=paste0(output_dir, '/FE.', funcat_name ,'.tsv'), sep='\t', row.names = F, quote = F)

        # Delete data
        rm(gene_annotation, keyMap_funct, singleKeys, metaG_values, metaT_values, feature_ids)

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
        qsave(FE_physeq,
              preset = "balanced",
              file = paste0(phyloseq_dir, '/FE.', funcat_name ,'.qs'))
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
        n_row <- nrow(data.table(unclass(tax_table(x)))[Funct != 'Unknown',])
        ga_fa_df <- rbind(ga_fa_df, c(funcat_name, n_row, 'Functions'))
    }
    if (!is.null(gt_ft_df)) {
        x <- metaT_physeq
        n_row <- nrow(data.table(unclass(tax_table(x)))[Funct != 'Unknown',])
        gt_ft_df <- rbind(gt_ft_df, c(funcat_name, n_row, 'Functions'))
    }
    if (!is.null(ge_fe_df)) {
        x <- FE_physeq
        n_row <- nrow(data.table(unclass(tax_table(x)))[Funct != 'Unknown',])
        ge_fe_df <- rbind(ge_fe_df, c(funcat_name, n_row, 'Functions'))
    }

    message(format(Sys.time(), digits=0), " - ", funcat_name, " finished!")

} else {
    # create .tsv, .qs, .pdf
    file.create(paste0(output_dir, '/FA.', funcat_name ,'.tsv'))
    file.create(paste0(output_dir, '/FT.', funcat_name ,'.tsv'))
    file.create(paste0(output_dir, '/FE.', funcat_name ,'.tsv'))
    file.create(paste0(phyloseq_dir, 'FA.', funcat_name ,'.qs'))
    file.create(paste0(phyloseq_dir, 'FT.', funcat_name ,'.qs'))
    file.create(paste0(phyloseq_dir, 'FE.', funcat_name ,'.qs'))
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
