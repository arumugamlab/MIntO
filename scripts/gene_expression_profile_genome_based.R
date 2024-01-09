#!/usr/bin/env Rscript

# '''
# Generate gene expression profile from genome-based mode gene profiles
#
# Authors: Carmen Saenz, Mani Arumugam
#
# '''

##########################  ** Load libraries **  ##########################
library(dplyr)
library(data.table)
library(phyloseq)
library(KEGGREST)
library(tibble)
library(stringr)
library(purrr)
library(optparse)

# Parse command line arguments
opt_list <- list(
                make_option("--threads", type="integer", default=4, help="number of threads [default: %default]"),
                make_option("--outdir", type="character", default=NULL, help="output directory to write normalized counts", metavar="directory"),
                make_option("--omics", type="character", default=NULL, help="which omics to summarize: metaG, metaT or metaG_metaT", metavar="string"),
                make_option("--funcat-names", type="character", default=NULL, help="comma-delimited list of functional categores to summarize: e.g., 'eggNOG_OGs,dbCAN.EC,KEGG_KO'", metavar="string"),
                make_option("--gene-profile", type="character", default=NULL, help="file containing gene profiles in metaG and/or metaT space", metavar="file"),
                make_option("--annotation", type="character", default=NULL, help="file with gene annotations", metavar="file"),
                make_option("--metadata", type="character", default=NULL, help="file with sample metadata", metavar="file"),
                make_option("--main-factor", type="character", default=NULL, help="main factor variable to use when plotting data", metavar="string")
                )
opts <- parse_args(OptionParser(option_list=opt_list))

threads_n <- opts$threads
output_dir <- opts$outdir
omics <- opts$omics #'metaG_metaT'
annot_file <- opts$annotation
metadata_file <- opts$metadata
input_file <- opts[['gene-profile']]
funcat_args  <- opts[['funcat-names']]
main_factor <- opts[['main-factor']]

funcat_names <- unlist(strsplit(funcat_args, split = "\\,")[[1]])
print(funcat_names)

setDTthreads(threads = threads_n)
set.seed(1234)

plot_PCA <- function(physeq, color, label) {

    library(phyloseq)
    library(ggplot2)
    library(ggrepel)

    ##### metadata:
    sample_data_df <- data.frame(sample_data(physeq), stringsAsFactors = F)

    ### Counts - non-zero/inf
    ## Prepare data - replace NA values by 0
    otu_table_df <-  unclass(otu_table(physeq)) %>%
                      as.data.frame(stringAsFactors = F) %>%
                      tibble::rownames_to_column('rownames') %>%
                      dplyr::mutate_all(function(x) ifelse(is.infinite(x) | is.na(x), 0, x)) %>%
                      filter(rowSums(across(where(is.numeric), ~.x != 0))>0) %>%
                      tibble::column_to_rownames('rownames')

    ##### Replace infinite values by log2(min value * 1e-2) for the PCA
    replace_val <- log((min(otu_table_df[otu_table_df >0 & !is.na(otu_table_df)]))*1e-2)

    ### Counts - Log transformed
    otu_table_df_log <- log(otu_table_df) %>%
                          tibble::rownames_to_column('ID') %>%
                          dplyr::mutate_all(function(x) ifelse(is.infinite(x), replace_val, x)) %>%
                          tibble::column_to_rownames('ID')

    #min(otu_table_df_log[!is.na(otu_table_df_log)])

    #Plotting scores of PC1 and PC2 with log transformation
    gene_expression_pca <- prcomp(t(otu_table_df_log), center = T, sca=T)

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

prepare_PCA <- function(physeq, type, color) {

    library(stringr)

    # Colors
    manual_plot_colors =c('#9D0208', '#264653','#e9c46a','#D8DCDE','#B6D0E0',
                          '#FFC87E','#F4A261','#E34F33','#E9C46A',
                          '#A786C9','#D4C0E2','#975773','#6699FF','#000066',
                          '#7AAFCA','#006699','#A9D181','#2F8475','#264445')

    # Title
    title_name <- paste0('PCA - gene ', type)

    # File names
    # abundance  -> GA
    # transcript -> GT
    # expression -> GE
    out_name <- paste0(visual_dir, 'G', str_sub(toupper(type), 1, 1), '.PCA.pdf')

    plot_PCA_out <- plot_PCA(physeq=physeq, color=color, label="sample_alias")
    pdf(out_name,width=8,height=8,paper="special" )
    print(plot_PCA_out  + scale_color_manual(values=manual_plot_colors, name=color) +
          #coord_fixed() +
          theme(legend.position="bottom")+ggtitle(title_name))
    print(plot_PCA_out  +
          facet_wrap(.~group)+
          scale_color_manual(values=manual_plot_colors, name=color) +
          #coord_fixed() +
          theme(legend.position="bottom")+
          ggtitle(title_name))
    dev.off()
}

#####
# Gene abundance profile
#####

profiles_tpm <- input_file

print('#################################### Paths ####################################')

# Generate output directories ####

dir.create(file.path(output_dir), showWarnings = FALSE)
visual_dir=paste0(output_dir, '/plots/')
dir.create(file.path(visual_dir), showWarnings = FALSE)
phyloseq_dir=paste0(output_dir, '/phyloseq_obj/')
dir.create(file.path(phyloseq_dir), showWarnings = FALSE)

print('#################################### gene profiles ####################################')

########################################################
## Read input files
########################################################

# Get gene annotation
# Replace , by ; in annotation IDs
gene_annot_df <- as.data.frame(fread(annot_file,header=T), stringsAsFactors = F, row.names = T) %>%
                        mutate(across(-ID, ~ gsub("\\,", "\\;", .))) %>%
                        select(ID, all_of(unique(funcat_names))) %>%
                        add_row(ID = 'Unknown')

# Read the TPM/MG normalized profile csv file for this omics
gene_profile_df <- fread(profiles_tpm, header=T) %>%
                    as.data.frame(stringsAsFactors = F)

########################################################
## Get abundance and gene information into separate DFs
########################################################


bed_colnames <- c("chr","start","stop","name","score","strand","source","feature","frame","info")
if ( sum(bed_colnames %in% colnames(gene_profile_df)) == length(bed_colnames)) { # MAG-genes or refgenome-genes mode

    # Get useful non-abundance information columns (BED columns) in a separate DF
    # Modify some fields in the info DF
    # Keep also coord, which will be used to merge info and abundance later
    bed_info_only_df <- gene_profile_df %>%
                            mutate(info = ifelse(name=='.', coord, info)) %>%
                            select(coord, all_of(bed_colnames)) %>%
                            mutate(ID = info) %>%
                            mutate(name    = stringr::str_split(ID, '\\|') %>% map_chr(., 2)) # Get 2nd field

    # Get abundance information in a separate DF
    # And remove all rows with 0 sum
    abundance_only_df <- gene_profile_df %>%
                            select(-gene_length, -all_of(bed_colnames)) %>%
                            distinct(coord, .keep_all = TRUE) %>%
                            mutate(across(-coord, ~ as.numeric(replace(., . == '' | is.na(.), 0)))) %>%
                            tibble::column_to_rownames('coord') %>%
                            filter(rowSums(across(where(is.numeric), ~.x != 0))>0)

    # Merge bed-info DF and annotation DF
    names_list <- names(gene_annot_df)
    gene_info_annotation_df <- merge(gene_annot_df, bed_info_only_df, by='ID', all.y=T) %>%
                                    tibble::column_to_rownames('coord') %>%
                                    select(name, all_of(names_list))

    rm(bed_info_only_df)
} else { # db-genes mode
    # Profile file only has ID column other than abundance columns
    # Remove all rows with 0 sum
    abundance_only_df <- gene_profile_df %>%
                            mutate(across(-ID, ~ as.numeric(replace(., . == '' | is.na(.), 0)))) %>%
                            tibble::column_to_rownames('ID') %>%
                            filter(rowSums(across(where(is.numeric), ~.x != 0))>0)
    gene_info_annotation_df <- merge(gene_annot_df, gene_profile_df, by='ID', all.y=T) %>%
                                    mutate(ID2 = ID) %>%
                                    tibble::column_to_rownames('ID2')
}

# Write annotations into csv file
fwrite(gene_annot_df, file=paste0(output_dir, '/Annotations.csv'), row.names = F, quote = F)

# Free memory
rm(gene_profile_df, gene_annot_df)
gc()

########################################################
# Split df into metaG and metaT
# Keep only common samples
########################################################

metaG_tpm_profile <- NULL
metaT_tpm_profile <- NULL

# Subset profile file by metaG and metaT
if (omics %like% 'metaG') {
    metaG_tpm_profile = abundance_only_df %>%
                            select(starts_with("metaG")) %>%
                            rename_with(~ gsub("metaG.", "", .x, fixed=TRUE))
}
if (omics %like% 'metaT') {
    metaT_tpm_profile = abundance_only_df %>%
                            select(starts_with("metaT")) %>%
                            rename_with(~ gsub("metaT.", "", .x, fixed=TRUE))
}

# Free memory
rm(abundance_only_df)
gc()

########################################################
# Get samples in common metaG and metaT
########################################################

samples_intersect <- NULL
if (is.null(metaG_tpm_profile)) {
    samples_intersect <- colnames(metaT_tpm_profile)
} else if (is.null(metaT_tpm_profile)) {
    samples_intersect <- colnames(metaG_tpm_profile)
} else {
    samples_intersect <- intersect(colnames(metaG_tpm_profile), colnames(metaT_tpm_profile))
}

########################################################
# Write profile data to files
########################################################

# If there is metadata file, use it instead of the dummy metadata created later
metadata_df <- NULL
if (metadata_file != 'None'){
    metadata_df <- as.data.frame(fread(metadata_file,  header = T), stringsAsFactors = F)
}

# Write GA profile data for metaG
if (omics == 'metaG') {
    # Get profile
    metaG_tpm_profile_sub <- metaG_tpm_profile %>%
                                select(all_of(samples_intersect)) %>%
                                tibble::rownames_to_column('ID') %>%
                                dplyr::select(ID, everything())

    # Write csv file
    fwrite(metaG_tpm_profile_sub, file=paste0(output_dir, '/GA.csv'), row.names = F, quote = F)
    metaG_tpm_profile_sub <- metaG_tpm_profile_sub %>%
                                tibble::column_to_rownames('ID')

    # Write phyloseq
    if (is.null(metadata_df)) {
        metadata_df <- data.frame(sample=names(metaG_tpm_profile_sub), group='control', sample_alias=names(metaG_tpm_profile_sub))
    }
    samp <- sample_data(metadata_df)
    rownames(samp) <- samp$sample_alias
    metaG_physeq <- phyloseq(otu_table(as.matrix(metaG_tpm_profile_sub), taxa_are_rows = T),
                            tax_table(as.matrix(gene_info_annotation_df)), samp)
    saveRDS(metaG_physeq, file = paste0(phyloseq_dir, '/GA.rds'))

    # Make PCA plot
    prepare_PCA(physeq=metaG_physeq, type="abundance",  color=main_factor)

    # Free memory
    rm(metaG_tpm_profile, metaG_physeq)
    gc()
}

# Write GT data for metaT
if (omics == 'metaT') {
    # Get profile
    metaT_tpm_profile_sub <- metaT_tpm_profile %>%
                                select(all_of(samples_intersect)) %>%
                                tibble::rownames_to_column('ID') %>%
                                dplyr::select(ID, everything())

    # Write csv file
    fwrite(metaT_tpm_profile_sub, file=paste0(output_dir, '/GT.csv'), row.names = F, quote = F)
    metaT_tpm_profile_sub <- metaT_tpm_profile_sub %>%
                                tibble::column_to_rownames('ID')

    # Write phyloseq
    if (is.null(metadata_df)) {
        metadata_df <- data.frame(sample=names(metaT_tpm_profile_sub), group='control', sample_alias=names(metaT_tpm_profile_sub))
    }
    samp <- sample_data(metadata_df)
    rownames(samp) <- samp$sample_alias
    metaT_physeq <- phyloseq(otu_table(as.matrix(metaT_tpm_profile_sub), taxa_are_rows = T),
                            tax_table(as.matrix(gene_info_annotation_df)), samp)
    saveRDS(metaT_physeq, file = paste0(phyloseq_dir, '/GT.rds'))

    # Make PCA plot

    prepare_PCA(physeq=metaT_physeq, type="transcript", color=main_factor)

    # Free memory
    rm(metaT_tpm_profile, metaT_physeq)
    gc()
}

if (omics == 'metaG_metaT'){

    ################################## GENE EXPRESSION  ##################################
    print('#################################### GENE EXPRESSION  ####################################')

    # Get profiles
    metaG_tpm_profile_sub <- metaG_tpm_profile %>%
                                select(all_of(samples_intersect))

    metaT_tpm_profile_sub <- metaT_tpm_profile %>%
                                select(all_of(samples_intersect))


    # # Normalization
    # When metaG is 0, just get metaT value (replace metaG 0 values by 1) ####
    # We don't do this anymore. These will become NaN

    # Gene Expression
    print(paste("metaG-metaT match:", identical(names(metaG_tpm_profile_sub), names(metaT_tpm_profile_sub))))
    gene_expression_tpm=(metaT_tpm_profile_sub)/(metaG_tpm_profile_sub)
    gene_expression_tpm_noinf <- gene_expression_tpm %>%
                                               tibble::rownames_to_column('ID') %>%
                                               dplyr::mutate(across(-ID, function(x) ifelse(is.infinite(x) | is.nan(x), NA, x))) %>%
                                               dplyr::select(ID, everything())

    # Write csv file
    fwrite(gene_expression_tpm_noinf, file=paste0(output_dir, '/GE.csv'), row.names = F, quote = F)
    gene_expression_tpm_noinf <- gene_expression_tpm_noinf %>%
                                               tibble::column_to_rownames('ID')

    # Write phyloseq
    if (is.null(metadata_df)) {
        metadata_df <- data.frame(sample=names(gene_expression_tpm_noinf), group='control', sample_alias=names(gene_expression_tpm_noinf))
    }
    samp <- sample_data(metadata_df)
    rownames(samp) <- samp$sample_alias
    GE_physeq <- phyloseq(otu_table(as.matrix(gene_expression_tpm_noinf), taxa_are_rows = T),
                                                     tax_table(as.matrix(gene_info_annotation_df)), samp)
    saveRDS(GE_physeq, file = paste0(phyloseq_dir, '/GE.rds'))

    # Make PCA plot
    prepare_PCA(physeq=GE_physeq, type="expression", color=main_factor)

    # Delete data
    rm(GE_physeq, gene_expression_tpm)
}

rm(metadata_df, gene_info_annotation_df)

print('Done!')
