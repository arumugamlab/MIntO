#!/usr/bin/env Rscript

# '''
# Generate gene expression profile from genome-based mode gene profiles
#
# Authors: Carmen Saenz, Mani Arumugam
#
# '''

##########################  ** Load libraries **  ##########################
library(data.table)
library(KEGGREST)
library(optparse)
library(phyloseq)
library(stringr)

logmsg <- function(...) {
    message("# ", format(Sys.time(), digits=0), " - ", paste(list(...)), collapse=" ")
}

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

funcat_names <- unique(unlist(strsplit(funcat_args, split = "\\,")[[1]]))
logmsg("Including annotations for: ", paste(funcat_names, collapse=","))

setDTthreads(threads = threads_n)
set.seed(1234)

#####
# Gene abundance profile
#####

profiles_tpm <- input_file

logmsg('#################################### Making output directory ####################################')

# Generate output directories ####

dir.create(file.path(output_dir), showWarnings = FALSE)
visual_dir=paste0(output_dir, '/plots')
dir.create(file.path(visual_dir), showWarnings = FALSE)
phyloseq_dir=paste0(output_dir, '/phyloseq_obj')
dir.create(file.path(phyloseq_dir), showWarnings = FALSE)

########################################################
## Read input files
########################################################

logmsg('#################################### Reading input files ####################################')

# Get gene annotations
logmsg("Reading annotations")
gene_annot_dt <- fread(annot_file, header=T, select=c("ID", funcat_names))
logmsg("  done")
logmsg("Indexing")
setkey(gene_annot_dt, ID)
logmsg("  done")
logmsg("  annot: ", dim(gene_annot_dt))

# Read the TPM/MG normalized profile csv file for this omics
logmsg("Reading profile")
gene_profile_dt <- fread(profiles_tpm, header=T)
logmsg("  done")
logmsg("  profile: ", dim(gene_profile_dt))

# Subset relevant columns from the profile table
logmsg("Subsetting profile columns")

# Get all column names in profile
all_column_names <- names(gene_profile_dt)

# Rename the key column info --> ID in refgenome/MAG mode
if (!"ID" %in% all_column_names && "info" %in% all_column_names) {
    setnames(gene_profile_dt, c("info"), c("ID"))
}

# If single omic mode but multi-omic input file, subset columns to relevant samples only
sample_names <- grep(paste0("^", omics, "\\."), all_column_names, fixed=FALSE, perl=TRUE, value=TRUE)
if (omics == 'metaG_metaT') {
    sample_names <- grep(paste0("^meta[GT]\\."), all_column_names, fixed=FALSE, perl=TRUE, value=TRUE)
}
gene_profile_dt <- gene_profile_dt[, c("ID", sample_names), with=FALSE]
logmsg("  done")

# Index
logmsg("Indexing")
setkey(gene_profile_dt, ID)
logmsg("  done")

# Make a list of sample columns for easy lookup later
sample_cols <- setdiff(colnames(gene_profile_dt), c("ID"))

########################################################
## Get abundance and gene information into separate DFs
########################################################

logmsg("#################################### Processing data ####################################")

# Merge tables
logmsg("Combining profiles and annotations")
gene_combined_dt <- gene_annot_dt[gene_profile_dt, ] # outer join - profiles should be kept for unannotated genes
logmsg("  done")

# Free memory
logmsg("Freeing memory")
rm(gene_profile_dt, gene_annot_dt)
gc()
logmsg("  done")

# For metaG_metaT, the real profile is the ratio. So don't do any filtering now.
# For others, do the filtering

if (omics != 'metaG_metaT') {

    # Estimate rowSum but negate it so that setkey will sort by desc(rowSum)
    logmsg("Estimating rowSum")
    gene_combined_dt <- gene_combined_dt[, negRowSum := rowSums(.SD, na.rm = TRUE), .SDcols = sample_cols]
    logmsg("  done")

    # Set key
    logmsg("Indexing")
    setkey(gene_combined_dt, negRowSum)
    logmsg("  done")

    # Remove zerosum rows

    logmsg("Removing zero-sum rows")
    logmsg("  Before: ", dim(gene_combined_dt))
    gene_combined_dt <- gene_combined_dt[!J(0)] # This is the so-called 'not join' idiom
    logmsg("  After : ", dim(gene_combined_dt))
    logmsg("  done")

    # NOTE: By now, profile table is free of zerosum rows and sorted(desc) by rowSum
}

# Prepare separate tables
# 1. Annotations
#    Get ID and the relevant functional categories
# 2. Profile
#    Get the ID and sample columns
logmsg("Preparing separate annotation and profile tables")
gene_annot_dt <- gene_combined_dt[, c("ID", funcat_names), with=FALSE]
gene_combined_dt <- gene_combined_dt[, c("ID", sample_cols), with=FALSE]
logmsg("  done")

########################################################
# For metaG_metaT, split profiles into metaG and metaT
########################################################

metaG_profile <- NULL
metaT_profile <- NULL

all_column_names <- names(gene_combined_dt)

logmsg("Synchronizing metaG and metaT if applicable")
# Subset profile file by metaG and metaT
if (omics %like% 'metaG') {
    old_names     <- grep("^metaG\\.", all_column_names, fixed=FALSE, perl=TRUE, value=TRUE)
    new_names     <- gsub("^metaG\\.", "", old_names, fixed=FALSE, perl=TRUE)
    metaG_profile <- gene_combined_dt[, c("ID", old_names), with=FALSE]
    setnames(metaG_profile, old_names, new_names)
}
if (omics %like% 'metaT') {
    old_names     <- grep("^metaT\\.", all_column_names, fixed=FALSE, perl=TRUE, value=TRUE)
    new_names     <- gsub("^metaT\\.", "", old_names, fixed=FALSE, perl=TRUE)
    metaT_profile <- gene_combined_dt[, c("ID", old_names), with=FALSE]
    setnames(metaT_profile, old_names, new_names)
}
logmsg("  done")

# Free memory
logmsg("Freeing memory")
rm(gene_combined_dt)
gc()
logmsg("  done")

########################################################
# For metaG_metaT, subset profiles for common samples
########################################################

# Columns shared b/w metaG and metaT. This includes 'ID', which is also needed in final table.
if (omics == 'metaG_metaT') {
    intersect_columns <- intersect(colnames(metaG_profile), colnames(metaT_profile))
    metaG_profile <- metaG_profile[, intersect_columns, with=FALSE]
    metaT_profile <- metaT_profile[, intersect_columns, with=FALSE]
}

########################################################
# Make metadata df
########################################################

# If there is metadata file, use it instead of the dummy metadata created later
metadata_df <- NULL
if (metadata_file != 'None'){
    metadata_df <- as.data.frame(fread(metadata_file,  header = T), stringsAsFactors = F)
}
if (is.null(metadata_df)) {
    if (omics == 'metaT') {
        metadata_df <- data.frame(sample=names(metaT_profile), group='control', sample_alias=names(metaT_profile))
    } else {
        metadata_df <- data.frame(sample=names(metaG_profile), group='control', sample_alias=names(metaG_profile))
    }
}
sample_metadata <- sample_data(metadata_df)
rownames(sample_metadata) <- sample_metadata$sample_alias

########################################################
# Write output files
########################################################

logmsg("#################################### Writing output files ####################################")

# Common function to make output files for each type
make_output_files <- function(profile=NULL, type=NULL, label=NULL) {

    library(qs)

    logmsg("Processing ", type, " data")

    # Write csv file
    logmsg(" Writing CSV file")
    fwrite(profile, file=paste0(output_dir, '/', label, '.csv'), row.names = F, quote = F)

    # Prepare phyloseq
    logmsg(" Making phyloseq")
    physeq <- phyloseq(otu_table(as.matrix(profile, rownames="ID"), taxa_are_rows = T),
                       tax_table(as.matrix(gene_annot_dt, rownames="ID")),
                       sample_metadata)

    # Write phyloseq
    logmsg(" Writing phyloseq")
    qsave(physeq,
          file = paste0(phyloseq_dir, '/', label, '.qs'),
          preset = "balanced",
          nthreads = threads_n,
         )
}

# Requirements:
#  'profile' should only have numeric columns - IDs should have been removed
#  'profile' should have removed all zerosum rows
#  'profile' should be sorted (desc) by rowSum
plot_PCA <- function(profile, label, color, metadata) {

    library(ggplot2)
    library(ggrepel)

    ##### metadata:
    sample_data_df <- data.frame(metadata, stringsAsFactors = F)

    ### Counts - non-zero/inf
    ## Prepare data - replace NA values by 0
    ## Update 2024.10.05 - NA/Inf should not be there in the tables anymore.

    # Get the lowest value in table and its log() to replace 0s by log2(min value * 1e-2) for the PCA
    profile_copy <- (
                     profile
                     [, lapply(.SD, function(x) min(x[x>0], na.rm = TRUE))]
                     [, minv := min(.SD, na.rm = TRUE)]
                    )
    min_value <- (min(profile_copy[, .(minv)]))
    logmsg("  Min value in table: ", min_value)
    replace_val <- log(min_value*1e-2)

    # Log transform values
    pca_data <- profile[, lapply(.SD, function(x) ifelse(x==0, replace_val, log(x)))]
    gc()
    pca_data <- data.table::transpose(pca_data)
    gc()


    #Plotting scores of PC1 and PC2 with log transformation
    pca_results <- prcomp(pca_data, center = T, sca=T)

    # ***************** ggplot way - w coloring *****************
    dtp <- data.frame('sample' = sample_data_df[[label]],
                      'group' = sample_data_df[[color]],
                       pca_results$x[,1:2]) # the first two componets are selected
    total_variance <- sum(pca_results$sdev)
    axis_names <- head(colnames(data.table(pca_results$x)), 2)
    percentage <- round(head(pca_results$sdev, 2) / total_variance * 100, 2)
    percentage <- paste0(axis_names, " (", percentage, "%)")

    PCA_Sample_site_abundance <- (ggplot(data = dtp, aes(x = PC1, y = PC2, color = group)) +
                                  geom_text_repel(aes(label = sample),nudge_x = 0.04, size = 3.5, segment.alpha = 0.5) +
                                  geom_point(size = 2, shape = 16)+
                                  xlab(percentage[1]) + ylab(percentage[2]) +
                                  #labs(title = title) +
                                  theme_bw()+
                                  theme(plot.title = element_text(size=10), legend.position="bottom"))#+ stat_ellipse(type = "norm", linetype = 2))
    return(PCA_Sample_site_abundance)
}

prepare_PCA <- function(profile, label, color, metadata, title) {

    logmsg(" Making PCA plots")

    # Colors
    manual_plot_colors =c('#9D0208', '#264653','#e9c46a','#D8DCDE','#B6D0E0',
                          '#FFC87E','#F4A261','#E34F33','#E9C46A',
                          '#A786C9','#D4C0E2','#975773','#6699FF','#000066',
                          '#7AAFCA','#006699','#A9D181','#2F8475','#264445')

    # Title
    title_name <- paste0('PCA - ', title)

    # File names
    # abundance  -> GA
    # transcript -> GT
    # expression -> GE
    out_name <- paste0(visual_dir, '/', label, '.PCA.pdf')

    plot_PCA_out <- plot_PCA(profile=profile, color=color, label="sample_alias", metadata=metadata)
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

# For PCA, if the table is too big, we will limit ourselves to maxN rows
maxN = 500000

# Write GA profile data for metaG
if (omics == 'metaG') {
    make_output_files(profile=metaG_profile, type='abundance', label='GA')

    # Get first maxN rows and remove ID
    metaG_profile <- metaG_profile[, first(.SD, maxN)][, ID := NULL]
    gc()

    # Make PCA plot
    prepare_PCA(profile=metaG_profile, title='gene abundance', label='GA', color=main_factor, metadata=sample_metadata)
}

# Write GT data for metaT
if (omics == 'metaT') {
    make_output_files(profile=metaT_profile, type='transcript', label='GT')

    # Get first maxN rows and remove ID
    metaT_profile <- metaT_profile[, head(.SD, maxN)][, ID := NULL]
    gc()

    # Make PCA plot
    prepare_PCA(profile=metaT_profile, title='gene transcript', label='GT', color=main_factor, metadata=sample_metadata)
}

if (omics == 'metaG_metaT') {

    # # Normalization
    # When metaG is 0, just get metaT value (replace metaG 0 values by 1) ####
    # We don't do this anymore. These will become NaN

    # Gene Expression
    logmsg(paste("metaG-metaT match:", identical(names(metaG_profile), names(metaT_profile))))
    sample_cols <- setdiff(colnames(metaG_profile), c("ID"))

    # Since ID is non-numerical, we calculate the ratios after removing ID column
    IDs <- metaG_profile[, .(ID)]
    gene_expression <- metaT_profile[, ID := NULL]/metaG_profile[, ID := NULL]

    # Now let us add ID back
    gene_expression <- gene_expression[, lapply(.SD, function(x) ifelse(is.infinite(x) | is.nan(x), NA, x))][, ID := IDs]
    # Make ID the 1st column
    setcolorder(gene_expression, c('ID'))

    # Free up memory
    logmsg("Freeing memory")
    rm(metaG_profile, metaT_profile, IDs)
    gc()

    # Estimate rowSum but negate it so that setkey will sort by desc(rowSum)
    logmsg("Estimating rowSum")
    gene_expression <- gene_expression[, negRowSum := -rowSums(.SD, na.rm = TRUE), .SDcols = sample_cols]
    logmsg("  done")

    # Set key
    logmsg("Indexing")
    setkey(gene_expression, negRowSum)
    logmsg("  done")

    # Remove zerosum rows
    logmsg("Removing zero-sum rows")
    logmsg("  Before: ", dim(gene_expression))
    gene_expression <- gene_expression[!J(0)][, negRowSum := NULL]
    logmsg("  After : ", dim(gene_expression))
    logmsg("  done")

    # NOTE: By now, profile table is free of zerosum rows and sorted(desc) by rowSum

    make_output_files(profile=gene_expression, type='expression', label='GE')

    # Prepare data for PCA

    # Get first maxN rows
    gene_expression <- (
                        gene_expression
                        [, head(.SD, maxN)]
                        [, ID := NULL]
                        [, lapply(.SD, function(x) ifelse(is.na(x), 0, x))]
                       )
    gc()

    # Make PCA plot
    prepare_PCA(profile=gene_expression, title='gene expression', label='GE', color=main_factor, metadata=sample_metadata)
    rm(gene_expression)
}


# Write annotations into csv file
#    Replace , by ; in annotations, since output will be csv
#    Replace NA by "-" to catch later
#    Concat all annotations
#    Remove empty annotations (meaning '-' in every column
#    Delete annotation column
logmsg("Removing empty annotations")
logmsg("  Before  : ", dim(gene_annot_dt))
gene_annot_dt <- (
                  gene_annot_dt
                  [, lapply(.SD, function(x) str_replace_all(str_replace_na(x, replacement="-"), pattern=",", replacement=";"))]
                  [, annot := do.call(paste, c(.SD, sep = "")), .SDcols = funcat_names]
                 )
setkey(gene_annot_dt, annot)
empty_annotation <- strrep("-", length(funcat_names))
gene_annot_dt <- gene_annot_dt[!J(empty_annotation)][, annot := NULL]
logmsg("  After   : ", dim(gene_annot_dt))
logmsg("  done")
logmsg("Writing out annotations")
fwrite(gene_annot_dt,
       file = paste0(output_dir, '/Annotations.csv'),
       row.names = FALSE,
       quote = FALSE)
logmsg("  done")

# Freeing memory, but only to check what was the peak memory usage.
logmsg("Freeing memory")
gc()
logmsg("Done!")
