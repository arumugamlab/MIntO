#!/usr/bin/env Rscript

# '''
# Generate gene expression profile from genome-based mode gene profiles
#
# Authors: Carmen Saenz, Mani Arumugam
#
# '''

##########################  ** Load libraries **  ##########################
library(data.table)
library(optparse)
library(phyloseq)
library(stringr)
library(this.path)

##########################  ** Load functions **  ##########################

# Include common utility functions
source(this.path::here('include', 'utils.R'))

# Include PCA functions
source(this.path::here('include', 'plots_PCA.R'))

# Parse command line arguments
opt_list <- list(
                make_option("--threads", type="integer", default=4, help="number of threads [default: %default]"),
                make_option("--outdir", type="character", default=NULL, help="output directory to write normalized counts", metavar="directory"),
                make_option("--omics", type="character", default=NULL, help="which omics to summarize: metaG, metaT or metaG_metaT", metavar="string"),
                make_option("--funcat-names", type="character", default=NULL, help="comma-delimited list of functional categores to summarize: e.g., 'eggNOG.OGs,dbCAN.EC,kofam.KEGG_KO'", metavar="string"),
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


#####
# Process arguments
#####

funcat_names <- unique(unlist(strsplit(funcat_args, split = "\\,")[[1]]))
logmsg("Including annotations for: ", paste(funcat_names, collapse=","))

setDTthreads(threads = threads_n)
set.seed(1234)

# Input gene abundance profile

profiles_tpm <- input_file

#####
# Common function to make output files for each type
#####
make_output_files <- function(profile=NULL, metadata=NULL, type=NULL, label=NULL) {

    library(qs)

    logmsg("Processing ", type, " data")

    # Write tsv file
    logmsg(" Writing TSV file")
    fwrite(profile, file=paste0(output_dir, '/', label, '.tsv'), sep = "\t", row.names = F, quote = F)

    # Prepare phyloseq
    logmsg(" Making phyloseq")
    physeq <- phyloseq(otu_table(as.matrix(profile, rownames="ID"), taxa_are_rows = T),
                       tax_table(as.matrix(gene_annot_dt, rownames="ID")),
                       metadata)

    # Write phyloseq
    logmsg(" Writing phyloseq")
    qsave(physeq,
          file = paste0(phyloseq_dir, '/', label, '.qs'),
          preset = "high",
          nthreads = threads_n,
         )
}

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

# Read the TPM/MG normalized profile tsv file for this omics
logmsg("Reading profile")
gene_profile_dt <- fread(profiles_tpm, header=T)
logmsg("  done")
logmsg("  profile: ", dim(gene_profile_dt))

# Subset relevant columns from the profile table
logmsg("Subsetting profile columns")

# Get all column names in profile
all_column_names <- names(gene_profile_dt)

# If single omic mode but multi-omic input file, subset columns to relevant samples only
sample_names <- grep(paste0("^", omics, "\\."), all_column_names, fixed=FALSE, perl=TRUE, value=TRUE)
if (omics == 'metaG_metaT') {
    sample_names <- grep(paste0("^meta[GT]\\."), all_column_names, fixed=FALSE, perl=TRUE, value=TRUE)
}
gene_profile_dt <- gene_profile_dt[, c("ID", sample_names), with=FALSE]
logmsg("  ", length(all_column_names), " --> ", length(sample_names))
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
    gene_combined_dt <- gene_combined_dt[, negRowSum := -rowSums(.SD, na.rm = TRUE), .SDcols = sample_cols]
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
    if (length(metadata_df$sample_alias) != length(unique(metadata_df$sample_alias))) {
        stop(paste0("In sample metadata file '", metadata_file, "', column 'sample_alias' does not have unique values per row!"))
    }

    # Reorder metadata by same order in profile file, so that it can be matched properly when ID is removed
    # Easiest way is to do an inner join, which will also nicely fail if a sample doesnt have metadata

    col_names = names(metadata_df)
    if (omics == 'metaT') {
        anchor_df = metaT_profile %>%
                        dplyr::select(-ID) %>%
                        head(2) %>%
                        t() %>%
                        as.data.frame() %>%
                        tibble::rownames_to_column("sample_alias")
    } else {
        anchor_df = metaG_profile %>%
                        dplyr::select(-ID) %>%
                        head(2) %>%
                        t() %>%
                        as.data.frame() %>%
                        tibble::rownames_to_column("sample_alias")
    }
    metadata_df = dplyr::inner_join(anchor_df, metadata_df) %>%
                    dplyr::select(any_of(col_names))
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

# For PCA, if the table is too big, we will limit ourselves to maxN rows
maxN = 500000

# Write GA profile data for metaG
if (omics == 'metaG') {
    make_output_files(profile=metaG_profile, metadata=sample_metadata, type='abundance', label='GA')

    # Get first maxN rows and remove ID
    metaG_profile <- metaG_profile[, first(.SD, maxN)][, ID := NULL]
    gc()

    # Make PCA plot
    prepare_PCA(profile=metaG_profile, title='gene abundance', label='GA', color=main_factor, metadata=sample_metadata)
}

# Write GT data for metaT
if (omics == 'metaT') {
    make_output_files(profile=metaT_profile, metadata=sample_metadata, type='transcript', label='GT')

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

    make_output_files(profile=gene_expression, metadata=sample_metadata, type='expression', label='GE')

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


# Write annotations into tsv file
#    Replace NA by "-" to catch later
#    Concat all annotations
#    Remove empty annotations (meaning '-' in every column
#    Delete annotation column
logmsg("Removing empty annotations")
logmsg("  Before  : ", dim(gene_annot_dt))
gene_annot_dt <- (
                  gene_annot_dt
                  [, lapply(.SD, function(x) str_replace_na(x, replacement="-"))]
                  [, annot := do.call(paste, c(.SD, sep = "")), .SDcols = funcat_names]
                 )
setkey(gene_annot_dt, annot)
empty_annotation <- strrep("-", length(funcat_names))
gene_annot_dt <- gene_annot_dt[!J(empty_annotation)][, annot := NULL]
logmsg("  After   : ", dim(gene_annot_dt))
logmsg("  done")
logmsg("Writing out annotations")
fwrite(gene_annot_dt,
       file = paste0(output_dir, '/Annotations.tsv'),
       sep = '\t',
       row.names = FALSE,
       quote = FALSE)
logmsg("  done")

# Freeing memory, but only to check what was the peak memory usage.
logmsg("Freeing memory")
gc()
logmsg("Done!")
