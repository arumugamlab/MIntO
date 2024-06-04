#!/usr/bin/env Rscript

# '''
# Generate gene expression profile from genome-based mode gene profiles 
#
# Authors: Mani Arumugam
#
# '''

library(optparse)
opt_list <- list(
  make_option("--metadata", type="character", default=NULL, help="metadata file", metavar="file"),
  make_option("--from", type="character", default="sample", help="metadata column connected to sample files", metavar="string"),
  make_option("--to", type="character", default="sample_alias", help="metadata column connected to unique sample alias", metavar="string"),
  make_option("--input", type="character", default=NULL, help="input gene profile file", metavar="file"),
  make_option("--output", type="character", default=NULL, help="output relabelled gene profile file", metavar="file"),
  make_option("--prefix", type="character", default=NULL, help="prefix to add both for <from> string and <to> string", metavar="string"),
  make_option("--threads", type="integer", default=4, help="number of threads [default: %default]")
)
  
opt <- parse_args(OptionParser(option_list=opt_list))

profile_file <- opt$input
metadata_file <- opt$metadata
out_file <- opt$output
map_from <- opt$from
map_to <- opt$to
prefix <- opt$prefix

library(dplyr)
library(data.table)
library(tibble)

setDTthreads(threads = opt$threads)

# Read metadata
metadata_df <- fread(metadata_file, header=T) %>%
                        as.data.frame(stringsAsFactors = F)

# Read input file
profile_df <- fread(profile_file, header=T) %>%
                        as.data.frame(stringsAsFactors = F)

# Create the mapping
if (!is.null(prefix)) {
    mapping <- metadata_df %>%
                    select(all_of(c(map_from, map_to))) %>%
                    mutate_all( ~ paste0(prefix, .)) %>%
                    column_to_rownames(map_from)
} else {
    mapping <- metadata_df %>%
                    select(all_of(c(map_from, map_to))) %>%
                    column_to_rownames(map_from)
}

# Which columns not to relabel
# ID is common for all. gene_length comes from MAG/refgenome mode.
do_not_relabel <- c("ID", "gene_length")

# Which columns to relabel
columns_to_rename <- names(profile_df)
columns_to_rename <- columns_to_rename[!columns_to_rename %in% do_not_relabel]

# Do the relabeling
profile_df <- profile_df %>%
                rename_with(.fn = function(x) mapping[x, ], .cols = all_of(columns_to_rename))

# Write the file out
fwrite(profile_df, file=out_file, sep="\t", row.names=FALSE, quote=FALSE)
