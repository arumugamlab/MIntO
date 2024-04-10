#!/usr/bin/env Rscript

# Normalize gene abundances and write profiles with '<prefix><sample>' as sample-id.

# Parse command line arguments
library(optparse)
opt_list <- list(
                make_option("--normalize", type="character", default=NULL, help="normalization type: MG or TPM"),
                make_option("--bed", type="character", default=NULL, help="gene mapped-read-count bed file", metavar="file"),
                make_option(c("--out"), type="character", default=NULL, help="output file to write normalized counts", metavar="file"),
                make_option("--sample-prefix", type="character", default="", help="prefix for sample-id in the output table: typically, 'metaG.' or 'metaT.' [default: none]", metavar="sample-prefix"),
                make_option("--min-read-count", type="integer", default=2, help="minimum mapped read-count to consider gene is present [default: %default]"),
                make_option("--MG", type="character", default=NULL, help="marker gene table", metavar="file"),
                make_option("--threads", type="integer", default=4, help="number of threads [default: %default]"),
                make_option("--memory", type="integer", default=10, help="maximum memory to be used")
                )
opt <- parse_args(OptionParser(option_list=opt_list))

normalize <- opt$normalize
gene_abund_bed <- opt$bed
gene_norm_csv <- opt$out
threads_n <- opt$threads
memory_lim <- opt$memory
fetchMG_table <- opt$MG
prefix <- opt[['sample-prefix']]
read_n <- opt[['min-read-count']]

##########################  ** Load libraries **  ##########################

library(data.table)
library(purrr)

setDTthreads(threads = threads_n)

## Bed file colnames - not to mutate when applying threshold
gene_info <- c("chr","start","stop","name","score","strand","source","feature","frame","info")

## GENE ABUNDANCES per SAMPLE
# Load data - raw counts
# Filter number of mapped reads bellow the threshold
# Calculate gene length
# Add 'coord' field
# If name is '.', then his is refgenome mode, so set header=coord.
# fetchMG cannot handle '.' in gene names, so we had replaced '.' with '-' in gene IDs.
# To match that, we do a gsub from '.' to '-' here
gene_rpk_bed_DT <- (
                     unique(fread(gene_abund_bed, header=T, data.table=TRUE)[, coord := paste0(chr, '_', start, '_', stop)],
                            by='coord')
                     [, `:=`(
                             header =  gsub('\\.', '-', info),
                             gene_length = abs(stop-start) + 1
                            )
                     ]
                     [, ID_MAG := stringr::str_split(header, "\\|") %>% map_chr(., 1)] # Get the first field
                     [name=='.', header := coord]
                   )

# Make a list of sample columns for easy lookup later
sample_cols <- setdiff(colnames(gene_rpk_bed_DT), c(gene_info, "coord", "header", "gene_length", "ID_MAG"))

# Filter gene abundances df by gene coordinates
# Calculate RPK from gene abundances (gene abundances normalized by gene length)
if (normalize == 'MG') {
    # Define 10 marker genes as subset of COGs
    marker_COGs <- c('COG0012', 'COG0016', 'COG0018', 'COG0172', 'COG0215',
                     'COG0495', 'COG0525', 'COG0533', 'COG0541', 'COG0552')

    ## Get marker genes in this genome-set
    MGs_in_MAGs <- (
                     fread(fetchMG_table,
                                header=T,
                                data.table=TRUE,
                                col.names = c('header', 'HMM_score', 'COG', 'taxid.projectid')
                               )
                     [COG %in% marker_COGs, c('header', 'COG')]
                     [, ID_MAG := stringr::str_split(header, '\\|') %>% map_chr(., 1)] # Get the first field
                   )

    # Filter genes by min_read and then length-normalize

    cols_to_remove <- c("coord", "gene_length", gene_info)
    gene_rpk_only_DT <- (
                          gene_rpk_bed_DT
                          [, c(sample_cols) := lapply(.SD, function(x) ifelse(x <= read_n, 0, x / get("gene_length"))), .SDcols = sample_cols]
                          [, -cols_to_remove, with = FALSE]
                        )

    # Calculate median abundance of the 10 single-copy marker genes per MAG in each sample
    median_MG_rpk_per_MAG_DT <- (
                                   merge(MGs_in_MAGs, gene_rpk_only_DT, by=c('ID_MAG', 'header'), sort = FALSE)
                                   [, lapply(.SD, median) , by = ID_MAG, .SDcols = sample_cols]
                                   [, header := "MG"]
                                )

    # Merge median abundance of the 10 single-copy marker genes and transcript abundance
    gene_rpk_plus_median_MG_rpk <- rbind(median_MG_rpk_per_MAG_DT, gene_rpk_only_DT)

    # normalize with the median abundance of the 10 single-copy marker genes
    norm_rpk_per_MG <- function(x){
      if (is.character(x[1])){
        return(x)
      }else{
        x = x/x[1]
        return(x)
      }
    }

    rpk_final_norm <- (
                        gene_rpk_plus_median_MG_rpk
                        [, lapply(.SD, norm_rpk_per_MG), by = ID_MAG]
                        [! header == 'MG']
                        [, lapply(.SD, function(x) ifelse(is.infinite(x) | is.nan(x), 0, x))]
                        [, -c('ID_MAG')]
                      )

} else {
    # Remove extra columns
    cols_to_remove <- c("coord", gene_info)
    gene_rpk_only_DT <- (
                          gene_rpk_bed_DT
                          [, -cols_to_remove, with = FALSE]
                        )

    # Calculate colSums before min-read filtering
    sample_rpk_sums_DT <- (
                            gene_rpk_only_DT
                            [, lapply(.SD, function(x) x / get("gene_length")), .SDcols = sample_cols]
                            [, lapply(.SD, sum)]
                          )

    #Apply min-read filtering and then length-normalize
    sample_rpk_DT <- (
                        gene_rpk_only_DT
                        [, c(sample_cols) := lapply(.SD, function(x) ifelse(x <= read_n, 0, x / get("gene_length"))), .SDcols = sample_cols]
                        [,-c("ID_MAG", "header", "gene_length")]
                     )

    # Normalize with backed-up sum to get TPM, but it may not sum to 1M due to min-read
    sample_rpk_norm <- mapply('/', sample_rpk_DT*1e6 , sample_rpk_sums_DT)

    # Get ID_MAG and header back
    rpk_final_norm <- cbind(gene_rpk_only_DT[, c('ID_MAG', 'header')], sample_rpk_norm)
}

# For the sample columns, add 'metaG.' or 'metaT.' prefix from --sample-prefix
renamed_sample_cols <- paste0(prefix, sample_cols)
setnames(rpk_final_norm, sample_cols, renamed_sample_cols)

# Merge normalized rpk with gene_info and pick cols to print
cols_to_print <- c('coord', gene_info, 'gene_length', renamed_sample_cols)
final_table <- (
                 merge(gene_rpk_bed_DT, rpk_final_norm, by=c('header'), sort = FALSE)
                 [, cols_to_print, with = FALSE]
               )
print(dim(final_table))

# Write the table
fwrite(final_table, file = gene_norm_csv, row.names = F, col.names = T, sep = "\t", quote = F)
