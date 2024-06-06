#!/usr/bin/env Rscript

# Normalize gene abundances and write profiles.

# Parse command line arguments
library(optparse)
opt_list <- list(
                make_option("--normalize", type="character", default=NULL, help="normalization type: MG or TPM"),
                make_option("--bed", type="character", default=NULL, help="gene mapped-read-count bed file", metavar="file"),
                make_option(c("--out"), type="character", default=NULL, help="output file to write normalized counts", metavar="file"),
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
read_n <- opt[['min-read-count']]

##########################  ** Load libraries **  ##########################

library(data.table)

setDTthreads(threads = threads_n)

## GENE ABUNDANCES per SAMPLE
# Load data - raw counts
# Calculate gene length
# Remove useless columns
gene_abundance <- fread(gene_abund_bed, header=T, data.table=TRUE)
setkey(gene_abundance, ID)

# Make a list of sample columns for easy lookup later
sample_cols <- setdiff(colnames(gene_abundance), c("ID", "gene_length", "ID_MAG"))

# Normalize gene abundances by gene length and MG/TPM
if (normalize == 'MG') {
    # Define 10 marker genes as subset of COGs
    marker_COGs <- c('COG0012', 'COG0016', 'COG0018', 'COG0172', 'COG0215',
                     'COG0495', 'COG0525', 'COG0533', 'COG0541', 'COG0552')

    ## Get marker genes in this genome-set
    MGs_in_MAGs <- fread(fetchMG_table,
                         header=T,
                         data.table=TRUE,
                         select=c('#protein_sequence_id', 'COG')
                        )
    setnames(MGs_in_MAGs, c('ID', 'COG'))
    setkey(MGs_in_MAGs, COG)
    MGs_in_MAGs <- MGs_in_MAGs[J(marker_COGs), ]
    setkey(MGs_in_MAGs, ID)

    # Filter genes by min_read and then length-normalize
    gene_abundance <- (
                       gene_abundance
                       [, ID_MAG := sub('\\|.*', '', ID)] # Retain the first pipe-delimited field
                       [, c(sample_cols) := lapply(.SD, function(x) ifelse(x <= read_n, 0, x / get("gene_length"))), .SDcols = sample_cols]
                       [, gene_length := NULL]
                      )

    # Calculate median abundance of the 10 single-copy marker genes per MAG in each sample
    #
    # Sometimes there are multiple copies of the same COG in a genome.
    # This could be biological but also spurios from genome assembly.
    # To avoid issues, we first take max value for a MAG-COG pair, and then take median.
    # This way, each COG contributes at most once for estimating median COG abundance.
    median_MG_per_MAG <- (
                          merge(MGs_in_MAGs, gene_abundance, by='ID')
                          [, lapply(.SD, max) , by = c('COG', 'ID_MAG'), .SDcols = sample_cols]
                          [, lapply(.SD, median) , by = ID_MAG, .SDcols = sample_cols]
                          [, ID := "MG"]
                         )

    # Merge median abundance of the 10 single-copy marker genes and transcript abundance
    gene_abundance <- rbind(median_MG_per_MAG, gene_abundance)

    ########################################################################################################
    # MG NORMALIZATION FUNCTION:
    # --------------------------
    # normalize with the median abundance of the 10 single-copy marker genes.
    # Each time, x is a vector of RPK values for a single MAG within a single sample.
    # Since we rbind'ed median-MG-rpk and then gene-level-RPK, x[1] is always median-MG-rpk.
    # x[2:] are the actual genes.
    # When we divide each vector by its first element, we essentially normalize all genes by median-MG-rpk.
    ########################################################################################################
    norm_rpk_per_MG <- function(x){
      if (is.character(x[1])){
        return(x)
      }else{
        x = x/x[1]
        return(x)
      }
    }

    gene_abundance <- (
                       gene_abundance
                       [, lapply(.SD, norm_rpk_per_MG), by = ID_MAG]
                       [, ID_MAG := NULL]
                       [! ID == 'MG']
                       [, lapply(.SD, function(x) ifelse(is.infinite(x) | is.nan(x), 0, x))]
                      )
} else {

    # Calculate colSums before min-read filtering; and scale to million for TPM normalization
    sample_rpk_sums <- (
                        gene_abundance
                        [, lapply(.SD, function(x) x / get("gene_length")), .SDcols = sample_cols]
                        [, lapply(.SD, sum)]
                        [, lapply(.SD, function(x) x / 1e6)]
                        [, ID := 'SUM']
                       )

    #Apply min-read filtering, length-normalize
    gene_abundance <- (
                       gene_abundance
                       [, c(sample_cols) := lapply(.SD, function(x) ifelse(x <= read_n, 0, x / get("gene_length"))), .SDcols = sample_cols]
                       [, gene_length := NULL]
                      )

    # Merge median abundance of the 10 single-copy marker genes and transcript abundance
    gene_abundance <- rbind(sample_rpk_sums, gene_abundance)

    ########################################################################################################
    # TPM NORMALIZATION FUNCTION:
    # --------------------------
    # normalize with the total length-normalized read-counts.
    # Each time, x is a vector of RPK values for a single sample.
    # Since we rbind'ed SUM and then gene-level-RPK, x[1] is always SUM.
    # x[2:] are the actual genes.
    # When we divide each vector by its first element, we essentially normalize all genes by SUM.
    ########################################################################################################
    norm_rpk_TPM <- function(x){
      if (is.character(x[1])){
        return(x)
      }else{
        x = x/x[1]
        return(x)
      }
    }

    # Normalize with backed-up 'sum / 1M' to get TPM, but it may not sum to 1M due to min-read
    gene_abundance <- (
                       gene_abundance
                       [, lapply(.SD, norm_rpk_TPM)]
                       [! ID == 'SUM', ]
                      )
}

setcolorder(gene_abundance, "ID")

print(dim(gene_abundance))

# Write the table
fwrite(gene_abundance, file = gene_norm_csv, row.names = F, col.names = T, sep = "\t", quote = F)

# Call gc() to report peak memory usage
gc()
