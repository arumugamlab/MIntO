#!/usr/bin/env Rscript

# Parse command line arguments
library(optparse)
opt_list <- list(
                make_option("--normalize", type="character", default=NULL, help="normalization type: MG or TPM"),
                make_option("--bed", type="character", default=NULL, help="gene mapped-read-count bed file", metavar="file"),
                make_option(c("--out"), type="character", default=NULL, help="output file to write normalized counts", metavar="file"),
                make_option("--omics", type="character", default=NULL, help="which omics: metaG or metaT", metavar="omic-type"),
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
omics <- opt$omics
read_n <- opt[['min-read-count']]

##########################  ** Load libraries **  ##########################
library(dplyr)
library(data.table)
library(tibble)
library(purrr)

setDTthreads(threads = threads_n)


## Bed file colnames - not to mutate when applying threshold
gene_info <- c("chr","start","stop","name","score","strand","source","feature","frame","info")

## GENE ABUNDANCES per SAMPLE
# Load data - raw counts
# Filter number of mapped reads bellow the threshold
# Calculate gene length
# Add 'coord' field
gene_abund_bed_coord_df <- as.data.frame(fread(gene_abund_bed, header=T), stringsAsFactors = F) %>%
                                mutate(coord = paste0(chr, '_', start, '_', stop)) %>%
                                mutate(coord = gsub('\\-','\\_', coord)) %>%
                                mutate(gene_length = abs(stop-start) + 1)

# Filter gene abundances df by gene coordinates
# Calculate RPK from gene abundances (gene abundances normalized by gene length)
gene_rpk_bed_df <- gene_abund_bed_coord_df %>%
                                distinct(coord, .keep_all = TRUE) %>%
                                mutate(coord = paste0(chr, '_', start, '_', stop)) %>%
                                mutate(header =  gsub('[-.]', '_', info)) %>%
                                mutate(ID_MAG = stringr::str_split(header, "\\|") %>% map_chr(., 1)) %>% # Get the first field
                                mutate(header = ifelse(name=='.', coord, header)) %>%
                                select(ID_MAG, header, gene_length, everything())

gene_rpk_only_df <- gene_rpk_bed_df %>%
                                select(-all_of(c(gene_info, "coord")))

if (normalize == 'MG') {
    ## Get marker genes
    # replace '-' by '.'
    # Subset COGs
    marker_COGs <- c('COG0012', 'COG0016', 'COG0018', 'COG0172', 'COG0215',
                     'COG0495', 'COG0525', 'COG0533', 'COG0541', 'COG0552')

    MGs_in_MAGs <- fread(fetchMG_table,
                                header=T,
                                col.names = c('header', 'HMM_score', 'COG', 'taxid.projectid')
                               ) %>%
                         as.data.frame(stringsAsFactors=FALSE) %>%
                         mutate(header = gsub('-', '_', header)) %>%
                         select(c('header', 'COG')) %>%
                         filter(COG %in% marker_COGs) %>%
                         mutate(ID_MAG = stringr::str_split(header, '\\|') %>% map_chr(., 1)) # Get the first field

    # Filter genes by min_read and then length-normalize
    gene_rpk_only_df <- gene_rpk_only_df %>%
                            mutate(across(-c(ID_MAG, header, gene_length), \(x) {ifelse(x <= read_n, 0, x)})) %>%
                            mutate(across(-c(ID_MAG, header, gene_length), ~ . / gene_length)) %>%
                            select(-gene_length)

    # Calculate median abundance of the 10 single-copy marker genes per MAG
    median_MG_rpk_per_MAG <- merge(MGs_in_MAGs, gene_rpk_only_df, by=c('ID_MAG','header')) %>%
                                    select(-header, -COG) %>%
                                    dplyr::group_by(ID_MAG) %>%
                                    dplyr::summarise(across(everything(), median)) %>%
                                    mutate(header = 'MG') %>%
                                    select(ID_MAG, header, everything())

    # Merge median abundance of the 10 single-copy marker genes and transcript abundance
    gene_rpk_plus_median_MG_rpk <- rbind(median_MG_rpk_per_MAG, gene_rpk_only_df)

    # melt to make one gene rpk per row for operational convenience
    rpk_melt <- reshape2::melt(gene_rpk_plus_median_MG_rpk, id.vars = c("ID_MAG", "header"))

    # Normalize
    # Remove MG rows
    # Replace NaN and infinite values by 0
    rpk_melt_percell <- rpk_melt %>%
                             group_by(ID_MAG, variable) %>%
                             mutate(value = value/value[header == 'MG']) %>%
                             filter(header != 'MG') %>%
                             mutate(value = ifelse(is.infinite(value) | is.nan(value), 0, value))

    # unmelt to bring back table
    rpk_final_norm <- reshape2::dcast(rpk_melt_percell, ID_MAG + header ~ variable,  value.var="value") #fun.aggregate = list,

} else {
    # Calculate colSums before min-read filtering
    sample_rpk_sums <- gene_rpk_only_df %>%
                            mutate(across(-c(ID_MAG, header, gene_length), ~ . / gene_length)) %>%
                            select(-c(ID_MAG, header, gene_length)) %>%
                            colSums()

    # Apply min-read filtering and then length-normalize
    sample_rpk_df <- gene_rpk_only_df %>%
                            select(-c(ID_MAG, header)) %>%
                            mutate(across(-gene_length, \(x) {ifelse(x <= read_n, 0, x)})) %>%
                            mutate(across(-gene_length, ~ . / gene_length)) %>%
                            select(-gene_length)

    # Normalize with backed-up sum to get TPM, but it may not sum to 1M due to min-read
    sample_rpk_norm <- mapply('/', sample_rpk_df*1e6 , sample_rpk_sums)

    # Get ID_MAG and header back
    rpk_final_norm <- cbind(select(gene_rpk_only_df, ID_MAG, header), sample_rpk_norm)
}

# Get the right columns
rpk_final_norm <- rpk_final_norm %>%
                         select(-ID_MAG) %>%
                         rename_at(vars(!matches(c("header"))), ~ paste0(omics,".", .))

# Merge normalized rpk with gene_info
final_table <- gene_rpk_bed_df %>%
                  select(c('header','coord', all_of(gene_info), 'gene_length')) %>%
                  merge(., rpk_final_norm, by='header') %>%
                  select(-header)# %>%
                  #mutate_if(is.numeric, round, digits = 6)
print(dim(final_table))

write.table(format(final_table, digits=6), gene_norm_csv, row.names = F, col.names = T, sep = "\t", quote = F)
