#!/usr/bin/env Rscript

# ---
# title: "Cumulative read length plot"
# author: "Carmen Saenz"
# ---

##  ****** Pre-processing of metaG and metaT data step - quality and read length filtering ******
# Calculation of MINLEN parameter in Trimmomatic: 
# Check the distribution of the cumulative reads per base position on all the trimmed by quality paired-read fastq files

# cumulative reads distribution per bp = ((total reads - cumulative sum of read per bp)/ total reads)x 1e6

# Keep the % of forward-reverse reads specified in the config file
# If the estimated read length cutoff is below 50bp, trimmomatic will use 50bp as the minimum sequence length for MINLEN (Trimmomatic).

# Parse command line arguments
library(optparse)
option_list = list(
                make_option(c("--input"),      type="character", default=NULL, help="input file", metavar="character"),
                make_option(c("--frac"),       type="double",    default=NULL, help="minimum fraction of reads to keep", metavar="numeric"),
                make_option(c("--out_plot"),   type="character", default=NULL, help="output file with plots", metavar="character"),
                make_option(c("--out_cutoff"), type="character", default=NULL, help="output file with length cutfoff value", metavar="character")
                )
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (any(is.null(c(opt$input, opt$frac, opt$out_plot, opt$out_cutoff)))) {
  print_help(opt_parser)
  stop("Missing required arguments\n", call.=FALSE)
}


# Load libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

# Arguments
input_file = opt$input
fraction_remain = as.numeric(opt$frac)
output_plot = opt$out_plot
output_file = opt$out_cutoff

read_len_df <- as.data.frame(fread(file=input_file, sep=" "), header=FALSE, stringsAsFactors = F)
colnames(read_len_df) <- c('n_reads', 'len_reads', 'sample')

# Total reads per sample
read_len_df_total_reads = data.frame(read_len_df %>% 
                                       dplyr::group_by(sample) %>% 
                                       dplyr::summarise(total_reads = sum(n_reads)))

# Calculate cumulative reads distribution per bp per sample
read_len_df_cumsum = data.frame(read_len_df %>% 
                                  dplyr::group_by(sample) %>% 
                                  dplyr::arrange(len_reads) %>% # order by read length
                                  dplyr::mutate(cumsum = cumsum(n_reads)))

# Combine files to calculate total number of reads
read_len_df_cumsum_total_reads = merge(read_len_df_cumsum, read_len_df_total_reads, by = 'sample')
read_len_df_cumsum_total_reads$cumsum_total=as.numeric(read_len_df_cumsum_total_reads$total_reads -read_len_df_cumsum_total_reads$cumsum)

read_len_df_total_reads$cumsum_total <- read_len_df_total_reads$total_reads
read_len_df_total_reads$n_reads <- read_len_df_total_reads$len_reads <- read_len_df_total_reads$cumsum <- 0

read_len_df_plot <- merge(read_len_df_cumsum_total_reads, read_len_df_total_reads, all=TRUE)
read_len_df_plot$cumsum_total_perc=as.numeric((100*read_len_df_plot$cumsum_total)/read_len_df_plot$total_reads)

# https://stackoverflow.com/questions/31150028/insert-missing-time-rows-into-a-dataframe
read_len_df_plot <- merge(expand.grid(sample=unique(read_len_df_plot$sample), 
                          len_reads=min(read_len_df_plot$len_reads):max(read_len_df_plot$len_reads)),
              read_len_df_plot, all=TRUE)
read_len_df_plot <- read_len_df_plot[,(names(read_len_df_plot) %in% c('sample', 'len_reads', 'cumsum_total_perc'))]

read_len_df_plot <- read_len_df_plot[with(read_len_df_plot, order(sample, len_reads)),]
# https://stackoverflow.com/questions/23340150/replace-missing-values-na-with-most-recent-non-na-by-group
read_len_df_plot <- read_len_df_plot %>% dplyr::group_by(sample) %>% tidyr::fill(cumsum_total_perc)

# Prepare data to plot
read_len_df_plot$sample_name = gsub('(.*)\\_(.*)',"\\1",read_len_df_plot$sample)
read_len_df_plot$sample_pair  = gsub('(.*)\\_(.*)',"\\2",read_len_df_plot$sample)
read_len_df_plot$fwd_rv <- 'Forward reads'
read_len_df_plot$fwd_rv[read_len_df_plot$sample_pair == "2"]<- 'Reverse reads'

# Plot cumulative reads distribution per bp per sample
title = "Read length distribution after trimming"
read_len_plot = (ggplot(data=read_len_df_plot, aes(x=len_reads, y=cumsum_total_perc, group=len_reads)) +
                    geom_boxplot(outlier.colour="grey", outlier.shape=16, outlier.size=0.5, color="grey") +
                    scale_y_continuous(minor_breaks = seq(80, 100, 5), breaks = seq(0, 100, 20)) +
                    ylim(40, 100) +
                    theme_bw() +
                    theme(axis.text = element_text(size = 9.5), legend.position = 'none') +
                    labs(x ="read trim length (bp)", y="fraction remaining after trimming (%)", title = title) +
                    facet_grid(fwd_rv~.))

# Calculate MINLEN parameter in Trimmomatic to keep the % of forward-reverse reads specified in the config file
read_len_df_plot$cumsum_total_perc <- read_len_df_plot$cumsum_total_perc

read_len_df_plot_filter <- read_len_df_plot %>% 
mutate_at(vars(cumsum_total_perc), round, 2) %>% 
filter(cumsum_total_perc >= fraction_remain) #%>% 
#filter(cumsum_total_perc < (fraction_remain+1))

read_len_df_plot_filter2 <- as.data.frame(read_len_df_plot_filter %>% 
                                            dplyr::group_by(sample) %>% 
                                            dplyr::filter(cumsum_total_perc %in% min(cumsum_total_perc)) %>% 
                                            dplyr::filter(len_reads %in% max(len_reads)))

# MINLEN parameter in Trimmomatic

min_len_read = min(as.numeric(read_len_df_plot_filter2$len_reads))
print(paste0('Estimated trimming length cutoff to keep ', fraction_remain, '% of the reads: ',min_len_read, 'bp'))
if (min_len_read < 50){
  min_len_read = 50
}

print(paste0('Recommended trimming length cutoff to keep ', fraction_remain, '% of the reads: ',min_len_read, 'bp'))
cat(min_len_read, file=output_file, sep="\n")

# Output plot
read_len_plot_out <- (read_len_plot + 
                        geom_hline(yintercept = fraction_remain, color = "black", linetype = "dashed") + 
                        geom_vline(xintercept = min_len_read, color = "black", linetype = "dashed"))

pdf(output_plot, width=10,height=7,paper="special")
print(read_len_plot_out)
dev.off()

print('Done!')
