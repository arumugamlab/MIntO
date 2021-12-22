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


# Load libraries
library(data.table)
library(dplyr)
library(ggplot2)


# Arguments
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
fraction_remain = as.numeric(args[2])
output_plot = args[3]
output_file = args[4]

read_len_df <- as.data.frame(fread(file=input_file, sep=" "), header=FALSE, stringsAsFactors = F)
colnames(read_len_df) <- c('n_reads', 'len_reads', 'sample')

# Total reads per sample
read_len_df_total_reads = data.frame(read_len_df %>% 
                                       dplyr::group_by(sample) %>% 
                                       dplyr::summarise(total_reads = sum(n_reads)))

# Calculate cumulative reads distribution per bp per sample
read_len_df_cumsum = data.frame(read_len_df %>% 
                                  dplyr::group_by(sample) %>% 
                                  dplyr::arrange(len_reads) %>% # order by read lenght
                                  dplyr::mutate(cumsum = cumsum(n_reads)))

# Combine files to calculate total number of reads
read_len_df_cumsum_total_reads = merge(read_len_df_cumsum, read_len_df_total_reads, by = 'sample')
read_len_df_cumsum_total_reads$cumsum_total=as.numeric(read_len_df_cumsum_total_reads$total_reads -read_len_df_cumsum_total_reads$cumsum)

read_len_df_total_reads$cumsum_total <- read_len_df_total_reads$total_reads
read_len_df_total_reads$n_reads <- read_len_df_total_reads$len_reads <- read_len_df_total_reads$cumsum <- 0

read_len_df_plot <- merge(read_len_df_cumsum_total_reads, read_len_df_total_reads, all=TRUE)
read_len_df_plot$cumsum_total_perc=as.numeric((100*read_len_df_plot$cumsum_total)/read_len_df_plot$total_reads)

# Prepare data to plot
read_len_df_plot$sample_name = gsub('(.*)\\_(.*)',"\\1",read_len_df_plot$sample)
read_len_df_plot$sample_pair  = gsub('(.*)\\_(.*)',"\\2",read_len_df_plot$sample)
read_len_df_plot$fwd_rv <- 'Fwd'
read_len_df_plot$fwd_rv[read_len_df_plot$sample_pair == "2"]<- 'Rev'

# Plot cumulative reads distribution per bp per sample
title = "Read lengths after trimming"
read_len_plot= (ggplot(data=read_len_df_plot, aes(x=len_reads, y=cumsum_total_perc, group = sample_name, color=sample_name)) + 
geom_line()  + 
theme_bw() + 
theme(axis.text = element_text(size = 9.5),legend.position = 'none') + 
labs(x ="read bp (length cutoff)", y="fraction remaining (%)", title = title) + 
facet_grid(fwd_rv~.))

# Calculate MINLEN parameter in Trimmomatic to keep the % of forward-reverse reads specified in the config file
read_len_df_plot$cumsum_total_perc <- read_len_df_plot$cumsum_total_perc

read_len_df_plot_filter <- read_len_df_plot %>% 
mutate_at(vars(cumsum_total_perc), round, 2) %>% 
filter(cumsum_total_perc >= fraction_remain) %>% 
filter(cumsum_total_perc < (fraction_remain+1))

min_len_read = min(read_len_df_plot_filter$len_reads)
if (min_len_read < 50){
  min_len_read = 50
}

# MINLEN parameter in Trimmomatic
print(paste0('Recommended trimming lenght cutoff to keep ', fraction_remain, '% of the reads: ',min_len_read, 'bp'))
cat(paste0("TRIMMOMATIC_minlen: ",min_len_read),file=output_file,sep="\n")

# Output plot
read_len_plot_out <- (read_len_plot + 
geom_hline(yintercept = fraction_remain, color = "black", linetype = "dashed") + 
geom_vline(xintercept = min_len_read, color = "black", linetype = "dashed"))

pdf(output_plot, width=10,height=7,paper="special")
print(read_len_plot_out)
dev.off()

print('Done!')