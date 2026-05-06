#!/usr/bin/env Rscript

library(optparse)
parser = OptionParser()
parser = add_option(parser, c("-o", "--output"), type="character", default=NULL, help="output file")
parser = add_option(parser, c("-i", "--input"), type="character", default=NULL, help="input file")
parser = add_option(parser, c("-r", "--fragcount"), type="character", default=NULL, help="fragment counts file")
parser = add_option(parser, c("-m", "--metadata"), type="character", default=NULL, help="metadata file")
parser = add_option(parser, c("-f", "--factor"), type="character", default=NULL, help="variable to color")

opt = parse_args(parser)

if (any(is.null(c(opt$input, opt$output, opt$metadata)))) {
  print_help(opt_parser)
  stop("Missing required arguments\n", call.=FALSE)
}

library(dplyr)
library(ggplot2)
library(ggrepel)

df <- read.delim(opt$input, header = F, col.names = c("sample", "maprate"))
m_df <- read.delim(opt$metadata, header = T)
f_df <- read.delim(opt$fragcount, header = T)

df <- left_join(df, m_df, by = "sample")
df <- left_join(df, f_df, by = "sample")

pdf(opt$output, width=7, height=7, paper="special")

maprate_plot <- ggplot(df, aes(fragment_count_clean, maprate)) +
  geom_point(aes(color = as.factor(get(opt$factor)))) +
  geom_text_repel(aes(label = sample), 
                  size = 3.0, 
                  segment.alpha = 0.5, 
                  max.overlaps = 5) +
  xlab("Fragment count") +
  ylab("Mapping rate") +
  guides(color=guide_legend(title="Color"))

print(maprate_plot)
dev.off()
