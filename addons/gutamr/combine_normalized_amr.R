#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(magrittr)
library(tidyr)
library(optparse)
opt_list <- list(
                make_option("--folder", type="character", default=NULL, help="input folder"),
                make_option("--type", type="character", default=NULL, help="analysis type"),
                make_option("--output", type="character", default=NULL, help="output file to write normalized counts", metavar="file"),
                make_option("--min_frag_count", type="integer", default=10, help="minimum mapped fragment-count [default: %default]"),
                make_option("--min_cov", type="integer", default=60, help="minimum template coverage [default: %default]"),
                make_option("--max_id_diff", type="integer", default=5, help="maximum difference between template coverage and identity [default: %default]")
                )
opt <- parse_args(OptionParser(option_list=opt_list))

read_input = F
if ( opt$type == "reads") {
    read_input = T
}

if (! endsWith(opt$folder, "/")) {
  opt$folder <- paste0(opt$folder, "/")
}

# file read loop func
read_tsv_loop <- function(folder_path="", glob_pattern="", skip = 0, usecolnames = T){
  file.names <- dir(folder_path, pattern = glob_pattern)
  dflist_obj <- vector("list", length(file.names))
  
  for(i in 1:length(file.names)){
  	file_id <- tools::file_path_sans_ext(file.names[i])
  	df_tmp <- read_tsv(paste0(folder_path,file.names[i]), skip = skip, col_names = usecolnames, show_col_types = FALSE)
  	df_tmp$file_id <- file_id
    dflist_obj[[i]] <- df_tmp
  }
  tsv_df <- do.call("rbind", dflist_obj)
  return(tsv_df)
}

read_mapstat_loop <- function(folder_path="", glob_pattern=""){
  tsv_df <- read_tsv_loop(folder_path, glob_pattern, skip = 6)
  file.names <- dir(folder_path, pattern = glob_pattern)
  dflist_obj <- vector("list", length(file.names))
  
  for(i in 1:length(file.names)){
    file_id <- tools::file_path_sans_ext(file.names[i])
    df_tmp <- read_tsv(paste0(folder_path,file.names[i]), n_max = 4, col_names = c("X1", "totalFragmentCount"), show_col_types = FALSE) %>% select(-X1)
    df_tmp$file_id <- file_id
    dflist_obj[[i]] <- df_tmp[4,]
  }
  md_df <- do.call("rbind", dflist_obj)
  md_df %<>% mutate(totalFragmentCount=as.numeric(totalFragmentCount)) 
  tsv_df %<>% left_join(md_df, by = "file_id")
  return(tsv_df)
}

##### read input files

# res files for both types
amr_res <- read_tsv_loop(opt$folder, ".*\\.res")
amr_res %<>% rename(Template="#Template")
amr_res %<>% tidyr::extract(Template, into = "ARO_accession"  , regex = "(ARO:\\d+)", remove = F)

amr_res_filter <- amr_res %>% 
  filter((Template_Coverage >= opt$min_cov) & (Template_Coverage <= (200 - opt$min_cov))) %>% filter(Template_Coverage - Template_Identity <= opt$max_id_diff)

# mapstat only for reads
if (read_input) {
  amr_mst <- read_mapstat_loop(opt$folder, ".*\\.mapstat")
  amr_mst %<>% rename(Template="# refSequence")
  amr_mst %<>% tidyr::extract(Template, into = "ARO_accession"  , regex = "(ARO:\\d+)", remove = F)
  
  amr_mst_filter <- left_join(amr_res %>% select(file_id, ARO_accession), amr_mst)
  amr_mst_filter %<>% filter(fragmentCountAln > opt$min_frag_count)
  amr_mst_filter %<>% mutate(fpkm=fragmentCountAln/(refCoveredPositions/1000) * 1/(totalFragmentCount/10**6))
  
  amr_mst_filter_wide <- amr_mst_filter %>% select(file_id, ARO_accession, fpkm) %>% 
  pivot_wider(names_from = ARO_accession, values_from = fpkm, values_fill = 0)
  
  write_tsv(amr_mst_filter_wide, opt$output)
  write_tsv(amr_res_filter, paste0(tools::file_path_sans_ext(opt$output), ".filter.res"))
} else {
  # write filtered res file to disk
  write_tsv(amr_res_filter, opt$output)
}


