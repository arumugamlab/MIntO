#!/usr/bin/env Rscript

# '''
# Generate gene expression profile from genome-based mode gene profiles
#
# Authors: Carmen Saenz, Mani Arumugam
#
# '''

args = commandArgs(trailingOnly=TRUE)

##########################  ** Load arguments **  ##########################
threads_n <- args[1]
integration_dir <- args[2]
omics <- args[3] #'metaG_metaT'
annot_file <- args[4]
metadata_file <- args[5]
input_file <- args[6]
annot_arg  <- args[7]
annot_names <- unlist(strsplit(annot_arg, split = "\\,")[[1]])
main_factor <- args[8]
print(annot_names)

##########################  ** Load libraries **  ##########################
library(dplyr)
library(data.table)
library(phyloseq)
library(KEGGREST)
library(tibble)
library(stringr)
library(purrr)

#library(unix)

setDTthreads(threads = threads_n)
set.seed(1234)

plot_PCA <- function(physeq, color, label) {

    library(phyloseq)
    library(ggplot2)
    library(ggrepel)

    ##### metadata:
    sample_data_df <- data.frame(sample_data(physeq), stringsAsFactors = F)

    ### Counts - non-zero/inf
    ## Prepare data - replace NA values by 0
    otu_table_df <-  unclass(otu_table(physeq)) %>%
                      as.data.frame(stringAsFactors = F) %>%
                      tibble::rownames_to_column('rownames') %>%
                      dplyr::mutate_all(function(x) ifelse(is.infinite(x) | is.na(x), 0, x)) %>%
                      filter(rowSums(across(where(is.numeric), ~.x != 0))>0) %>%
                      tibble::column_to_rownames('rownames')

    ##### Replace infinite values by log2(min value * 1e-2) for the PCA
    replace_val <- log((min(otu_table_df[otu_table_df >0 & !is.na(otu_table_df)]))*1e-2)

    ### Counts - Log transformed
    otu_table_df_log <- log(otu_table_df) %>%
                          tibble::rownames_to_column('ID') %>%
                          dplyr::mutate_all(function(x) ifelse(is.infinite(x), replace_val, x)) %>%
                          tibble::column_to_rownames('ID')

    #min(otu_table_df_log[!is.na(otu_table_df_log)])

    #Plotting scores of PC1 and PC2 with log transformation
    gene_expression_pca <- prcomp(t(otu_table_df_log), center = T, sca=T)

    # ***************** ggplot way - w coloring *****************
    dtp <- data.frame('sample' = sample_data_df[[label]],
                      'group' = sample_data_df[[color]],
                       gene_expression_pca$x[,1:2]) # the first two componets are selected (NB: you can also select 3 for 3D plottings or 3+)
    df_out <- as.data.frame(gene_expression_pca$x)
    percentage <- round(gene_expression_pca$sdev / sum(gene_expression_pca$sdev) * 100, 2)
    percentage <- paste0( colnames(df_out), " (", paste0(as.character(percentage), "%", ")"))

    PCA_Sample_site_abundance <- (ggplot(data = dtp, aes(x = PC1, y = PC2, color = group)) +
                                  geom_text_repel(aes(label = sample),nudge_x = 0.04, size = 3.5, segment.alpha = 0.5) +
                                  geom_point(size = 2, shape = 16)+
                                  xlab(percentage[1]) + ylab(percentage[2]) +
                                  #labs(title = title) +
                                  theme_bw()+
                                  theme(plot.title = element_text(size=10), legend.position="bottom"))#+ stat_ellipse(type = "norm", linetype = 2))
    return(PCA_Sample_site_abundance)
}

prepare_PCA <- function(physeq, type, color) {

    library(stringr)

    # Colors
    manual_plot_colors =c('#9D0208', '#264653','#e9c46a','#D8DCDE','#B6D0E0',
                          '#FFC87E','#F4A261','#E34F33','#E9C46A',
                          '#A786C9','#D4C0E2','#975773','#6699FF','#000066',
                          '#7AAFCA','#006699','#A9D181','#2F8475','#264445')

    # Title
    title_name <- paste0('PCA - gene ', type)

    # File names
    # abundance  -> FA
    # transcript -> FT
    # expression -> FE
    out_name <- paste0(visual_dir, 'G', str_sub(toupper(type), 1, 1), '.PCA.pdf')

    plot_PCA_out <- plot_PCA(physeq=physeq, color=color, label="sample_alias")
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

#####
# Gene abundance profile
profiles_tpm <- input_file

print('#################################### Paths ####################################')
# Generate output directories ####
dir.create(file.path(integration_dir), showWarnings = FALSE)
visual_dir=paste0(integration_dir,'/plots/')
dir.create(file.path(visual_dir), showWarnings = FALSE)
phyloseq_dir=paste0(integration_dir,'/phyloseq_obj/')
dir.create(file.path(phyloseq_dir), showWarnings = FALSE)

bed_colnames <- c("chr","start","stop","name","score","strand","source","feature","frame","info")

if (omics == 'metaG_metaT'){
    print('#################################### DB profiles ####################################')
    tpm_profile_df <- fread(profiles_tpm, header=T) %>%
                        as.data.frame(stringsAsFactors = F) %>%
                        mutate(info = ifelse(name=='.', coord, info))
    tpm_profile_bed_df <- tpm_profile_df %>%
                            select(coord, all_of(bed_colnames))

    # tpm_profile will remove all rows with 0 sum
    tpm_profile <- tpm_profile_df %>%
                            select(-gene_length, -all_of(bed_colnames)) %>%
                            distinct(coord, .keep_all = TRUE) %>%
                            mutate(across(-coord, ~ as.numeric(replace(., . == '' | is.na(.), 0)))) %>%
                            tibble::column_to_rownames('coord') %>%
                            filter(rowSums(across(where(is.numeric), ~.x != 0))>0)
    rm(tpm_profile_df)
    #rownames(tpm_profile) <- NULL

    ## Bed file ####
    ## Gene abundances
    tpm_profile_bed_df <- tpm_profile_bed_df %>%
                            mutate(source2 = stringr::str_split(source,  '\\,') %>% map_chr(., 1)) %>% # Get the first field
                            mutate(name2   = stringr::str_split(name,    '\\;') %>% map_chr(., 1)) %>%
                            mutate(genome  = stringr::str_split(chr,     '\\|') %>% map_chr(., 1)) %>%
                            mutate(contig  = stringr::str_split(chr,     '\\|') %>% map_chr(., 2)) %>%
                            mutate(ID_gene = info) %>%
                            mutate(name3   = stringr::str_split(ID_gene, '\\|') %>% map_chr(., 2)) %>%
                            mutate(name    = name3)

    # Annotation ####
    ## Replace , by ; in annotation IDs
    profiles_annot_df <- as.data.frame(fread(annot_file,header=T), stringsAsFactors = F, row.names = T) %>%
                            mutate(across(-ID, ~ gsub("\\,", "\\;", .))) %>%
                            rename(ID_gene = ID) %>%
                            select(ID_gene, all_of(unique(annot_names)))

    #### Bed file + Annotations ####

    ## Subset df by colname
    gene_info <- c("coord","chr","start","stop", "name", "source","feature", "ID_gene")
    ## Gene abundances
    tpm_profile_bed_sub_df <- tpm_profile_bed_df %>%
                                select(all_of(gene_info))
    names_list <- names(profiles_annot_df)
    profiles_annot_tpm_df_sub <- merge(profiles_annot_df, tpm_profile_bed_sub_df, by='ID_gene', all.y=T) %>%
                                    tibble::column_to_rownames('coord') %>%
                                    select(name, all_of(names_list))

    # Delete data
    rm(profiles_annot_df, tpm_profile_bed_sub_df, tpm_profile_bed_df)

    # Split df into metaG and metaT ####

    # Subset metaG and metaT - keep only common samples

    ## Subset profile file by metaG and metaT
    metaG_tpm_profile = tpm_profile %>%
                            select(starts_with("metaG")) %>%
                            rename_with(~ gsub("metaG.", "", .x, fixed=TRUE))
    metaT_tpm_profile = tpm_profile %>%
                            select(starts_with("metaT")) %>%
                            rename_with(~ gsub("metaT.", "", .x, fixed=TRUE))
    rm(tpm_profile)

    ## Keep only samples in common meta G and metaT
    metaG_samples <- colnames(metaG_tpm_profile)
    metaT_samples <- colnames(metaT_tpm_profile)
    samples_intersect <- intersect(metaG_samples, metaT_samples)

    # Save GA and GT data ####
    ## Subset metaG and metaT - keep only common samples
    ### Replace column names with the same extention

    # metaG
    metaG_tpm_profile_sub <- metaG_tpm_profile %>%
                                select(all_of(samples_intersect)) %>%
                                tibble::rownames_to_column('info') %>%
                                dplyr::select(info, everything())
    fwrite(metaG_tpm_profile_sub, file=paste0(integration_dir, '/GA.csv'), row.names = F, quote = F)
    metaG_tpm_profile_sub <- metaG_tpm_profile_sub %>%
                                tibble::column_to_rownames('info')

    # metaT
    metaT_tpm_profile_sub <- metaT_tpm_profile %>%
                                select(all_of(samples_intersect)) %>%
                                tibble::rownames_to_column('info') %>%
                                dplyr::select(info, everything())
    fwrite(metaT_tpm_profile_sub, file=paste0(integration_dir, '/GT.csv'), row.names = F, quote = F)
    metaT_tpm_profile_sub <- metaT_tpm_profile_sub %>%
                                tibble::column_to_rownames('info')

    # Delete data
    rm(metaG_tpm_profile, metaT_tpm_profile)


    ################################## GENE EXPRESSION  ##################################
    print('#################################### GENE EXPRESSION  ####################################')

    # Annotations
    profiles_annot_sub <- profiles_annot_tpm_df_sub[rownames(profiles_annot_tpm_df_sub) %in% rownames(metaG_tpm_profile_sub),]
    profiles_annot_sub <- profiles_annot_sub[rownames(metaG_tpm_profile_sub),]

    fwrite(profiles_annot_sub, file=paste0(integration_dir, '/Annotations.csv'), row.names = F, quote = F)


    #### phyloseq object ####
    # Metadata
    if (metadata_file=='None'){
        metadata_df <- data.frame(sample=names(metaG_tpm_profile_sub), group='control', sample_alias=names(metaG_tpm_profile_sub))
    } else{
        metadata_df <- as.data.frame(fread(metadata_file,  header = T), stringsAsFactors = F)
    }
    samp <- sample_data(metadata_df)
    rownames(samp) <- samp$sample_alias

    ## metaG
    metaG_physeq <- phyloseq(otu_table(as.matrix(metaG_tpm_profile_sub), taxa_are_rows = T),
                            tax_table(as.matrix(profiles_annot_sub)), samp)
    saveRDS(metaG_physeq, file = paste0(phyloseq_dir, '/GA.rds'))

    ## metaT
    metaT_physeq <- phyloseq(otu_table(as.matrix(metaT_tpm_profile_sub), taxa_are_rows = T),
                            tax_table(as.matrix(profiles_annot_sub)), samp)
    saveRDS(metaT_physeq, file = paste0(phyloseq_dir, '/GT.rds'))

    # # Normalization
    # When metaG is 0, just get metaT value (replace metaG 0 values by 1) ####
    # We don't do this anymore. These will become NaN

    # Gene Expression
    print(paste("metaG-metaT match:", identical(names(metaG_tpm_profile_sub), names(metaT_tpm_profile_sub))))
    gene_expression_tpm=(metaT_tpm_profile_sub)/(metaG_tpm_profile_sub)
    gene_expression_tpm_noinf <- gene_expression_tpm %>%
                                               tibble::rownames_to_column('info') %>%
                                               dplyr::mutate(across(-info, function(x) ifelse(is.infinite(x) | is.nan(x), NA, x))) %>%
                                               dplyr::select(info, everything())


    # Save GE data ####
    fwrite(gene_expression_tpm_noinf, file=paste0(integration_dir, '/GE.csv'), row.names = F, quote = F)
    gene_expression_tpm_noinf <- gene_expression_tpm_noinf %>% tibble::column_to_rownames('info')

    #### phyloseq object
    GE_physeq <- phyloseq(otu_table(as.matrix(gene_expression_tpm_noinf), taxa_are_rows = T),
                                                     tax_table(as.matrix(profiles_annot_sub)), samp)
    saveRDS(GE_physeq, file = paste0(phyloseq_dir, '/GE.rds'))

    # PCA - metaG ####

    prepare_PCA(physeq=metaG_physeq, type="abundance",  color=main_factor)
    prepare_PCA(physeq=metaT_physeq, type="transcript", color=main_factor)
    prepare_PCA(physeq=GE_physeq,    type="expression", color=main_factor)

    print('Done!')

    # Delete data
    rm(GE_physeq, metaG_physeq, metaT_physeq,
       gene_expression_tpm, metadata_df, profiles_annot_sub)
    rm(metaG_tpm_profile_sub, metaT_tpm_profile_sub, profiles_annot_tpm_df_sub)
    # Free memory
    gc()

} else{
  ######################################################
  # DUP_START
  ######################################################

  print('#################################### DB profiles ####################################')
  tpm_profile_df = as.data.frame(fread(profiles_tpm, header=T), stringsAsFactors = F)
  tpm_profile_df$info[tpm_profile_df$name=='.'] <- tpm_profile_df$coord[tpm_profile_df$name=='.']
  tpm_profile_bed_df <- tpm_profile_df[colnames(tpm_profile_df) %in% c('coord',bed_colnames)]
  tpm_profile_sub_df <- tpm_profile_df[!colnames(tpm_profile_df) %in% c(bed_colnames, 'gene_length')]
  tpm_profile_sub2 <- tpm_profile_sub_df[!duplicated(tpm_profile_sub_df[,'coord']),]
  rownames(tpm_profile_sub2) <- NULL
  tpm_profile <- tpm_profile_sub2 %>%
    mutate(across(-1,  ~ as.numeric(replace(., . == '', 0))))
  tpm_profile[is.na(tpm_profile)] <- 0
  rm(tpm_profile_sub_df)
  rm(tpm_profile_df)
  rm(tpm_profile_sub2)

  ## Bed file ####
  ## Subset df by colname
  gene_info <- c("coord","chr","start","stop", "name", "source","feature")
  ## Gene abundances
  tpm_profile_bed_df$coord <- paste0(tpm_profile_bed_df$chr, '_',tpm_profile_bed_df$start, '_', tpm_profile_bed_df$stop)
  #tpm_profile_bed_df$coord <- gsub('-', '_', tpm_profile_bed_df$coord)
  tpm_profile_bed_df$source2 <- unlist(lapply(strsplit(as.character(tpm_profile_bed_df$source),'\\,'), `[`, 1))
  tpm_profile_bed_df$name2 <- unlist(lapply(strsplit(as.character(tpm_profile_bed_df$name),'\\;'), `[`, 1))
  tpm_profile_bed_df$genome <- unlist(lapply(strsplit(as.character(tpm_profile_bed_df$chr),'\\|'), `[`, 1))
  tpm_profile_bed_df$contig <- unlist(lapply(strsplit(as.character(tpm_profile_bed_df$chr),'\\|'), `[`, 2))

  tpm_profile_bed_df$ID_gene <- tpm_profile_bed_df$info
  tpm_profile_bed_df$name3 <- unlist(lapply(stringr::str_split(pattern = '\\|', string = tpm_profile_bed_df$ID_gene), `[`, 3))
  tpm_profile_bed_df$name3[is.na(tpm_profile_bed_df$name3)] <- unlist(lapply(stringr::str_split(pattern = '\\|', string = tpm_profile_bed_df$ID_gene[is.na(tpm_profile_bed_df$name3)]), `[`, 2))
  tpm_profile_bed_df$name <- tpm_profile_bed_df$name3

  # Annotation ####
  profiles_annot_df <- as.data.frame(fread(annot_file,header=T), stringsAsFactors = F, row.names = T)

  ## Replace , by ; in annotation IDs
  profiles_annot_df[,2:ncol(profiles_annot_df)] <- data.frame(lapply(profiles_annot_df[,2:ncol(profiles_annot_df)], function(x) {
    gsub("\\,", "\\;", x)
  }))

  profiles_annot_df$ID_gene <- profiles_annot_df$ID
  profiles_annot_df$ID <- NULL
  profiles_annot_df <- subset(profiles_annot_df, select=unique(c("ID_gene", c(annot_names))))

  #### Bed file + Annotations ####
  ## Subset df by colname
  gene_info <- c("coord","chr","start","stop", "name", "source","feature", "ID_gene")
  ## Gene abundances
  tpm_profile_bed_sub_df <- tpm_profile_bed_df[colnames(tpm_profile_bed_df) %in% gene_info]
  profiles_annot_tpm_df <- merge(profiles_annot_df, tpm_profile_bed_sub_df, by='ID_gene', all.y=T)
  rownames(profiles_annot_tpm_df) <- profiles_annot_tpm_df$coord
  profiles_annot_tpm_df$coord <- NULL

  names_list <- names(profiles_annot_df)
  profiles_annot_tpm_df_sub <- subset(profiles_annot_tpm_df, select= c('name',names_list))#"KEGG_ko","kofam_KO"
  mappings_ncol <- ncol(profiles_annot_tpm_df_sub)
  # Delete data
  rm(profiles_annot_df, tpm_profile_bed_sub_df, tpm_profile_bed_df)

  ######################################################
  # DUP_END
  ######################################################

  # Subset data - exclude gene abundances = 0 for all samples####
  rownames(tpm_profile) <- tpm_profile$coord
  tpm_profile$coord <- NULL
  omics_tpm_profile = tpm_profile[grepl(omics, colnames(tpm_profile))]
  colnames(omics_tpm_profile) <- gsub(x = colnames(omics_tpm_profile), pattern = paste0(omics,"\\."), replacement = "")
  # Delete data
  rm(tpm_profile)

  ## Delete rows that sum 0 ####
  omics_tpm_profile_rowsum_0 = c(which(rowSums(omics_tpm_profile)==0))
  omics_tpm_profile_rowsum_0_list <- rownames(omics_tpm_profile[omics_tpm_profile_rowsum_0,])

  if (length(omics_tpm_profile_rowsum_0_list) != 0) {
    tpm_profile_rowsum_0 <- c(omics_tpm_profile_rowsum_0_list)
    omics_tpm_profile_rowsum_not_0 <- omics_tpm_profile[! rownames(omics_tpm_profile) %in% tpm_profile_rowsum_0,]
  }else{
    omics_tpm_profile_rowsum_not_0 <- omics_tpm_profile
  }
  # Save gene abundance data ####
  omics_tpm_profile_rowsum_not_0$info <- rownames(omics_tpm_profile_rowsum_not_0)
  omics_tpm_profile_rowsum_not_0 <- as.data.frame(omics_tpm_profile_rowsum_not_0 %>%
                                                    dplyr::select(info, everything()))
  if (omics == 'metaG'){
    omics_label <-'GA'
    omics_label_plot <- 'gene abundance'
  }else if (omics == 'metaT'){
    omics_label <-'GT'
    omics_label_plot <- 'gene transcript'
  }
  write.csv(omics_tpm_profile_rowsum_not_0,paste0(integration_dir, '/', omics_label, '.csv'), row.names = F, quote = F)
  omics_tpm_profile_rowsum_not_0$info <- NULL

  #### Annotations ####
  profiles_annot_sub <- profiles_annot_tpm_df_sub[rownames(profiles_annot_tpm_df_sub) %in% rownames(omics_tpm_profile_rowsum_not_0),]
  # Abundances
  all(rownames(omics_tpm_profile_rowsum_not_0) %in% rownames(profiles_annot_sub))
  ## [1] TRUE
  all(rownames(omics_tpm_profile_rowsum_not_0) == rownames(profiles_annot_sub))
  ## [1] TRUE
  profiles_annot_sub <- profiles_annot_sub[rownames(omics_tpm_profile_rowsum_not_0),]
  all(rownames(omics_tpm_profile_rowsum_not_0) == rownames(profiles_annot_sub))
  ## [1] TRUE

  write.csv(profiles_annot_sub,paste0(integration_dir, '/Annotations.csv'), row.names = F, quote = F)

  # Metadata
  if (metadata_file=='None'){
    metadata_df <- data.frame(sample=names(omics_tpm_profile_rowsum_not_0), group='control', sample_alias=names(omics_tpm_profile_rowsum_not_0))
  } else{
    metadata_df <- as.data.frame(fread(metadata_file,  header = T), stringsAsFactors = F)
  }
  samp <- sample_data(metadata_df)
  rownames(samp) <- samp[["sample_alias"]]

  #### phyloseq object ####
  omics_tpm_profile_rowsum_not_0_phyloseq <- phyloseq(otu_table(as.matrix(omics_tpm_profile_rowsum_not_0), taxa_are_rows = T),
                                                      tax_table(as.matrix(profiles_annot_sub)), samp)
  saveRDS(omics_tpm_profile_rowsum_not_0_phyloseq, file = paste0(phyloseq_dir, '/', omics_label, '.rds'))

  print('#################################### PCA  ####################################')#####
  manual_plot_colors =c('#9D0208', '#264653','#e9c46a','#D8DCDE','#B6D0E0',
                        '#FFC87E','#F4A261','#E34F33','#E9C46A',
                        '#A786C9','#D4C0E2','#975773','#6699FF','#000066',
                        '#7AAFCA','#006699','#A9D181','#2F8475','#264445')
  #distance_lab='bray'
  # PCA ####
  ## Prepare data - replace NA values by 0
  omics_tpm_profile_rowsum_noinf <- as.data.frame(omics_tpm_profile_rowsum_not_0 %>%
                                                    tibble::rownames_to_column('rownames') %>%
                                                    dplyr::mutate_all(function(x) ifelse(is.infinite(x), 0, x)) %>%
                                                    tibble::column_to_rownames('rownames'))
  omics_tpm_profile_rowsum_noinf <- as.data.frame(omics_tpm_profile_rowsum_noinf %>%
                                                    tibble::rownames_to_column('rownames') %>%
                                                    dplyr::mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
                                                    tibble::column_to_rownames('rownames'))

  omics_tpm_profile_rowsum_noinf_plot <- phyloseq(otu_table(as.matrix(omics_tpm_profile_rowsum_noinf), taxa_are_rows = T), tax_table(as.matrix(profiles_annot_sub)), samp)

  # title_name <- paste0('PCoA - ', omics_label_plot)
  # out_name <- paste0(visual_dir, omics_label, '.PCoA-',distance_lab,'.pdf')
  # plot_PCoA_out <- plot_PCoA(distance_lab, omics_tpm_profile_rowsum_noinf_plot, "condition", "sample_alias")
  #
  # pdf(out_name,width=8,height=8,paper="special" )
  # print(plot_PCoA_out  +
  #         scale_color_manual(values=manual_plot_colors, name=main_factor) +
  #         #coord_fixed() +
  #         theme(legend.position="bottom")+
  #         ggtitle(title_name, 'Bray Curtis'))
  # print(plot_PCoA_out  +
  #         facet_wrap(.~condition)+
  #         scale_color_manual(values=manual_plot_colors, name=main_factor) +
  #         #coord_fixed() +
  #         theme(legend.position="bottom")+
  #         ggtitle(title_name, 'Bray Curtis'))
  # dev.off()

  title_name <- paste0('PCA -  ', omics_label_plot)
  out_name <- paste0(visual_dir, omics_label, '.PCA.pdf')
  plot_PCA_out <- plot_PCA(omics_tpm_profile_rowsum_noinf_plot, "condition", "sample_alias")
  pdf(out_name,width=8,height=8,paper="special" )
  print(plot_PCA_out  + scale_color_manual(values=manual_plot_colors, name=main_factor) +
          #coord_fixed() +
          theme(legend.position="bottom")+
          ggtitle(title_name))
  print(plot_PCA_out  +
          facet_wrap(.~condition)+
          scale_color_manual(values=manual_plot_colors, name=main_factor) +
          #coord_fixed() +
          theme(legend.position="bottom")+
          ggtitle(title_name))
  dev.off()

  print('Done!')
}
