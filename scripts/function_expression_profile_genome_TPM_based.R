#!/usr/bin/env Rscript

# '''
# Generate function expression profile from genome-based mode gene profiles (TPM normalized)
#
# Authors: Carmen Saenz
#
# '''
args = commandArgs(trailingOnly=TRUE)

##########################  ** Load arguments **  ########################## 
threads_n <- args[1]
memory_lim <- args[2]
wd <- args[3]
omics <- args[4] #'metaG_metaT'
map_reference <- args[5]
normalization <- args[6]
identity <- args[7]
annot_file <- args[8]
metadata_file <- args[9]
#input_file <- args[10]
minto_dir <- args[10]
annot_arg  <- args[11]
annot_names <- unlist(strsplit(annot_arg, split = "\\,")[[1]])
print(annot_names)

dir_DB <- paste0(wd, "/output/data_integration/", map_reference)
input_file <- paste0(dir_DB, "/",omics,".genes_abundances.p", identity,".bed")

##########################  ** Load libraries **  ########################## 
library(dplyr)
library(data.table)
library(phyloseq)
library(KEGGREST)
library(ggplot2)
library(ggrepel)
library(PFAM.db)

#library(unix)

setDTthreads(threads = threads_n)
set.seed(1234)
##########################  ** Load functions **  ########################## 
plot_PCoA <- function(distance_lab, data_phyloseq, color, label){ #out_name,title_name,
  label2 <- noquote(label)
  #### PCoA
  ord<- ordinate(data_phyloseq, method = "PCoA", distance = distance_lab)
  PcoA_Sample_site_abundance <- plot_ordination(data_phyloseq, ord, color = color, shape=NULL, label=label) #type = type,
  PcoA_Sample_site_abundance <- PcoA_Sample_site_abundance +
    ggtitle(title_name, distance_lab)  +
    geom_point(size = 2.5) +
    theme_bw() +
    theme(legend.position = "top")
  PcoA_Sample_site_abundance$layers[[1]] <- NULL
  PcoA_Sample_site_abundance$layers[[2]] <- NULL
  PcoA_Sample_site_abundance$layers[[3]] <- NULL
  PcoA_Sample_site_abundance_2<- PcoA_Sample_site_abundance +
    geom_text_repel(aes(label = sample_alias), size = 3.0, segment.alpha = 0.5, max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
    geom_point(size = 2, shape = 16)
  PcoA_Sample_site_abundance_2$layers[[1]] <- NULL
  return(PcoA_Sample_site_abundance_2)
}

# # *******************************************************************************************
# # **adonis/adonis2, Permutational Multivariate Analysis of Variance Using Distance Matrix**:
# library(vegan)
# #library(ztable)
# metadata_adonis <- as(sample_data(metaG_tpm_profile_rowsum_not_0_phyloseq), "data.frame")
# set.seed(2)
# adonis_Status<-vegan::adonis(phyloseq::distance(metaG_tpm_profile_rowsum_not_0_phyloseq, method= distance_lab) ~ Condition+time, data = metadata_adonis)
# #ztable(data.frame(adonis_Status$aov.tab),digits = 3)
# adonis_Status_df <- data.frame(adonis_Status$aov.tab)
# adonis_Status_df$Coefficients <- rownames(adonis_Status_df)
# write.table(adonis_Status_df, paste0(dir_diff_analysis_plot, db , '.FE.',filename3,'_',compartment,"_no0in_metaG.PCoA-",distance_lab,"_PERMANOVA.csv"), row.names = F, col.names = T, sep = "\t", quote = F)
# # *******************************************************************************************

# *************************** Use GA and GT raw counts to generate FE profile *************************** ####
# # GA and GT raw counts will be normalized by gene length
# # The genes will be clustered by function IDs and TPM normalized
# # FE profile will be generated: (FT_TPM) / (FA_TPM)
#Raw gene abundances #####

profiles_tpm <- input_file

print('#################################### Paths ####################################')
filename=paste0(omics,".genes_abundances.p", identity,".", normalization)
# Generate output directories ####
integration_dir=paste0(dir_DB, '/',filename)
dir.create(file.path(integration_dir), showWarnings = FALSE)
visual_dir=paste0(integration_dir,'/plots/')
dir.create(file.path(visual_dir), showWarnings = FALSE)
phyloseq_dir=paste0(integration_dir,'/phyloseq_obj/')
dir.create(file.path(phyloseq_dir), showWarnings = FALSE)

if (omics == 'metaG_metaT'){
  print('#################################### DB profiles ####################################')
  tpm_profile_df = as.data.frame(fread(profiles_tpm, header=T), stringsAsFactors = F)
  tpm_profile_df$coord <- paste0(tpm_profile_df$chr, '_',tpm_profile_df$start, '_', tpm_profile_df$stop)
  tpm_profile_df$coord <- gsub('\\-','\\_', tpm_profile_df$coord)
  tpm_profile_df$gene_lenght <- abs(tpm_profile_df$stop- tpm_profile_df$start) +1
  tpm_profile_df <- as.data.frame(tpm_profile_df %>%
                                    dplyr::select(coord, chr,start,stop,name,score,strand,source,feature,frame,info,gene_lenght , everything()))
  tpm_profile_df$info[tpm_profile_df$name=='.'] <- tpm_profile_df$coord[tpm_profile_df$name=='.']
  bed_colnames <- c("chr","start","stop","name","score","strand","source","feature","frame","info")
  tpm_profile_bed_df <- tpm_profile_df[colnames(tpm_profile_df) %in% c('coord',bed_colnames)]
  tpm_profile_sub_df <- tpm_profile_df[!colnames(tpm_profile_df) %in% c(bed_colnames)]
  tpm_profile <- tpm_profile_sub_df[!duplicated(tpm_profile_sub_df[,'coord']),]
  rownames(tpm_profile) <- tpm_profile$coord
  tpm_profile$coord <- NULL
  ## Calculate RPK from gene abundances
  gene_rpk_profile <- tpm_profile[2:ncol(tpm_profile)]/tpm_profile$gene_lenght
  ## RPK sum per sample
  rpk_sum <- as.data.frame(colSums(gene_rpk_profile))
  names(rpk_sum) <- 'total_sum'
  
  gene_rpk_profile$coord <- rownames(gene_rpk_profile)
  tpm_profile <- as.data.frame(gene_rpk_profile %>%
                                 dplyr::select(coord, everything()))
  rm(tpm_profile_sub_df)
  rm(tpm_profile_df)
  rm(gene_rpk_profile)
  rownames(tpm_profile) <- NULL
  
  ## Bed file ####
  ## Subset df by colname
  gene_info <- c("coord","chr","start","stop", "name", "source","feature")
  ## Gene abundances
  tpm_profile_bed_df$coord <- paste0(tpm_profile_bed_df$chr, '_',tpm_profile_bed_df$start, '_', tpm_profile_bed_df$stop)
  tpm_profile_bed_df$coord <- gsub('-', '_', tpm_profile_bed_df$coord)
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
  profiles_annot_df$ID_gene <- profiles_annot_df$ID
  profiles_annot_df$ID <- NULL
  profiles_annot_df <- subset(profiles_annot_df, select=unique(c("ID_gene", annot_names)))
  
  #### Bed file + Annotations ####
  ## Subset df by colname
  gene_info <- c("coord","chr","start","stop", "name", "source","feature", "ID_gene")
  ## Gene counts
  tpm_profile_bed_sub_df <- tpm_profile_bed_df[colnames(tpm_profile_bed_df) %in% gene_info]
  profiles_annot_tpm_df <- merge(profiles_annot_df, tpm_profile_bed_sub_df, by='ID_gene', all.y=T)
  rownames(profiles_annot_tpm_df) <- profiles_annot_tpm_df$coord
  profiles_annot_tpm_df$coord <- NULL
  
  names_list <- names(profiles_annot_df)
  profiles_annot_tpm_df_sub <- subset(profiles_annot_tpm_df, select= c('name',names_list))#"KEGG_ko","kofam_KO"
  # Delete data
  rm(profiles_annot_df, tpm_profile_bed_sub_df, tpm_profile_bed_df, profiles_annot_tpm_df)
  
  # Split df into metaG and metaT ####
  rownames(tpm_profile) <- tpm_profile$coord
  tpm_profile$coord <- NULL
  
  # Subset metaG and metaT - keep only common samples
  ## Subset profile file by metaG and metaT
  metaG_tpm_profile = tpm_profile[grepl("metaG", colnames(tpm_profile))]
  colnames(metaG_tpm_profile) <- gsub(x = colnames(metaG_tpm_profile), pattern = "metaG\\.", replacement = "")
  metaT_tpm_profile = tpm_profile[grepl("metaT", colnames(tpm_profile))]
  colnames(metaT_tpm_profile) <- gsub(x = colnames(metaT_tpm_profile), pattern = "metaT\\.", replacement = "")
  ## Keep only samples in common metaG and metaT
  metaG_samples <- colnames(metaG_tpm_profile)
  metaT_samples <- colnames(metaT_tpm_profile)
  samples_intersect <- intersect(metaG_samples, metaT_samples)
  ## Subset metaG and metaT - keep only common samples
  ### Replace column names with the same extention
  metaG_tpm_profile_sub <- metaG_tpm_profile[, colnames(metaG_tpm_profile) %in% samples_intersect]
  metaT_tpm_profile_sub <- metaT_tpm_profile[, colnames(metaT_tpm_profile) %in% samples_intersect]
  # Delete data
  rm(metaG_tpm_profile, metaT_tpm_profile, tpm_profile)
  
  ## Select gene IDs that are in the GE profile ####
  # GE profile - phyloseq object
  ge_tpm <- readRDS(paste0(phyloseq_dir, '/GE.rds'))
  # Taxa (Expression)
  ge_tpm_taxa <- as.data.frame(unclass(tax_table(ge_tpm)), stringsAsFactors = F)
  ge_tpm_taxa_genes <- rownames(ge_tpm_taxa)
  
  # Metadata
  metadata_df <- data.frame(sample_data(ge_tpm), stringsAsFactors = F)
  samp <- sample_data(metadata_df)
  rownames(samp) <- samp$sample
  
  metaG_tpm_profile_rowsum_not_0 <- metaG_tpm_profile_sub[rownames(metaG_tpm_profile_sub) %in% ge_tpm_taxa_genes,]
  metaT_tpm_profile_rowsum_not_0 <- metaT_tpm_profile_sub[rownames(metaT_tpm_profile_sub) %in% ge_tpm_taxa_genes,]
  
  # Create df to plot number of genes and functions
  ge_fe_df <- data.frame(DB='genes', feature_n=length(ge_tpm_taxa_genes), feature = 'Genes')
  
  #################################### FUNCTION EXPRESSION  ####################################
  print('#################################### FUNCTION EXPRESSION  ####################################')
  manual_plot_colors =c('#9D0208', '#264653','#e9c46a','#D8DCDE','#B6D0E0',
                        '#FFC87E','#F4A261','#E34F33','#E9C46A',
                        '#A786C9','#D4C0E2','#975773','#6699FF','#000066',
                        '#7AAFCA','#006699','#A9D181','#2F8475','#264445') 
  metaG_tpm_profile_rowsum_not_0$coord <- rownames(metaG_tpm_profile_rowsum_not_0)
  rownames(metaG_tpm_profile_rowsum_not_0) <- NULL
  metaT_tpm_profile_rowsum_not_0$coord <- rownames(metaT_tpm_profile_rowsum_not_0)
  rownames(metaT_tpm_profile_rowsum_not_0) <- NULL
  
  profiles_annot_tpm_df_sub$coord <- rownames(profiles_annot_tpm_df_sub)
  profiles_annot_tpm_df_sub <- as.data.frame(profiles_annot_tpm_df_sub %>%
                                               dplyr::select(coord, ID_gene, name, everything()))
  profiles_annot_tpm_df_sub <- profiles_annot_tpm_df_sub[!colnames(profiles_annot_tpm_df_sub) %in% c("KEGG_ko", "kofam_KO")]
  #profiles_annot_tpm_df_sub <- profiles_annot_tpm_df_sub[colnames(profiles_annot_tpm_df_sub) %in% c("coord","ID_gene" ,"dbCAN.enzclass")]
  
  profiles_annot_tpm_df_sub <- profiles_annot_tpm_df_sub[profiles_annot_tpm_df_sub$coord %in% metaG_tpm_profile_rowsum_not_0$coord,]
  mappings_ncol <- ncol(profiles_annot_tpm_df_sub)
  for(db in colnames(profiles_annot_tpm_df_sub)[4:mappings_ncol]){#change it deppending on the number of selectedMappings columns
    print(db)
    message(Sys.time(), ": Running ", db)
    dbMap <- profiles_annot_tpm_df_sub[,c("coord", "ID_gene","name",db)]
    dbMap <- dbMap[!is.na(dbMap[[db]]),]
    dbMap <- dbMap[dbMap[[db]] != "",]
    dbMap <- dbMap[dbMap[[db]] != "-",]
    dbMap <- dbMap[rowSums(is.na(dbMap)) != ncol(dbMap),]
    dbkeys <- unique(unlist(strsplit(profiles_annot_tpm_df_sub[,db], split = "\\; ")))
    dbkeys <- unique(unlist(strsplit(dbkeys, split = "\\,")))
    singleKeys <- strsplit(dbMap[,4], split = "\\; |\\;|\\, |\\,")
    db_name <- gsub(' ', '.', db)
    if (length(singleKeys) > 0){
      keyMap <- data.frame(
        coord = rep(dbMap$coord, sapply(singleKeys, length)),
        ID_gene = rep(dbMap$ID_gene, sapply(singleKeys, length)),
        name = rep(dbMap$name, sapply(singleKeys, length)),
        Funct = unlist(singleKeys),
        stringsAsFactors = FALSE
      )
      if (db == "KEGG_Pathway"){
        dim(keyMap)
        keyMap$Funct<- gsub('ko', '',keyMap$Funct)
        keyMap$Funct<- gsub('map', '',keyMap$Funct)
        keyMap <- keyMap[!duplicated(keyMap),]
        keyMap$Funct <- paste0('map', keyMap$Funct)
        dim(keyMap)
      }
      keyMap_u <- keyMap[!duplicated(keyMap[c(1,2,3,4)]),]
      message(Sys.time(), ": ", paste(db, "processed!"))
      keyMap_funct <- data.frame(Funct = c(unique(keyMap_u$Funct)))
      
      ### Annotation descriptions ####
      # library(purrr)
      # library(magrittr)
      if (db_name == 'eggNOG.OGs'){
        cog_df <- as.data.frame(fread(paste0(minto_dir,'/data/cog-20.def.tab'), sep ="\t", header = F, stringsAsFactors = F))
        names(cog_df) <- c('ID_COG', 'funct_category_COG', 'name_COG', 'name_gene', 'funct_pathway', 'ID_PubMed', 'ID_PDB')
        #cog_df <- subset(cog_df, select=c('ID_COG', 'name_COG', 'funct_pathway'))
        cog_df <- subset(cog_df, select=c('ID_COG', 'name_COG'))
        names(cog_df) <- c('Funct', 'Description')
        annot_df <- merge(keyMap_funct, cog_df, by = 'Funct', all.x = T)
      } else if  (db_name == 'merged_KO'){
        annot_df <- data.frame(Description=unlist(lapply(keyMap_funct$Funct, function (x) keggFind('ko', x))),stringsAsFactors = F)
        annot_df$Funct <- rownames(annot_df)
        annot_df$Funct <- gsub('ko:', '', annot_df$Funct)
      }else if  (db_name == 'KEGG_Pathway'){
        annot_df <- data.frame(Description=unlist(lapply(keyMap_funct$Funct, function (x) keggFind('pathway', x))),stringsAsFactors = F)
        annot_df$Funct <- rownames(annot_df)
        annot_df$Funct <- gsub('path:', '', annot_df$Funct)
      }else if  (db_name == 'KEGG_Module'){
        modules_list <- as.data.frame(fread(paste0(minto_dir,'/data/KEGG_Modules_20171212.csv'), header=T), stringsAsFactors = F)
        modules_list_sub <- subset(modules_list, select=c('Module', 'Definition'))
        names(modules_list_sub) <- c('Funct', 'Description')
        annot_df <- merge(keyMap_funct, modules_list_sub, by='Funct', all=T)
      }else if  (db_name %in% c('PFAMs')){
        x <- PFAMID
        mapped_keys <- mappedkeys(x)
        xx <- as.list(x[mapped_keys])
        test <- as.data.frame(do.call(rbind, xx))
        names(test) <- 'Pfam_family'
        test$Pfam_id <- rownames(test)
        x <- PFAMDE
        mapped_keys <- mappedkeys(x)
        xx <- as.list(x[mapped_keys])
        test2 <- as.data.frame(do.call(rbind, xx))
        names(test2) <- 'Description'
        test2$Pfam_id <- rownames(test2)
        test3 <- merge(test, test2, by='Pfam_id')
        annot_df <- merge(keyMap_funct, test3, by.x='Funct' , by.y='Pfam_family', all.x=T)
        annot_df$Pfam_id <- NULL
        rm(test, test2, test3, x, xx)
      } else if  (db_name %in% c('dbCAN.mod', 'dbCAN.enzclass', 'CAZy')){
        x <- PFAMCAZY
        mapped_keys <- mappedkeys(x)
        xx <- as.list(x[mapped_keys])
        test <- as.data.frame(do.call(rbind, xx))
        #data.frame(matrix((xx), nrow=length(xx), byrow=F))
        #names(test) <- 'CAZY_id'
        test$Pfam_id <- rownames(test)
        test_cazy <- as.data.frame(test %>% group_by(Pfam_id) %>% 
                                     mutate(CAZY_id = paste(unique(c(V1, V2, V3, V4)), collapse = ',')))
        test_cazy$V1 <- test_cazy$V2 <- test_cazy$V3 <- test_cazy$V4 <- NULL
        test_cazy.singleKeys <- strsplit(test_cazy$CAZY_id, split = "\\  |\\,")
        if (length(test_cazy.singleKeys) > 0){
          test_cazy.keyMap <- data.frame(
            Pfam_id = rep(test_cazy$Pfam_id, sapply(test_cazy.singleKeys, length)),
            Funct = unlist(test_cazy.singleKeys),
            stringsAsFactors = FALSE
          )
        }
        x <- PFAMDE
        mapped_keys <- mappedkeys(x)
        xx <- as.list(x[mapped_keys])
        test2 <- as.data.frame(do.call(rbind, xx))
        names(test2) <- 'Description'
        test2$Pfam_id <- rownames(test2)
        test3 <- merge(test_cazy.keyMap, test2, by='Pfam_id')
        annot_df <- merge(keyMap_funct, test3, by='Funct', all.x=T)
        annot_df$Pfam_id <- NULL
        rm(test, test2, test3, test_cazy, test_cazy.keyMap, test_cazy.singleKeys, x, xx)
      } else {
        annot_df <- keyMap_funct
        annot_df$Description <- '-'
      }
      
      #### metaG ####
      metaG_tpm_profile_rowsum_not_0_annot <- merge(keyMap_u, metaG_tpm_profile_rowsum_not_0, by='coord', all=T)
      metaG_tpm_profile_rowsum_not_0_annot$coord <- NULL
      metaG_tpm_profile_rowsum_not_0_annot$ID_gene <- NULL
      metaG_tpm_profile_rowsum_not_0_annot$name <- NULL
      
      metaG_funct_profile_rowsum_not_0_annot_dplyr <- as.data.frame(metaG_tpm_profile_rowsum_not_0_annot %>%
                                                                      dplyr::group_by(Funct) %>%
                                                                      dplyr::summarise(across(everything(),sum)))
      metaG_funct_profile_rowsum_not_0_annot_dplyr$Funct[which(is.na(metaG_funct_profile_rowsum_not_0_annot_dplyr$Funct))] <- 'Unknown'
      names(metaG_funct_profile_rowsum_not_0_annot_dplyr) <- gsub('X', '', names(metaG_funct_profile_rowsum_not_0_annot_dplyr))
      
      # Calculate TPM from Function abundances (normalized to gene length) ####
      rownames(metaG_funct_profile_rowsum_not_0_annot_dplyr) <- metaG_funct_profile_rowsum_not_0_annot_dplyr$Funct
      metaG_funct_profile_rowsum_not_0_annot_dplyr$Funct <- NULL
      metaG_funct_profile_rowsum_not_0_annot_dplyr[is.na(metaG_funct_profile_rowsum_not_0_annot_dplyr)] <- 0
      
      funct_rpk_samples_u_t_df <- as.data.frame(t(metaG_funct_profile_rowsum_not_0_annot_dplyr))
      ## RPK sum per sample
      rpk_sum <- as.data.frame(colSums(metaG_funct_profile_rowsum_not_0_annot_dplyr))
      names(rpk_sum) <- 'total_sum'
      
      funct_rpk_samples_total_reads_t_df <- cbind(funct_rpk_samples_u_t_df, rpk_sum)
      # ## Calculate TPM from RPK and SCALING FACTOR
      # gene_tpm_samples_u_t_df <- gene_rpk_samples_total_reads_t_df[1:ncol(gene_rpk_samples_total_reads_t_df)-1]/gene_rpk_samples_total_reads_t_df$scaling_factor
      
      ## Calculate TPM from RPK and RPK sum per sample ####
      funct_tpm_samples_u_t_df <- funct_rpk_samples_total_reads_t_df[1:ncol(funct_rpk_samples_total_reads_t_df)-1]/funct_rpk_samples_total_reads_t_df$total_sum
      funct_tpm_samples_u_t_df <- funct_tpm_samples_u_t_df*1e6
      
      metaG_funct_profile_rowsum_not_0_annot_dplyr <- as.data.frame(t(funct_tpm_samples_u_t_df))
      #metaG_funct_profile_rowsum_not_0_annot_dplyr[is.na(metaG_funct_profile_rowsum_not_0_annot_dplyr)] <- 0
      
      colSums(metaG_funct_profile_rowsum_not_0_annot_dplyr)
      
      metaG_funct_profile_rowsum_not_0_annot_dplyr$Funct <- rownames(metaG_funct_profile_rowsum_not_0_annot_dplyr)
      metaG_funct_profile_rowsum_not_0_annot_dplyr <- as.data.frame(metaG_funct_profile_rowsum_not_0_annot_dplyr %>%
                                                                      dplyr::select(Funct, everything()))
      metaG_funct_profile_rowsum_not_0_annot_desc <- merge(annot_df, metaG_funct_profile_rowsum_not_0_annot_dplyr, by='Funct', all.y=T)
      write.csv(metaG_funct_profile_rowsum_not_0_annot_desc,paste0(integration_dir, '/FA.', db_name ,'.csv'), row.names = F, quote = F)
      
      #### phyloseq object ####
      rownames(metaG_funct_profile_rowsum_not_0_annot_dplyr) <- metaG_funct_profile_rowsum_not_0_annot_dplyr$Funct
      metaG_funct_profile_rowsum_not_0_annot_dplyr$Funct <- NULL
      
      ## metadata
      #metadata_df <- data.frame(sample=names(metaG_funct_profile_rowsum_not_0_annot_dplyr), sample_color='sample')
      samp <- sample_data(metadata_df)
      rownames(samp) <- samp$sample
      
      # taxa
      metaG_annot_desc <- subset(metaG_funct_profile_rowsum_not_0_annot_desc, select=c("Funct","Description"))
      rownames(metaG_annot_desc) <- metaG_annot_desc$Funct
      
      metaG_function_expression_phyloseq <- phyloseq(otu_table(as.matrix(metaG_funct_profile_rowsum_not_0_annot_dplyr), 
                                                               taxa_are_rows = T), tax_table(as.matrix(metaG_annot_desc)), samp)
      head(taxa_names(metaG_function_expression_phyloseq))
      saveRDS(metaG_function_expression_phyloseq, file = paste0(phyloseq_dir, '/FA.', db_name ,'.rds'))
      
      #### metaT ####
      metaT_tpm_profile_rowsum_not_0_annot <- merge(keyMap_u, metaT_tpm_profile_rowsum_not_0, by='coord', all=T)
      metaT_tpm_profile_rowsum_not_0_annot$coord <- NULL
      metaT_tpm_profile_rowsum_not_0_annot$ID_gene <- NULL
      metaT_tpm_profile_rowsum_not_0_annot$name <- NULL
      
      metaT_funct_profile_rowsum_not_0_annot_dplyr <- as.data.frame(metaT_tpm_profile_rowsum_not_0_annot %>%
                                                                      dplyr::group_by(Funct) %>%
                                                                      dplyr::summarise(across(everything(),sum)))
      metaT_funct_profile_rowsum_not_0_annot_dplyr$Funct[which(is.na(metaT_funct_profile_rowsum_not_0_annot_dplyr$Funct))] <- 'Unknown'
      
      names(metaT_funct_profile_rowsum_not_0_annot_dplyr) <- gsub('X', '', names(metaT_funct_profile_rowsum_not_0_annot_dplyr))
      
      # Calculate TPM from Function abundances (normalized to gene length) ####
      rownames(metaT_funct_profile_rowsum_not_0_annot_dplyr) <- metaT_funct_profile_rowsum_not_0_annot_dplyr$Funct
      metaT_funct_profile_rowsum_not_0_annot_dplyr$Funct <- NULL
      metaT_funct_profile_rowsum_not_0_annot_dplyr[is.na(metaT_funct_profile_rowsum_not_0_annot_dplyr)] <- 0
      
      funct_rpk_samples_u_t_df <- as.data.frame(t(metaT_funct_profile_rowsum_not_0_annot_dplyr))
      ## RPK sum per sample
      rpk_sum <- as.data.frame(colSums(metaT_funct_profile_rowsum_not_0_annot_dplyr))
      names(rpk_sum) <- 'total_sum'
      
      funct_rpk_samples_total_reads_t_df <- cbind(funct_rpk_samples_u_t_df, rpk_sum)
      # ## Calculate TPM from RPK and SCALING FACTOR
      # gene_tpm_samples_u_t_df <- gene_rpk_samples_total_reads_t_df[1:ncol(gene_rpk_samples_total_reads_t_df)-1]/gene_rpk_samples_total_reads_t_df$scaling_factor
      
      ## Calculate TPM from RPK and RPK sum per sample ####
      funct_tpm_samples_u_t_df <- funct_rpk_samples_total_reads_t_df[1:ncol(funct_rpk_samples_total_reads_t_df)-1]/funct_rpk_samples_total_reads_t_df$total_sum
      funct_tpm_samples_u_t_df <- funct_tpm_samples_u_t_df*1e6
      
      metaT_funct_profile_rowsum_not_0_annot_dplyr <- as.data.frame(t(funct_tpm_samples_u_t_df))
      #metaT_funct_profile_rowsum_not_0_annot_dplyr[is.na(metaT_funct_profile_rowsum_not_0_annot_dplyr)] <- 0
      
      colSums(metaT_funct_profile_rowsum_not_0_annot_dplyr)
      
      metaT_funct_profile_rowsum_not_0_annot_dplyr$Funct <- rownames(metaT_funct_profile_rowsum_not_0_annot_dplyr)
      metaT_funct_profile_rowsum_not_0_annot_dplyr <- as.data.frame(metaT_funct_profile_rowsum_not_0_annot_dplyr %>%
                                                                      dplyr::select(Funct, everything()))
      metaT_funct_profile_rowsum_not_0_annot_desc <- merge(annot_df, metaT_funct_profile_rowsum_not_0_annot_dplyr, by='Funct', all.y=T)
      write.csv(metaT_funct_profile_rowsum_not_0_annot_desc,paste0(integration_dir, '/FT.', db_name ,'.csv'), row.names = F, quote = F)
      
      #### phyloseq object ####
      rownames(metaT_funct_profile_rowsum_not_0_annot_dplyr) <- metaT_funct_profile_rowsum_not_0_annot_dplyr$Funct
      metaT_funct_profile_rowsum_not_0_annot_dplyr$Funct <- NULL
      
      ## metadata
      #metadata_df <- data.frame(sample=names(metaT_funct_profile_rowsum_not_0_annot_dplyr), sample_color='sample')
      samp <- sample_data(metadata_df)
      rownames(samp) <- samp$sample
      
      # taxa
      metaT_annot_desc <- subset(metaT_funct_profile_rowsum_not_0_annot_desc, select=c("Funct","Description"))
      rownames(metaT_annot_desc) <- metaT_annot_desc$Funct
      
      metaT_function_expression_phyloseq <- phyloseq(otu_table(as.matrix(metaT_funct_profile_rowsum_not_0_annot_dplyr), taxa_are_rows = T), tax_table(as.matrix(metaT_annot_desc)), samp)
      head(taxa_names(metaT_function_expression_phyloseq))
      saveRDS(metaT_function_expression_phyloseq, file = paste0(phyloseq_dir, '/FT.', db_name ,'.rds'))
      
      # Expression - When metaG is 0, just get metaT value (replace metaG 0 values by 1) ####
      metaG_funct_profile_rowsum_not_0_1 <- metaG_funct_profile_rowsum_not_0_annot_dplyr
      metaG_funct_profile_rowsum_not_0_1[metaG_funct_profile_rowsum_not_0_1==0] <- 1
      
      metaG_funct_profile_pseud <- metaG_funct_profile_rowsum_not_0_1
      metaT_funct_profile_pseud <- metaT_funct_profile_rowsum_not_0_annot_dplyr
      
      # # metaG and metaT matrix same function order
      function_ids = rownames(metaG_funct_profile_pseud)
      relat_num_transc_cell_funct = metaT_funct_profile_pseud[match(rownames(metaG_funct_profile_pseud), rownames(metaT_funct_profile_pseud)),]
      relat_num_functions_cell_funct = metaG_funct_profile_pseud
      # Function Expression
      function_expression_tpm=(relat_num_transc_cell_funct)/(relat_num_functions_cell_funct)
      function_expression_tpm <- as.data.frame(function_expression_tpm %>%
                                                 tibble::rownames_to_column('rownames') %>%
                                                 dplyr::mutate_all(function(x) ifelse(is.infinite(x), NA, x)) %>%
                                                 tibble::column_to_rownames('rownames'))
      function_expression_tpm_rowsum_0 = c(c(which(rowSums(function_expression_tpm)==0)), c(which(is.na(rowSums(function_expression_tpm)))))
      function_expression_tpm_rowsum_0_list <- rownames(function_expression_tpm[function_expression_tpm_rowsum_0,])
      function_expression_tpm_noNA <- function_expression_tpm[! rownames(function_expression_tpm) %in% function_expression_tpm_rowsum_0_list,]
      rm(function_expression_tpm)
      function_expression_tpm <- function_expression_tpm_noNA
      
      function_expression_tpm$Funct <- rownames(function_expression_tpm)
      function_expression_tpm <- as.data.frame(function_expression_tpm %>%
                                                 dplyr::select(Funct, everything()))
      
      function_expression_tpm_desc <- merge(annot_df, function_expression_tpm, by='Funct', all.y=T)
      write.csv(function_expression_tpm_desc,paste0(integration_dir, '/FE.', db_name ,'.csv'), row.names = F, quote = F)
      
      #### phyloseq object - metaT relative to metaG -####
      rownames(function_expression_tpm) <- function_expression_tpm$Funct
      function_expression_tpm$Funct <- NULL
      #min(function_expression_tpm[function_expression_tpm >0 & !is.na(function_expression_tpm)])
      
      ## metadata
      #metadata_df <- data.frame(sample=names(function_expression_tpm), sample_color='sample')
      samp <- sample_data(metadata_df)
      rownames(samp) <- samp$sample
      
      # taxa
      FE_annot_desc <- subset(function_expression_tpm_desc, select=c("Funct","Description"))
      rownames(FE_annot_desc) <- FE_annot_desc$Funct
      
      function_expression_phyloseq <- phyloseq(otu_table(as.matrix(function_expression_tpm), taxa_are_rows = T), tax_table(as.matrix(FE_annot_desc)), samp)
      head(taxa_names(function_expression_phyloseq))
      saveRDS(function_expression_phyloseq, file = paste0(phyloseq_dir, '/FE.', db_name ,'.rds'))
      
      # Create df to plot number of genes and functions ####
      ge_fe_df <- rbind(ge_fe_df, c(db_name, nrow(FE_annot_desc[FE_annot_desc$Funct!='Unknown',]), 'Functions'))
      
      print('#################################### PCoA  ####################################')#####
      distance_lab='bray'
      # PCoA - metaG ####
      ## Prepare data - replace NA values by 0
      metaG_function_expression_noinf <- as.data.frame(metaG_funct_profile_rowsum_not_0_annot_dplyr %>%
                                                         tibble::rownames_to_column('rownames') %>%
                                                         dplyr::mutate_all(function(x) ifelse(is.infinite(x), 0, x)) %>%
                                                         tibble::column_to_rownames('rownames'))
      metaG_function_expression_noinf <- as.data.frame(metaG_function_expression_noinf %>%
                                                         tibble::rownames_to_column('rownames') %>%
                                                         dplyr::mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
                                                         tibble::column_to_rownames('rownames'))
      
      metaG_function_expression_noinf_plot <- phyloseq(otu_table(as.matrix(metaG_function_expression_noinf), taxa_are_rows = T), tax_table(as.matrix(metaG_annot_desc)), samp)
      
      title_name <- paste0('PCoA - function abundance - ', db_name)
      out_name <- paste0(visual_dir, 'FA.', db_name ,'.PCoA-',distance_lab,'.pdf')
      plot_PCoA_out <- plot_PCoA(distance_lab, metaG_function_expression_noinf_plot, "condition", "sample_alias")
      pdf(out_name,width=8,height=8,paper="special" )
      print(plot_PCoA_out  + scale_color_manual(values=manual_plot_colors, name='Condition') +
              #coord_fixed() +
              theme(legend.position="bottom")+ggtitle(title_name, 'Bray Curtis'))
      print(plot_PCoA_out  + 
              facet_wrap(.~condition)+
              scale_color_manual(values=manual_plot_colors, name='Condition') +
              #coord_fixed() +
              theme(legend.position="bottom")+
              ggtitle(title_name, 'Bray Curtis'))
      dev.off()
      
      # PCoA - metaT ####
      ## Prepare data - replace NA values by 0
      metaT_function_expression_noinf <- as.data.frame(metaT_funct_profile_rowsum_not_0_annot_dplyr %>%
                                                         tibble::rownames_to_column('rownames') %>%
                                                         dplyr::mutate_all(function(x) ifelse(is.infinite(x), 0, x)) %>%
                                                         tibble::column_to_rownames('rownames'))
      metaT_function_expression_noinf <- as.data.frame(metaT_function_expression_noinf %>%
                                                         tibble::rownames_to_column('rownames') %>%
                                                         dplyr::mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
                                                         tibble::column_to_rownames('rownames'))
      
      metaT_function_expression_noinf_plot <- phyloseq(otu_table(as.matrix(metaT_function_expression_noinf), taxa_are_rows = T), tax_table(as.matrix(metaT_annot_desc)), samp)
      
      title_name <- paste0('PCoA - function transcript - ', db_name)
      out_name <- paste0(visual_dir, 'FT.', db_name ,'.PCoA-',distance_lab,'.pdf')
      plot_PCoA_out <- plot_PCoA(distance_lab, metaT_function_expression_noinf_plot, "condition", "sample_alias")
      pdf(out_name,width=8,height=8,paper="special" )
      print(plot_PCoA_out  + scale_color_manual(values=manual_plot_colors, name='Condition') +
              #coord_fixed() +
              theme(legend.position="bottom")+ggtitle(title_name, 'Bray Curtis'))
      print(plot_PCoA_out  + 
              facet_wrap(.~condition)+
              scale_color_manual(values=manual_plot_colors, name='Condition') +
              #coord_fixed() +
              theme(legend.position="bottom")+
              ggtitle(title_name, 'Bray Curtis'))
      dev.off()
      
      # PCoA - GE ####
      ## Prepare data - replace NA values by 0
      function_expression_noinf <- as.data.frame(function_expression_tpm %>%
                                                   tibble::rownames_to_column('rownames') %>%
                                                   dplyr::mutate_all(function(x) ifelse(is.infinite(x), 0, x)) %>%
                                                   tibble::column_to_rownames('rownames'))
      function_expression_noinf <- as.data.frame(function_expression_noinf %>%
                                                   tibble::rownames_to_column('rownames') %>%
                                                   dplyr::mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
                                                   tibble::column_to_rownames('rownames'))
      
      function_expression_noinf_plot <- phyloseq(otu_table(as.matrix(function_expression_noinf), taxa_are_rows = T), tax_table(as.matrix(FE_annot_desc)), samp)
      title_name <- paste0('PCoA - function expression - ', db_name)
      out_name <- paste0(visual_dir, 'FE.', db_name ,'.PCoA-',distance_lab,'.pdf')
      plot_PCoA_out <- plot_PCoA(distance_lab, function_expression_noinf_plot, "condition", "sample_alias")
      pdf(out_name,width=8,height=8,paper="special" )
      print(plot_PCoA_out  + scale_color_manual(values=manual_plot_colors, name='Condition') +
              #coord_fixed() +
              theme(legend.position="bottom")+ggtitle(title_name, 'Bray Curtis'))
      print(plot_PCoA_out  + 
              facet_wrap(.~condition)+
              scale_color_manual(values=manual_plot_colors, name='Condition') +
              #coord_fixed() +
              theme(legend.position="bottom")+
              ggtitle(title_name, 'Bray Curtis'))
      dev.off()
    } else{
      # create .csv, .rds, .pdf
      distance_lab='bray'
      file.create(paste0(integration_dir, '/FA.', db_name ,'.csv'))
      file.create(paste0(integration_dir, '/FT.', db_name ,'.csv'))
      file.create(paste0(integration_dir, '/FE.', db_name ,'.csv'))
      file.create(paste0(phyloseq_dir, 'FA.', db_name ,'.rds'))
      file.create(paste0(phyloseq_dir, 'FT.', db_name ,'.rds'))
      file.create(paste0(phyloseq_dir, 'FE.', db_name ,'.rds'))
      file.create(paste0(visual_dir, 'FA.', db_name ,'.PCoA-',distance_lab,'.pdf'))
      file.create(paste0(visual_dir, 'FT.', db_name ,'.PCoA-',distance_lab,'.pdf'))
      file.create(paste0(visual_dir, 'FE.', db_name ,'.PCoA-',distance_lab,'.pdf'))
      ge_fe_df <- rbind(ge_fe_df, c(db_name, 0, 'Functions'))
    }
  }
  write.csv(ge_fe_df,paste0(integration_dir, '/GE_FE_features.csv'), row.names = F, quote = F)
  
  ge_fe_df<- as.data.frame(fread(paste0(integration_dir,'/GE_FE_features.csv'), header=T), stringsAsFactors = F, row.names = T)
  ge_fe_df$feature_n <- as.numeric(ge_fe_df$feature_n)
  feature_order <- ge_fe_df$DB[order(ge_fe_df$feature_n)]
  ge_fe_df$DB <- factor(ge_fe_df$DB, levels=feature_order)
  
  pdf(paste0(visual_dir,'/GE_FE_features.pdf'),width=6,height=5,paper="special" )
  print(ggplot(data=ge_fe_df, aes(x=DB, y=feature_n)) +
          geom_bar(stat="identity")+ theme_minimal()+
          facet_wrap(feature~., scales= "free") + labs(y='Number of features', x='')+
          theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))+
          geom_text(aes(label=feature_n), position=position_dodge(width=0.9), vjust=-0.25, size = 3)
        
  )
  dev.off()
} else {
  print('#################################### DB profiles ####################################')
  tpm_profile_df = as.data.frame(fread(profiles_tpm, header=T), stringsAsFactors = F)
  tpm_profile_df$coord <- paste0(tpm_profile_df$chr, '_',tpm_profile_df$start, '_', tpm_profile_df$stop)
  tpm_profile_df$coord <- gsub('\\-','\\_', tpm_profile_df$coord)
  tpm_profile_df$gene_lenght <- abs(tpm_profile_df$stop- tpm_profile_df$start) +1
  tpm_profile_df <- as.data.frame(tpm_profile_df %>%
                                    dplyr::select(coord, chr,start,stop,name,score,strand,source,feature,frame,info,gene_lenght , everything()))
  tpm_profile_df$info[tpm_profile_df$name=='.'] <- tpm_profile_df$coord[tpm_profile_df$name=='.']
  bed_colnames <- c("chr","start","stop","name","score","strand","source","feature","frame","info")
  tpm_profile_bed_df <- tpm_profile_df[colnames(tpm_profile_df) %in% c('coord',bed_colnames)]
  tpm_profile_sub_df <- tpm_profile_df[!colnames(tpm_profile_df) %in% c(bed_colnames)]
  tpm_profile <- tpm_profile_sub_df[!duplicated(tpm_profile_sub_df[,'coord']),]
  rownames(tpm_profile) <- tpm_profile$coord
  tpm_profile$coord <- NULL
  ## Calculate RPK from gene abundances
  gene_rpk_profile <- tpm_profile[2:ncol(tpm_profile)]/tpm_profile$gene_lenght
  ## RPK sum per sample
  rpk_sum <- as.data.frame(colSums(gene_rpk_profile))
  names(rpk_sum) <- 'total_sum'
  
  gene_rpk_profile$coord <- rownames(gene_rpk_profile)
  tpm_profile <- as.data.frame(gene_rpk_profile %>%
                                 dplyr::select(coord, everything()))
  rm(tpm_profile_sub_df)
  rm(tpm_profile_df)
  rm(gene_rpk_profile)
  rownames(tpm_profile) <- NULL
  
  ## Bed file ####
  ## Subset df by colname
  gene_info <- c("coord","chr","start","stop", "name", "source","feature")
  ## Gene abundances
  tpm_profile_bed_df$coord <- paste0(tpm_profile_bed_df$chr, '_',tpm_profile_bed_df$start, '_', tpm_profile_bed_df$stop)
  tpm_profile_bed_df$coord <- gsub('-', '_', tpm_profile_bed_df$coord)
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
  profiles_annot_df$ID_gene <- profiles_annot_df$ID
  profiles_annot_df$ID <- NULL
  profiles_annot_df <- subset(profiles_annot_df, select=unique(c("ID_gene", annot_names)))
  
  #### Bed file + Annotations ####
  ## Subset df by colname
  gene_info <- c("coord","chr","start","stop", "name", "source","feature", "ID_gene")
  ## Gene counts
  tpm_profile_bed_sub_df <- tpm_profile_bed_df[colnames(tpm_profile_bed_df) %in% gene_info]
  profiles_annot_tpm_df <- merge(profiles_annot_df, tpm_profile_bed_sub_df, by='ID_gene', all.y=T)
  rownames(profiles_annot_tpm_df) <- profiles_annot_tpm_df$coord
  profiles_annot_tpm_df$coord <- NULL
  
  names_list <- names(profiles_annot_df)
  profiles_annot_tpm_df_sub <- subset(profiles_annot_tpm_df, select= c('name',names_list))#"KEGG_ko","kofam_KO"
  # Delete data
  rm(profiles_annot_df, tpm_profile_bed_sub_df, tpm_profile_bed_df, profiles_annot_tpm_df)
  
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
  # Gene abundance data ####
  if (omics == 'metaG'){
    omics_label <-'FA'
    omics_label2 <- 'GA'
    omics_label_plot <- 'function abundance'
  }else if (omics == 'metaT'){
    omics_label <-'FT'
    omics_label2 <- 'GT'
    omics_label_plot <- 'function transcript'
  }
  # omics
  ga_norm <- readRDS(paste0(phyloseq_dir, '/', omics_label2, '.rds'))
  # Taxa (omics)
  ga_norm_taxa <- as.data.frame(unclass(tax_table(ga_norm)), stringsAsFactors = F)
  ga_norm_atxa_genes <- rownames(ga_norm_taxa)
  
  # Metadata
  # Metadata
  metadata_df <- data.frame(sample_data(ga_norm), stringsAsFactors = F)
  samp <- sample_data(metadata_df)
  rownames(samp) <- samp$sample
  
  ge_tpm_taxa_genes <- rownames(omics_tpm_profile_rowsum_not_0)
  ge_fe_df <- data.frame(DB='genes', feature_n=length(ge_tpm_taxa_genes), feature = 'Genes')
  
  #################################### FUNCTION EXPRESSION  ####################################
  print('#################################### FUNCTION EXPRESSION  ####################################')
  manual_plot_colors =c('#9D0208', '#264653','#e9c46a','#D8DCDE','#B6D0E0',
                        '#FFC87E','#F4A261','#E34F33','#E9C46A',
                        '#A786C9','#D4C0E2','#975773','#6699FF','#000066',
                        '#7AAFCA','#006699','#A9D181','#2F8475','#264445') 
  omics_tpm_profile_rowsum_not_0$coord <- rownames(omics_tpm_profile_rowsum_not_0)
  rownames(omics_tpm_profile_rowsum_not_0) <- NULL

  profiles_annot_tpm_df_sub$coord <- rownames(profiles_annot_tpm_df_sub)
  profiles_annot_tpm_df_sub <- as.data.frame(profiles_annot_tpm_df_sub %>%
                                               dplyr::select(coord, ID_gene, name, everything()))
  profiles_annot_tpm_df_sub <- profiles_annot_tpm_df_sub[!colnames(profiles_annot_tpm_df_sub) %in% c("KEGG_ko", "kofam_KO")]
  
  profiles_annot_tpm_df_sub <- profiles_annot_tpm_df_sub[profiles_annot_tpm_df_sub$coord %in% omics_tpm_profile_rowsum_not_0$coord,]
  mappings_ncol <- ncol(profiles_annot_tpm_df_sub)
  for(db in colnames(profiles_annot_tpm_df_sub)[4:mappings_ncol]){#change it deppending on the number of selectedMappings columns
    print(db)
    message(Sys.time(), ": Running ", db)
    dbMap <- profiles_annot_tpm_df_sub[,c("coord", "ID_gene","name",db)]
    dbMap <- dbMap[!is.na(dbMap[[db]]),]
    dbMap <- dbMap[dbMap[[db]] != "",]
    dbMap <- dbMap[dbMap[[db]] != "-",]
    dbMap <- dbMap[rowSums(is.na(dbMap)) != ncol(dbMap),]
    dbkeys <- unique(unlist(strsplit(profiles_annot_tpm_df_sub[,db], split = "\\; ")))
    dbkeys <- unique(unlist(strsplit(dbkeys, split = "\\,")))
    singleKeys <- strsplit(dbMap[,4], split = "\\; |\\;|\\, |\\,")
    db_name <- gsub(' ', '.', db)
    if (length(singleKeys) > 0){
      keyMap <- data.frame(
        coord = rep(dbMap$coord, sapply(singleKeys, length)),
        ID_gene = rep(dbMap$ID_gene, sapply(singleKeys, length)),
        name = rep(dbMap$name, sapply(singleKeys, length)),
        Funct = unlist(singleKeys),
        stringsAsFactors = FALSE
      )
      if (db == "KEGG_Pathway"){
        dim(keyMap)
        keyMap$Funct<- gsub('ko', '',keyMap$Funct)
        keyMap$Funct<- gsub('map', '',keyMap$Funct)
        keyMap <- keyMap[!duplicated(keyMap),]
        keyMap$Funct <- paste0('map', keyMap$Funct)
        dim(keyMap)
      }
      keyMap_u <- keyMap[!duplicated(keyMap[c(1,2,3,4)]),]
      message(Sys.time(), ": ", paste(db, "processed!"))
      keyMap_funct <- data.frame(Funct = c(unique(keyMap_u$Funct)))
      
      ### Annotation descriptions ####
      # library(purrr)
      # library(magrittr)
      if (db_name == 'eggNOG.OGs'){
        cog_df <- as.data.frame(fread(paste0(minto_dir,'/data/cog-20.def.tab'), sep ="\t", header = F, stringsAsFactors = F))
        names(cog_df) <- c('ID_COG', 'funct_category_COG', 'name_COG', 'name_gene', 'funct_pathway', 'ID_PubMed', 'ID_PDB')
        #cog_df <- subset(cog_df, select=c('ID_COG', 'name_COG', 'funct_pathway'))
        cog_df <- subset(cog_df, select=c('ID_COG', 'name_COG'))
        names(cog_df) <- c('Funct', 'Description')
        annot_df <- merge(keyMap_funct, cog_df, by = 'Funct', all.x = T)
      } else if  (db_name == 'merged_KO'){
        annot_df <- data.frame(Description=unlist(lapply(keyMap_funct$Funct, function (x) keggFind('ko', x))),stringsAsFactors = F)
        annot_df$Funct <- rownames(annot_df)
        annot_df$Funct <- gsub('ko:', '', annot_df$Funct)
      }else if  (db_name == 'KEGG_Pathway'){
        annot_df <- data.frame(Description=unlist(lapply(keyMap_funct$Funct, function (x) keggFind('pathway', x))),stringsAsFactors = F)
        annot_df$Funct <- rownames(annot_df)
        annot_df$Funct <- gsub('path:', '', annot_df$Funct)
      }else if  (db_name == 'KEGG_Module'){
        modules_list <- as.data.frame(fread(paste0(minto_dir,'/data/KEGG_Modules_20171212.csv'), header=T), stringsAsFactors = F)
        modules_list_sub <- subset(modules_list, select=c('Module', 'Definition'))
        names(modules_list_sub) <- c('Funct', 'Description')
        annot_df <- merge(keyMap_funct, modules_list_sub, by='Funct', all=T)
      }else if  (db_name %in% c('PFAMs')){
        x <- PFAMID
        mapped_keys <- mappedkeys(x)
        xx <- as.list(x[mapped_keys])
        test <- as.data.frame(do.call(rbind, xx))
        names(test) <- 'Pfam_family'
        test$Pfam_id <- rownames(test)
        x <- PFAMDE
        mapped_keys <- mappedkeys(x)
        xx <- as.list(x[mapped_keys])
        test2 <- as.data.frame(do.call(rbind, xx))
        names(test2) <- 'Description'
        test2$Pfam_id <- rownames(test2)
        test3 <- merge(test, test2, by='Pfam_id')
        annot_df <- merge(keyMap_funct, test3, by.x='Funct' , by.y='Pfam_family', all.x=T)
        annot_df$Pfam_id <- NULL
        rm(test, test2, test3, x, xx)
      } else if  (db_name %in% c('dbCAN.mod', 'dbCAN.enzclass', 'CAZy')){
        x <- PFAMCAZY
        mapped_keys <- mappedkeys(x)
        xx <- as.list(x[mapped_keys])
        test <- as.data.frame(do.call(rbind, xx))
        #data.frame(matrix((xx), nrow=length(xx), byrow=F))
        #names(test) <- 'CAZY_id'
        test$Pfam_id <- rownames(test)
        test_cazy <- as.data.frame(test %>% group_by(Pfam_id) %>% 
                                     mutate(CAZY_id = paste(unique(c(V1, V2, V3, V4)), collapse = ',')))
        test_cazy$V1 <- test_cazy$V2 <- test_cazy$V3 <- test_cazy$V4 <- NULL
        test_cazy.singleKeys <- strsplit(test_cazy$CAZY_id, split = "\\  |\\,")
        if (length(test_cazy.singleKeys) > 0){
          test_cazy.keyMap <- data.frame(
            Pfam_id = rep(test_cazy$Pfam_id, sapply(test_cazy.singleKeys, length)),
            Funct = unlist(test_cazy.singleKeys),
            stringsAsFactors = FALSE
          )
        }
        x <- PFAMDE
        mapped_keys <- mappedkeys(x)
        xx <- as.list(x[mapped_keys])
        test2 <- as.data.frame(do.call(rbind, xx))
        names(test2) <- 'Description'
        test2$Pfam_id <- rownames(test2)
        test3 <- merge(test_cazy.keyMap, test2, by='Pfam_id')
        annot_df <- merge(keyMap_funct, test3, by='Funct', all.x=T)
        annot_df$Pfam_id <- NULL
        rm(test, test2, test3, test_cazy, test_cazy.keyMap, test_cazy.singleKeys, x, xx)
      } else {
        annot_df <- keyMap_funct
        annot_df$Description <- '-'
      }
      
      #### omics ####
      omics_tpm_profile_rowsum_not_0_annot <- merge(keyMap_u, omics_tpm_profile_rowsum_not_0, by='coord', all=T)
      omics_tpm_profile_rowsum_not_0_annot$coord <- NULL
      omics_tpm_profile_rowsum_not_0_annot$ID_gene <- NULL
      omics_tpm_profile_rowsum_not_0_annot$name <- NULL
      
      omics_tpm_profile_rowsum_not_0_annot_dplyr <- as.data.frame(omics_tpm_profile_rowsum_not_0_annot %>%
                                                                    dplyr::group_by(Funct) %>%
                                                                    dplyr::summarise(across(everything(),sum)))
      omics_tpm_profile_rowsum_not_0_annot_dplyr$Funct[which(is.na(omics_tpm_profile_rowsum_not_0_annot_dplyr$Funct))] <- 'Unknown'
      names(omics_tpm_profile_rowsum_not_0_annot_dplyr) <- gsub('X', '', names(omics_tpm_profile_rowsum_not_0_annot_dplyr))
      
      # Calculate TPM from Function abundances (normalized to gene length) ####
      rownames(omics_tpm_profile_rowsum_not_0_annot_dplyr) <- omics_tpm_profile_rowsum_not_0_annot_dplyr$Funct
      omics_tpm_profile_rowsum_not_0_annot_dplyr$Funct <- NULL
      omics_tpm_profile_rowsum_not_0_annot_dplyr[is.na(omics_tpm_profile_rowsum_not_0_annot_dplyr)] <- 0
      
      funct_rpk_samples_u_t_df <- as.data.frame(t(omics_tpm_profile_rowsum_not_0_annot_dplyr))
      ## RPK sum per sample
      rpk_sum <- as.data.frame(colSums(omics_tpm_profile_rowsum_not_0_annot_dplyr))
      names(rpk_sum) <- 'total_sum'
      
      funct_rpk_samples_total_reads_t_df <- cbind(funct_rpk_samples_u_t_df, rpk_sum)
      # ## Calculate TPM from RPK and SCALING FACTOR
      # gene_tpm_samples_u_t_df <- gene_rpk_samples_total_reads_t_df[1:ncol(gene_rpk_samples_total_reads_t_df)-1]/gene_rpk_samples_total_reads_t_df$scaling_factor
      
      ## Calculate TPM from RPK and RPK sum per sample ####
      funct_tpm_samples_u_t_df <- funct_rpk_samples_total_reads_t_df[1:ncol(funct_rpk_samples_total_reads_t_df)-1]/funct_rpk_samples_total_reads_t_df$total_sum
      funct_tpm_samples_u_t_df <- funct_tpm_samples_u_t_df*1e6
      
      omics_funct_profile_rowsum_not_0_annot_dplyr <- as.data.frame(t(funct_tpm_samples_u_t_df))
      
      colSums(omics_funct_profile_rowsum_not_0_annot_dplyr)
      
      omics_funct_profile_rowsum_not_0_annot_dplyr$Funct <- rownames(omics_funct_profile_rowsum_not_0_annot_dplyr)
      omics_funct_profile_rowsum_not_0_annot_dplyr <- as.data.frame(omics_funct_profile_rowsum_not_0_annot_dplyr %>%
                                                                      dplyr::select(Funct, everything()))
      omics_funct_profile_rowsum_not_0_annot_desc <- merge(annot_df, omics_funct_profile_rowsum_not_0_annot_dplyr, by='Funct', all.y=T)
      write.csv(omics_funct_profile_rowsum_not_0_annot_desc,paste0(integration_dir, '/', omics_label, '.', db_name ,'.csv'), row.names = F, quote = F)
      
      #### phyloseq object ####
      rownames(omics_funct_profile_rowsum_not_0_annot_dplyr) <- omics_funct_profile_rowsum_not_0_annot_dplyr$Funct
      omics_funct_profile_rowsum_not_0_annot_dplyr$Funct <- NULL
      
      ## metadata
      #metadata_df <- data.frame(sample=names(omics_funct_profile_rowsum_not_0_annot_dplyr), sample_color='sample')
      samp <- sample_data(metadata_df)
      rownames(samp) <- samp$sample
      
      # taxa
      omics_annot_desc <- subset(omics_funct_profile_rowsum_not_0_annot_desc, select=c("Funct","Description"))
      rownames(omics_annot_desc) <- omics_annot_desc$Funct
      
      omics_function_expression_phyloseq <- phyloseq(otu_table(as.matrix(omics_funct_profile_rowsum_not_0_annot_dplyr), taxa_are_rows = T), 
                                                     tax_table(as.matrix(omics_annot_desc)), samp)
      saveRDS(omics_function_expression_phyloseq, file = paste0(phyloseq_dir, '/', omics_label, '.', db_name ,'.rds'))
      
      # Create df to plot number of genes and functions ####
      ge_fe_df <- rbind(ge_fe_df, c(db_name, nrow(omics_annot_desc[omics_annot_desc$Funct!='Unknown',]), 'Functions'))
      
      print('#################################### PCoA  ####################################')#####
      distance_lab='bray'
      # PCoA - metaG ####
      ## Prepare data - replace NA values by 0
      omics_function_expression_noinf <- as.data.frame(omics_funct_profile_rowsum_not_0_annot_dplyr %>%
                                                         tibble::rownames_to_column('rownames') %>%
                                                         dplyr::mutate_all(function(x) ifelse(is.infinite(x), 0, x)) %>%
                                                         tibble::column_to_rownames('rownames'))
      omics_function_expression_noinf <- as.data.frame(omics_function_expression_noinf %>%
                                                         tibble::rownames_to_column('rownames') %>%
                                                         dplyr::mutate_all(function(x) ifelse(is.na(x), 0, x)) %>%
                                                         tibble::column_to_rownames('rownames'))
      
      omics_function_expression_noinf_plot <- phyloseq(otu_table(as.matrix(omics_function_expression_noinf), taxa_are_rows = T), 
                                                       tax_table(as.matrix(omics_annot_desc)), samp)
      
      title_name <- paste0('PCoA - ', omics_label_plot, ' - ', db_name)
      out_name <- paste0(visual_dir, omics_label, '.', db_name ,'.PCoA-',distance_lab,'.pdf')
      plot_PCoA_out <- plot_PCoA(distance_lab, omics_function_expression_noinf_plot, "condition", "sample_alias")
      pdf(out_name,width=8,height=8,paper="special" )
      print(plot_PCoA_out  + scale_color_manual(values=manual_plot_colors, name='Condition') +
              #coord_fixed() +
              theme(legend.position="bottom")+ggtitle(title_name, 'Bray Curtis'))
      print(plot_PCoA_out  + 
              facet_wrap(.~condition)+
              scale_color_manual(values=manual_plot_colors, name='Condition') +
              #coord_fixed() +
              theme(legend.position="bottom")+
              ggtitle(title_name, 'Bray Curtis'))
      dev.off()
    } else{
      # create .csv, .rds, .pdf
      distance_lab='bray'
      file.create(paste0(integration_dir, '/FA.', db_name ,'.csv'))
      file.create(paste0(integration_dir, '/FT.', db_name ,'.csv'))
      file.create(paste0(integration_dir, '/FE.', db_name ,'.csv'))
      file.create(paste0(phyloseq_dir, 'FA.', db_name ,'.rds'))
      file.create(paste0(phyloseq_dir, 'FT.', db_name ,'.rds'))
      file.create(paste0(phyloseq_dir, 'FE.', db_name ,'.rds'))
      file.create(paste0(visual_dir, 'FA.', db_name ,'.PCoA-',distance_lab,'.pdf'))
      file.create(paste0(visual_dir, 'FT.', db_name ,'.PCoA-',distance_lab,'.pdf'))
      file.create(paste0(visual_dir, 'FE.', db_name ,'.PCoA-',distance_lab,'.pdf'))
      ge_fe_df <- rbind(ge_fe_df, c(db_name, 0, 'Functions'))
    }
  }
  
  features_file <- paste0(omics_label2,'_', omics_label,'_features')
  write.csv(ge_fe_df,paste0(integration_dir, '/', features_file, '.csv'), row.names = F, quote = F)
  
  ge_fe_df<- as.data.frame(fread(paste0(integration_dir, '/', features_file, '.csv'), header=T), stringsAsFactors = F, row.names = T)
  ge_fe_df$feature_n <- as.numeric(ge_fe_df$feature_n)
  feature_order <- ge_fe_df$DB[order(ge_fe_df$feature_n)]
  ge_fe_df$DB <- factor(ge_fe_df$DB, levels=feature_order)
  
  pdf(paste0(visual_dir,'/', features_file, '.pdf'),width=6,height=5,paper="special" )
  print(ggplot(data=ge_fe_df, aes(x=DB, y=feature_n)) +
          geom_bar(stat="identity")+ theme_minimal()+
          facet_wrap(feature~., scales= "free") + labs(y='Number of features', x='')+
          theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))+
          geom_text(aes(label=feature_n), position=position_dodge(width=0.9), vjust=-0.25, size = 3)
        
  )
  dev.off()
}

print('done!')


