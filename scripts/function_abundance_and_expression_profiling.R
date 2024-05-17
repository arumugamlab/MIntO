#!/usr/bin/env Rscript

# '''
# Generate function expression profile from genome-based mode gene profiles (MG normalized)
#
# Authors: Carmen Saenz, Mani Arumugam
#
# '''


##########################  ** Load libraries **  ##########################
library(data.table)
library(phyloseq)
library(tidyr)
library(stringr)
library(purrr)
library(optparse)
library(qs)
library(dplyr)

logmsg <- function(...) {
    message("# ", format(Sys.time(), digits=0), " - ", paste(list(...)), collapse=" ")
}

# Parse command line arguments
opt_list <- list(
                make_option("--threads", type="integer", default=4, help="number of threads [default: %default]"),
                make_option("--outdir", type="character", default=NULL, help="output directory to write normalized counts", metavar="directory"),
                make_option("--omics", type="character", default=NULL, help="which omics to summarize: metaG, metaT or metaG_metaT", metavar="string"),
                make_option("--normalization", type="character", default=NULL, help="normalization type: MG or TPM"),
                make_option("--funcat-name", type="character", default=NULL, help="name of functional category to summarize: e.g., eggNOG_OGs, dbCAN.EC", metavar="string"),
                make_option("--funcat-desc", type="character", default=NULL, help="file containing descriptions of functions", metavar="file"),
                make_option("--genome-weights-metaG", type="character", default=NULL, help="file with metaG profiles", metavar="file"),
                make_option("--genome-weights-metaT", type="character", default=NULL, help="file with metaT profiles", metavar="file"),
                make_option("--main-factor", type="character", default=NULL, help="main factor variable to use when plotting data", metavar="string")
                )
opts <- parse_args(OptionParser(option_list=opt_list))


##########################  ** Load arguments **  ##########################
threads_n <- opts$threads
output_dir <- opts$outdir
omics <- opts$omics #'metaG_metaT'
normalization <- opts$normalization
funcat_name  <- opts[['funcat-name']]
funcat_desc_file <- opts[['funcat-desc']]
main_factor <- opts[['main-factor']]
metaG_profile_file <- opts[['genome-weights-metaG']]
metaT_profile_file <- opts[['genome-weights-metaT']]

setDTthreads(threads = as.numeric(threads_n))
set.seed(1234)

# genome profiles to use as weights

metaG_genome_weights <- NULL
metaT_genome_weights <- NULL

if (normalization == 'MG') {
    if (omics %in% c('metaG', 'metaG_metaT')) {
        if (is.null(metaG_profile_file)) {
            stop("--genome-weights-metaG must be provided with --omics=metaG")
        } else {
            metaG_genome_weights <- (
                                     fread(metaG_profile_file, header=T)
                                     [, Funct := 'REL']
                                    )
            setnames(metaG_genome_weights, "ID", "MAG")
            setnames(metaG_genome_weights, colnames(metaG_genome_weights),
                     gsub(x = colnames(metaG_genome_weights), pattern = "metaG\\.", replacement = ""))
        }
    }

    if (omics %in% c('metaT', 'metaG_metaT')) {
        if (is.null(metaT_profile_file)) {
            stop("--genome-weights-metaT must be provided with --omics=metaT")
        } else {
            metaT_genome_weights <- (
                                     fread(metaT_profile_file, header=T)
                                     [, Funct := 'REL']
                                    )
            setnames(metaT_genome_weights, "ID", "MAG")
            setnames(metaT_genome_weights, colnames(metaT_genome_weights),
                     gsub(x = colnames(metaT_genome_weights), pattern = "metaT\\.", replacement = ""))
        }
    }
}



##########################  ** Load functions **  ##########################

# Requirements:
#  'profile' should only have numeric columns - IDs should have been removed
#  'profile' should have removed all zerosum rows
#  'profile' should be sorted (desc) by rowSum
plot_PCA <- function(profile, label, color, metadata) {

    library(ggplot2)
    library(ggrepel)

    ##### metadata:
    sample_data_df <- data.frame(metadata, stringsAsFactors = F)

    ### Counts - non-zero/inf
    ## Prepare data - replace NA values by 0
    ## Update 2024.10.05 - NA/Inf should not be there in the tables anymore.

    # Get the lowest value in table and its log() to replace 0s by log2(min value * 1e-2) for the PCA
    profile_copy <- (
                     profile
                     [, lapply(.SD, function(x) min(x[x>0], na.rm = TRUE))]
                     [, minv := min(.SD, na.rm = TRUE)]
                    )
    min_value <- (min(profile_copy[, .(minv)]))
    logmsg("  Min value in table: ", min_value)
    replace_val <- log(min_value*1e-2)

    # Log transform values
    pca_data <- profile[, lapply(.SD, function(x) ifelse(x==0, replace_val, log(x)))]
    gc()
    pca_data <- data.table::transpose(pca_data)
    gc()


    #Plotting scores of PC1 and PC2 with log transformation
    pca_results <- prcomp(pca_data, center = T, sca=T)

    # ***************** ggplot way - w coloring *****************
    dtp <- data.frame('sample' = sample_data_df[[label]],
                      'group' = sample_data_df[[color]],
                       pca_results$x[,1:2]) # the first two componets are selected
    total_variance <- sum(pca_results$sdev)
    axis_names <- head(colnames(data.table(pca_results$x)), 2)
    percentage <- round(head(pca_results$sdev, 2) / total_variance * 100, 2)
    percentage <- paste0(axis_names, " (", percentage, "%)")

    PCA_Sample_site_abundance <- (ggplot(data = dtp, aes(x = PC1, y = PC2, color = group)) +
                                  geom_text_repel(aes(label = sample),nudge_x = 0.04, size = 3.5, segment.alpha = 0.5) +
                                  geom_point(size = 2, shape = 16)+
                                  xlab(percentage[1]) + ylab(percentage[2]) +
                                  #labs(title = title) +
                                  theme_bw()+
                                  theme(plot.title = element_text(size=10), legend.position="bottom"))#+ stat_ellipse(type = "norm", linetype = 2))
    return(PCA_Sample_site_abundance)
}

prepare_PCA <- function(profile, label, color, metadata, title) {

    logmsg(" Making PCA plots")

    # Colors
    manual_plot_colors =c('#9D0208', '#264653','#e9c46a','#D8DCDE','#B6D0E0',
                          '#FFC87E','#F4A261','#E34F33','#E9C46A',
                          '#A786C9','#D4C0E2','#975773','#6699FF','#000066',
                          '#7AAFCA','#006699','#A9D181','#2F8475','#264445')

    # Title
    title_name <- paste0('PCA - ', title)

    # File names
    # abundance  -> GA
    # transcript -> GT
    # expression -> GE
    out_name <- paste0(visual_dir, '/', label, '.PCA.pdf')

    plot_PCA_out <- plot_PCA(profile=profile, color=color, label="sample_alias", metadata=metadata)
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

get_annotation_descriptions <- function(db_name, functionListDT) {
    ### Annotation descriptions ####
    # library(purrr)
    # library(magrittr)
    if  (db_name %in% c('eggNOG_OGs', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_KO', 'kofam_Pathway', 'kofam_Module', 'kofam_KO', 'merged_KO', 'dbCAN.EC')){
        annotations <- (fread(funcat_desc_file, header=T)
                        [, .(Funct, Description)]
                        [, Description := iconv(Description, from = "ISO-8859-1", to = "UTF-8")]
                       )
        setkey(annotations, Funct)
        annotations <- unique(merge(functionListDT, annotations, by='Funct', all=T)
                              [, Description := paste(sort(unique(Description)), collapse=';'), by=Funct]
                             )
    }
    else if  (db_name %in% c('PFAMs')){

        library(PFAM.db)

        # Get PFAM IDs
        pfam_names <- as.list(PFAMID[mappedkeys(PFAMID)])
        setDT(pfam_names, keep.rownames = TRUE)
        pfam_names <- data.table::transpose(pfam_names, keep.names = "Pfam_id")
        setnames(pfam_names, "V1", "Funct")

        # Get PFAM descriptions
        pfam_desc <- as.list(PFAMDE[mappedkeys(PFAMDE)])
        setDT(pfam_desc, keep.rownames = TRUE)
        pfam_desc <- data.table::transpose(pfam_desc, keep.names = "Pfam_id")
        setnames(pfam_desc, "V1", "Description")

        # Merge based on PFAM_ID and delete that join key
        pfam_table <- merge(pfam_names, pfam_desc, by='Pfam_id')[, Pfam_id := NULL]

        # Merge descriptions for the relevant PFAMs
        annotations <- unique(merge(functionListDT, pfam_table, by='Funct', all.x=T)
                                    [,Description := paste(sort(unique(Description)), collapse=';'), by='Funct']
                             )
        rm(pfam_names, pfam_desc, pfam_table, x, xx, keyMap_funct_desc,keyMap_funct_desc2)
    }
    else if  (db_name %in% c('dbCAN.module', 'dbCAN.enzclass', 'CAZy')){

        library(PFAM.db)

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
          test_cazy.keyMap <- data.table(
            Pfam_id = rep(test_cazy$Pfam_id, sapply(test_cazy.singleKeys, length)),
            Funct = unlist(test_cazy.singleKeys)
          )
        }

        # Get PFAM descriptions
        pfam_desc <- as.list(PFAMDE[mappedkeys(PFAMDE)])
        setDT(pfam_desc, keep.rownames = TRUE)
        pfam_desc <- data.table::transpose(pfam_desc, keep.names = "Pfam_id")
        setnames(pfam_desc, "V1", "Description")

        # Merge based on PFAM_ID and delete that join key
        cazy_table <- merge(test_cazy.keyMap, pfam_desc, by='Pfam_id')[, Pfam_id := NULL]

        # Merge descriptions for the relevant PFAMs
        annotations <- unique(merge(functionListDT, cazy_table, by='Funct', all.x=T)
                              [,Description := paste(sort(unique(Description)), collapse=';'), by='Funct']
                             )
    }
    else {
        annotations <- unique(functionListDT[, Description := '-'])
    }
    return(annotations)
}

make_profile_files <- function(keys, profile, file_label, database, annotations, weights) {

      # Make a list of sample-names
      sample_cols <- setdiff(colnames(profile), c('gene_id'))

      # Get abundance information by gene and select only relevant columns.
      logmsg("total genes ", nrow(profile))
      counts <- (merge.data.table(keys, profile, on="Funct", all.x=FALSE, all.y=TRUE)
                 [, MAG := stringr::str_split(gene_id, "\\|") %>% map_chr(., 1)] # Get the first field
                 [, gene_id := NULL]
                )
      logmsg("clean genes ", nrow(counts))

      # From count-per-gene, make count-per-function
      # If necessary, weight the functions by REL of genomes
      if (!is.null(weights)) {

          # Summarize by MAG,Funct
          # And prepend with weights
          counts <- rbind(weights,
                          counts[, lapply(.SD, sum, na.rm=TRUE), by = c('MAG', 'Funct')]
                         )
          logmsg("MAGxfunctions ", nrow(counts))

          # weigh by the given weights per MAG
          weight_MAG_by_proportions <- function(x){
            if (is.character(x[1])){
              return(x)
            }else{
              x = x*x[1]
              return(x)
            }
          }

          # Normalize by weighting by 'weights'
          # Remove REL rows
          # Summarize by Funct
          counts <- (
                      counts
                      [, lapply(.SD, weight_MAG_by_proportions), by = MAG]
                      [! Funct == 'REL']
                      [, MAG := NULL]
                      [, lapply(.SD, sum, na.rm=TRUE), .SDcols = sample_cols, by = Funct]
                    )
      } else {

          # Summarize by Funct
          counts <- (
                      counts
                      [, MAG := NULL]
                      [, lapply(.SD, sum, na.rm=TRUE), .SDcols = sample_cols, by = Funct]
                    )
      }

      # Estimate 'Unknown' functions - genes with no annotation
      counts <- counts[is.na(Funct), Funct := "Unknown"]
#print(counts[Funct == '28HG4',])

      logmsg("functions ", nrow(counts))

      # Add annotations
      counts <- merge(annotations, counts, by='Funct', all.x=FALSE, all.y=TRUE)

      # Write annotated functional profile as tsv
      fwrite(counts,
             file=paste0(output_dir, '/', file_label, '.', database ,'.tsv'),
             sep='\t',
             row.names = F,
             quote = F)

      #### Create phyloseq object ####

      my_otu_table <- as.matrix(counts[, -c('Funct', 'Description')])
      rownames(my_otu_table) <- counts$Funct

      # taxa
      my_tax_table <- as.matrix(counts)
      rownames(my_tax_table) <- counts$Funct

      physeq <- phyloseq(otu_table(my_otu_table, taxa_are_rows = T),
                         tax_table(my_tax_table),
                         sample_metadata)

      # Write phyloseq object
      qsave(physeq,
            preset = "balanced",
            file = paste0(phyloseq_dir, '/', file_label, '.', database ,'.qs'))

      # Index on Funct
      setkey(counts, Funct)

      # Remove Unknown functions
      # It has been written out in tsv files for metaG and metaT
      # And it should NOT be calculated for FE
      counts <- counts[!J("Unknown")]

      return(counts)
}

process_phyloseq_obj <- function(physeq_file) {

    logmsg("began reading phyloseq")

    physeq <- qread(physeq_file, nthreads=threads_n)

    # Get otu_table and set rownames as 'gene_id'
    # Convert NA's to 0's
    # Reorder features by decreasing sum across samples, by setting negative rowSum as key
    # Remove features with zero sum across all samples
    columns_wanted <- sample_names(physeq)
    counts <- (
                data.table(unclass(otu_table(physeq)), keep.rownames=TRUE)
                [, lapply(.SD, function(x) ifelse(is.na(x), 0, x))]
                [, negRowSum := -rowSums(.SD, na.rm = FALSE), .SDcols = columns_wanted]
              )
    setkey(counts, negRowSum)
    counts <- counts[!J(0),][, negRowSum := NULL]
    setnames(counts, "rn", "gene_id")
    setkey(counts, gene_id)

    # Get tax_table and set rownames as 'gene_id'
    taxa <- data.table(unclass(tax_table(physeq)), keep.rownames=TRUE)
    setnames(taxa, "rn", "gene_id")

    # Get metadata table
    # Note: Don't skip the as.data.frame step - it causes trouble otherwise
    metadata <- data.table(as.data.frame(unclass(sample_data(physeq))), keep.rownames=FALSE)

    logmsg("done  reading phyloseq")

    return(list(counts=counts, taxa=taxa, metadata=metadata))
}

# *************************** Use GA and GT raw counts to generate FE profile *************************** ####
# The genes will be clustered by function IDs.
# TPM normalization - sum all TPMs for the member genes
#                     (there will be redundancy and total might be over 1 million)
# MG normalization  - weight individual member genes by the rel. abundance of the genome and sum.
# FE profile will be generated: (FT_norm) / (FA_norm)
#Raw gene abundances #####

print('#################################### Paths ####################################')

# Generate output directories ####
dir.create(file.path(output_dir), showWarnings = FALSE)
visual_dir=paste0(output_dir,'/plots/')
dir.create(file.path(visual_dir), showWarnings = FALSE)
phyloseq_dir=paste0(output_dir,'/phyloseq_obj/')
dir.create(file.path(phyloseq_dir), showWarnings = FALSE)


print('#################################### DB profiles ####################################')

# Initialize some variables
metadata_df <- NULL
gene_annotation <- NULL
ga_fa_df <- NULL
gt_ft_df <- NULL
ge_fe_df <- NULL

gene_annotation <- NULL

# metaG
if (omics %like% "metaG") {
    res <- process_phyloseq_obj(paste0(phyloseq_dir, '/GA.qs'))

    ga_fa_df <- data.frame(DB='genes', feature_n=nrow(res$taxa), feature = 'Genes')

    # gene profile
    metaG_profile <- res$counts

    # gene annotations
    # Get entries for this funcat_name
    gene_annotation <- res$taxa[, c('gene_id', funcat_name), with = FALSE]
    setkeyv(gene_annotation, funcat_name)
    gene_annotation <- gene_annotation[!J(c(NA, "", "-"))]

    # metadata
    metadata_df <- res$metadata
}

# metaT
if (omics %like% "metaT") {
    res <- process_phyloseq_obj(paste0(phyloseq_dir, '/GT.qs'))

    gt_ft_df <- data.frame(DB='genes', feature_n=nrow(res$taxa), feature = 'Genes')

    metaT_profile <- res$counts

    # gene annotations
    metaT_gene_annotation <- res$taxa[, c('gene_id', funcat_name), with = FALSE]
    setkeyv(metaT_gene_annotation, funcat_name)
    metaT_gene_annotation <- metaT_gene_annotation[!J(c(NA, "", "-"))]

    # Only need entries that are not already present in gene_annotation
    if (omics == 'metaT') {
        gene_annotation <- metaT_gene_annotation
        metadata_df <- res$metadata
    } else { # metaG_metaT

        # Register genes common in metaG/metaT
        # Get union of metaG and metaT annotations

        mG_genes <- nrow(gene_annotation)
        mT_genes <- nrow(metaT_gene_annotation)

        gene_annotation <- funion(gene_annotation, metaT_gene_annotation, all=FALSE)
        total_genes <- nrow(gene_annotation)

        n_common <- mG_genes + mT_genes - total_genes

        ge_fe_df <- data.frame(DB='genes', feature_n=n_common, feature = 'Genes')
    }

    rm(metaT_gene_annotation)
}

# Delete data
rm(res)

# By now, we will have:
#   metadata_df
#   gene_annotation
#   ge_fe_df

## metadata for phyloseq
sample_metadata <- sample_data(metadata_df)
rownames(sample_metadata) <- sample_metadata[["sample_alias"]]


#################################### FUNCTION EXPRESSION  ####################################
print('#################################### FUNCTION EXPRESSION  ####################################')

# Filter tax.table based on features in metaG, if metaG_metaT. metaG/metaT have been filtered earlier.

#if (omics == 'metaG_metaT') {
#    gene_annotation <- gene_annotation[gene_id %in% metaG_profile$gene_id, ]
#}


if (nrow(gene_annotation) > 0) {

    logmsg(funcat_name, " started")

    ################################################
    # Prepare annotations and their descriptions
    ################################################

    # Replicate the gene annotation table by having N rows for gene with N annotations.

    singleKeys <- strsplit(gene_annotation[[funcat_name]], split = "\\; |\\;|\\, |\\,")
    Gene2FuncMap <- unique(data.table(gene_id = rep(gene_annotation$gene_id, sapply(singleKeys, length)),
                                Funct = unlist(singleKeys)))
    setkey(Gene2FuncMap, gene_id)
    rm(gene_annotation)

    # Since KEGG annotations can have map12345 and ko12345 for the same pathway,
    # remove such redundancy.
    if (funcat_name %in% c("KEGG_Pathway", "kofam_Pathway")) {
        Gene2FuncMap <- unique(Gene2FuncMap[, Funct := gsub('ko', 'map', Funct)])
    }

    myFunctions <- data.table(Funct = c(unique(Gene2FuncMap$Funct)))
    setkey(myFunctions, Funct)

    # Get annotation descriptions from external sources

    annot_dt <- get_annotation_descriptions(funcat_name, myFunctions)
    setkey(annot_dt, Funct)
    rm(myFunctions)

    # Free up memory
    logmsg("Freeing memory")
    gc()

    #####################################################
    # Make profiles and write them to 'output' directory
    #####################################################

    # metaG or metaG_metaT
    if (omics %like% "metaG") {
        metaG_counts <- make_profile_files(keys=Gene2FuncMap,
                                           profile=metaG_profile,
                                           file_label="FA",
                                           database=funcat_name,
                                           annotations=annot_dt,
                                           weights=metaG_genome_weights)
    }

    # metaT or metaG_metaT
    if (omics %like% "metaT") {
        metaT_counts <- make_profile_files(keys=Gene2FuncMap,
                                           profile=metaT_profile,
                                           file_label="FT",
                                           database=funcat_name,
                                           annotations=annot_dt,
                                           weights=metaT_genome_weights)
    }

    rm(singleKeys, Gene2FuncMap)

    # metaG_metaT
    if (omics == 'metaG_metaT') {

        sample_cols <- setdiff(colnames(metaG_counts), c("Funct", "Description"))

        # Get common features
        feature_ids <- intersect(metaG_counts[['Funct']], metaT_counts[['Funct']])
        metaG_counts <- metaG_counts[J(feature_ids),]
        metaT_counts <- metaT_counts[J(feature_ids),]

        # Backup annotations for future use
        FE_annotations <- data.table(Funct       = metaG_counts[['Funct']],
                                     Description = metaG_counts[['Description']])

        # Remove annotations
        metaG_counts <- metaG_counts[, `:=` (Funct = NULL, Description = NULL) ]
        metaT_counts <- metaT_counts[, `:=` (Funct = NULL, Description = NULL) ]

        # Function Expression
        function_expression <- (metaT_counts)/(metaG_counts)

        # Convert Inf and NAN to NA
        function_expression <- (
                                function_expression
                                [, lapply(.SD, function(x) ifelse(is.infinite(x) | is.nan(x), NA, x))]
                               )

        # cbind annotations back for writing tsv file
        function_expression <- cbind(FE_annotations, function_expression)

        # Estimate rowSum but negate it so that setkey will sort by desc(rowSum)
        logmsg("Estimating rowSum")
        function_expression <- function_expression[, negRowSum := -rowSums(.SD, na.rm = TRUE), .SDcols = sample_cols]
        logmsg("  done")

        # Set key
        logmsg("Indexing")
        setkey(function_expression, negRowSum)
        logmsg("  done")

        # Remove zerosum rows
        logmsg("Removing zero-sum rows")
        logmsg("  Before: ", dim(function_expression))
        function_expression <- function_expression[!J(0)][, negRowSum := NULL]
        logmsg("  After : ", dim(function_expression))
        logmsg("  done")

        # Write txt file
        fwrite(function_expression, file=paste0(output_dir, '/FE.', funcat_name ,'.tsv'), sep='\t', row.names = F, quote = F)

        #### phyloseq object - metaT relative to metaG -####

        function_expression <- function_expression[, Description := NULL]
        FE_physeq <- phyloseq(otu_table(as.matrix(function_expression, rownames="Funct"), taxa_are_rows = T),
                              tax_table(as.matrix(FE_annotations, rownames="Funct")),
                              sample_metadata)
        qsave(FE_physeq,
              preset = "balanced",
              file = paste0(phyloseq_dir, '/FE.', funcat_name ,'.qs'))

        # Free up memory
        logmsg("Freeing memory")
        rm(feature_ids, FE_annotations, FE_phyloseq)
        gc()
    }


    ###########
    # Plot PCA
    ###########

    # PCA - metaG ####
    if (omics %like% 'metaG') {
        # Record feature count
        n_row <- nrow(metaG_counts)
        ga_fa_df <- rbind(ga_fa_df, c(funcat_name, n_row, 'Functions'))

        # Remove Funct column
        metaG_counts <- metaG_counts[, `:=` (Funct = NULL, Description = NULL) ]

        # Make PCA
        prepare_PCA(profile=metaG_counts, title=paste0("function abundance - ", funcat_name), label=paste0('FA.', funcat_name), color=main_factor, metadata=sample_metadata)
    }

    # PCA - metaT ####
    if (omics %like% 'metaT') {
        # Record feature count
        n_row <- nrow(metaT_counts)
        gt_ft_df <- rbind(gt_ft_df, c(funcat_name, n_row, 'Functions'))

        # Remove Funct column
        metaT_counts <- metaT_counts[, `:=` (Funct = NULL, Description = NULL) ]

        # Make PCA
        prepare_PCA(profile=metaT_counts, title=paste0("function transcript - ", funcat_name), label=paste0('FT.', funcat_name), color=main_factor, metadata=sample_metadata)
    }

    # PCA - GE ####
    if (omics == 'metaG_metaT') {
        # Record feature count
        n_row <- nrow(function_expression)
        ge_fe_df <- rbind(ge_fe_df, c(funcat_name, n_row, 'Functions'))

        # Remove Funct column
        function_expression <- function_expression[, Funct := NULL][, lapply(.SD, function(x) ifelse(is.na(x), 0, x))]

        # Make PCA
        prepare_PCA(profile=function_expression, title=paste0("function expression - ", funcat_name), label=paste0('FE', funcat_name), color=main_factor, metadata=sample_metadata)
    }

    logmsg(funcat_name, " finished!")

} else {
    # create .tsv, .qs, .pdf
    file.create(paste0(output_dir, '/FA.', funcat_name ,'.tsv'))
    file.create(paste0(output_dir, '/FT.', funcat_name ,'.tsv'))
    file.create(paste0(output_dir, '/FE.', funcat_name ,'.tsv'))
    file.create(paste0(phyloseq_dir, 'FA.', funcat_name ,'.qs'))
    file.create(paste0(phyloseq_dir, 'FT.', funcat_name ,'.qs'))
    file.create(paste0(phyloseq_dir, 'FE.', funcat_name ,'.qs'))
    file.create(paste0(visual_dir, 'FA.', funcat_name ,'.PCA.pdf'))
    file.create(paste0(visual_dir, 'FT.', funcat_name ,'.PCA.pdf'))
    file.create(paste0(visual_dir, 'FE.', funcat_name ,'.PCA.pdf'))
    ge_fe_df <- rbind(ge_fe_df, c(funcat_name, 0, 'Functions'))
}

# Free memory
rm(metaG_counts, metaT_counts, function_expression)
gc()

#####################################
# Write out feature count statistics
#####################################

print('#################################### FEATURE COUNTS  ####################################')#####

# Write counts
if (!is.null(ga_fa_df)) {
    fwrite(ga_fa_df, file=paste0(output_dir, '/GA_FA_features.', funcat_name, '.tsv'), sep='\t', row.names = F, quote = F)
}

if (!is.null(gt_ft_df)) {
    fwrite(gt_ft_df, file=paste0(output_dir, '/GT_FT_features.', funcat_name, '.tsv'), sep='\t', row.names = F, quote = F)
}

if (!is.null(ge_fe_df)) {
    fwrite(ge_fe_df, file=paste0(output_dir, '/GE_FE_features.', funcat_name, '.tsv'), sep='\t', row.names = F, quote = F)
}

print('done!')
