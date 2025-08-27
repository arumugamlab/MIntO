#!/usr/bin/env Rscript

# '''
# Generate function expression profile from genome-based mode gene profiles (MG or TPM normalized)
# Depending on --omics, it creates a set of files:
#   metaG => FA.*
#   metaT => FT.*
#   metaG_metaT => FE.*
#
# Authors: Carmen Saenz, Mani Arumugam
#
# '''


##########################  ** Load libraries **  ##########################
library(phyloseq)
library(optparse)
library(qs)
library(data.table)
library(this.path)
library(stringr) # not for strings, but for '%>%'

##########################  ** Load functions **  ##########################

# Include common utility functions
source(this.path::here('include', 'utils.R'))

# Include PCA functions
source(this.path::here('include', 'plots_PCA.R'))

# Parse command line arguments
opt_list <- list(
                make_option("--gene-profile-metaG", type="character", default=NULL, help="file containing metaG gene abundances"),
                make_option("--gene-profile-metaT", type="character", default=NULL, help="file containing metaT gene abundances"),
                make_option("--gene-annotation-file", type="character", default=NULL, help="file containing gene annotations"),
                make_option("--metadata-file", type="character", default=NULL, help="file containing sample metadata"),
                make_option("--min-completeness", type="double", default=75, help="minimum module completeness [default: %default]"),
                make_option("--completeness-file", type="character", default=NULL, help="file containing module completeness"),
                make_option("--threads", type="integer", default=4, help="number of threads [default: %default]"),
                make_option("--outdir", type="character", default=NULL, help="output directory to write normalized counts", metavar="directory"),
                make_option("--omics", type="character", default=NULL, help="which omics to summarize: metaG, metaT or metaG_metaT", metavar="string"),
                make_option("--normalization", type="character", default=NULL, help="normalization type: MG or TPM"),
                make_option("--funcat-name", type="character", default=NULL, help="name of functional category to summarize: e.g., eggNOG.OGs, dbCAN.EC", metavar="string"),
                make_option("--funcat-desc", type="character", default=NULL, help="file containing descriptions of functions", metavar="file"),
                make_option("--genome-weights-metaG", type="character", default=NULL, help="file with metaG profiles", metavar="file"),
                make_option("--genome-weights-metaT", type="character", default=NULL, help="file with metaT profiles", metavar="file"),
                make_option("--color-factor", type="character", default=NULL, help="factor used for color when plotting data", metavar="string"),
                make_option("--shape-factor", type="character", default=NULL, help="factor used for shape when plotting data", metavar="string"),
                make_option("--label-factor", type="character", default="sample_alias", help="factor used to label points when plotting data", metavar="string")
                )
opts <- parse_args(OptionParser(option_list=opt_list))


##########################  ** Load arguments **  ##########################
threads_n <- opts$threads
metaG_gene_profile  <- opts[['gene-profile-metaG']]
metaT_gene_profile  <- opts[['gene-profile-metaT']]
gene_annotation_file  <- opts[['gene-annotation-file']]
metadata_file  <- opts[['metadata-file']]
output_dir <- opts$outdir
omics <- opts$omics #'metaG_metaT'
normalization <- opts$normalization
funcat_name  <- opts[['funcat-name']]
funcat_desc_file <- opts[['funcat-desc']]
color_factor <- opts[['color-factor']]
shape_factor <- opts[['shape-factor']]
label_factor <- opts[['label-factor']]
metaG_genome_weights_file <- opts[['genome-weights-metaG']]
metaT_genome_weights_file <- opts[['genome-weights-metaT']]
min_completeness <- opts[['min-completeness']]
completeness_file <- opts[['completeness-file']]

setDTthreads(threads = as.numeric(threads_n))
set.seed(1234)

# genome profiles to use as weights

metaG_genome_weights <- NULL
metaT_genome_weights <- NULL

if (normalization == 'MG') {
    if (omics %in% c('metaG', 'metaG_metaT')) {
        if (is.null(metaG_genome_weights_file)) {
            stop(paste0("--genome-weights-metaG must be provided with --omics=", omics))
        } else {
            metaG_genome_weights <- fread(metaG_genome_weights_file, header=T)
            setnames(metaG_genome_weights, "ID", "MAG")
            setnames(metaG_genome_weights, colnames(metaG_genome_weights),
                     gsub(x = colnames(metaG_genome_weights), pattern = "metaG\\.", replacement = ""))
        }
    }

    if (omics %in% c('metaT', 'metaG_metaT')) {
        if (is.null(metaT_genome_weights_file)) {
            stop(paste0("--genome-weights-metaT must be provided with --omics=", omics))
        } else {
            metaT_genome_weights <- fread(metaT_genome_weights_file, header=T)
            setnames(metaT_genome_weights, "ID", "MAG")
            setnames(metaT_genome_weights, colnames(metaT_genome_weights),
                     gsub(x = colnames(metaT_genome_weights), pattern = "metaT\\.", replacement = ""))
        }
    }
}



##########################  ** Define functions **  ##########################

get_annotation_descriptions <- function(db_name, functionListDT) {
    ### Annotation descriptions ####
    if  (db_name %in% c('eggNOG.OGs', 'eggNOG.KEGG_Pathway', 'eggNOG.KEGG_Module', 'eggNOG.KEGG_KO', 'kofam.KEGG_Pathway', 'kofam.KEGG_Module', 'kofam.KEGG_KO', 'merged.KEGG_KO', 'dbCAN.EC')){
        annotations <- (fread(funcat_desc_file,
                              header=TRUE,
                              data.table=TRUE,
                              select=c('Funct', 'Description'))
                        [, Description := iconv(Description, from = "ISO-8859-1", to = "UTF-8")]
                       )
        setkey(annotations, Funct)
        annotations <- unique(merge(functionListDT, annotations, by='Funct', all=T)
                              [, Description := paste(sort(unique(Description)), collapse=';'), by=Funct]
                             )
    }
    else if  (db_name %in% c('eggNOG.PFAMs')){

        suppressMessages(library(PFAM.db))

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
        rm(pfam_names, pfam_desc, pfam_table)
    }
    else if  (db_name %in% c('dbCAN.module', 'dbCAN.enzclass')){

        suppressMessages(library(PFAM.db))

        # Get the PFAM identifiers that are mapped to a CAZY
        mapped_keys <- mappedkeys(PFAMCAZY)
        # Convert to a list
        map_list <- as.list(PFAMCAZY[mapped_keys])
        map_table <- data.table(do.call(rbind, map_list), keep.rownames=TRUE)
        setnames(map_table, 'rn', 'Pfam_id')

        # Get unique dbCAN entities for each PFAM key
        map_table <- (
                      map_table
                      [, CAZY_id := paste(unique(unlist(.SD)), collapse = ','), by='Pfam_id']
                      [, c("V1", "V2", "V3", "V4") := NULL]
                     )

        cazy_singleKeys <- strsplit(map_table[['CAZY_id']], split = "\\  |\\,")
        if (length(cazy_singleKeys) > 0) {
            cazy_keyMap <- data.table(
                                      Pfam_id = rep(map_table$Pfam_id, sapply(cazy_singleKeys, length)),
                                      Funct = unlist(cazy_singleKeys))
        }

        # Get PFAM descriptions
        pfam_desc <- as.list(PFAMDE[mappedkeys(PFAMDE)])
        setDT(pfam_desc, keep.rownames = TRUE)
        pfam_desc <- data.table::transpose(pfam_desc, keep.names = "Pfam_id")
        setnames(pfam_desc, "V1", "Description")

        # Merge based on PFAM_ID and delete that join key
        cazy_table <- merge(cazy_keyMap, pfam_desc, by='Pfam_id')[, Pfam_id := NULL]

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

# Create function profiles and return a data.table.
# Write profiles as tsv and phyloseq files, if write_output==TRUE.
# Assumptions:
#   keys and profile are already indexed by gene_id

make_profile_files_batched <- function(keys, profile, file_label, database, annotations, weights, write_output = TRUE, batch_size = 500) {
  
  logmsg("Making ", file_label, " profile files for ", database, decorate=TRUE)
  
  # Make a list of sample-names
  sample_cols <- setdiff(colnames(profile), c('gene_id'))
  
  # Create slices for too large profiles
  mag_start_ind <- which(str_detect(profile$gene_id, "^\\w{10}_0000\\d"))
  if (length(mag_start_ind) != 0){
    tempdt <- data.table(row_id=mag_start_ind, gene_id=profile$gene_id[mag_start_ind], MAG=sub('_.*', '', profile$gene_id[mag_start_ind]))
    mag_start_ind <- {unique(tempdt, by = c("MAG"))}$row_id
    mag_start_ind[[1]] <- 1
    rm(tempdt)
  } else {
    mag_start_ind <- seq(1,length(profile$gene_id), 2000)
  }
  
  row_batches <- mag_start_ind[seq(1, length(mag_start_ind), batch_size)]
  counts_list <- vector("list", length(row_batches))
  
  batch_next = 1
  for (row_i in row_batches) {
    # profile row slice
    batch_next = batch_next + 1
    if (batch_next <= length(row_batches)){
      row_end  <- row_batches[[batch_next]] - 1
    } else {
      row_end  <- dim(profile)[[1]]
    }
    logmsg("start: ", row_i,", end:", row_end)
    
    # Merge abundance and annotation information by gene
    logmsg("Mapping genes to ", database, " items")
    counts <- merge(keys, profile[row_i:row_end], all.x=FALSE, all.y=TRUE) # right join, implicitly by 'gene_id' that is keyed in both
    logmsg("  done")
    logmsg("Profiled genes      ", nrow(profile[row_i:row_end]))
    logmsg("Gene-Function pairs ", nrow(counts))
    
    # From count-per-gene, make count-per-function
    # If necessary, weight the functions by REL of genomes
    if (!is.null(weights)) {
      
      # Summarize by MAG,Funct
      logmsg("Summarizing by MAG-Function pairs")
      counts <- (
        counts
        [, MAG := sub('_.*', '', gene_id)] # Retain the first underscore-delimited field
        [, gene_id := NULL]
        [, lapply(.SD, sum, na.rm=TRUE), by = c('MAG', 'Funct')]
      )
      logmsg("  done")
      
      # Get completeness
      if (database == "kofam.KEGG_Module" & !is.null(completeness_file)) {
        mod_comp <- fread(completeness_file,
                          header=TRUE,
                          sep="\t",
                          data.table=TRUE,
                          select=c('ID', 'module_accession', 'completeness'))
        setkey(mod_comp, completeness)
        mod_comp = mod_comp[completeness >= min_completeness, ]
        setnames(mod_comp, c('ID', 'module_accession'), c('MAG', 'Funct'))
        counts <- merge(mod_comp, counts, by=c('MAG', 'Funct'))[, completeness := NULL]
      }
      
      # Mark the weight-rows using Funct = 'REL'
      weights[, Funct := "REL"]
      
      # Ensure identical column orders before rbind()
      setcolorder(weights, colnames(counts))
      
      # Prepend with weights to prepare for weighing by MAG abundance
      counts <- rbind(weights, counts)
      logmsg("MAG-Function  pairs ", nrow(counts))
      
      # weigh by the given weights per MAG
      weight_MAG_by_proportions <- function(x){
        if (is.character(x[1])){
          return(x)
        }else{
          return(x*x[1])
        }
      }
      
      # Normalize by weighting by 'weights'
      # Remove REL rows
      logmsg("Weighting functions by MAG abundance, and summarizing")
      # Summarize by Funct
      counts <- (
        counts
        [, lapply(.SD, weight_MAG_by_proportions), by = MAG]
        [Funct != 'REL', ]
        [, MAG := NULL]
        [, lapply(.SD, sum, na.rm=TRUE), .SDcols = sample_cols, by = Funct]
      )
      logmsg("  done")
    } else {
      
      # Summarize by Funct
      logmsg("Summarizing functions")
      counts <- (
        counts
        [, gene_id := NULL]
        [, lapply(.SD, sum, na.rm=TRUE), .SDcols = sample_cols, by = Funct]
      )
      logmsg("  done")
    }
    
    # Add to list
    counts_list[[batch_next - 1]] <- counts
    rm(counts)
  }
  
  # Put together sliced counts
  counts <- do.call("rbind", counts_list)
  rm(counts_list)
  gc()
  
  # Summarize by Funct
  logmsg("Summarizing batches functions")
  counts <- (
    counts
    [, lapply(.SD, sum, na.rm=TRUE), .SDcols = sample_cols, by = Funct]
  )
  logmsg("  done")
  
  setkey(counts, 'Funct')
  
  # Estimate 'Unknown' functions - genes with no annotation
  counts <- counts[.(NA), Funct := "Unknown"]
  
  logmsg("Profiled functions  ", nrow(counts))
  
  # Add annotations
  logmsg("Merging function profiles and annotations")
  counts <- merge(annotations, counts, by='Funct', all.x=FALSE, all.y=TRUE) # Right join
  logmsg("  done")
  
  if (write_output) {
    # Write annotated functional profile as tsv
    logmsg("Writing function profile into text file")
    fwrite(counts,
           file=paste0(output_dir, '/', file_label, '.', database ,'.tsv'),
           sep='\t',
           row.names = F,
           quote = F)
    logmsg("  done")
    
    #### Create phyloseq object ####
    
    logmsg("Creating phyloseq object")
    
    my_otu_table <- as.matrix(counts[, -c('Funct', 'Description')])
    rownames(my_otu_table) <- counts$Funct
    
    # taxa
    my_tax_table <- as.matrix(counts)
    rownames(my_tax_table) <- counts$Funct
    
    physeq <- phyloseq(otu_table(my_otu_table, taxa_are_rows = T),
                       tax_table(my_tax_table),
                       sample_metadata)
    logmsg("  done")
    
    # Write phyloseq object
    logmsg("Writing phyloseq object")
    qsave(physeq,
          preset = "high",
          file = paste0(phyloseq_dir, '/', file_label, '.', database ,'.qs'))
    logmsg("  done")
  }
  
  # Index on Funct
  setkey(counts, Funct)
  
  # Remove Unknown functions
  # It has been written out in tsv files for metaG and metaT
  # And it should NOT be calculated for FE
  counts <- counts[!J("Unknown"),]
  
  return(counts)
}

get_gene_count <- function(filename) {

    # Set key column as 'gene_id'
    # Convert NA's to 0's - not needed since no NA in input (12.03.2025)
    # Remove features with zero sum across all samples - not needed since this is done in input files already (12.03.2025)
    logmsg("Reading profile")
    counts = fread(filename, sep="\t")
    setnames(counts, "ID", "gene_id")
    setkey(counts, gene_id)
    logmsg("  done")

    return(counts)
}

get_annotations <- function(filename) {

    # Get annotations for only this funcat_name and set rownames as 'gene_id'
    # Remove empty annotations
    logmsg("Reading annotations")
    annotations = fread(filename, sep="\t", select=c("ID", funcat_name))
    setnames(annotations, "ID", "gene_id")
    setkeyv(annotations, funcat_name)
    annotations = annotations[!J(c(NA, "", "-"))]
    logmsg("  done")

    return(annotations)
}

get_metadata <- function(filename) {

    # Get metadata table
    metadata = fread(filename, sep="\t")

    return(metadata)
}

# *************************** Use GA and GT raw counts to generate FE profile *************************** ####
# The genes will be clustered by function IDs.
# TPM normalization - sum all TPMs for the member genes
#                     (there will be redundancy and total might be over 1 million)
# MG normalization  - weight individual member genes by the rel. abundance of the genome and sum.
# FE profile will be generated: (FT_norm) / (FA_norm)
#Raw gene abundances #####

# Generate output directories ####
dir.create(file.path(output_dir), showWarnings = FALSE)
visual_dir=paste0(output_dir,'/plots/')
dir.create(file.path(visual_dir), showWarnings = FALSE)
phyloseq_dir=paste0(output_dir,'/phyloseq_obj/')
dir.create(file.path(phyloseq_dir), showWarnings = FALSE)


logmsg("Reading input files", decorate=TRUE)

# Initialize some variables
metadata_df <- NULL
gene_annotation <- NULL
ga_fa_df <- NULL
gt_ft_df <- NULL
ge_fe_df <- NULL

gene_annotation <- NULL

# metaG
if (omics %like% "metaG") {

    # gene profile
    metaG_profile <- get_gene_count(metaG_gene_profile)
    ga_fa_df <- data.frame(DB='genes', feature_n=nrow(metaG_profile), feature = 'Genes')

    # gene annotations
    gene_annotation <- get_annotations(gene_annotation_file)

    # metadata
    metadata_df <- get_metadata(metadata_file)

    if (normalization == 'MG') {
        # subset columns in weights
        intersect_columns    = union(c("MAG"), intersect(colnames(metaG_profile), colnames(metaG_genome_weights)))
        metaG_genome_weights = metaG_genome_weights[, intersect_columns, with=FALSE]
    }
}


# metaT
if (omics %like% "metaT") {


    # gene profile
    metaT_profile <- get_gene_count(metaT_gene_profile)
    gt_ft_df <- data.frame(DB='genes', feature_n=nrow(metaT_profile), feature = 'Genes')

    # gene annotations
    metaT_gene_annotation <- get_annotations(gene_annotation_file)

    # Only need entries that are not already present in gene_annotation
    if (omics == 'metaT') {
        gene_annotation <- metaT_gene_annotation
        metadata_df <- get_metadata(metadata_file)
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

    if (normalization == 'MG') {
        # subset columns in weights
        intersect_columns    = union(c("MAG"), intersect(colnames(metaT_profile), colnames(metaT_genome_weights)))
        metaT_genome_weights = metaT_genome_weights[, intersect_columns, with=FALSE]
    }

    rm(metaT_gene_annotation)
}

########################################################
# For metaG_metaT, subset profiles for common samples
########################################################

if (omics == 'metaG_metaT') {
    logmsg("began getting common samples")
    # Columns shared b/w metaG and metaT. This includes 'gene_id', which is also needed in final table.
    intersect_columns = intersect(colnames(metaG_profile), colnames(metaT_profile))
    metaG_profile = metaG_profile[, intersect_columns, with=FALSE]
    metaT_profile = metaT_profile[, intersect_columns, with=FALSE]

    if (normalization == 'MG') {
        # Filter weights to get the same columns, but skip 'gene_id'
        intersect_columns    = union(c("MAG"), intersect_columns[intersect_columns != 'gene_id'])
        metaG_genome_weights = metaG_genome_weights[, intersect_columns, with=FALSE]
        metaT_genome_weights = metaT_genome_weights[, intersect_columns, with=FALSE]
    }
    logmsg("done  getting common samples")
}

# By now, we will have:
#   metadata_df
#   gene_annotation
#   ga_fa_df, gt_ft_df, ge_fe_df

# Reorder metadata by same order in profile file, so that it can be matched properly when ID is removed
# Easiest way is to do an inner join, which will also nicely fail if a sample doesnt have metadata

logmsg("began reordering metadata rows")
col_names = names(metadata_df)
if (omics == 'metaT') {
    anchor_df = metaT_profile %>%
                    dplyr::select(-gene_id) %>%
                    head(2) %>%
                    t() %>%
                    as.data.frame() %>%
                    tibble::rownames_to_column("sample_alias")
} else {
    anchor_df = metaG_profile %>%
                    dplyr::select(-gene_id) %>%
                    head(2) %>%
                    t() %>%
                    as.data.frame() %>%
                    tibble::rownames_to_column("sample_alias")
}
metadata_df = dplyr::inner_join(anchor_df, metadata_df) %>%
                    dplyr::select(any_of(col_names))
logmsg("done  reordering metadata rows")

## metadata for phyloseq
sample_metadata <- sample_data(metadata_df)
rownames(sample_metadata) <- sample_metadata[["sample_alias"]]


#################################### FUNCTION PROFILING  ####################################

# Filter tax.table based on features in metaG, if metaG_metaT. metaG/metaT have been filtered earlier.

#if (omics == 'metaG_metaT') {
#    gene_annotation <- gene_annotation[gene_id %in% metaG_profile$gene_id, ]
#}


if (nrow(gene_annotation) > 0) {

    ################################################
    # Prepare annotations and their descriptions
    ################################################

    # Replicate the gene annotation table by having N rows for gene with N annotations.

    singleKeys <- strsplit(gene_annotation[[funcat_name]], split = "\\; |\\;|\\, |\\,")
    Gene2FuncMap <- data.table(gene_id = rep(gene_annotation$gene_id, sapply(singleKeys, length)),
                               Funct = unlist(singleKeys))

    # Since KEGG annotations can have map12345 and ko12345 for the same pathway,
    # remove such redundancy.
    if (funcat_name %in% c("eggNOG.KEGG_Pathway", "kofam.KEGG_Pathway")) {
        Gene2FuncMap <- Gene2FuncMap[, Funct := gsub('ko', 'map', Funct)]
    }

    # Get unique annotations
    Gene2FuncMap <- unique(Gene2FuncMap)
    setkey(Gene2FuncMap, gene_id)

    rm(singleKeys, gene_annotation)

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
        metaG_counts <- make_profile_files_batched(keys=Gene2FuncMap,
                                           profile=metaG_profile,
                                           file_label="FA",
                                           database=funcat_name,
                                           annotations=annot_dt,
                                           weights=metaG_genome_weights,
                                           write_output=(omics=="metaG"))
        # Record feature count
        ga_fa_df <- rbind(ga_fa_df, c(funcat_name, nrow(metaG_counts), 'Functions'))

    }

    # metaT or metaG_metaT
    if (omics %like% "metaT") {
        metaT_counts <- make_profile_files_batched(keys=Gene2FuncMap,
                                           profile=metaT_profile,
                                           file_label="FT",
                                           database=funcat_name,
                                           annotations=annot_dt,
                                           weights=metaT_genome_weights,
                                           write_output=(omics=="metaT"))
        # Record feature count
        gt_ft_df <- rbind(gt_ft_df, c(funcat_name, nrow(metaT_counts), 'Functions'))

    }

    rm(Gene2FuncMap)

    # metaG_metaT
    if (omics == 'metaG_metaT') {

        logmsg("Making FE profile files for ", funcat_name, decorate=TRUE)
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
              preset = "high",
              file = paste0(phyloseq_dir, '/FE.', funcat_name ,'.qs'))

        # Free up memory
        logmsg("Freeing memory")
        rm(feature_ids, FE_annotations, FE_physeq)
        gc()
    }


    ###########
    # Plot PCA
    ###########

    logmsg("Drawing PCoA plots for FE in ", funcat_name, decorate=TRUE)

    # PCA - metaG ####
    if (omics == 'metaG') {
        # Remove Funct column
        metaG_counts <- metaG_counts[, `:=` (Funct = NULL, Description = NULL) ]

        # Make PCA
        prepare_PCA(profile=metaG_counts, title=paste0("function abundance - ", funcat_name), datatype=paste0('FA.', funcat_name), color=color_factor, shape=shape_factor, label=label_factor, metadata=sample_metadata)

        # Write feature count
        fwrite(ga_fa_df, file=paste0(output_dir, '/GA_FA_features.', funcat_name, '.tsv'), sep='\t', row.names = F, quote = F)

        rm(metaG_counts)
    }

    # PCA - metaT ####
    if (omics == 'metaT') {
        # Remove Funct column
        metaT_counts <- metaT_counts[, `:=` (Funct = NULL, Description = NULL) ]

        # Make PCA
        prepare_PCA(profile=metaT_counts, title=paste0("function transcript - ", funcat_name), datatype=paste0('FT.', funcat_name), color=color_factor, shape=shape_factor, label=label_factor, metadata=sample_metadata)

        # Write feature count
        fwrite(gt_ft_df, file=paste0(output_dir, '/GT_FT_features.', funcat_name, '.tsv'), sep='\t', row.names = F, quote = F)

        rm(metaT_counts)
    }

    # PCA - GE ####
    if (omics == 'metaG_metaT') {
        # Record feature count
        n_row <- nrow(function_expression)
        ge_fe_df <- rbind(ge_fe_df, c(funcat_name, n_row, 'Functions'))

        # Remove Funct column
        function_expression <- function_expression[, Funct := NULL][, lapply(.SD, function(x) ifelse(is.na(x), 0, x))]

        # Make PCA
        prepare_PCA(profile=function_expression, title=paste0("function expression - ", funcat_name), datatype=paste0('FE.', funcat_name), color=color_factor, shape=shape_factor, label=label_factor, metadata=sample_metadata)

        # Write feature count
        fwrite(ge_fe_df, file=paste0(output_dir, '/GE_FE_features.', funcat_name, '.tsv'), sep='\t', row.names = F, quote = F)

        rm(function_expression)
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

# Garbage-collect
gc()

logmsg("  DONE  ", decorate=TRUE)
