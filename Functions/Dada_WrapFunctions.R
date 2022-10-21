#######################################
### FUNCTION: Dada2_wrap
#######################################

# Function that runs the dada2 pipeline up to the sequence table (seqtab)
## Input
# path: The path to the folder containing the sample folders with the fastq files. 
# F_pattern: a regular expression to find the fastq files with the forward reads in the sample folders
# R_pattern: a regular expression to find the fastq files with the reverse reads in the sample folders
# path2: default = NULL, the function creates three folders Dada_Plots, Dada_Data and Dada_FilteredFastqs, it creates them by default in the path folder, when path2 is given 
# the folders will instead be generated in path2
# trimLeft: = trimLeft from the fastqPairedFilter command (how many nucleotids will be clipped from the beginning of the reads)
# truncLen: = truncLen from the fastqPairedFilter command (where the reads will be truncated)
# maxN: = maxN from the fastqPairedFilter command, After truncation, sequences with more than maxN Ns will be discarded, NB: Dada2 requires no Ns!!
# maxEE: = maxEE from the fastqPairedFilter command, the maximum expected error allowed to let a read pass the filter
# trunQ: = truncQ from the fastqPairedFilter command, Truncate reads at the first instance of a quality score less than or equal to truncQ, ? actually not sure if these reads are then discarded?
#    # NB: currently run before trimLeft/truncLen which does the size check. So unless truncated after the truncLen
# the size check will kick the affected reads out, so makes little sense, see: https://github.com/benjjneb/dada2/issues/140 
# nreadsLearn = 1e+06, # the number of reads (distributed over far less unique reads) used to learn the error matrixes, i.e. nreads in dada2:::learnErrors 
# Should account for a number of samples that in total have about 1 milliion filtered reads (number will be given in the log file)
# err_F: the error matrix for the dada command for the forward reads. Default NULL then the error matrix will be estimated (using SelfConsit = TRUE)
# err_R: analogous to err_F for the reverse reads
# minOverlap: = minOverlap from the mergePairs command
# maxMismatch: = maxMismatch from the mergePairs command
# F_QualityStats and R_QualityStats: NULL by default, if given e.g. from Dada_QualityCheck the quality stats collection part is jumped over
# filtFs and filtRs: NULL by default, if given and FilteredFolder with files exist, Filtering can be jumped over.
# pool: pool option of dada command
## Output
# PLOTS:
# Several Plots saved as pdfs in generated Dada_Plots folder # NB: the plots are currently assuming 250 nt
# DATA AND LOG FILE:
# Data and the file DadaWrapper.log are saved in the generated folder: Dada_Data
# DenoisedData.RData contains: seqtab (final sequence table), mergers, mergers.nochim, bimFs, bimRs, 
# ReadSummary (data frame summarising the reads and amplicons at the different stages),
# QualityStats.RData contains PackageVersions, F_QualityStats, R_QualityStats
# FILTERED FASTQ files
# are stored in the generated folder Dada_FilteredFastqs
# NB the filtered fastq files are currently named "SampleName"_F_Filtered.fastq.gz and "SampleName"_R_Filtered.fastq.gz

Dada2_wrap <- function(path, F_pattern, R_pattern, path2 = NULL,
                       trimLeft = c(10,10), truncLen = c(220, 160), 
                       maxN = 0, maxEE = 2, truncQ = 2,
                       nreadsLearn = 1e+06,
                       err_F = NULL,
                       err_R = NULL,
                       minOverlap = 20,
                       maxMismatch = 0,
                       F_QualityStats = NULL,
                       R_QualityStats = NULL,
                       filtFs = NULL,
                       filtRs = NULL,
                       pool = FALSE,
                       justConcatenate=FALSE) {
  
  ##############################
  ### call the required packages
  ##############################
  
  ## dada2:
  source("https://bioconductor.org/biocLite.R")
  try(library(dada2), biocLite("dada2"))
  
  ## Short Read
  try(library(ShortRead), biocLite("ShortRead"))
  
  ## ggplot2
  try(library(ggplot2), install.packages("ggplot2"))
  
  ## dplyr
  try(library(dplyr), install.packages("dplyr"))
  
  ## dplyr
  try(library(tidyr), install.packages("tidyr"))
  
  
  message("*********************** All packages loaded ***********************
          ********************************************************************")
  
  ##############################
  ### save the R and package Versions
  ##############################
  # NB: outputs an error and stops function if a Package is not installed
  RVer <- R.Version()
  RVer <- RVer$version.string
  PackageVersions <- data.frame(Package = c("R", "dada2", "ShortRead", "ggplot2", "dplyr", "tidyr"),
                                Version = c(RVer,
                                            as.character(packageVersion("dada2")),
                                            as.character(packageVersion("ShortRead")),
                                            as.character(packageVersion("ggplot2")),
                                            as.character(packageVersion("dplyr")),
                                            as.character(packageVersion("tidyr"))))
  
  message(paste(PackageVersions$Package, ": ", PackageVersions$Version, "; ", sep= ""))
  
  
  ##############################
  ### Construct character vectors to the FW and RV fastq files
  ##############################
  
  if(is.null(path2)){
    
    path2 = path
  }
  
  
  folders <- list.dirs(path, recursive = FALSE, full.names = FALSE)
  
  if(length(folders) != 0) {
    
    SampleNames <- folders
    
    if(sum(grepl("^Dada", SampleNames)) != 0){
      warning("**There are folders starting with \"Dada\" in your path folder, maybe you have run the function on this path folder before\n.
              The folders starting with Dada will not be considered sample folders!!\n Files within Dada_Data Dada_FilteredFastqs and Dada_Plots will be overwritten**")
    }
    
    # exclude Folders starting with "Dada" from the folders considered as sample fodlers
    if(length(grep("^Dada", SampleNames))!=0) {
      SampleNames <- SampleNames[-grep("^Dada", SampleNames)]
    }
    
    F_fastq <- character(length = length(SampleNames))
    R_fastq <- character(length = length(SampleNames))
    
    for (i in 1:length(SampleNames)) {
      CurrentPath <- file.path(path, SampleNames[i])
      
      if(sum(grepl(F_pattern, list.files(CurrentPath))) == 0) {
        stop(paste("F_pattern fits no file in ", CurrentPath))
      }
      if(sum(grepl(F_pattern, list.files(CurrentPath))) > 1) {
        stop(paste("F_pattern fits several files in ", CurrentPath))
      }
      if(sum(grepl(R_pattern, list.files(CurrentPath))) == 0) {
        stop(paste("R_pattern fits no file in ", CurrentPath))
      }
      if(sum(grepl(R_pattern, list.files(CurrentPath))) > 1) {
        stop(paste("R_pattern fits several files in ", CurrentPath))
      }
      F_fastq[i] <- file.path(CurrentPath, list.files(CurrentPath)[grepl(F_pattern, list.files(CurrentPath))])
      R_fastq[i] <- file.path(CurrentPath, list.files(CurrentPath)[grepl(R_pattern, list.files(CurrentPath))])
      
    }
    
    } else {
      
      stop("No sample folders were found in the given path! Currently the Dada2_wrap function can only handle the situation where the fastq files are in separate folders for each sample.
           These folders have to be in the \"path\" folder. Other situations have to be added.")
    }
  
  ##############################
  ### Start the log file
  ##############################
  
  DataFolder <- file.path(path2, "Dada_Data")
  
  # ATTENTION: FUNCTION DELETES ALL FILES ALREADY PRESRENT in DataFolder
  if(file.exists(DataFolder)){
    file.remove(list.files(DataFolder, full.names = TRUE))
  }
  
  
  dir.create(DataFolder, showWarnings = FALSE)
  
  if(!file.exists(DataFolder)){
    stop("DataFolder could not be created! something wrong with path2, maybe permission denied.")
  }
  
  ptm <- proc.time()
  
  LogFile <- file.path(DataFolder, "DadaWrapper.log")
  cat("Time after package installation: ", file = LogFile)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  cat("Package Versions: ", file = LogFile, append = TRUE)
  cat(paste(PackageVersions$Package, ": ", PackageVersions$Version, "; ", sep= ""), file = LogFile, append = TRUE)
  
  
  # Collect inputs in a data frame for saving it later
  Input <- list(path = path, F_pattern = F_pattern, R_pattern = R_pattern, path2 = path2,
                trimLeft = trimLeft, truncLen = truncLen, 
                maxN = maxN, maxEE = maxEE, truncQ = truncQ,
                nreadsLearn = nreadsLearn, err_F = err_F, err_R = err_R,
                minOverlap = minOverlap, maxMismatch = maxMismatch, F_QualityStats = F_QualityStats,
                R_QualityStats = R_QualityStats, filtFs = filtFs, filtRs = filtRs)
  
  
  ##############################
  ### Determine the quality scores and save the stats in Data folder
  ##############################
  
  # Collect quality Score data of the fastQ files and store a data.frame for each sample in the lists FW_QualityStats and RV_QualityStats ####################
  
  if (is.null(F_QualityStats) || is.null(R_QualityStats)) {
    
    F_QualityStats <- list()
    R_QualityStats <- list()
    
    for (i in seq_along(F_fastq)) {
      
      Current_FWfq <- F_fastq[i]
      Current_dfFW <- qa(Current_FWfq, n = 1e06)[["perCycle"]]$quality
      # df is a data frame containing for each cycle (nt) the distribution of Quality scores
      # e.g. Cycle 1 had 7 different quality scores, then 7 rows of cycle one, for each score the count says how many reads had this score
      Current_dfFW <- dplyr::group_by(Current_dfFW, Cycle)
      
      Current_dfQStatsFW <- dplyr::summarise(
        Current_dfFW,
        NoReads = sum(Count),
        Mean_QS = sum(Count*Score)/sum(Count),
        SD_QS = sqrt(sum(Count*((Score-Mean_QS)^2))/(NoReads-1)),
        Median_QS = Score[which(cumsum(Count)/sum(Count) >= .5)][[1]],
        q25_QS = Score[which(cumsum(Count)/sum(Count) >= .25)][[1]],
        q75_QS = Score[which(cumsum(Count)/sum(Count) >= .75)][[1]])
      
      # Check that all reads are of same length
      # x <- range(Current_dfQStatsFW$NoReads)
      # if (!all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5)) {
      #         stop("Not all reads of same length in file", F_fastq[i])
      # }
      
      F_QualityStats[[i]] <- as.data.frame(Current_dfQStatsFW)
      # as.data.frame to un-dplyr the data.frame
      
      # collect the same stats for the RV FastQ files
      Current_RVfq <- R_fastq[i]
      Current_dfRV <- qa(Current_RVfq, n = 1e06)[["perCycle"]]$quality
      
      Current_dfRV <- dplyr::group_by(Current_dfRV, Cycle)
      
      Current_dfQStatsRV <- dplyr::summarise(
        Current_dfRV,
        NoReads = sum(Count),
        Mean_QS = sum(Count*Score)/sum(Count),
        SD_QS = sqrt(sum(Count*((Score-Mean_QS)^2))/(NoReads-1)),
        Median_QS = Score[which(cumsum(Count)/sum(Count) >= .5)][[1]],
        q25_QS = Score[which(cumsum(Count)/sum(Count) >= .25)][[1]],
        q75_QS = Score[which(cumsum(Count)/sum(Count) >= .75)][[1]])
      
      ## Check that all reads are of same length
      # x <- range(Current_dfQStatsRV$NoReads)
      #  if (!all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5)) {
      #          stop("Not all reads of same length in file", R_fastq[i])
      #  }
      # 
      R_QualityStats[[i]] <- as.data.frame(Current_dfQStatsRV)
      
      rm(Current_FWfq, Current_dfFW, Current_dfQStatsFW,
         Current_RVfq, Current_dfRV, Current_dfQStatsRV, x)
      
    }
    
    # add the sample names as names to the lists:
    names(F_QualityStats) <- SampleNames
    names(R_QualityStats) <- SampleNames
    
    message("*********************** Quality Stats Collected ***********************
            ********************************************************************")
    cat("\n*** Quality Stats Collected ***", file = LogFile, append = TRUE)
    
  } else {
    
    if((names(F_QualityStats) != SampleNames) || (names(R_QualityStats) != SampleNames)) {
      
      stop("The given F_QualityStats or R_QualityStats do not fit to the SampleNames in path!")
    }
    
    message("*********************** Quality Stats given ***********************
            ********************************************************************")
    cat("\n*** Quality Stats given***", file = LogFile, append = TRUE)
    
  }
  
  save(PackageVersions, F_QualityStats, R_QualityStats, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  TimePassed <- proc.time()-ptm
  cat("\nTime after Quality Stats collection: ", file = LogFile, append = TRUE)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
  cat("\n*** Start generating Quality Plots ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Generate and save some quality plots
  ##############################
  
  PlotFolder <- file.path(path2, "Dada_Plots")
  
  # ATTENTION: FUNCTION DELETES ALL FILES ALREADY PRESRENT in PlotFolder
  if(file.exists(PlotFolder)){
    file.remove(list.files(PlotFolder, full.names = TRUE))
  }
  
  dir.create(PlotFolder, showWarnings = FALSE)
  
  
  pdf(file = file.path(PlotFolder, "MedianQScore_F_AllNucleotides.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(F_QualityStats, SampleNames))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_F_150to240.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(F_QualityStats, SampleNames, xlim_low = 150, xlim_high = 240))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_R_AllNucleotides.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(R_QualityStats, SampleNames))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_R_150to240.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(R_QualityStats, SampleNames, xlim_low = 150, xlim_high = 240))
  dev.off()
  
  message("*********************** Plots generated start filtering ***********************
          ********************************************************************")
  
  cat("\n*** Quality Plots generated ***", file = LogFile, append = TRUE)
  cat("\n*** Start Filtering ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Filtering
  ##############################
  
  FilteredFolder <- file.path(path2, "Dada_FilteredFastqs")
  
  
  if(is.null(filtFs) || is.null(filtRs)){
    
    ## Not sure if this is wanted, deleting all files that are already in the DataFolder folder
    if(file.exists(FilteredFolder)){
      file.remove(list.files(FilteredFolder, full.names = TRUE))
    }
    
    dir.create(FilteredFolder, showWarnings = FALSE)
    
    filtFs <- file.path(FilteredFolder, paste0(SampleNames, "_F_Filtered.fastq.gz"))
    filtRs <- file.path(FilteredFolder, paste0(SampleNames, "_R_Filtered.fastq.gz"))
    names(filtFs) <- SampleNames
    names(filtRs) <- SampleNames
    
    
    for(i in seq_along(F_fastq)) {
      
      message(paste("**Filtering sample No:", i, "called:", SampleNames[i]), " **")
      cat(paste("\n**Filtering sample No:", i, "called:", SampleNames[i]), file = LogFile, append = TRUE)
      fastqPairedFilter(c(F_fastq[i], R_fastq[i]), c(filtFs[i], filtRs[i]), trimLeft = trimLeft, 
                        truncLen = truncLen, maxN = maxN, maxEE = maxEE, truncQ = truncQ, 
                        compress=TRUE, verbose=TRUE)
    }
    
    # check if files have been created
    if(!all(file.exists(filtFs)) || !all(file.exists(filtRs))) {
      cat("\n*** ERROR: Not all filtered files were created, maybe trimming impossible***", file = LogFile, append = TRUE)
      #stop("** Not all filtered files were created, maybe trimming impossible")
      
    }
    
    message("*********************** Filtering Done ***********************
            ********************************************************************")
    cat("\n*** Filtering Done ***", file = LogFile, append = TRUE)
    
  } else {
    
    if((names(filtFs) != SampleNames) || (names(filtRs) != SampleNames)) {
      
      cat("\n*** ERROR: The given filtFs or filtRs do not have sample names!***", file = LogFile, append = TRUE)
      stop("The given filtFs or filtRs do not have sample names!")
    }
    
    if(!all(file.exists(filtFs)) || !all(file.exists(filtRs))) {
      cat("\n*** ERROR: Not all files in the given filtFs or filtRs existed***", file = LogFile, append = TRUE)
      stop("** Not all files in the given filtFs or filtRs existed")
      
    }
    
    message("*********************** Filtered files were given ***********************
            ********************************************************************")
    cat("\n*** Filtered files were given ***", file = LogFile, append = TRUE)
    
    
  }
  
  
  save(PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  cat("\nTime after filtering step: ", file = LogFile, append = TRUE)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  TimePassed <- proc.time()-ptm
  cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
  cat("\n*** Start estimating err_F if not given ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Estimate err_F matrix
  ##############################
  
  if(is.null(err_F)){
    
    # the new learnErrors from dada2 allows setting nreads (usually 1 million), it then reads in as many samples until
    # the number of unique sequences (drp$uniques) reaches nreads.
    # learnErrorsAdj is basically dada2::learnErrors, just adds the NREADS from the function to the output list and also adds
    # dreplicated (drps) and denoised (dds) datasets if all samples were used for error matrix estimation. If not drps and dds
    # remain NULL
    errorsFW <- learnErrorsAdj(fls = filtFs, nreads = nreadsLearn, multithread = TRUE)
    
    message(paste("**", errorsFW$NREADS, "reads have been used for err_F estimation **"))
    
    err_F <- errorsFW$err_out
    
    pdf(file = file.path(PlotFolder, "errorRates_F.pdf"), width = 10, height = 10)
    print(plotErrors(errorsFW, nominalQ=TRUE))
    dev.off()
    
    save(err_F, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
    
    message("*********************** err_F has been estimated ***********************
            ********************************************************************")
    cat("\n*** err_F has been estimated ***", file = LogFile, append = TRUE)
    cat(paste("\nReads used for err_F estimation: ", errorsFW$NREADS), file = LogFile, append = TRUE)
    TimePassed <- proc.time()-ptm
    cat("\nTime after err_F estimation: ", file = LogFile, append = TRUE)
    cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
    cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
    cat("\n*** Start estimating err_R if not given ***", file = LogFile, append = TRUE)
    
  }
  
  # in case both err_F and err_R had been given, they would not have been saved if not:
  save(err_F, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  ##############################
  ### Estimate err_R matrix
  ##############################
  
  if(is.null(err_R)){
    
    
    errorsRV <- learnErrorsAdj(fls = filtRs, nreads = nreadsLearn, multithread = TRUE)
    
    message(paste("**", errorsRV$NREADS, "reads have been used for err_R estimation **"))
    
    err_R <- errorsRV$err_out
    
    pdf(file = file.path(PlotFolder, "errorRates_R.pdf"), width = 10, height = 10)
    print(plotErrors(errorsRV, nominalQ=TRUE))
    dev.off()
    
    save(err_F, err_R, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
    
    message("*********************** err_R has been estimated ***********************
            ********************************************************************")
    cat("\n*** err_R has been estimated ***", file = LogFile, append = TRUE)
    cat(paste("\nReads used for err_R estimation: ", errorsRV$NREADS), file = LogFile, append = TRUE)
    TimePassed <- proc.time()-ptm
    cat("\nTime after err_R estimation: ", file = LogFile, append = TRUE)
    cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
    cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
    cat("\n*** Start denoising data, bimera detection, and merging of reads into amplicons***", file = LogFile, append = TRUE)
    
  }
  
  save(err_F, err_R, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  ##############################
  ### Denoising (dada command for all samples)
  ##############################
  
  message("*********************** Start denoising and bimera detection ***********************
          ********************************************************************")
  
  if (pool) { # very time and memory consuming but might help identifying rare SVs
    
    drp_F <- derepFastq(filtFs, verbose = T)
    dd_F <- dada(drp_F, err = err_F, pool = pool, multithread = TRUE)
    drp_R <- derepFastq(filtRs, verbose = TRUE)
    dd_R <- dada(drp_R, err = err_R, pool = pool, multithread = TRUE)
    
    
    # NB: the following demands that the samples were denoised in the given order, I therefore removed tha randomize option from learnErrorsAdj
    names(dd_F) <- SampleNames
    names(drp_F) <- SampleNames
    names(dd_R) <- SampleNames
    names(drp_R) <- SampleNames
    
    # only for the ReadSummary later
    Uniques_F <- sapply(drp_F, function(x) length(x$uniques))
    Uniques_R <- sapply(drp_R, function(x) length(x$uniques))
    Denoised_F <- sapply(dd_F, function(x) length(x$denoised))
    Denoised_R <- sapply(dd_R, function(x) length(x$denoised))
    NoFilteredReads <- sapply(drp_F, function(x) sum(x$uniques)) # would be the same for drp_R, or sum dd_F$denoised
    
    mergers <- mergePairs(dd_F, drp_F, dd_R, drp_R, minOverlap = minOverlap, maxMismatch = maxMismatch) #!justConcatenate=TRUE
    
    rm(drp_F, drp_R, dd_F, dd_R)
    
    
  } else if (exists("errorsFW", inherits = FALSE) && exists("errorsRV", inherits = FALSE) && # do not do dereplication and denoising again if all samples had been used for determining the error matrixes
             !is.null(errorsFW$dds) && !is.null(errorsRV$dds)) {
    
    dd_F <- errorsFW$dds
    drp_F <- errorsFW$drps
    dd_R <- errorsRV$dds
    drp_R <- errorsRV$drps
    
    # just to not save the huge files:
    errorsFW$dds <- NULL
    errorsFW$drps <- NULL
    errorsRV$dds <- NULL
    errorsRV$drps <- NULL
    
    # NB: the following demands that the samples were denoised in the given order, I therefore removed tha randomize option from learnErrorsAdj
    names(dd_F) <- SampleNames
    names(drp_F) <- SampleNames
    names(dd_R) <- SampleNames
    names(drp_R) <- SampleNames
    
    # only for the ReadSummary later
    Uniques_F <- sapply(drp_F, function(x) length(x$uniques))
    Uniques_R <- sapply(drp_R, function(x) length(x$uniques))
    Denoised_F <- sapply(dd_F, function(x) length(x$denoised))
    Denoised_R <- sapply(dd_R, function(x) length(x$denoised))
    NoFilteredReads <- sapply(drp_F, function(x) sum(x$uniques)) # would be the same for drp_R, or sum dd_F$denoised
    
    mergers <- mergePairs(dd_F, drp_F, dd_R, drp_R, minOverlap = minOverlap, maxMismatch = maxMismatch)
    
    rm(drp_F, drp_R, dd_F, dd_R)
    
  } else {
    
    if (exists("errorsFW", inherits = FALSE)) {
      errorsFW$dds <- NULL
      errorsFW$drps <- NULL
    }
    if (exists("errorsRV", inherits = FALSE)) {
      errorsRV$dds <- NULL
      errorsRV$drps <- NULL
    }
    
    mergers <- vector("list", length(SampleNames))
    names(mergers) <- SampleNames
    NoFilteredReads <- vector("numeric", length(SampleNames))
    names(NoFilteredReads) <- SampleNames
    Uniques_F <- vector("numeric", length(SampleNames))
    names(Uniques_F) <- SampleNames
    Uniques_R <- vector("numeric", length(SampleNames))
    names(Uniques_R) <- SampleNames
    Denoised_F <- vector("numeric", length(SampleNames))
    names(Denoised_F) <- SampleNames
    Denoised_R <- vector("numeric", length(SampleNames))
    names(Denoised_R) <- SampleNames
    
    for(sam in SampleNames) {
      cat("Processing:", sam, "\n")
      cat(paste("\nDenoising sample:", sam), file = LogFile, append = TRUE)
      derepF <- derepFastq(filtFs[[sam]])
      NoFilteredReads[sam] <- sum(derepF$uniques)
      ddF <- dada(derepF, err=err_F, multithread=TRUE) 
      Uniques_F[sam] <- length(derepF$uniques)
      Denoised_F[sam] <- length(ddF$denoised)
      
      derepR <- derepFastq(filtRs[[sam]])
      ddR <- dada(derepR, err=err_R, multithread=TRUE)
      Uniques_R[sam] <- length(derepR$uniques)
      Denoised_R[sam] <- length(ddR$denoised)
      merger <- mergePairs(ddF, derepF, ddR, derepR, minOverlap = minOverlap, maxMismatch = maxMismatch)
      mergers[[sam]] <- merger
      
      rm(derepF, derepR, ddF, ddR, merger)
    }
    
  }
  
  
  
  
  if(!all((sapply(1:length(mergers), function(i) dim(mergers[[i]])[1]))!=0)){
    message("**Not in all samples merged amplicons could be found! maybe minOverlap too strict??**")
    cat("\n*** Not in all samples merged amplicons could be found! maybe minOverlap too strict?? ***", file = LogFile, append = TRUE)
  }
  
  if(all((sapply(1:length(mergers), function(i) dim(mergers[[i]])[1]))==0)){
    cat("\n*** ERROR: In none of the samples merged amplicons could be found! maybe minOverlap too strict?? ***", file = LogFile, append = TRUE)
    stop("**In none of the samples merged amplicons could be found! maybe minOverlap too strict??**")
  }
  
  message("*********************** all samples denoised mergerd amplicons generated ***********************
          ********************************************************************")
  cat("\n*** all samples denoised and mergerd amplicons generated ***", file = LogFile, append = TRUE)
  TimePassed <- proc.time()-ptm
  cat("\nTime after denoising: ", file = LogFile, append = TRUE)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
  cat("\n*** Start generating Sequence table and removing bimeras ***", file = LogFile, append = TRUE)
  
  
  ##############################
  ### Generating sequence table
  ##############################
  
  seqtab <- makeSequenceTable(mergers)
  
  
  ##############################
  ### Bimera removal
  ##############################
  
  # removes denoised sequence FW RV combi that was paired in merger for which the FW or the RV sequence was identified as bimera
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  # 1-(dim(seqtab.nochim)[2]/dim(seqtab)[2]) # The fraction chimeras made up of detected sequence variants
  # 1-(sum(seqtab.nochim)/sum(seqtab)) # The fraction of total reads that were chimeras (over all samples)
  
  message("*********************** Seq table generated and Bimeras removed ***********************
          ********************************************************************")
  cat("\n*** seqtable generated, bimeras removed ***", file = LogFile, append = TRUE)
  cat("\n*** Start saving the data ***", file = LogFile, append = TRUE)
  
  
  # generate also a read summary data frame
  ReadSummary <- data.frame(Sample = SampleNames, NoReads = sapply(F_QualityStats, function(df){df$NoReads[1]}))
  # NB, No total reads would be the same from R_QualityStats, could be used for sanity check
  ReadSummary$FilteredReads <- NoFilteredReads
  ReadSummary$MergedReads <- rowSums(seqtab)
  ReadSummary$MergedReadsWOBimera <- rowSums(seqtab.nochim)
  ReadSummary$UniqueSequences_F <- Uniques_F
  ReadSummary$DenoisedSequences_F <- Denoised_F
  ReadSummary$UniqueSequences_R <- Uniques_R
  ReadSummary$DenoisedSequences_R <- Denoised_R
  ReadSummary$Unique_Amplicons <- rowSums(seqtab != 0)
  ReadSummary$UniqueAmpliconsWOBimera <- rowSums(seqtab.nochim != 0)
  rownames(ReadSummary) <- NULL
  
  # I now also save errorsFW and RV to be able to recapitulate the error plots in later analyses
  if (!exists("errorsFW", inherits = FALSE)) {
    errorsFW <- NULL
  }
  if (!exists("errorsRV", inherits = FALSE)) {
    errorsRV <- NULL
  }
  
  save(errorsFW, errorsRV, seqtab.nochim, seqtab, mergers, ReadSummary, file = file.path(DataFolder, "DenoisedData.RData"))
  
  message("*********************** Sequence table generated, Data saved ***********************
          ********************************************************************")
  cat("\n*** Sequence table generated, Data saved ***", file = LogFile, append = TRUE)
  cat("\n*** Start generating Summary Plots ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Summary Plots
  ##############################
  
  # the width of the plots probably needs adjustment
  width = 5 + 0.5*(length(SampleNames)/10)
  
  # Plot the number of reads at the different steps for each sample (see also ReadSummary)
  pdf(file = file.path(PlotFolder, "NoReads_AllSamples.pdf"), width = width, height = 6)
  print(NoReads_StepsSimple(ReadSummary = ReadSummary, SampleNames = SampleNames, sort = TRUE))
  dev.off()
  
  # PLot the total number of amplicons against the number of unique amplicons for each sample
  FinalNumbers <- dplyr::select(ReadSummary, Sample, UniqueAmplicons = UniqueAmpliconsWOBimera, NoAmplicons = MergedReadsWOBimera)
  # FinalNumbers <- data.frame(Sample = rownames(seqtab.nochim), UniqueAmplicons = rowSums(seqtab.nochim != 0), NoAmplicons = rowSums(seqtab.nochim))
  
  Tr <- TotalandUniqueAmplicons(FinalNumbers = FinalNumbers, seqtab = seqtab.nochim)
  pdf(file = file.path(PlotFolder, "TotalvsUniqueAmplicons.pdf"), width = width, height = 6)
  print(Tr)
  dev.off()
  
  # Plot a linear regression line to illustrate the possible association between the total number of amplicons and the number of 
  # unique amplicons
  Tr2 <- ggplot(FinalNumbers, aes(x = NoAmplicons, y = UniqueAmplicons)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    xlab("Total Amplicons") +
    ylab("Unique Amplicons") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = "#999999", size = .15),
          panel.grid.major.x = element_line(color = "#999999", size = .15))
  
  
  pdf(file = file.path(PlotFolder, "AssociationTotaltoUniqueAmplicons.pdf"), width = 7, height = 6)
  print(Tr2)
  dev.off()
  
  # Plot a histogram illustrating in how many samples the amplicons are present
  if(length(SampleNames) < 10) {binwidth = 1}
  if(length(SampleNames) > 10 && length(SampleNames) < 100) {binwidth = 2}
  if(length(SampleNames) > 100) {binwidth = 3}
  
  FinalNumbers2 <- data.frame(InNoSamples = colSums(seqtab.nochim != 0))
  
  Trr <- ggplot(data = FinalNumbers2, aes(x = InNoSamples))  + 
    geom_histogram(binwidth = binwidth, col = "black", fill = "#E69F00")+
    geom_rug() +
    xlab("Present in No Samples") + 
    ylab("Count") +
    theme_bw() + 
    ggtitle(paste("Total No of unique amplicons:", dim(seqtab.nochim)[2], "Only in 1 sample:", sum(FinalNumbers2 ==1))) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(color = "#999999", size = .15)) +
    coord_cartesian(ylim = c(0,150))
  
  pdf(file = file.path(PlotFolder, "HistogramAmpliconsinNoSamples.pdf"), width = 7, height = 6)
  print(Trr)
  dev.off()
  cat("\n*** Summary Plots generated, Function done ***", file = LogFile, append = TRUE)
  message("*********************** Dada_wrap done ***********************
          ********************************************************************")
}


  
Dada2_wrap_ITS <- function(
  path = "P:/Fungome_Project/ITS2_paper/",
  F_pattern = "*1.fastq" ,
  R_pattern = "*2.fastq",
  path2 = NULL ,
  trimLeft = c(10,10),
  truncLen = 0,
  maxEE = c(2,5),
  maxN = 0,
  truncQ = 2,
  nreadsLearn = 1.2e+06 ,
  err_F = NULL,
  err_R = NULL,
  minOverlap = 20,
  maxMismatch = 0,
  F_QualityStats = NULL,
  R_QualityStats = NULL,
  filtFs = NULL,
  filtRs = NULL,
  justConcatenate=FALSE,
  pool = FALSE){
  
  
  ##############################
  ### call the required packages
  ##############################
  
  ## dada2:
  # source("https://bioconductor.org/biocLite.R")
  try(library(dada2), biocLite("dada2"))
  
  ## Short Read
  try(library(ShortRead), biocLite("ShortRead"))
  
  ## ggplot2
  try(library(ggplot2), install.packages("ggplot2"))
  
  ## dplyr
  try(library(dplyr), install.packages("dplyr"))
  
  ## dplyr
  try(library(tidyr), install.packages("tidyr"))
  
  message("*********************** All packages loaded ***********************
          ********************************************************************")
  
  ##############################
  ### save the R and package Versions
  ##############################
  # NB: outputs an error and stops function if a Package is not installed
  RVer <- R.Version()
  RVer <- RVer$version.string
  PackageVersions <- data.frame(Package = c("R", "dada2", "ShortRead", "ggplot2", "dplyr", "tidyr"),
                                Version = c(RVer,
                                            as.character(packageVersion("dada2")),
                                            as.character(packageVersion("ShortRead")),
                                            as.character(packageVersion("ggplot2")),
                                            as.character(packageVersion("dplyr")),
                                            as.character(packageVersion("tidyr"))))
  
  message(paste(PackageVersions$Package, ": ", PackageVersions$Version, "; ", sep= ""))
  
  ##############################
  ### Construct character vectors to the FW and RV fastq files
  ##############################
  
  if(is.null(path2)){
    
    path2 = path
  }
  
  
  folders <- list.dirs(path, recursive = FALSE, full.names = FALSE)
  
  if(length(folders) != 0) {
    
    SampleNames <- folders
    
    if(sum(grepl("^Dada", SampleNames)) != 0){
      warning("**There are folders starting with \"Dada\" in your path folder, maybe you have run the function on this path folder before\n.
              The folders starting with Dada will not be considered sample folders!!\n Files within Dada_Data Dada_FilteredFastqs and Dada_Plots will be overwritten**")
    }
    
    # exclude Folders starting with "Dada" from the folders considered as sample fodlers
    if(length(grep("^Dada", SampleNames))!=0) {
      SampleNames <- SampleNames[-grep("^Dada", SampleNames)]
    }
    
    F_fastq <- character(length = length(SampleNames))
    R_fastq <- character(length = length(SampleNames))
    
    for (i in 1:length(SampleNames)) {
      CurrentPath <- file.path(path, SampleNames[i])
      
      if(sum(grepl(F_pattern, list.files(CurrentPath))) == 0) {
        stop(paste("F_pattern fits no file in ", CurrentPath))
      }
      if(sum(grepl(F_pattern, list.files(CurrentPath))) > 1) {
        stop(paste("F_pattern fits several files in ", CurrentPath))
      }
      if(sum(grepl(R_pattern, list.files(CurrentPath))) == 0) {
        stop(paste("R_pattern fits no file in ", CurrentPath))
      }
      if(sum(grepl(R_pattern, list.files(CurrentPath))) > 1) {
        stop(paste("R_pattern fits several files in ", CurrentPath))
      }
      F_fastq[i] <- file.path(CurrentPath, list.files(CurrentPath)[grepl(F_pattern, list.files(CurrentPath))])
      R_fastq[i] <- file.path(CurrentPath, list.files(CurrentPath)[grepl(R_pattern, list.files(CurrentPath))])
      
    }
    
    } else {
      
      stop("No sample folders were found in the given path! Currently the Dada2_wrap function can only handle the situation where the fastq files are in separate folders for each sample.
           These folders have to be in the \"path\" folder. Other situations have to be added.")
    }
  
  ##############################
  ### Start the log file
  ##############################
  
  DataFolder <- file.path(path2, "Dada_Data")
  
  # ATTENTION: FUNCTION DELETES ALL FILES ALREADY PRESRENT in DataFolder
  if(file.exists(DataFolder)){
    file.remove(list.files(DataFolder, full.names = TRUE))
  }
  
  
  dir.create(DataFolder, showWarnings = FALSE)
  
  if(!file.exists(DataFolder)){
    stop("DataFolder could not be created! something wrong with path2, maybe permission denied.")
  }
  
  ptm <- proc.time()
  
  LogFile <- file.path(DataFolder, "DadaWrapper.log")
  cat("Time after package installation: ", file = LogFile)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  cat("Package Versions: ", file = LogFile, append = TRUE)
  cat(paste(PackageVersions$Package, ": ", PackageVersions$Version, "; ", sep= ""), file = LogFile, append = TRUE)
  
  
  # Collect inputs in a data frame for saving it later
  Input <- list(path = path, F_pattern = F_pattern, R_pattern = R_pattern, path2 = path2,
                trimLeft = trimLeft, truncLen = truncLen, 
                maxN = maxN, maxEE = maxEE, truncQ = truncQ,
                nreadsLearn = nreadsLearn, err_F = err_F, err_R = err_R,
                minOverlap = minOverlap, maxMismatch = maxMismatch, F_QualityStats = F_QualityStats,
                R_QualityStats = R_QualityStats, filtFs = filtFs, filtRs = filtRs, justConcatenate=justConcatenate)
  
  
  ##############################
  ### Determine the quality scores and save the stats in Data folder
  ##############################
  
  # Collect quality Score data of the fastQ files and store a data.frame for each sample in the lists FW_QualityStats and RV_QualityStats ####################
  
  if (is.null(F_QualityStats) || is.null(R_QualityStats)) {
    
    F_QualityStats <- list()
    R_QualityStats <- list()
    
    for (i in seq_along(F_fastq)) {
      
      Current_FWfq <- F_fastq[i]
      Current_dfFW <- qa(Current_FWfq, n = 1e06)[["perCycle"]]$quality
      # df is a data frame containing for each cycle (nt) the distribution of Quality scores
      # e.g. Cycle 1 had 7 different quality scores, then 7 rows of cycle one, for each score the count says how many reads had this score
      Current_dfFW <- dplyr::group_by(Current_dfFW, Cycle)
      
      Current_dfQStatsFW <- dplyr::summarise(
        Current_dfFW,
        NoReads = sum(Count),
        Mean_QS = sum(Count*Score)/sum(Count),
        SD_QS = sqrt(sum(Count*((Score-Mean_QS)^2))/(NoReads-1)),
        Median_QS = Score[which(cumsum(Count)/sum(Count) >= .5)][[1]],
        q25_QS = Score[which(cumsum(Count)/sum(Count) >= .25)][[1]],
        q75_QS = Score[which(cumsum(Count)/sum(Count) >= .75)][[1]])
      
      
      F_QualityStats[[i]] <- as.data.frame(Current_dfQStatsFW)
      # as.data.frame to un-dplyr the data.frame
      
      # collect the same stats for the RV FastQ files
      Current_RVfq <- R_fastq[i]
      Current_dfRV <- qa(Current_RVfq, n = 1e06)[["perCycle"]]$quality
      
      Current_dfRV <- dplyr::group_by(Current_dfRV, Cycle)
      
      Current_dfQStatsRV <- dplyr::summarise(
        Current_dfRV,
        NoReads = sum(Count),
        Mean_QS = sum(Count*Score)/sum(Count),
        SD_QS = sqrt(sum(Count*((Score-Mean_QS)^2))/(NoReads-1)),
        Median_QS = Score[which(cumsum(Count)/sum(Count) >= .5)][[1]],
        q25_QS = Score[which(cumsum(Count)/sum(Count) >= .25)][[1]],
        q75_QS = Score[which(cumsum(Count)/sum(Count) >= .75)][[1]])
      
      
      R_QualityStats[[i]] <- as.data.frame(Current_dfQStatsRV)
      
      rm(Current_FWfq, Current_dfFW, Current_dfQStatsFW,
         Current_RVfq, Current_dfRV, Current_dfQStatsRV)
      
    }
    
    # add the sample names as names to the lists:
    names(F_QualityStats) <- SampleNames
    names(R_QualityStats) <- SampleNames
    
    message("*********************** Quality Stats Collected ***********************
            ********************************************************************")
    cat("\n*** Quality Stats Collected ***", file = LogFile, append = TRUE)
    
  } else {
    
    if((names(F_QualityStats) != SampleNames) || (names(R_QualityStats) != SampleNames)) {
      
      stop("The given F_QualityStats or R_QualityStats do not fit to the SampleNames in path!")
    }
    
    message("*********************** Quality Stats given ***********************
            ********************************************************************")
    cat("\n*** Quality Stats given***", file = LogFile, append = TRUE)
    
  }
  
  save(PackageVersions, F_QualityStats, R_QualityStats, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  TimePassed <- proc.time()-ptm
  cat("\nTime after Quality Stats collection: ", file = LogFile, append = TRUE)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
  cat("\n*** Start generating Quality Plots ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Generate and save some quality plots
  ##############################
  
  PlotFolder <- file.path(path2, "Dada_Plots")
  
  # ATTENTION: FUNCTION DELETES ALL FILES ALREADY PRESRENT in PlotFolder
  if(file.exists(PlotFolder)){
    file.remove(list.files(PlotFolder, full.names = TRUE))
  }
  
  dir.create(PlotFolder, showWarnings = FALSE)
  
  
  pdf(file = file.path(PlotFolder, "MedianQScore_F_AllNucleotides.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(F_QualityStats, SampleNames,ylim_low = 0, ylim_high = 42))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_R_AllNucleotides.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(R_QualityStats, SampleNames,ylim_low = 0, ylim_high = 42))
  dev.off()
  
  
  message("*********************** Plots generated start filtering ***********************
          ********************************************************************")
  
  cat("\n*** Quality Plots generated ***", file = LogFile, append = TRUE)
  cat("\n*** Start Filtering ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Filtering
  ##############################
  
  FilteredFolder <- file.path(path2, "Dada_FilteredFastqs")
  
  
  if(is.null(filtFs) || is.null(filtRs)){
    
    ## Not sure if this is wanted, deleting all files that are already in the DataFolder folder
    if(file.exists(FilteredFolder)){
      file.remove(list.files(FilteredFolder, full.names = TRUE))
    }
    
    dir.create(FilteredFolder, showWarnings = FALSE)
    
    filtFs <- file.path(FilteredFolder, paste0(SampleNames, "_F_Filtered.fastq.gz"))
    filtRs <- file.path(FilteredFolder, paste0(SampleNames, "_R_Filtered.fastq.gz"))
    names(filtFs) <- SampleNames
    names(filtRs) <- SampleNames
    
    
    for(i in seq_along(F_fastq)) {
      
      message(paste("**Filtering sample No:", i, "called:", SampleNames[i]), " **")
      cat(paste("\n**Filtering sample No:", i, "called:", SampleNames[i]), file = LogFile, append = TRUE)
      fastqPairedFilter(c(F_fastq[i], R_fastq[i]), c(filtFs[i], filtRs[i]), trimLeft = trimLeft, 
                        truncLen = truncLen, maxN = maxN, maxEE = maxEE, truncQ = truncQ, 
                        compress=TRUE, verbose=TRUE)
    }
    
    # check if files have been created
    if(!all(file.exists(filtFs)) || !all(file.exists(filtRs))) {
      
      cat("\n*** ERROR: Not all filtered files were created, maybe trimming impossible***", file = LogFile, append = TRUE)
      
      remove<-as.vector(which(!file.exists(filtFs)))
      names(filtFs[remove])
      message("Not all filtered files were created. Samples", remove, "were removed")
      filtFs<-filtFs[-c(remove)]
      filtRs<-filtRs[-c(remove)]
      SampleNames<-SampleNames[-c(remove)]
      F_QualityStats<-F_QualityStats[-c(remove)]
      R_QualityStats<-R_QualityStats[-c(remove)]
    }
    
    message("*********************** Filtering Done ***********************
            ********************************************************************")
    
    cat("\n*** Filtering Done ***", file = LogFile, append = TRUE)
    
  } else {
    
    if((names(filtFs) != SampleNames) || (names(filtRs) != SampleNames)) {
      
      cat("\n*** ERROR: The given filtFs or filtRs do not have sample names!***", file = LogFile, append = TRUE)
      stop("The given filtFs or filtRs do not have sample names!")
    }
    
    if(!all(file.exists(filtFs)) || !all(file.exists(filtRs))) {
      cat("\n*** ERROR: Not all files in the given filtFs or filtRs existed***", file = LogFile, append = TRUE)
      message("Not all filtered files were created")
      
    }
    
    message("*********************** Filtered files were given ***********************
            ********************************************************************")
    cat("\n*** Filtered files were given ***", file = LogFile, append = TRUE)
    
    
  }
  
  
  save(PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  cat("\nTime after filtering step: ", file = LogFile, append = TRUE)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  TimePassed <- proc.time()-ptm
  cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
  cat("\n*** Start estimating err_F if not given ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Estimate err_F matrix
  ##############################
  
  if(is.null(err_F)){
    
    # the new learnErrors from dada2 allows setting nreads (usually 1 million), it then reads in as many samples until
    # the number of unique sequences (drp$uniques) reaches nreads.
    # learnErrorsAdj is basically dada2::learnErrors, just adds the NREADS from the function to the output list and also adds
    # dreplicated (drps) and denoised (dds) datasets if all samples were used for error matrix estimation. If not drps and dds
    # remain NULL
    errorsFW <- learnErrorsAdj(fls = filtFs, nreads = nreadsLearn, multithread = TRUE)
    
    message(paste("**", errorsFW$NREADS, "reads have been used for err_F estimation **"))
    
    err_F <- errorsFW$err_out
    
    pdf(file = file.path(PlotFolder, "errorRates_F.pdf"), width = 10, height = 10)
    print(plotErrors(errorsFW, nominalQ=TRUE))
    dev.off()
    
    save(err_F, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
    
    message("*********************** err_F has been estimated ***********************
            ********************************************************************")
    cat("\n*** err_F has been estimated ***", file = LogFile, append = TRUE)
    cat(paste("\nReads used for err_F estimation: ", errorsFW$NREADS), file = LogFile, append = TRUE)
    TimePassed <- proc.time()-ptm
    cat("\nTime after err_F estimation: ", file = LogFile, append = TRUE)
    cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
    cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
    cat("\n*** Start estimating err_R if not given ***", file = LogFile, append = TRUE)
    
  }
  
  # in case both err_F and err_R had been given, they would not have been saved if not:
  save(err_F, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  ##############################
  ### Estimate err_R matrix
  ##############################
  
  if(is.null(err_R)){
    
    
    errorsRV <- learnErrorsAdj(fls = filtRs, nreads = nreadsLearn, multithread = TRUE)
    
    message(paste("**", errorsRV$NREADS, "reads have been used for err_R estimation **"))
    
    err_R <- errorsRV$err_out
    
    pdf(file = file.path(PlotFolder, "errorRates_R.pdf"), width = 10, height = 10)
    print(plotErrors(errorsRV, nominalQ=TRUE))
    dev.off()
    
    save(err_F, err_R, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
    
    message("*********************** err_R has been estimated ***********************
            ********************************************************************")
    cat("\n*** err_R has been estimated ***", file = LogFile, append = TRUE)
    cat(paste("\nReads used for err_R estimation: ", errorsRV$NREADS), file = LogFile, append = TRUE)
    TimePassed <- proc.time()-ptm
    cat("\nTime after err_R estimation: ", file = LogFile, append = TRUE)
    cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
    cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
    cat("\n*** Start denoising data, bimera detection, and merging of reads into amplicons***", file = LogFile, append = TRUE)
    
  }
  
  save(err_F, err_R, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  ##############################
  ### Denoising (dada command for all samples)
  ##############################
  
  message("*********************** Start denoising and bimera detection ***********************
          ********************************************************************")
  
  if (pool) { # very time and memory consuming but might help identifying rare SVs
    
    drp_F <- derepFastq(filtFs, verbose = T)
    dd_F <- dada(drp_F, err = err_F, pool = pool, multithread = TRUE)
    drp_R <- derepFastq(filtRs, verbose = TRUE)
    dd_R <- dada(drp_R, err = err_R, pool = pool, multithread = TRUE)
    
    
    # NB: the following demands that the samples were denoised in the given order, I therefore removed tha randomize option from learnErrorsAdj
    names(dd_F) <- SampleNames
    names(drp_F) <- SampleNames
    names(dd_R) <- SampleNames
    names(drp_R) <- SampleNames
    
    # only for the ReadSummary later
    Uniques_F <- sapply(drp_F, function(x) length(x$uniques))
    Uniques_R <- sapply(drp_R, function(x) length(x$uniques))
    Denoised_F <- sapply(dd_F, function(x) length(x$denoised))
    Denoised_R <- sapply(dd_R, function(x) length(x$denoised))
    NoFilteredReads <- sapply(drp_F, function(x) sum(x$uniques)) # would be the same for drp_R, or sum dd_F$denoised
    
    mergers <- mergePairs(dd_F, drp_F, dd_R, drp_R, minOverlap = minOverlap, maxMismatch = maxMismatch, justConcatenate = justConcatenate) #!justConcatenate=TRUE
    
    rm(drp_F, drp_R, dd_F, dd_R)
    
    
  } else if (exists("errorsFW", inherits = FALSE) && exists("errorsRV", inherits = FALSE) && # do not do dereplication and denoising again if all samples had been used for determining the error matrixes
             !is.null(errorsFW$dds) && !is.null(errorsRV$dds)) {
    
    dd_F <- errorsFW$dds
    drp_F <- errorsFW$drps
    dd_R <- errorsRV$dds
    drp_R <- errorsRV$drps
    
    # just to not save the huge files:
    errorsFW$dds <- NULL
    errorsFW$drps <- NULL
    errorsRV$dds <- NULL
    errorsRV$drps <- NULL
    
    # NB: the following demands that the samples were denoised in the given order, I therefore removed tha randomize option from learnErrorsAdj
    names(dd_F) <- SampleNames
    names(drp_F) <- SampleNames
    names(dd_R) <- SampleNames
    names(drp_R) <- SampleNames
    
    # only for the ReadSummary later
    Uniques_F <- sapply(drp_F, function(x) length(x$uniques))
    Uniques_R <- sapply(drp_R, function(x) length(x$uniques))
    Denoised_F <- sapply(dd_F, function(x) length(x$denoised))
    Denoised_R <- sapply(dd_R, function(x) length(x$denoised))
    NoFilteredReads <- sapply(drp_F, function(x) sum(x$uniques)) # would be the same for drp_R, or sum dd_F$denoised
    
    mergers <- mergePairs(dd_F, drp_F, dd_R, drp_R, minOverlap = minOverlap, maxMismatch = maxMismatch,  justConcatenate = justConcatenate)
    
    rm(drp_F, drp_R, dd_F, dd_R)
    
  } else {
    
    if (exists("errorsFW", inherits = FALSE)) {
      errorsFW$dds <- NULL
      errorsFW$drps <- NULL
    }
    if (exists("errorsRV", inherits = FALSE)) {
      errorsRV$dds <- NULL
      errorsRV$drps <- NULL
    }
    
    mergers <- vector("list", length(SampleNames))
    names(mergers) <- SampleNames
    NoFilteredReads <- vector("numeric", length(SampleNames))
    names(NoFilteredReads) <- SampleNames
    Uniques_F <- vector("numeric", length(SampleNames))
    names(Uniques_F) <- SampleNames
    Uniques_R <- vector("numeric", length(SampleNames))
    names(Uniques_R) <- SampleNames
    Denoised_F <- vector("numeric", length(SampleNames))
    names(Denoised_F) <- SampleNames
    Denoised_R <- vector("numeric", length(SampleNames))
    names(Denoised_R) <- SampleNames
    
    for(sam in SampleNames) {
      cat("Processing:", sam, "\n")
      cat(paste("\nDenoising sample:", sam), file = LogFile, append = TRUE)
      derepF <- derepFastq(filtFs[[sam]])
      NoFilteredReads[sam] <- sum(derepF$uniques)
      ddF <- dada(derepF, err=err_F, multithread=TRUE) 
      Uniques_F[sam] <- length(derepF$uniques)
      Denoised_F[sam] <- length(ddF$denoised)
      
      derepR <- derepFastq(filtRs[[sam]])
      ddR <- dada(derepR, err=err_R, multithread=TRUE)
      Uniques_R[sam] <- length(derepR$uniques)
      Denoised_R[sam] <- length(ddR$denoised)
      merger <- mergePairs(ddF, derepF, ddR, derepR, minOverlap = minOverlap, maxMismatch = maxMismatch, justConcatenate = justConcatenate)
      mergers[[sam]] <- merger
      
      rm(derepF, derepR, ddF, ddR, merger)
    }
    
  }
  
  
  if(!all((sapply(1:length(mergers), function(i) dim(mergers[[i]])[1]))!=0)){
    message("**Not in all samples merged amplicons could be found! maybe minOverlap too strict??**")
    cat("\n*** Not in all samples merged amplicons could be found! maybe minOverlap too strict?? ***", file = LogFile, append = TRUE)
  }
  
  if(all((sapply(1:length(mergers), function(i) dim(mergers[[i]])[1]))==0)){
    cat("\n*** ERROR: In none of the samples merged amplicons could be found! maybe minOverlap too strict?? ***", file = LogFile, append = TRUE)
    stop("**In none of the samples merged amplicons could be found! maybe minOverlap too strict??**")
  }
  
  message("*********************** all samples denoised mergerd amplicons generated ***********************
          ********************************************************************")
  cat("\n*** all samples denoised and mergerd amplicons generated ***", file = LogFile, append = TRUE)
  TimePassed <- proc.time()-ptm
  cat("\nTime after denoising: ", file = LogFile, append = TRUE)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
  cat("\n*** Start generating Sequence table and removing bimeras ***", file = LogFile, append = TRUE)
  
  
  ##############################
  ### Generating sequence table
  ##############################
  
  seqtab <- makeSequenceTable(mergers)
  
  
  ##############################
  ### Bimera removal
  ##############################
  
  # removes denoised sequence FW RV combi that was paired in merger for which the FW or the RV sequence was identified as bimera
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  # 1-(dim(seqtab.nochim)[2]/dim(seqtab)[2]) # The fraction chimeras made up of detected sequence variants
  # 1-(sum(seqtab.nochim)/sum(seqtab)) # The fraction of total reads that were chimeras (over all samples)
  
  message("*********************** Seq table generated and Bimeras removed ***********************
          ********************************************************************")
  cat("\n*** seqtable generated, bimeras removed ***", file = LogFile, append = TRUE)
  cat("\n*** Start saving the data ***", file = LogFile, append = TRUE)
  
  
  # generate also a read summary data frame
  ReadSummary <- data.frame(Sample = SampleNames, NoReads = sapply(F_QualityStats, function(df){df$NoReads[1]}))
  # NB, No total reads would be the same from R_QualityStats, could be used for sanity check
  ReadSummary$FilteredReads <- NoFilteredReads
  ReadSummary$MergedReads <- rowSums(seqtab)
  ReadSummary$MergedReadsWOBimera <- rowSums(seqtab.nochim)
  ReadSummary$UniqueSequences_F <- Uniques_F
  ReadSummary$DenoisedSequences_F <- Denoised_F
  ReadSummary$UniqueSequences_R <- Uniques_R
  ReadSummary$DenoisedSequences_R <- Denoised_R
  ReadSummary$Unique_Amplicons <- rowSums(seqtab != 0)
  ReadSummary$UniqueAmpliconsWOBimera <- rowSums(seqtab.nochim != 0)
  rownames(ReadSummary) <- NULL
  
  # I now also save errorsFW and RV to be able to recapitulate the error plots in later analyses
  if (!exists("errorsFW", inherits = FALSE)) {
    errorsFW <- NULL
  }
  if (!exists("errorsRV", inherits = FALSE)) {
    errorsRV <- NULL
  }
  
  save(errorsFW, errorsRV, seqtab.nochim, seqtab, mergers, ReadSummary, file = file.path(DataFolder, "DenoisedData.RData"))
  
  message("*********************** Sequence table generated, Data saved ***********************
          ********************************************************************")
  cat("\n*** Sequence table generated, Data saved ***", file = LogFile, append = TRUE)
  cat("\n*** Start generating Summary Plots ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Summary Plots
  ##############################
  
  # the width of the plots probably needs adjustment
  width = 5 + 0.5*(length(SampleNames)/10)
  
  # Plot the number of reads at the different steps for each sample (see also ReadSummary)
  pdf(file = file.path(PlotFolder, "NoReads_AllSamples.pdf"), width = width, height = 6)
  print(NoReads_StepsSimple(ReadSummary = ReadSummary, SampleNames = SampleNames, sort = TRUE))
  dev.off()
  
  # PLot the total number of amplicons against the number of unique amplicons for each sample
  FinalNumbers <- dplyr::select(ReadSummary, Sample, UniqueAmplicons = UniqueAmpliconsWOBimera, NoAmplicons = MergedReadsWOBimera)
  # FinalNumbers <- data.frame(Sample = rownames(seqtab.nochim), UniqueAmplicons = rowSums(seqtab.nochim != 0), NoAmplicons = rowSums(seqtab.nochim))
  
  Tr <- TotalandUniqueAmplicons(FinalNumbers = FinalNumbers, seqtab = seqtab.nochim)
  pdf(file = file.path(PlotFolder, "TotalvsUniqueAmplicons.pdf"), width = width, height = 6)
  print(Tr)
  dev.off()
  
  # Plot a linear regression line to illustrate the possible association between the total number of amplicons and the number of 
  # unique amplicons
  Tr2 <- ggplot(FinalNumbers, aes(x = NoAmplicons, y = UniqueAmplicons)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    xlab("Total Amplicons") +
    ylab("Unique Amplicons") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = "#999999", size = .15),
          panel.grid.major.x = element_line(color = "#999999", size = .15))
  
  
  pdf(file = file.path(PlotFolder, "AssociationTotaltoUniqueAmplicons.pdf"), width = 7, height = 6)
  print(Tr2)
  dev.off()
  
  # Plot a histogram illustrating in how many samples the amplicons are present
  if(length(SampleNames) < 10) {binwidth = 1}
  if(length(SampleNames) > 10 && length(SampleNames) < 100) {binwidth = 2}
  if(length(SampleNames) > 100) {binwidth = 3}
  
  FinalNumbers2 <- data.frame(InNoSamples = colSums(seqtab.nochim != 0))
  
  Trr <- ggplot(data = FinalNumbers2, aes(x = InNoSamples))  + 
    geom_histogram(binwidth = binwidth, col = "black", fill = "#E69F00")+
    geom_rug() +
    xlab("Present in No Samples") + 
    ylab("Count") +
    theme_bw() + 
    ggtitle(paste("Total No of unique amplicons:", dim(seqtab.nochim)[2], "Only in 1 sample:", sum(FinalNumbers2 ==1))) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(color = "#999999", size = .15)) +
    coord_cartesian(ylim = c(0,150))
  
  pdf(file = file.path(PlotFolder, "HistogramAmpliconsinNoSamples.pdf"), width = 7, height = 6)
  print(Trr)
  dev.off()
  cat("\n*** Summary Plots generated, Function done ***", file = LogFile, append = TRUE)
  message("*********************** Dada_wrap done ***********************
          ********************************************************************")
}

################
#Forward reds
####

Dada2_wrap_2 <- function(path, F_pattern, R_pattern, path2 = NULL,
                       trimLeft = c(10,10), truncLen = c(220, 160), 
                       maxN = 0, maxEE = 2, truncQ = 2,
                       nreadsLearn = 1e+06,
                       err_F = NULL,
                       err_R = NULL,
                       minOverlap = 20,
                       maxMismatch = 0,
                       F_QualityStats = NULL,
                       R_QualityStats = NULL,
                       filtFs = NULL,
                       filtRs = NULL,
                       pool = FALSE,
                       justConcatenate=FALSE
                       
                       ) {
  
  ##############################
  ### call the required packages
  ##############################
  
  ## dada2:
  # source("https://bioconductor.org/biocLite.R")
  try(library(dada2), biocLite("dada2"))
  
  ## Short Read
  try(library(ShortRead), biocLite("ShortRead"))
  
  ## ggplot2
  try(library(ggplot2), install.packages("ggplot2"))
  
  ## dplyr
  try(library(dplyr), install.packages("dplyr"))
  
  ## dplyr
  try(library(tidyr), install.packages("tidyr"))
  
  
  message("*********************** All packages loaded ***********************
          ********************************************************************")
  
  ##############################
  ### save the R and package Versions
  ##############################
  # NB: outputs an error and stops function if a Package is not installed
  RVer <- R.Version()
  RVer <- RVer$version.string
  PackageVersions <- data.frame(Package = c("R", "dada2", "ShortRead", "ggplot2", "dplyr", "tidyr"),
                                Version = c(RVer,
                                            as.character(packageVersion("dada2")),
                                            as.character(packageVersion("ShortRead")),
                                            as.character(packageVersion("ggplot2")),
                                            as.character(packageVersion("dplyr")),
                                            as.character(packageVersion("tidyr"))))
  
  message(paste(PackageVersions$Package, ": ", PackageVersions$Version, "; ", sep= ""))
  
  
  ##############################
  ### Construct character vectors to the FW and RV fastq files
  ##############################
  
  if(is.null(path2)){
    
    path2 = path
  }
  
  
  folders <- list.dirs(path, recursive = FALSE, full.names = FALSE)
  
  if(length(folders) != 0) {
    
    SampleNames <- folders
    
    if(sum(grepl("^Dada", SampleNames)) != 0){
      warning("**There are folders starting with \"Dada\" in your path folder, maybe you have run the function on this path folder before\n.
              The folders starting with Dada will not be considered sample folders!!\n Files within Dada_Data Dada_FilteredFastqs and Dada_Plots will be overwritten**")
    }
    
    # exclude Folders starting with "Dada" from the folders considered as sample fodlers
    if(length(grep("^Dada", SampleNames))!=0) {
      SampleNames <- SampleNames[-grep("^Dada", SampleNames)]
    }
    
    F_fastq <- character(length = length(SampleNames))
    R_fastq <- character(length = length(SampleNames))
    
    for (i in 1:length(SampleNames)) {
      CurrentPath <- file.path(path, SampleNames[i])
      
      if(sum(grepl(F_pattern, list.files(CurrentPath))) == 0) {
        stop(paste("F_pattern fits no file in ", CurrentPath))
      }
      if(sum(grepl(F_pattern, list.files(CurrentPath))) > 1) {
        stop(paste("F_pattern fits several files in ", CurrentPath))
      }
      if(sum(grepl(R_pattern, list.files(CurrentPath))) == 0) {
        stop(paste("R_pattern fits no file in ", CurrentPath))
      }
      if(sum(grepl(R_pattern, list.files(CurrentPath))) > 1) {
        stop(paste("R_pattern fits several files in ", CurrentPath))
      }
      F_fastq[i] <- file.path(CurrentPath, list.files(CurrentPath)[grepl(F_pattern, list.files(CurrentPath))])
      R_fastq[i] <- file.path(CurrentPath, list.files(CurrentPath)[grepl(R_pattern, list.files(CurrentPath))])
      
    }
    
    } else {
      
      stop("No sample folders were found in the given path! Currently the Dada2_wrap function can only handle the situation where the fastq files are in separate folders for each sample.
           These folders have to be in the \"path\" folder. Other situations have to be added.")
    }
  
  ##############################
  ### Start the log file
  ##############################
  
  DataFolder <- file.path(path2, "Dada_Data")
  
  # ATTENTION: FUNCTION DELETES ALL FILES ALREADY PRESRENT in DataFolder
  if(file.exists(DataFolder)){
    file.remove(list.files(DataFolder, full.names = TRUE))
  }
  
  
  dir.create(DataFolder, showWarnings = FALSE)
  
  if(!file.exists(DataFolder)){
    stop("DataFolder could not be created! something wrong with path2, maybe permission denied.")
  }
  
  ptm <- proc.time()
  
  LogFile <- file.path(DataFolder, "DadaWrapper.log")
  cat("Time after package installation: ", file = LogFile)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  cat("Package Versions: ", file = LogFile, append = TRUE)
  cat(paste(PackageVersions$Package, ": ", PackageVersions$Version, "; ", sep= ""), file = LogFile, append = TRUE)
  
  
  # Collect inputs in a data frame for saving it later
  Input <- list(path = path, F_pattern = F_pattern, R_pattern = R_pattern, path2 = path2,
                trimLeft = trimLeft, truncLen = truncLen, 
                maxN = maxN, maxEE = maxEE, truncQ = truncQ,
                nreadsLearn = nreadsLearn, err_F = err_F, err_R = err_R,
                minOverlap = minOverlap, maxMismatch = maxMismatch, F_QualityStats = F_QualityStats,
                R_QualityStats = R_QualityStats, filtFs = filtFs, filtRs = filtRs, justConcatenate=justConcatenate)
  
  
  ##############################
  ### Determine the quality scores and save the stats in Data folder
  ##############################
  
  # Collect quality Score data of the fastQ files and store a data.frame for each sample in the lists FW_QualityStats and RV_QualityStats ####################
  
  if (is.null(F_QualityStats) || is.null(R_QualityStats)) {
    
    F_QualityStats <- list()
    R_QualityStats <- list()
    
    for (i in seq_along(F_fastq)) {
      
      Current_FWfq <- F_fastq[i]
      Current_dfFW <- qa(Current_FWfq, n = 1e06)[["perCycle"]]$quality
      # df is a data frame containing for each cycle (nt) the distribution of Quality scores
      # e.g. Cycle 1 had 7 different quality scores, then 7 rows of cycle one, for each score the count says how many reads had this score
      Current_dfFW <- dplyr::group_by(Current_dfFW, Cycle)
      
      Current_dfQStatsFW <- dplyr::summarise(
        Current_dfFW,
        NoReads = sum(Count),
        Mean_QS = sum(Count*Score)/sum(Count),
        SD_QS = sqrt(sum(Count*((Score-Mean_QS)^2))/(NoReads-1)),
        Median_QS = Score[which(cumsum(Count)/sum(Count) >= .5)][[1]],
        q25_QS = Score[which(cumsum(Count)/sum(Count) >= .25)][[1]],
        q75_QS = Score[which(cumsum(Count)/sum(Count) >= .75)][[1]])
      
      # Check that all reads are of same length
      # x <- range(Current_dfQStatsFW$NoReads)
      # if (!all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5)) {
      #         stop("Not all reads of same length in file", F_fastq[i])
      # }
      
      F_QualityStats[[i]] <- as.data.frame(Current_dfQStatsFW)
      # as.data.frame to un-dplyr the data.frame
      
      # collect the same stats for the RV FastQ files
      Current_RVfq <- R_fastq[i]
      Current_dfRV <- qa(Current_RVfq, n = 1e06)[["perCycle"]]$quality
      
      Current_dfRV <- dplyr::group_by(Current_dfRV, Cycle)
      
      Current_dfQStatsRV <- dplyr::summarise(
        Current_dfRV,
        NoReads = sum(Count),
        Mean_QS = sum(Count*Score)/sum(Count),
        SD_QS = sqrt(sum(Count*((Score-Mean_QS)^2))/(NoReads-1)),
        Median_QS = Score[which(cumsum(Count)/sum(Count) >= .5)][[1]],
        q25_QS = Score[which(cumsum(Count)/sum(Count) >= .25)][[1]],
        q75_QS = Score[which(cumsum(Count)/sum(Count) >= .75)][[1]])
      
      ## Check that all reads are of same length
      # x <- range(Current_dfQStatsRV$NoReads)
      # if (!all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5)) {
      #         stop("Not all reads of same length in file", R_fastq[i])
      # }
      # 
      R_QualityStats[[i]] <- as.data.frame(Current_dfQStatsRV)
      
      rm(Current_FWfq, Current_dfFW, Current_dfQStatsFW,
         Current_RVfq, Current_dfRV, Current_dfQStatsRV, x)
      
    }
    
    # add the sample names as names to the lists:
    names(F_QualityStats) <- SampleNames
    names(R_QualityStats) <- SampleNames
    
    message("*********************** Quality Stats Collected ***********************
            ********************************************************************")
    cat("\n*** Quality Stats Collected ***", file = LogFile, append = TRUE)
    
  } else {
    
    if((names(F_QualityStats) != SampleNames) || (names(R_QualityStats) != SampleNames)) {
      
      stop("The given F_QualityStats or R_QualityStats do not fit to the SampleNames in path!")
    }
    
    message("*********************** Quality Stats given ***********************
            ********************************************************************")
    cat("\n*** Quality Stats given***", file = LogFile, append = TRUE)
    
  }
  
  save(PackageVersions, F_QualityStats, R_QualityStats, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  TimePassed <- proc.time()-ptm
  cat("\nTime after Quality Stats collection: ", file = LogFile, append = TRUE)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
  cat("\n*** Start generating Quality Plots ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Generate and save some quality plots
  ##############################
  
  PlotFolder <- file.path(path2, "Dada_Plots")
  
  # ATTENTION: FUNCTION DELETES ALL FILES ALREADY PRESRENT in PlotFolder
  if(file.exists(PlotFolder)){
    file.remove(list.files(PlotFolder, full.names = TRUE))
  }
  
  dir.create(PlotFolder, showWarnings = FALSE)
  
  
  pdf(file = file.path(PlotFolder, "MedianQScore_F_AllNucleotides.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(F_QualityStats, SampleNames))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_F_150to240.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(F_QualityStats, SampleNames, xlim_low = 150, xlim_high = 240))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_R_AllNucleotides.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(R_QualityStats, SampleNames))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_R_150to240.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(R_QualityStats, SampleNames, xlim_low = 150, xlim_high = 240))
  dev.off()
  
  message("*********************** Plots generated start filtering ***********************
          ********************************************************************")
  
  cat("\n*** Quality Plots generated ***", file = LogFile, append = TRUE)
  cat("\n*** Start Filtering ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Filtering
  ##############################
  
  FilteredFolder <- file.path(path2, "Dada_FilteredFastqs")
  
  
  if(is.null(filtFs) || is.null(filtRs)){
    
    ## Not sure if this is wanted, deleting all files that are already in the DataFolder folder
    if(file.exists(FilteredFolder)){
      file.remove(list.files(FilteredFolder, full.names = TRUE))
    }
    
    dir.create(FilteredFolder, showWarnings = FALSE)
    
    filtFs <- file.path(FilteredFolder, paste0(SampleNames, "_F_Filtered.fastq.gz"))
    filtRs <- file.path(FilteredFolder, paste0(SampleNames, "_R_Filtered.fastq.gz"))
    names(filtFs) <- SampleNames
    names(filtRs) <- SampleNames
    
    
    for(i in seq_along(F_fastq)) {
      
      message(paste("**Filtering sample No:", i, "called:", SampleNames[i]), " **")
      cat(paste("\n**Filtering sample No:", i, "called:", SampleNames[i]), file = LogFile, append = TRUE)
      fastqPairedFilter(c(F_fastq[i], R_fastq[i]), c(filtFs[i], filtRs[i]), trimLeft = trimLeft, 
                        truncLen = truncLen, maxN = maxN, maxEE = maxEE, truncQ = truncQ, 
                        compress=TRUE, verbose=TRUE)
    }
    
    # check if files have been created
    if(!all(file.exists(filtFs)) || !all(file.exists(filtRs))) {
      cat("\n*** ERROR: Not all filtered files were created, maybe trimming impossible***", file = LogFile, append = TRUE)
      #stop("** Not all filtered files were created, maybe trimming impossible")
      
    }
    
    message("*********************** Filtering Done ***********************
            ********************************************************************")
    cat("\n*** Filtering Done ***", file = LogFile, append = TRUE)
    
  } else {
    
    if((names(filtFs) != SampleNames) || (names(filtRs) != SampleNames)) {
      
      cat("\n*** ERROR: The given filtFs or filtRs do not have sample names!***", file = LogFile, append = TRUE)
      stop("The given filtFs or filtRs do not have sample names!")
    }
    
    if(!all(file.exists(filtFs)) || !all(file.exists(filtRs))) {
      cat("\n*** ERROR: Not all files in the given filtFs or filtRs existed***", file = LogFile, append = TRUE)
      stop("** Not all files in the given filtFs or filtRs existed")
      
    }
    
    message("*********************** Filtered files were given ***********************
            ********************************************************************")
    cat("\n*** Filtered files were given ***", file = LogFile, append = TRUE)
    
    
  }
  
  
  save(PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  cat("\nTime after filtering step: ", file = LogFile, append = TRUE)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  TimePassed <- proc.time()-ptm
  cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
  cat("\n*** Start estimating err_F if not given ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Estimate err_F matrix
  ##############################
  
  if(is.null(err_F)){
    
    # the new learnErrors from dada2 allows setting nreads (usually 1 million), it then reads in as many samples until
    # the number of unique sequences (drp$uniques) reaches nreads.
    # learnErrorsAdj is basically dada2::learnErrors, just adds the NREADS from the function to the output list and also adds
    # dreplicated (drps) and denoised (dds) datasets if all samples were used for error matrix estimation. If not drps and dds
    # remain NULL
    errorsFW <- learnErrorsAdj(fls = filtFs, nreads = nreadsLearn, multithread = TRUE)
    
    message(paste("**", errorsFW$NREADS, "reads have been used for err_F estimation **"))
    
    err_F <- errorsFW$err_out
    
    pdf(file = file.path(PlotFolder, "errorRates_F.pdf"), width = 10, height = 10)
    print(plotErrors(errorsFW, nominalQ=TRUE))
    dev.off()
    
    save(err_F, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
    
    message("*********************** err_F has been estimated ***********************
            ********************************************************************")
    cat("\n*** err_F has been estimated ***", file = LogFile, append = TRUE)
    cat(paste("\nReads used for err_F estimation: ", errorsFW$NREADS), file = LogFile, append = TRUE)
    TimePassed <- proc.time()-ptm
    cat("\nTime after err_F estimation: ", file = LogFile, append = TRUE)
    cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
    cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
    cat("\n*** Start estimating err_R if not given ***", file = LogFile, append = TRUE)
    
  }
  
  # in case both err_F and err_R had been given, they would not have been saved if not:
  save(err_F, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  ##############################
  ### Estimate err_R matrix
  ##############################
  
  if(is.null(err_R)){
    
    
    errorsRV <- learnErrorsAdj(fls = filtRs, nreads = nreadsLearn, multithread = TRUE)
    
    message(paste("**", errorsRV$NREADS, "reads have been used for err_R estimation **"))
    
    err_R <- errorsRV$err_out
    
    pdf(file = file.path(PlotFolder, "errorRates_R.pdf"), width = 10, height = 10)
    print(plotErrors(errorsRV, nominalQ=TRUE))
    dev.off()
    
    save(err_F, err_R, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
    
    message("*********************** err_R has been estimated ***********************
            ********************************************************************")
    cat("\n*** err_R has been estimated ***", file = LogFile, append = TRUE)
    cat(paste("\nReads used for err_R estimation: ", errorsRV$NREADS), file = LogFile, append = TRUE)
    TimePassed <- proc.time()-ptm
    cat("\nTime after err_R estimation: ", file = LogFile, append = TRUE)
    cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
    cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
    cat("\n*** Start denoising data, bimera detection, and merging of reads into amplicons***", file = LogFile, append = TRUE)
    
  }
  
  save(err_F, err_R, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  ##############################
  ### Denoising (dada command for all samples)
  ##############################
  
  message("*********************** Start denoising and bimera detection ***********************
          ********************************************************************")
  
  if (pool) { # very time and memory consuming but might help identifying rare SVs
    
    drp_F <- derepFastq(filtFs, verbose = T)
    dd_F <- dada(drp_F, err = err_F, pool = pool, multithread = TRUE)
    drp_R <- derepFastq(filtRs, verbose = TRUE)
    dd_R <- dada(drp_R, err = err_R, pool = pool, multithread = TRUE)
    
    
    # NB: the following demands that the samples were denoised in the given order, I therefore removed tha randomize option from learnErrorsAdj
    names(dd_F) <- SampleNames
    names(drp_F) <- SampleNames
    names(dd_R) <- SampleNames
    names(drp_R) <- SampleNames
    
    # only for the ReadSummary later
    Uniques_F <- sapply(drp_F, function(x) length(x$uniques))
    Uniques_R <- sapply(drp_R, function(x) length(x$uniques))
    Denoised_F <- sapply(dd_F, function(x) length(x$denoised))
    Denoised_R <- sapply(dd_R, function(x) length(x$denoised))
    NoFilteredReads <- sapply(drp_F, function(x) sum(x$uniques)) # would be the same for drp_R, or sum dd_F$denoised
    
    mergers <- mergePairs(dd_F, drp_F, dd_R, drp_R, minOverlap = minOverlap, maxMismatch = maxMismatch, justConcatenate = justConcatenate) #!justConcatenate=TRUE
    
    rm(drp_F, drp_R, dd_F, dd_R)
    
    
  } else if (exists("errorsFW", inherits = FALSE) && exists("errorsRV", inherits = FALSE) && # do not do dereplication and denoising again if all samples had been used for determining the error matrixes
             !is.null(errorsFW$dds) && !is.null(errorsRV$dds)) {
    
    dd_F <- errorsFW$dds
    drp_F <- errorsFW$drps
    dd_R <- errorsRV$dds
    drp_R <- errorsRV$drps
    
    # just to not save the huge files:
    errorsFW$dds <- NULL
    errorsFW$drps <- NULL
    errorsRV$dds <- NULL
    errorsRV$drps <- NULL
    
    # NB: the following demands that the samples were denoised in the given order, I therefore removed tha randomize option from learnErrorsAdj
    names(dd_F) <- SampleNames
    names(drp_F) <- SampleNames
    names(dd_R) <- SampleNames
    names(drp_R) <- SampleNames
    
    # only for the ReadSummary later
    Uniques_F <- sapply(drp_F, function(x) length(x$uniques))
    Uniques_R <- sapply(drp_R, function(x) length(x$uniques))
    Denoised_F <- sapply(dd_F, function(x) length(x$denoised))
    Denoised_R <- sapply(dd_R, function(x) length(x$denoised))
    NoFilteredReads <- sapply(drp_F, function(x) sum(x$uniques)) # would be the same for drp_R, or sum dd_F$denoised
    
    mergers <- mergePairs(dd_F, drp_F, dd_R, drp_R, minOverlap = minOverlap, maxMismatch = maxMismatch,  justConcatenate = justConcatenate)
    
    rm(drp_F, drp_R, dd_F, dd_R)
    
  } else {
    
    if (exists("errorsFW", inherits = FALSE)) {
      errorsFW$dds <- NULL
      errorsFW$drps <- NULL
    }
    if (exists("errorsRV", inherits = FALSE)) {
      errorsRV$dds <- NULL
      errorsRV$drps <- NULL
    }
    
    mergers <- vector("list", length(SampleNames))
    names(mergers) <- SampleNames
    NoFilteredReads <- vector("numeric", length(SampleNames))
    names(NoFilteredReads) <- SampleNames
    Uniques_F <- vector("numeric", length(SampleNames))
    names(Uniques_F) <- SampleNames
    Uniques_R <- vector("numeric", length(SampleNames))
    names(Uniques_R) <- SampleNames
    Denoised_F <- vector("numeric", length(SampleNames))
    names(Denoised_F) <- SampleNames
    Denoised_R <- vector("numeric", length(SampleNames))
    names(Denoised_R) <- SampleNames
    
    for(sam in SampleNames) {
      cat("Processing:", sam, "\n")
      cat(paste("\nDenoising sample:", sam), file = LogFile, append = TRUE)
      derepF <- derepFastq(filtFs[[sam]])
      NoFilteredReads[sam] <- sum(derepF$uniques)
      ddF <- dada(derepF, err=err_F, multithread=TRUE) 
      Uniques_F[sam] <- length(derepF$uniques)
      Denoised_F[sam] <- length(ddF$denoised)
      
      derepR <- derepFastq(filtRs[[sam]])
      ddR <- dada(derepR, err=err_R, multithread=TRUE)
      Uniques_R[sam] <- length(derepR$uniques)
      Denoised_R[sam] <- length(ddR$denoised)
      merger <- mergePairs(ddF, derepF, ddR, derepR, minOverlap = minOverlap, maxMismatch = maxMismatch, justConcatenate = justConcatenate)
      mergers[[sam]] <- merger
      
      #rm(derepF, derepR, ddF, ddR, merger)
    }
    
  }
  
  
  
  
  if(!all((sapply(1:length(mergers), function(i) dim(mergers[[i]])[1]))!=0)){
    message("**Not in all samples merged amplicons could be found! maybe minOverlap too strict??**")
    cat("\n*** Not in all samples merged amplicons could be found! maybe minOverlap too strict?? ***", file = LogFile, append = TRUE)
  }
  
  if(all((sapply(1:length(mergers), function(i) dim(mergers[[i]])[1]))==0)){
    cat("\n*** ERROR: In none of the samples merged amplicons could be found! maybe minOverlap too strict?? ***", file = LogFile, append = TRUE)
    stop("**In none of the samples merged amplicons could be found! maybe minOverlap too strict??**")
  }
  
  message("*********************** all samples denoised mergerd amplicons generated ***********************
          ********************************************************************")
  cat("\n*** all samples denoised and mergerd amplicons generated ***", file = LogFile, append = TRUE)
  TimePassed <- proc.time()-ptm
  cat("\nTime after denoising: ", file = LogFile, append = TRUE)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
  cat("\n*** Start generating Sequence table and removing bimeras ***", file = LogFile, append = TRUE)
  
  
  ##############################
  ### Generating sequence table
  ##############################
  
  seqtab <- makeSequenceTable(mergers)
  seqtabForward<-makeSequenceTable(ddF)
  
  ##############################
  ### Bimera removal
  ##############################
  
  # removes denoised sequence FW RV combi that was paired in merger for which the FW or the RV sequence was identified as bimera
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  seqtab.nochim_forward<- removeBimeraDenovo(seqtabForward, method="consensus", multithread=TRUE, verbose=TRUE)
  # 1-(dim(seqtab.nochim)[2]/dim(seqtab)[2]) # The fraction chimeras made up of detected sequence variants
  # 1-(sum(seqtab.nochim)/sum(seqtab)) # The fraction of total reads that were chimeras (over all samples)
  
  message("*********************** Seq table generated and Bimeras removed ***********************
          ********************************************************************")
  cat("\n*** seqtable generated, bimeras removed ***", file = LogFile, append = TRUE)
  cat("\n*** Start saving the data ***", file = LogFile, append = TRUE)
  
  
  # generate also a read summary data frame
  ReadSummary <- data.frame(Sample = SampleNames, NoReads = sapply(F_QualityStats, function(df){df$NoReads[1]}))
  # NB, No total reads would be the same from R_QualityStats, could be used for sanity check
  ReadSummary$FilteredReads <- NoFilteredReads
  ReadSummary$MergedReads <- rowSums(seqtab)
  ReadSummary$MergedReadsWOBimera <- rowSums(seqtab.nochim)
  ReadSummary$MergedReadsWOBimera_forward <- rowSums(seqtab.nochim_forward)
    ReadSummary$UniqueSequences_F <- Uniques_F
  ReadSummary$DenoisedSequences_F <- Denoised_F
  ReadSummary$UniqueSequences_R <- Uniques_R
  ReadSummary$DenoisedSequences_R <- Denoised_R
  ReadSummary$Unique_Amplicons <- rowSums(seqtab != 0)
  ReadSummary$Unique_Amplicons_forward <- rowSums(seqtabForward != 0)
  ReadSummary$UniqueAmpliconsWOBimera <- rowSums(seqtab.nochim != 0)
  ReadSummary$UniqueAmpliconsWOBimera_forward <- rowSums(seqtab.nochim_forward != 0)
    ReadSummary$UniqueAmpliconsWOBimera_forward <- rowSums(seqtab.nochim_forward != 0)
  rownames(ReadSummary) <- NULL
  
  # I now also save errorsFW and RV to be able to recapitulate the error plots in later analyses
  if (!exists("errorsFW", inherits = FALSE)) {
    errorsFW <- NULL
  }
  if (!exists("errorsRV", inherits = FALSE)) {
    errorsRV <- NULL
  }
  
  save(errorsFW, errorsRV, seqtab.nochim, seqtab.nochim_forward, seqtab, seqtabForward, mergers, ReadSummary, file = file.path(DataFolder, "DenoisedData.RData"))
  
  message("*********************** Sequence table generated, Data saved ***********************
          ********************************************************************")
  cat("\n*** Sequence table generated, Data saved ***", file = LogFile, append = TRUE)
  cat("\n*** Start generating Summary Plots ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Summary Plots
  ##############################
  
  # the width of the plots probably needs adjustment
  width = 5 + 0.5*(length(SampleNames)/10)
  
  # Plot the number of reads at the different steps for each sample (see also ReadSummary)
  pdf(file = file.path(PlotFolder, "NoReads_AllSamples.pdf"), width = width, height = 6)
  print(NoReads_StepsSimple(ReadSummary = ReadSummary, SampleNames = SampleNames, sort = TRUE))
  dev.off()
  
  # PLot the total number of amplicons against the number of unique amplicons for each sample
  FinalNumbers <- dplyr::select(ReadSummary, Sample, UniqueAmplicons = UniqueAmpliconsWOBimera, NoAmplicons = MergedReadsWOBimera)
  FinalNumbers_forward<-dplyr::select(ReadSummary, Sample, UniqueAmplicons_forward = UniqueAmpliconsWOBimera_forward, NoAmplicons = MergedReadsWOBimera_forward)
  # FinalNumbers <- data.frame(Sample = rownames(seqtab.nochim), UniqueAmplicons = rowSums(seqtab.nochim != 0), NoAmplicons = rowSums(seqtab.nochim))
  
  Tr <- TotalandUniqueAmplicons(FinalNumbers = FinalNumbers, seqtab = seqtab.nochim)
  Tr_forward <- TotalandUniqueAmplicons(FinalNumbers = FinalNumbers_forward, seqtab = seqtab.nochim_forward)
  
  pdf(file = file.path(PlotFolder, "TotalvsUniqueAmplicons.pdf"), width = width, height = 6)
  print(Tr)
  print(Tr_forward)
  dev.off()
  
  # Plot a linear regression line to illustrate the possible association between the total number of amplicons and the number of 
  # unique amplicons
  Tr2 <- ggplot(FinalNumbers, aes(x = NoAmplicons, y = UniqueAmplicons)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    xlab("Total Amplicons") +
    ylab("Unique Amplicons") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = "#999999", size = .15),
          panel.grid.major.x = element_line(color = "#999999", size = .15))
  
  
  pdf(file = file.path(PlotFolder, "AssociationTotaltoUniqueAmplicons.pdf"), width = 7, height = 6)
  print(Tr2)
  dev.off()
  
  # Plot a histogram illustrating in how many samples the amplicons are present
  if(length(SampleNames) < 10) {binwidth = 1}
  if(length(SampleNames) > 10 && length(SampleNames) < 100) {binwidth = 2}
  if(length(SampleNames) > 100) {binwidth = 3}
  
  FinalNumbers2 <- data.frame(InNoSamples = colSums(seqtab.nochim != 0))
  
  Trr <- ggplot(data = FinalNumbers2, aes(x = InNoSamples))  + 
    geom_histogram(binwidth = binwidth, col = "black", fill = "#E69F00")+
    geom_rug() +
    xlab("Present in No Samples") + 
    ylab("Count") +
    theme_bw() + 
    ggtitle(paste("Total No of unique amplicons:", dim(seqtab.nochim)[2], "Only in 1 sample:", sum(FinalNumbers2 ==1))) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(color = "#999999", size = .15)) +
    coord_cartesian(ylim = c(0,150))
  
  pdf(file = file.path(PlotFolder, "HistogramAmpliconsinNoSamples.pdf"), width = 7, height = 6)
  print(Trr)
  dev.off()
  cat("\n*** Summary Plots generated, Function done ***", file = LogFile, append = TRUE)
  message("*********************** Dada_wrap done ***********************
          ********************************************************************")
}


#######################################
### FUNCTION: Dada2_QualityCheck
#######################################

# runs the first part of the Dada2_wrap function, creating the quality plots and data. Based on these one can decide on the filtering
# parameters to subsequently use the Dada2_wrap function
## Input
# path: The path to the folder containing the sample folders with the fastq files: 
# F_pattern: a regular expression to find the fastq files with the forward reads in the sample folders
# R_pattern: a regular expression to find the fastq files with the reverse reads in the sample folders
# path2: default = NULL, the function creates the folders Dada_Plots and Dada_Data, it creates them by default in the path folder, when path2 is given 
# the folders will instead be generated in path2
## Output
# PLOTS:
# Quality Plots saved as pdfs in the generated Dada_Plots folder # NB: the plots are currently assuming 250 nt
# DATA:
# saved in the Dada_Data folder:
# QualityStats.RData contains PackageVersions, F_QualityStats, and R_QualityStats

Dada2_QualityCheck <- function(path, F_pattern, R_pattern, path2 = NULL) {
  
  ##############################
  ### call the required packages
  ##############################
  
  ## dada2:
  # source("https://bioconductor.org/biocLite.R")
  try(library(dada2), biocLite("dada2"))
  
  ## Short Read
  try(library(ShortRead), biocLite("ShortRead"))
  
  ## ggplot2
  try(library(ggplot2), install.packages("ggplot2"))
  
  ## dplyr
  try(library(dplyr), install.packages("dplyr"))
  
  ## dplyr
  try(library(tidyr), install.packages("tidyr"))
  
  
  message("*********************** All packages loaded ***********************
          ********************************************************************")
  
  ##############################
  ### save the package Versions
  ##############################
  # NB: outputs an error and stops function if a Package is not installed
  PackageVersions <- data.frame(Package = c("dada2", "ShortRead", "ggplot2", "dplyr", "tidyr"),
                                Version = c(packageVersion("dada2"),
                                            packageVersion("ShortRead"),
                                            packageVersion("ggplot2"),
                                            packageVersion("dplyr"),
                                            packageVersion("tidyr")))
  
  message(paste(PackageVersions$Package, ": ", PackageVersions$Version, "; ", sep= ""))
  
  ##############################
  ### Construct character vectors to the FW and RV fastq files
  ##############################
  if(is.null(path2)){
    
    path2 = path
  }
  
  
  folders <- list.dirs(path, recursive = FALSE, full.names = FALSE)
  
  if(length(folders) != 0) {
    
    SampleNames <- folders
    
    if(sum(grepl("^Dada", SampleNames)) != 0){
      warning("**There are folders starting with \"Dada\" in your path folder, maybe you have run the function on this path folder before\n.
              The folders starting with Dada will not be considered sample folders!!\n Files within Dada_Data Dada_FilteredFastqs and Dada_Plots will be overwritten**")
    }
    
    # exclude Folders starting with "Dada" from the folders considered as sample fodlers
    if(length(grep("^Dada", SampleNames))!=0) {
      SampleNames <- SampleNames[-grep("^Dada", SampleNames)]
    }
    
    F_fastq <- character(length = length(SampleNames))
    R_fastq <- character(length = length(SampleNames))
    
    for (i in 1:length(SampleNames)) {
      CurrentPath <- file.path(path, SampleNames[i])
      
      if(sum(grepl(F_pattern, list.files(CurrentPath))) == 0) {
        stop(paste("F_pattern fits no file in ", CurrentPath))
      }
      if(sum(grepl(F_pattern, list.files(CurrentPath))) > 1) {
        stop(paste("F_pattern fits several files in ", CurrentPath))
      }
      if(sum(grepl(R_pattern, list.files(CurrentPath))) == 0) {
        stop(paste("R_pattern fits no file in ", CurrentPath))
      }
      if(sum(grepl(R_pattern, list.files(CurrentPath))) > 1) {
        stop(paste("R_pattern fits several files in ", CurrentPath))
      }
      F_fastq[i] <- file.path(CurrentPath, list.files(CurrentPath)[grepl(F_pattern, list.files(CurrentPath))])
      R_fastq[i] <- file.path(CurrentPath, list.files(CurrentPath)[grepl(R_pattern, list.files(CurrentPath))])
      
    } 
    
    } else {
      
      stop("No sample folders were found in the given path! Currently the Dada2_wrap function can only handle the situation where the fastq files are in separate folders for each sample.
           These folders have to be in the \"path\" folder. Other situations have to be added.")
    }
  
  
  ##############################
  ### Determine the quality scores and save the stats in Data folder
  ##############################
  
  DataFolder <- file.path(path2, "Dada_Data")
  
  # Not sure if this is wanted, deleting all files that are already in the DataFolder folder
  if(file.exists(DataFolder)){
    file.remove(list.files(DataFolder, full.names = TRUE))
  }
  
  dir.create(DataFolder, showWarnings = TRUE)
  
  if(!file.exists(DataFolder)){
    stop("DataFolder could not be created! something wrong with path2, maybe permission denied.")
  }
  
  # Collect quality Score data of the fastQ files and store a data.frame for each sample in the lists FW_QualityStats and RV_QualityStats ####################
  
  F_QualityStats <- list()
  R_QualityStats <- list()
  
  for (i in seq_along(F_fastq)) {
    
    Current_FWfq <- F_fastq[i]
    Current_dfFW <- qa(Current_FWfq, n = 1e06)[["perCycle"]]$quality
    # df is a data frame containing for each cycle (nt) the distribution of Quality scores,
    # e.g. Cycle 1 had 7 different quality scores, then 7 rows of cycle one, for each score the count says how many reads had this score
    Current_dfFW <- dplyr::group_by(Current_dfFW, Cycle)
    
    Current_dfQStatsFW <- dplyr::summarise(
      Current_dfFW,
      NoReads = sum(Count),
      Mean_QS = sum(Count*Score)/sum(Count),
      SD_QS = sqrt(sum(Count*((Score-Mean_QS)^2))/(NoReads-1)),
      Median_QS = Score[which(cumsum(Count)/sum(Count) >= .5)][[1]],
      q25_QS = Score[which(cumsum(Count)/sum(Count) >= .25)][[1]],
      q75_QS = Score[which(cumsum(Count)/sum(Count) >= .75)][[1]])
    
    ## Check that all reads are of same length
    # x <- range(Current_dfQStatsFW$NoReads)
    # if (!all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5)) {
    #         stop("Not all reads of same length in file", F_fastq[i])
    # }
    
    F_QualityStats[[i]] <- as.data.frame(Current_dfQStatsFW)
    # as.data.frame to un-dplyr the data.frame
    
    ####### collect the same stats for the RV FastQ files
    Current_RVfq <- R_fastq[i]
    Current_dfRV <- qa(Current_RVfq, n = 1e06)[["perCycle"]]$quality
    
    Current_dfRV <- dplyr::group_by(Current_dfRV, Cycle)
    
    Current_dfQStatsRV <- dplyr::summarise(
      Current_dfRV,
      NoReads = sum(Count),
      Mean_QS = sum(Count*Score)/sum(Count),
      SD_QS = sqrt(sum(Count*((Score-Mean_QS)^2))/(NoReads-1)),
      Median_QS = Score[which(cumsum(Count)/sum(Count) >= .5)][[1]],
      q25_QS = Score[which(cumsum(Count)/sum(Count) >= .25)][[1]],
      q75_QS = Score[which(cumsum(Count)/sum(Count) >= .75)][[1]])
    
    ## Check that all reads are of same length
    # x <- range(Current_dfQStatsRV$NoReads)
    # if (!all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5)) {
    #         stop("Not all reads of same length in file", R_fastq[i])
    # }
    
    R_QualityStats[[i]] <- as.data.frame(Current_dfQStatsRV)
    
    rm(Current_FWfq, Current_dfFW, Current_dfQStatsFW,
       Current_RVfq, Current_dfRV, Current_dfQStatsRV, x)
    
  }
  
  # add the sample names as names to the lists:
  names(F_QualityStats) <- SampleNames
  names(R_QualityStats) <- SampleNames
  
  save(PackageVersions, F_QualityStats, R_QualityStats, file = file.path(DataFolder, "QualityStats.RData"))
  
  message("*********************** Quality Stats Collected ***********************
          ********************************************************************")
  
  ##############################
  ### Generate and save some quality plots
  ##############################
  
  PlotFolder <- file.path(path2, "Dada_Plots")
  
  ## Not sure if this is wanted, deleting all files that are already in the PlotFolder folder
  if(file.exists(PlotFolder)){
    file.remove(list.files(PlotFolder, full.names = TRUE))
  }
  
  dir.create(PlotFolder, showWarnings = FALSE)
  
  
  pdf(file = file.path(PlotFolder, "MedianQScore_F_AllNucleotides.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(F_QualityStats, SampleNames))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_F_150to240.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(F_QualityStats, SampleNames, xlim_low = 150, xlim_high = 240))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_R_AllNucleotides.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(R_QualityStats, SampleNames))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_R_150to240.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(R_QualityStats, SampleNames, xlim_low = 150, xlim_high = 240))
  dev.off()
  
  message("*********************** Quality plots generated ***********************
          ********************************************************************")
}


#######################################
### FUNCTION: learnErrorsAdj
#######################################
# is basically the dada2:::learnErrors function unchanged, only that the NREADS used are added to the output list,
# in addition: if dada command has been run on all samples, dds is also saved to not have to run it again in 
# DadaWrapper.
# NB: this makes this function a bit dangerous if fls were already dereplicated samples, so keep fls to be file names!


learnErrorsAdj <- function (fls, nreads = 1e+06, errorEstimationFunction = loessErrfun, 
                            multithread = FALSE) #, randomize = FALSE 
{
  NREADS <- 0
  drps <- vector("list", length(fls))
  # if (randomize) {
  #         fls <- sample(fls)
  # }
  for (i in seq_along(fls)) {
    if (dada2:::is.list.of(fls, "derep")) {
      drps[[i]] <- fls[[i]]
    }
    else {
      drps[[i]] <- derepFastq(fls[[i]])
    }
    NREADS <- NREADS + sum(drps[[i]]$uniques)
    if (NREADS > nreads) {
      break
    }
  }
  drps <- drps[1:i]
  dds <- dada(drps, err = NULL, selfConsist = TRUE, multithread = multithread)
  # cat("Total reads used: ", NREADS, "\n")
  ErrorList <- getErrors(dds, detailed = TRUE)
  ErrorList$NREADS <- NREADS
  if (length(drps) == length(fls)){
    ErrorList$dds <- dds 
    ErrorList$drps <- drps
  } else {
    ErrorList$dds <- NULL
    ErrorList$drps <- NULL
  }
  return(ErrorList)
}
































##################################### alternatives ########################################

#######################################
### FUNCTION: Dada2_wrap_BimFWRV
#######################################

# basically the same as Dada2_wrap, just that the bimeras are calculated individually for FW and RV reads.
# That was done in the original dada2 paper but has been replaced by removeBimeraDenovo on the final seqtab
# in the dada2 tutorials. 
# Removal of bimeras individually on FW and RV reads as in this function usually removes slightly more amplicons as
# bimeras than does Dada2_wrap

Dada2_wrap_BimFWRV <- function(path, F_pattern, R_pattern, path2 = NULL,
                               trimLeft = c(10,10), truncLen = c(220, 160), 
                               maxN = 0, maxEE = 2, truncQ = 2,
                               nreadsLearn = 1e+06,
                               err_F = NULL,
                               err_R = NULL,
                               minOverlap = 20,
                               maxMismatch = 0,
                               F_QualityStats = NULL,
                               R_QualityStats = NULL,
                               filtFs = NULL,
                               filtRs = NULL) {
  
  ##############################
  ### call the required packages
  ##############################
  
  ## dada2:
  # source("https://bioconductor.org/biocLite.R")
  try(library(dada2), biocLite("dada2"))
  
  ## Short Read
  try(library(ShortRead), biocLite("ShortRead"))
  
  ## ggplot2
  try(library(ggplot2), install.packages("ggplot2"))
  
  ## dplyr
  try(library(dplyr), install.packages("dplyr"))
  
  ## dplyr
  try(library(tidyr), install.packages("tidyr"))
  
  
  message("*********************** All packages loaded ***********************
          ********************************************************************")
  
  ##############################
  ### save the package Versions
  ##############################
  # NB: outputs an error and stops function if a Package is not installed
  PackageVersions <- data.frame(Package = c("dada2", "ShortRead", "ggplot2", "dplyr", "tidyr"),
                                Version = c(packageVersion("dada2"),
                                            packageVersion("ShortRead"),
                                            packageVersion("ggplot2"),
                                            packageVersion("dplyr"),
                                            packageVersion("tidyr")))
  
  message(paste(PackageVersions$Package, ": ", PackageVersions$Version, "; ", sep= ""))
  
  
  ##############################
  ### Construct character vectors to the FW and RV fastq files
  ##############################
  
  if(is.null(path2)){
    
    path2 = path
  }
  
  
  folders <- list.dirs(path, recursive = FALSE, full.names = FALSE)
  
  if(length(folders) != 0) {
    
    SampleNames <- folders
    
    if(sum(grepl("^Dada", SampleNames)) != 0){
      warning("**There are folders starting with \"Dada\" in your path folder, maybe you have run the function on this path folder before\n.
              The folders starting with Dada will not be considered sample folders!!\n Files within Dada_Data Dada_FilteredFastqs and Dada_Plots will be overwritten**")
    }
    
    # exclude Folders starting with "Dada" from the folders considered as sample fodlers
    if(length(grep("^Dada", SampleNames))!=0) {
      SampleNames <- SampleNames[-grep("^Dada", SampleNames)]
    }
    
    F_fastq <- character(length = length(SampleNames))
    R_fastq <- character(length = length(SampleNames))
    
    for (i in 1:length(SampleNames)) {
      CurrentPath <- file.path(path, SampleNames[i])
      
      if(sum(grepl(F_pattern, list.files(CurrentPath))) == 0) {
        stop(paste("F_pattern fits no file in ", CurrentPath))
      }
      if(sum(grepl(F_pattern, list.files(CurrentPath))) > 1) {
        stop(paste("F_pattern fits several files in ", CurrentPath))
      }
      if(sum(grepl(R_pattern, list.files(CurrentPath))) == 0) {
        stop(paste("R_pattern fits no file in ", CurrentPath))
      }
      if(sum(grepl(R_pattern, list.files(CurrentPath))) > 1) {
        stop(paste("R_pattern fits several files in ", CurrentPath))
      }
      F_fastq[i] <- file.path(CurrentPath, list.files(CurrentPath)[grepl(F_pattern, list.files(CurrentPath))])
      R_fastq[i] <- file.path(CurrentPath, list.files(CurrentPath)[grepl(R_pattern, list.files(CurrentPath))])
      
    }
    
    } else {
      
      stop("No sample folders were found in the given path! Currently the Dada2_wrap function can only handle the situation where the fastq files are in separate folders for each sample.
           These folders have to be in the \"path\" folder. Other situations have to be added.")
    }
  
  ##############################
  ### Start the log file
  ##############################
  
  DataFolder <- file.path(path2, "Dada_Data")
  
  # ATTENTION: FUNCTION DELETES ALL FILES ALREADY PRESRENT in DataFolder
  if(file.exists(DataFolder)){
    file.remove(list.files(DataFolder, full.names = TRUE))
  }
  
  
  dir.create(DataFolder, showWarnings = FALSE)
  
  if(!file.exists(DataFolder)){
    stop("DataFolder could not be created! something wrong with path2, maybe permission denied.")
  }
  
  ptm <- proc.time()
  
  LogFile <- file.path(DataFolder, "DadaWrapper.log")
  cat("Time after package installation: ", file = LogFile)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  cat("Package Versions: ", file = LogFile, append = TRUE)
  cat(paste(PackageVersions$Package, ": ", PackageVersions$Version, "; ", sep= ""), file = LogFile, append = TRUE)
  
  
  # Collect inputs in a data frame for saving it later
  Input <- list(path = path, F_pattern = F_pattern, R_pattern = R_pattern, path2 = path2,
                trimLeft = trimLeft, truncLen = truncLen, 
                maxN = maxN, maxEE = maxEE, truncQ = truncQ,
                nreadsLearn = nreadsLearn, err_F = err_F, err_R = err_R,
                minOverlap = minOverlap, maxMismatch = maxMismatch, F_QualityStats = F_QualityStats,
                R_QualityStats = R_QualityStats, filtFs = filtFs, filtRs = filtRs)
  
  
  ##############################
  ### Determine the quality scores and save the stats in Data folder
  ##############################
  
  # Collect quality Score data of the fastQ files and store a data.frame for each sample in the lists FW_QualityStats and RV_QualityStats ####################
  
  if (is.null(F_QualityStats) || is.null(R_QualityStats)) {
    
    F_QualityStats <- list()
    R_QualityStats <- list()
    
    for (i in seq_along(F_fastq)) {
      
      Current_FWfq <- F_fastq[i]
      Current_dfFW <- qa(Current_FWfq, n = 1e06)[["perCycle"]]$quality
      # df is a data frame containing for each cycle (nt) the distribution of Quality scores
      # e.g. Cycle 1 had 7 different quality scores, then 7 rows of cycle one, for each score the count says how many reads had this score
      Current_dfFW <- dplyr::group_by(Current_dfFW, Cycle)
      
      Current_dfQStatsFW <- dplyr::summarise(
        Current_dfFW,
        NoReads = sum(Count),
        Mean_QS = sum(Count*Score)/sum(Count),
        SD_QS = sqrt(sum(Count*((Score-Mean_QS)^2))/(NoReads-1)),
        Median_QS = Score[which(cumsum(Count)/sum(Count) >= .5)][[1]],
        q25_QS = Score[which(cumsum(Count)/sum(Count) >= .25)][[1]],
        q75_QS = Score[which(cumsum(Count)/sum(Count) >= .75)][[1]])
      
      # Check that all reads are of same length
      # x <- range(Current_dfQStatsFW$NoReads)
      # if (!all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5)) {
      #         stop("Not all reads of same length in file", F_fastq[i])
      # }
      # 
      F_QualityStats[[i]] <- as.data.frame(Current_dfQStatsFW)
      # as.data.frame to un-dplyr the data.frame
      
      # collect the same stats for the RV FastQ files
      Current_RVfq <- R_fastq[i]
      Current_dfRV <- qa(Current_RVfq, n = 1e06)[["perCycle"]]$quality
      
      Current_dfRV <- dplyr::group_by(Current_dfRV, Cycle)
      
      Current_dfQStatsRV <- dplyr::summarise(
        Current_dfRV,
        NoReads = sum(Count),
        Mean_QS = sum(Count*Score)/sum(Count),
        SD_QS = sqrt(sum(Count*((Score-Mean_QS)^2))/(NoReads-1)),
        Median_QS = Score[which(cumsum(Count)/sum(Count) >= .5)][[1]],
        q25_QS = Score[which(cumsum(Count)/sum(Count) >= .25)][[1]],
        q75_QS = Score[which(cumsum(Count)/sum(Count) >= .75)][[1]])
      
      ## Check that all reads are of same length
      # x <- range(Current_dfQStatsRV$NoReads)
      # if (!all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5)) {
      #         stop("Not all reads of same length in file", R_fastq[i])
      # }
      # 
      R_QualityStats[[i]] <- as.data.frame(Current_dfQStatsRV)
      
      rm(Current_FWfq, Current_dfFW, Current_dfQStatsFW,
         Current_RVfq, Current_dfRV, Current_dfQStatsRV, x)
      
    }
    
    # add the sample names as names to the lists:
    names(F_QualityStats) <- SampleNames
    names(R_QualityStats) <- SampleNames
    
    message("*********************** Quality Stats Collected ***********************
            ********************************************************************")
    cat("\n*** Quality Stats Collected ***", file = LogFile, append = TRUE)
    
  } else {
    
    if((names(F_QualityStats) != SampleNames) || (names(R_QualityStats) != SampleNames)) {
      
      stop("The given F_QualityStats or R_QualityStats do not fit to the SampleNames in path!")
    }
    
    message("*********************** Quality Stats given ***********************
            ********************************************************************")
    cat("\n*** Quality Stats given***", file = LogFile, append = TRUE)
    
  }
  
  save(PackageVersions, F_QualityStats, R_QualityStats, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  TimePassed <- proc.time()-ptm
  cat("\nTime after Quality Stats collection: ", file = LogFile, append = TRUE)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
  cat("\n*** Start generating Quality Plots ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Generate and save some quality plots
  ##############################
  
  PlotFolder <- file.path(path2, "Dada_Plots")
  
  # ATTENTION: FUNCTION DELETES ALL FILES ALREADY PRESRENT in PlotFolder
  if(file.exists(PlotFolder)){
    file.remove(list.files(PlotFolder, full.names = TRUE))
  }
  
  dir.create(PlotFolder, showWarnings = FALSE)
  
  
  pdf(file = file.path(PlotFolder, "MedianQScore_F_AllNucleotides.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(F_QualityStats, SampleNames))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_F_150to240.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(F_QualityStats, SampleNames, xlim_low = 150, xlim_high = 240))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_R_AllNucleotides.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(R_QualityStats, SampleNames))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_R_150to240.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(R_QualityStats, SampleNames, xlim_low = 150, xlim_high = 240))
  dev.off()
  
  message("*********************** Plots generated start filtering ***********************
          ********************************************************************")
  
  cat("\n*** Quality Plots generated ***", file = LogFile, append = TRUE)
  cat("\n*** Start Filtering ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Filtering
  ##############################
  
  FilteredFolder <- file.path(path2, "Dada_FilteredFastqs")
  
  
  if(is.null(filtFs) || is.null(filtRs)){
    
    ## Not sure if this is wanted, deleting all files that are already in the DataFolder folder
    if(file.exists(FilteredFolder)){
      file.remove(list.files(FilteredFolder, full.names = TRUE))
    }
    
    dir.create(FilteredFolder, showWarnings = FALSE)
    
    filtFs <- file.path(FilteredFolder, paste0(SampleNames, "_F_Filtered.fastq.gz"))
    filtRs <- file.path(FilteredFolder, paste0(SampleNames, "_R_Filtered.fastq.gz"))
    names(filtFs) <- SampleNames
    names(filtRs) <- SampleNames
    
    if(!is.null(truncLen)){
    for(i in seq_along(F_fastq)) {
      
      message(paste("**Filtering sample No:", i, "called:", SampleNames[i]), " **")
      cat(paste("\n**Filtering sample No:", i, "called:", SampleNames[i]), file = LogFile, append = TRUE)
      fastqPairedFilter(c(F_fastq[i], R_fastq[i]), c(filtFs[i], filtRs[i]), trimLeft = trimLeft, 
                        truncLen = truncLen, maxN = maxN, maxEE = maxEE, truncQ = truncQ, 
                        compress=TRUE, verbose=TRUE)
    }}else {
      message(paste("** It is undesirable to truncate reads to a fixed length**"))
      message(paste("**Filtering sample No:", i, "called:", SampleNames[i]), " **")
      cat(paste("\n**Filtering sample No:", i, "called:", SampleNames[i]), file = LogFile, append = TRUE)
      fastqPairedFilter(c(F_fastq[i], R_fastq[i]), c(filtFs[i], filtRs[i]), trimLeft = trimLeft, 
                        maxN = maxN, maxEE = maxEE, truncQ = truncQ, 
                        compress=TRUE, verbose=TRUE)}
    
    # check if files have been created
    if(!all(file.exists(filtFs)) || !all(file.exists(filtRs))) {
      cat("\n*** ERROR: Not all filtered files were created, maybe trimming impossible***", file = LogFile, append = TRUE)
     # stop("** Not all filtered files were created, maybe trimming impossible")
      
    }
    
    message("*********************** Filtering Done ***********************
            ********************************************************************")
    cat("\n*** Filtering Done ***", file = LogFile, append = TRUE)
    
  } else {
    
    if((names(filtFs) != SampleNames) || (names(filtRs) != SampleNames)) {
      
      cat("\n*** ERROR: The given filtFs or filtRs do not have sample names!***", file = LogFile, append = TRUE)
      stop("The given filtFs or filtRs do not have sample names!")
    }
    
    if(!all(file.exists(filtFs)) || !all(file.exists(filtRs))) {
      cat("\n*** ERROR: Not all files in the given filtFs or filtRs existed***", file = LogFile, append = TRUE)
      stop("** Not all files in the given filtFs or filtRs existed")
      
    }
    
    message("*********************** Filtered files were given ***********************
            ********************************************************************")
    cat("\n*** Filtered files were given ***", file = LogFile, append = TRUE)
    
    
  }
  
  
  save(PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  cat("\nTime after filtering step: ", file = LogFile, append = TRUE)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  TimePassed <- proc.time()-ptm
  cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
  cat("\n*** Start estimating err_F if not given ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Estimate err_F matrix
  ##############################
  
  if(is.null(err_F)){
    
    # the new learnErrors from dada2 allows setting nreads (usually 1 million), it then reads in as many samples until
    # the number of unique sequences (drp$uniques) reaches nreads.
    # learnErrorsAdj is basically dada2::learnErrors, just adds the NREADS from the function to the output list and also adds
    # dreplicated (drps) and denoised (dds) datasets if all samples were used for error matrix estimation. If not drps and dds
    # remain NULL
    errorsFW <- learnErrorsAdj(fls = filtFs, nreads = nreadsLearn, multithread = TRUE)
    
    message(paste("**", errorsFW$NREADS, "reads have been used for err_F estimation **"))
    
    err_F <- errorsFW$err_out
    
    pdf(file = file.path(PlotFolder, "errorRates_F.pdf"), width = 10, height = 10)
    print(plotErrors(errorsFW, nominalQ=TRUE))
    dev.off()
    
    save(err_F, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
    
    message("*********************** err_F has been estimated ***********************
            ********************************************************************")
    cat("\n*** err_F has been estimated ***", file = LogFile, append = TRUE)
    cat(paste("\nReads used for err_F estimation: ", errorsFW$NREADS), file = LogFile, append = TRUE)
    TimePassed <- proc.time()-ptm
    cat("\nTime after err_F estimation: ", file = LogFile, append = TRUE)
    cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
    cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
    cat("\n*** Start estimating err_R if not given ***", file = LogFile, append = TRUE)
    
  }
  
  ##############################
  ### Estimate err_R matrix
  ##############################
  
  if(is.null(err_R)){
    
    
    errorsRV <- learnErrorsAdj(fls = filtRs, nreads = nreadsLearn, multithread = TRUE)
    
    message(paste("**", errorsRV$NREADS, "reads have been used for err_R estimation **"))
    
    err_R <- errorsRV$err_out
    
    pdf(file = file.path(PlotFolder, "errorRates_R.pdf"), width = 10, height = 10)
    print(plotErrors(errorsRV, nominalQ=TRUE))
    dev.off()
    
    save(err_F, err_R, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
    
    message("*********************** err_R has been estimated ***********************
            ********************************************************************")
    cat("\n*** err_R has been estimated ***", file = LogFile, append = TRUE)
    cat(paste("\nReads used for err_R estimation: ", errorsRV$NREADS), file = LogFile, append = TRUE)
    TimePassed <- proc.time()-ptm
    cat("\nTime after err_R estimation: ", file = LogFile, append = TRUE)
    cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
    cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
    cat("\n*** Start denoising data, bimera detection, and merging of reads into amplicons***", file = LogFile, append = TRUE)
    
  }
  
  
  ##############################
  ### Denoising (dada command for all samples) and bimeara identification
  ##############################
  
  message("*********************** Start denoising and bimera detection ***********************
          ********************************************************************")
  
  # do not do dereplication and denoising again if all samples had been used for determining the error matrixes
  if (exists("errorsFW", inherits = FALSE) && exists("errorsRV", inherits = FALSE) &&
      !is.null(errorsFW$dds) && !is.null(errorsRV$dds)) {
    
    dd_F <- errorsFW$dds
    drp_F <- errorsFW$drps
    dd_R <- errorsRV$dds
    drp_R <- errorsRV$drps
    
    rm(errorsRV, errorsFW)
    
    # NB: the following demands that the samples were denoised in the given order, I therefore removed tha randomize option from learnErrorsAdj
    names(dd_F) <- SampleNames
    names(drp_F) <- SampleNames
    names(dd_R) <- SampleNames
    names(drp_R) <- SampleNames
    
    # only for the ReadSummary later
    Uniques_F <- sapply(drp_F, function(x) length(x$uniques))
    Uniques_R <- sapply(drp_R, function(x) length(x$uniques))
    Denoised_F <- sapply(dd_F, function(x) length(x$denoised))
    Denoised_R <- sapply(dd_R, function(x) length(x$denoised))
    NoFilteredReads <- sapply(drp_F, function(x) sum(x$uniques)) # would be the same for drp_R, or sum dd_F$denoised
    
    mergers <- mergePairs(dd_F, drp_F, dd_R, drp_R, minOverlap = minOverlap, maxMismatch = maxMismatch)
    
    bimFs <- vector("list", length(SampleNames))
    names(bimFs) <- SampleNames
    bimRs <- vector("list", length(SampleNames))
    names(bimRs) <- SampleNames
    
    for(sam in SampleNames) {
      cat("Processing:", sam, "\n")
      #cat(paste("\nDenoising sample:", sam), file = LogFile, append = TRUE)
      ddF <- dd_F[[sam]] 
      bimFs[[sam]] <- isBimeraDenovo(ddF, verbose=TRUE)
      ddR <- dd_R[[sam]]
      bimRs[[sam]] <- isBimeraDenovo(ddR, verbose=TRUE)
      
      rm(ddF, ddR)
    }
    
    rm(drp_F, drp_R, dd_F, dd_R)
    
  } else {
    
    if (exists("errorsFW", inherits = FALSE)) {
      rm(errorsFW)
    }
    if (exists("errorsRV", inherits = FALSE)) {
      rm(errorsRV)
    }
    
    mergers <- vector("list", length(SampleNames))
    names(mergers) <- SampleNames
    NoFilteredReads <- vector("numeric", length(SampleNames))
    names(NoFilteredReads) <- SampleNames
    bimFs <- vector("list", length(SampleNames))
    names(bimFs) <- SampleNames
    bimRs <- vector("list", length(SampleNames))
    names(bimRs) <- SampleNames
    Uniques_F <- vector("numeric", length(SampleNames))
    names(Uniques_F) <- SampleNames
    Uniques_R <- vector("numeric", length(SampleNames))
    names(Uniques_R) <- SampleNames
    Denoised_F <- vector("numeric", length(SampleNames))
    names(Denoised_F) <- SampleNames
    Denoised_R <- vector("numeric", length(SampleNames))
    names(Denoised_R) <- SampleNames
    
    for(sam in SampleNames) {
      cat("Processing:", sam, "\n")
      cat(paste("\nDenoising sample:", sam), file = LogFile, append = TRUE)
      derepF <- derepFastq(filtFs[[sam]])
      NoFilteredReads[sam] <- sum(derepF$uniques)
      ddF <- dada(derepF, err=err_F, multithread=TRUE) 
      bimFs[[sam]] <- isBimeraDenovo(ddF, verbose=TRUE)
      Uniques_F[sam] <- length(derepF$uniques)
      Denoised_F[sam] <- length(ddF$denoised)
      
      derepR <- derepFastq(filtRs[[sam]])
      ddR <- dada(derepR, err=err_R, multithread=TRUE)
      bimRs[[sam]] <- isBimeraDenovo(ddR, verbose=TRUE)
      Uniques_R[sam] <- length(derepR$uniques)
      Denoised_R[sam] <- length(ddR$denoised)
      merger <- mergePairs(ddF, derepF, ddR, derepR, minOverlap = minOverlap, maxMismatch = maxMismatch)
      mergers[[sam]] <- merger
      
      rm(derepF, derepR, ddF, ddR, merger)
    }
    
  }
  
  if(!all((sapply(1:length(mergers), function(i) dim(mergers[[i]])[1]))!=0)){
    message("**Not in all samples merged amplicons could be found! maybe minOverlap too strict??**")
    cat("\n*** Not in all samples merged amplicons could be found! maybe minOverlap too strict?? ***", file = LogFile, append = TRUE)
  }
  
  if(all((sapply(1:length(mergers), function(i) dim(mergers[[i]])[1]))==0)){
    cat("\n*** ERROR: In none of the samples merged amplicons could be found! maybe minOverlap too strict?? ***", file = LogFile, append = TRUE)
    stop("**In none of the samples merged amplicons could be found! maybe minOverlap too strict??**")
  }
  
  message("*********************** all samples denoised, bimeras identified, mergerd amplicons generated ***********************
          ********************************************************************")
  cat("\n*** all samples denoised, bimeras identified, mergerd amplicons generated ***", file = LogFile, append = TRUE)
  TimePassed <- proc.time()-ptm
  cat("\nTime after denoising: ", file = LogFile, append = TRUE)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
  cat("\n*** Start removing bimera ***", file = LogFile, append = TRUE)
  
  
  ##############################
  ### Bimera removal
  ##############################
  
  # removes denoised sequence FW RV combi that was paired in merger for which the FW or the RV sequence was identified as bimera
  mergers.nochim <- mergers
  for (i in seq_along(mergers)) {
    mergers.nochim[[i]] <- mergers[[i]][!bimFs[[i]][mergers[[i]]$forward] & !bimRs[[i]][mergers[[i]]$reverse],]
  }
  
  message("*********************** Bimeras removed ***********************
          ********************************************************************")
  cat("\n*** bimeras removed ***", file = LogFile, append = TRUE)
  cat("\n*** Start generating Sequence table ***", file = LogFile, append = TRUE)
  
  
  ##############################
  ### Generating sequence table
  ##############################
  seqtab <- makeSequenceTable(mergers)
  seqtab.nochim <- makeSequenceTable(mergers.nochim)
  # 1-(dim(seqtab.nochimFWRW)[2]/dim(seqtab)[2]) # The fraction chimeras made up of detected sequence variants
  # 1-(sum(seqtab.nochimFWRW)/sum(seqtab)) # The fraction of total reads that were chimeras (over all samples)
  
  # generate also a read summary data frame
  ReadSummary <- data.frame(Sample = SampleNames, NoReads = sapply(F_QualityStats, function(df){df$NoReads[1]}))
  # NB, No total reads would be the same from R_QualityStats, could be used for sanity check
  ReadSummary$FilteredReads <- NoFilteredReads
  ReadSummary$MergedReads <- rowSums(seqtab)
  ReadSummary$MergedReadsWOBimera <- rowSums(seqtab.nochim)
  ReadSummary$UniqueSequences_F <- Uniques_F
  ReadSummary$DenoisedSequences_F <- Denoised_F
  ReadSummary$UniqueSequences_R <- Uniques_R
  ReadSummary$DenoisedSequences_R <- Denoised_R
  ReadSummary$Unique_Amplicons <- rowSums(seqtab != 0)
  ReadSummary$UniqueAmpliconsWOBimera <- rowSums(seqtab.nochim != 0)
  rownames(ReadSummary) <- NULL
  
  save(seqtab, seqtab.nochim, mergers, mergers.nochim, bimFs, bimRs, ReadSummary,
       file = file.path(DataFolder, "DenoisedData.RData"))
  
  message("*********************** Sequence table generated, Data saved ***********************
          ********************************************************************")
  cat("\n*** Sequence table generated, Data saved ***", file = LogFile, append = TRUE)
  cat("\n*** Start generating Summary Plots ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Summary Plots
  ##############################
  
  # the width of the plots probably needs adjustment
  width = 5 + 0.5*(length(SampleNames)/10)
  
  # Plot the number of reads at the different steps for each sample (see also ReadSummary)
  pdf(file = file.path(PlotFolder, "NoReads_AllSamples.pdf"), width = width, height = 6)
  print(NoReads_StepsSimple(ReadSummary = ReadSummary, SampleNames = SampleNames, sort = TRUE))
  dev.off()
  
  # PLot the total number of amplicons against the number of unique amplicons for each sample
  FinalNumbers <- dplyr::select(ReadSummary, Sample, UniqueAmplicons = UniqueAmpliconsWOBimera, NoAmplicons = MergedReadsWOBimera)
  # FinalNumbers <- data.frame(Sample = rownames(seqtab.nochim), UniqueAmplicons = rowSums(seqtab.nochim != 0), NoAmplicons = rowSums(seqtab.nochim))
  
  Tr <- TotalandUniqueAmplicons(FinalNumbers = FinalNumbers, seqtab = seqtab.nochim)
  pdf(file = file.path(PlotFolder, "TotalvsUniqueAmplicons.pdf"), width = width, height = 6)
  print(Tr)
  dev.off()
  
  # Plot a linear regression line to illustrate the possible association between the total number of amplicons and the number of 
  # unique amplicons
  Tr2 <- ggplot(FinalNumbers, aes(x = NoAmplicons, y = UniqueAmplicons)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    xlab("Total Amplicons") +
    ylab("Unique Amplicons") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = "#999999", size = .15),
          panel.grid.major.x = element_line(color = "#999999", size = .15))
  
  
  pdf(file = file.path(PlotFolder, "AssociationTotaltoUniqueAmplicons.pdf"), width = 7, height = 6)
  print(Tr2)
  dev.off()
  
  # Plot a histogram illustrating in how many samples the amplicons are present
  if(length(SampleNames) < 10) {binwidth = 1}
  if(length(SampleNames) > 10 && length(SampleNames) < 100) {binwidth = 2}
  if(length(SampleNames) > 100) {binwidth = 3}
  
  FinalNumbers2 <- data.frame(InNoSamples = colSums(seqtab.nochim != 0))
  
  Trr <- ggplot(data = FinalNumbers2, aes(x = InNoSamples))  + 
    geom_histogram(binwidth = binwidth, col = "black", fill = "#E69F00")+
    geom_rug() +
    xlab("Present in No Samples") + 
    ylab("Count") +
    theme_bw() + 
    ggtitle(paste("Total No of unique amplicons:", dim(seqtab.nochim)[2], "Only in 1 sample:", sum(FinalNumbers2 ==1))) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(color = "#999999", size = .15)) +
    coord_cartesian(ylim = c(0,150))
  
  pdf(file = file.path(PlotFolder, "HistogramAmpliconsinNoSamples.pdf"), width = 7, height = 6)
  print(Trr)
  dev.off()
  cat("\n*** Summary Plots generated, Function done ***", file = LogFile, append = TRUE)
}
#######################################
### FUNCTION: Dada2_wrap
#######################################

# Function that runs the dada2 pipeline up to the sequence table (seqtab)
## Input
# path: The path to the folder containing the sample folders with the fastq files. 
# F_pattern: a regular expression to find the fastq files with the forward reads in the sample folders
# R_pattern: a regular expression to find the fastq files with the reverse reads in the sample folders
# path2: default = NULL, the function creates three folders Dada_Plots, Dada_Data and Dada_FilteredFastqs, it creates them by default in the path folder, when path2 is given 
# the folders will instead be generated in path2
# trimLeft: = trimLeft from the fastqPairedFilter command (how many nucleotids will be clipped from the beginning of the reads)
# truncLen: = truncLen from the fastqPairedFilter command (where the reads will be truncated)
# maxN: = maxN from the fastqPairedFilter command, After truncation, sequences with more than maxN Ns will be discarded, NB: Dada2 requires no Ns!!
# maxEE: = maxEE from the fastqPairedFilter command, the maximum expected error allowed to let a read pass the filter
# trunQ: = truncQ from the fastqPairedFilter command, Truncate reads at the first instance of a quality score less than or equal to truncQ, ? actually not sure if these reads are then discarded?
#    # NB: currently run before trimLeft/truncLen which does the size check. So unless truncated after the truncLen
# the size check will kick the affected reads out, so makes little sense, see: https://github.com/benjjneb/dada2/issues/140 
# nreadsLearn = 1e+06, # the number of reads (distributed over far less unique reads) used to learn the error matrixes, i.e. nreads in dada2:::learnErrors 
# Should account for a number of samples that in total have about 1 milliion filtered reads (number will be given in the log file)
# err_F: the error matrix for the dada command for the forward reads. Default NULL then the error matrix will be estimated (using SelfConsit = TRUE)
# err_R: analogous to err_F for the reverse reads
# minOverlap: = minOverlap from the mergePairs command
# maxMismatch: = maxMismatch from the mergePairs command
# F_QualityStats and R_QualityStats: NULL by default, if given e.g. from Dada_QualityCheck the quality stats collection part is jumped over
# filtFs and filtRs: NULL by default, if given and FilteredFolder with files exist, Filtering can be jumped over.
# pool: pool option of dada command
## Output
# PLOTS:
# Several Plots saved as pdfs in generated Dada_Plots folder # NB: the plots are currently assuming 250 nt
# DATA AND LOG FILE:
# Data and the file DadaWrapper.log are saved in the generated folder: Dada_Data
# DenoisedData.RData contains: seqtab (final sequence table), mergers, mergers.nochim, bimFs, bimRs, 
# ReadSummary (data frame summarising the reads and amplicons at the different stages),
# QualityStats.RData contains PackageVersions, F_QualityStats, R_QualityStats
# FILTERED FASTQ files
# are stored in the generated folder Dada_FilteredFastqs
# NB the filtered fastq files are currently named "SampleName"_F_Filtered.fastq.gz and "SampleName"_R_Filtered.fastq.gz

Dada2_wrap <- function(path, F_pattern, R_pattern, path2 = NULL,
                       trimLeft = c(10,10), truncLen = c(220, 160), 
                       maxN = 0, maxEE = 2, truncQ = 2,
                       nreadsLearn = 1e+06,
                       err_F = NULL,
                       err_R = NULL,
                       minOverlap = 20,
                       maxMismatch = 0,
                       F_QualityStats = NULL,
                       R_QualityStats = NULL,
                       filtFs = NULL,
                       filtRs = NULL,
                       pool = FALSE) {
  
  ##############################
  ### call the required packages
  ##############################
  
  ## dada2:
  # source("https://bioconductor.org/biocLite.R")
  try(library(dada2), biocLite("dada2"))
  
  ## Short Read
  try(library(ShortRead), biocLite("ShortRead"))
  
  ## ggplot2
  try(library(ggplot2), install.packages("ggplot2"))
  
  ## dplyr
  try(library(dplyr), install.packages("dplyr"))
  
  ## dplyr
  try(library(tidyr), install.packages("tidyr"))
  
  
  message("*********************** All packages loaded ***********************
          ********************************************************************")
  
  ##############################
  ### save the R and package Versions
  ##############################
  # NB: outputs an error and stops function if a Package is not installed
  RVer <- R.Version()
  RVer <- RVer$version.string
  PackageVersions <- data.frame(Package = c("R", "dada2", "ShortRead", "ggplot2", "dplyr", "tidyr"),
                                Version = c(RVer,
                                            as.character(packageVersion("dada2")),
                                            as.character(packageVersion("ShortRead")),
                                            as.character(packageVersion("ggplot2")),
                                            as.character(packageVersion("dplyr")),
                                            as.character(packageVersion("tidyr"))))
  
  message(paste(PackageVersions$Package, ": ", PackageVersions$Version, "; ", sep= ""))
  
  
  ##############################
  ### Construct character vectors to the FW and RV fastq files
  ##############################
  
  if(is.null(path2)){
    
    path2 = path
  }
  
  
  folders <- list.dirs(path, recursive = FALSE, full.names = FALSE)
  
  if(length(folders) != 0) {
    
    SampleNames <- folders
    
    if(sum(grepl("^Dada", SampleNames)) != 0){
      warning("**There are folders starting with \"Dada\" in your path folder, maybe you have run the function on this path folder before\n.
              The folders starting with Dada will not be considered sample folders!!\n Files within Dada_Data Dada_FilteredFastqs and Dada_Plots will be overwritten**")
    }
    
    # exclude Folders starting with "Dada" from the folders considered as sample fodlers
    if(length(grep("^Dada", SampleNames))!=0) {
      SampleNames <- SampleNames[-grep("^Dada", SampleNames)]
    }
    
    F_fastq <- character(length = length(SampleNames))
    R_fastq <- character(length = length(SampleNames))
    
    for (i in 1:length(SampleNames)) {
      CurrentPath <- file.path(path, SampleNames[i])
      
      if(sum(grepl(F_pattern, list.files(CurrentPath))) == 0) {
        stop(paste("F_pattern fits no file in ", CurrentPath))
      }
      if(sum(grepl(F_pattern, list.files(CurrentPath))) > 1) {
        stop(paste("F_pattern fits several files in ", CurrentPath))
      }
      if(sum(grepl(R_pattern, list.files(CurrentPath))) == 0) {
        stop(paste("R_pattern fits no file in ", CurrentPath))
      }
      if(sum(grepl(R_pattern, list.files(CurrentPath))) > 1) {
        stop(paste("R_pattern fits several files in ", CurrentPath))
      }
      F_fastq[i] <- file.path(CurrentPath, list.files(CurrentPath)[grepl(F_pattern, list.files(CurrentPath))])
      R_fastq[i] <- file.path(CurrentPath, list.files(CurrentPath)[grepl(R_pattern, list.files(CurrentPath))])
      
    }
    
    } else {
      
      stop("No sample folders were found in the given path! Currently the Dada2_wrap function can only handle the situation where the fastq files are in separate folders for each sample.
           These folders have to be in the \"path\" folder. Other situations have to be added.")
    }
  
  ##############################
  ### Start the log file
  ##############################
  
  DataFolder <- file.path(path2, "Dada_Data")
  
  # ATTENTION: FUNCTION DELETES ALL FILES ALREADY PRESRENT in DataFolder
  if(file.exists(DataFolder)){
    file.remove(list.files(DataFolder, full.names = TRUE))
  }
  
  
  dir.create(DataFolder, showWarnings = FALSE)
  
  if(!file.exists(DataFolder)){
    stop("DataFolder could not be created! something wrong with path2, maybe permission denied.")
  }
  
  ptm <- proc.time()
  
  LogFile <- file.path(DataFolder, "DadaWrapper.log")
  cat("Time after package installation: ", file = LogFile)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  cat("Package Versions: ", file = LogFile, append = TRUE)
  cat(paste(PackageVersions$Package, ": ", PackageVersions$Version, "; ", sep= ""), file = LogFile, append = TRUE)
  
  
  # Collect inputs in a data frame for saving it later
  Input <- list(path = path, F_pattern = F_pattern, R_pattern = R_pattern, path2 = path2,
                trimLeft = trimLeft, truncLen = truncLen, 
                maxN = maxN, maxEE = maxEE, truncQ = truncQ,
                nreadsLearn = nreadsLearn, err_F = err_F, err_R = err_R,
                minOverlap = minOverlap, maxMismatch = maxMismatch, F_QualityStats = F_QualityStats,
                R_QualityStats = R_QualityStats, filtFs = filtFs, filtRs = filtRs)
  
  
  ##############################
  ### Determine the quality scores and save the stats in Data folder
  ##############################
  
  # Collect quality Score data of the fastQ files and store a data.frame for each sample in the lists FW_QualityStats and RV_QualityStats ####################
  
  if (is.null(F_QualityStats) || is.null(R_QualityStats)) {
    
    F_QualityStats <- list()
    R_QualityStats <- list()
    
    for (i in seq_along(F_fastq)) {
      
      Current_FWfq <- F_fastq[i]
      Current_dfFW <- qa(Current_FWfq, n = 1e06)[["perCycle"]]$quality
      # df is a data frame containing for each cycle (nt) the distribution of Quality scores
      # e.g. Cycle 1 had 7 different quality scores, then 7 rows of cycle one, for each score the count says how many reads had this score
      Current_dfFW <- dplyr::group_by(Current_dfFW, Cycle)
      
      Current_dfQStatsFW <- dplyr::summarise(
        Current_dfFW,
        NoReads = sum(Count),
        Mean_QS = sum(Count*Score)/sum(Count),
        SD_QS = sqrt(sum(Count*((Score-Mean_QS)^2))/(NoReads-1)),
        Median_QS = Score[which(cumsum(Count)/sum(Count) >= .5)][[1]],
        q25_QS = Score[which(cumsum(Count)/sum(Count) >= .25)][[1]],
        q75_QS = Score[which(cumsum(Count)/sum(Count) >= .75)][[1]])
      
      # Check that all reads are of same length
      # x <- range(Current_dfQStatsFW$NoReads)
      # if (!all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5)) {
      #         stop("Not all reads of same length in file", F_fastq[i])
      # }
      
      F_QualityStats[[i]] <- as.data.frame(Current_dfQStatsFW)
      # as.data.frame to un-dplyr the data.frame
      
      # collect the same stats for the RV FastQ files
      Current_RVfq <- R_fastq[i]
      Current_dfRV <- qa(Current_RVfq, n = 1e06)[["perCycle"]]$quality
      
      Current_dfRV <- dplyr::group_by(Current_dfRV, Cycle)
      
      Current_dfQStatsRV <- dplyr::summarise(
        Current_dfRV,
        NoReads = sum(Count),
        Mean_QS = sum(Count*Score)/sum(Count),
        SD_QS = sqrt(sum(Count*((Score-Mean_QS)^2))/(NoReads-1)),
        Median_QS = Score[which(cumsum(Count)/sum(Count) >= .5)][[1]],
        q25_QS = Score[which(cumsum(Count)/sum(Count) >= .25)][[1]],
        q75_QS = Score[which(cumsum(Count)/sum(Count) >= .75)][[1]])
      
      ## Check that all reads are of same length
      # x <- range(Current_dfQStatsRV$NoReads)
      # if (!all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5)) {
      #         stop("Not all reads of same length in file", R_fastq[i])
      # }
      # 
      R_QualityStats[[i]] <- as.data.frame(Current_dfQStatsRV)
      
      rm(Current_FWfq, Current_dfFW, Current_dfQStatsFW,
         Current_RVfq, Current_dfRV, Current_dfQStatsRV, x)
      
    }
    
    # add the sample names as names to the lists:
    names(F_QualityStats) <- SampleNames
    names(R_QualityStats) <- SampleNames
    
    message("*********************** Quality Stats Collected ***********************
            ********************************************************************")
    cat("\n*** Quality Stats Collected ***", file = LogFile, append = TRUE)
    
  } else {
    
    if((names(F_QualityStats) != SampleNames) || (names(R_QualityStats) != SampleNames)) {
      
      stop("The given F_QualityStats or R_QualityStats do not fit to the SampleNames in path!")
    }
    
    message("*********************** Quality Stats given ***********************
            ********************************************************************")
    cat("\n*** Quality Stats given***", file = LogFile, append = TRUE)
    
  }
  
  save(PackageVersions, F_QualityStats, R_QualityStats, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  TimePassed <- proc.time()-ptm
  cat("\nTime after Quality Stats collection: ", file = LogFile, append = TRUE)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
  cat("\n*** Start generating Quality Plots ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Generate and save some quality plots
  ##############################
  
  PlotFolder <- file.path(path2, "Dada_Plots")
  
  # ATTENTION: FUNCTION DELETES ALL FILES ALREADY PRESRENT in PlotFolder
  if(file.exists(PlotFolder)){
    file.remove(list.files(PlotFolder, full.names = TRUE))
  }
  
  dir.create(PlotFolder, showWarnings = FALSE)
  
  
  pdf(file = file.path(PlotFolder, "MedianQScore_F_AllNucleotides.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(F_QualityStats, SampleNames))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_F_150to240.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(F_QualityStats, SampleNames, xlim_low = 150, xlim_high = 240))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_R_AllNucleotides.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(R_QualityStats, SampleNames))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_R_150to240.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(R_QualityStats, SampleNames, xlim_low = 150, xlim_high = 240))
  dev.off()
  
  message("*********************** Plots generated start filtering ***********************
          ********************************************************************")
  
  cat("\n*** Quality Plots generated ***", file = LogFile, append = TRUE)
  cat("\n*** Start Filtering ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Filtering
  ##############################
  
  FilteredFolder <- file.path(path2, "Dada_FilteredFastqs")
  
  
  if(is.null(filtFs) || is.null(filtRs)){
    
    ## Not sure if this is wanted, deleting all files that are already in the DataFolder folder
    if(file.exists(FilteredFolder)){
      file.remove(list.files(FilteredFolder, full.names = TRUE))
    }
    
    dir.create(FilteredFolder, showWarnings = FALSE)
    
    filtFs <- file.path(FilteredFolder, paste0(SampleNames, "_F_Filtered.fastq.gz"))
    filtRs <- file.path(FilteredFolder, paste0(SampleNames, "_R_Filtered.fastq.gz"))
    names(filtFs) <- SampleNames
    names(filtRs) <- SampleNames
    
    
    for(i in seq_along(F_fastq)) {
      
      message(paste("**Filtering sample No:", i, "called:", SampleNames[i]), " **")
      cat(paste("\n**Filtering sample No:", i, "called:", SampleNames[i]), file = LogFile, append = TRUE)
      fastqPairedFilter(c(F_fastq[i], R_fastq[i]), c(filtFs[i], filtRs[i]), trimLeft = trimLeft, 
                        truncLen = truncLen, maxN = maxN, maxEE = maxEE, truncQ = truncQ, 
                        compress=TRUE, verbose=TRUE)
    }
    
    # check if files have been created
    if(!all(file.exists(filtFs)) || !all(file.exists(filtRs))) {
      cat("\n*** ERROR: Not all filtered files were created, maybe trimming impossible***", file = LogFile, append = TRUE)
      #stop("** Not all filtered files were created, maybe trimming impossible")
      
    }
    
    message("*********************** Filtering Done ***********************
            ********************************************************************")
    cat("\n*** Filtering Done ***", file = LogFile, append = TRUE)
    
  } else {
    
    if((names(filtFs) != SampleNames) || (names(filtRs) != SampleNames)) {
      
      cat("\n*** ERROR: The given filtFs or filtRs do not have sample names!***", file = LogFile, append = TRUE)
      stop("The given filtFs or filtRs do not have sample names!")
    }
    
    if(!all(file.exists(filtFs)) || !all(file.exists(filtRs))) {
      cat("\n*** ERROR: Not all files in the given filtFs or filtRs existed***", file = LogFile, append = TRUE)
      stop("** Not all files in the given filtFs or filtRs existed")
      
    }
    
    message("*********************** Filtered files were given ***********************
            ********************************************************************")
    cat("\n*** Filtered files were given ***", file = LogFile, append = TRUE)
    
    
  }
  
  
  save(PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  cat("\nTime after filtering step: ", file = LogFile, append = TRUE)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  TimePassed <- proc.time()-ptm
  cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
  cat("\n*** Start estimating err_F if not given ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Estimate err_F matrix
  ##############################
  
  if(is.null(err_F)){
    
    # the new learnErrors from dada2 allows setting nreads (usually 1 million), it then reads in as many samples until
    # the number of unique sequences (drp$uniques) reaches nreads.
    # learnErrorsAdj is basically dada2::learnErrors, just adds the NREADS from the function to the output list and also adds
    # dreplicated (drps) and denoised (dds) datasets if all samples were used for error matrix estimation. If not drps and dds
    # remain NULL
    errorsFW <- learnErrorsAdj(fls = filtFs, nreads = nreadsLearn, multithread = TRUE)
    
    message(paste("**", errorsFW$NREADS, "reads have been used for err_F estimation **"))
    
    err_F <- errorsFW$err_out
    
    pdf(file = file.path(PlotFolder, "errorRates_F.pdf"), width = 10, height = 10)
    print(plotErrors(errorsFW, nominalQ=TRUE))
    dev.off()
    
    save(err_F, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
    
    message("*********************** err_F has been estimated ***********************
            ********************************************************************")
    cat("\n*** err_F has been estimated ***", file = LogFile, append = TRUE)
    cat(paste("\nReads used for err_F estimation: ", errorsFW$NREADS), file = LogFile, append = TRUE)
    TimePassed <- proc.time()-ptm
    cat("\nTime after err_F estimation: ", file = LogFile, append = TRUE)
    cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
    cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
    cat("\n*** Start estimating err_R if not given ***", file = LogFile, append = TRUE)
    
  }
  
  # in case both err_F and err_R had been given, they would not have been saved if not:
  save(err_F, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  ##############################
  ### Estimate err_R matrix
  ##############################
  
  if(is.null(err_R)){
    
    
    errorsRV <- learnErrorsAdj(fls = filtRs, nreads = nreadsLearn, multithread = TRUE)
    
    message(paste("**", errorsRV$NREADS, "reads have been used for err_R estimation **"))
    
    err_R <- errorsRV$err_out
    
    pdf(file = file.path(PlotFolder, "errorRates_R.pdf"), width = 10, height = 10)
    print(plotErrors(errorsRV, nominalQ=TRUE))
    dev.off()
    
    save(err_F, err_R, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
    
    message("*********************** err_R has been estimated ***********************
            ********************************************************************")
    cat("\n*** err_R has been estimated ***", file = LogFile, append = TRUE)
    cat(paste("\nReads used for err_R estimation: ", errorsRV$NREADS), file = LogFile, append = TRUE)
    TimePassed <- proc.time()-ptm
    cat("\nTime after err_R estimation: ", file = LogFile, append = TRUE)
    cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
    cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
    cat("\n*** Start denoising data, bimera detection, and merging of reads into amplicons***", file = LogFile, append = TRUE)
    
  }
  
  save(err_F, err_R, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  ##############################
  ### Denoising (dada command for all samples)
  ##############################
  
  message("*********************** Start denoising and bimera detection ***********************
          ********************************************************************")
  
  if (pool) { # very time and memory consuming but might help identifying rare SVs
    
    drp_F <- derepFastq(filtFs, verbose = T)
    dd_F <- dada(drp_F, err = err_F, pool = pool, multithread = TRUE)
    drp_R <- derepFastq(filtRs, verbose = TRUE)
    dd_R <- dada(drp_R, err = err_R, pool = pool, multithread = TRUE)
    
    
    # NB: the following demands that the samples were denoised in the given order, I therefore removed tha randomize option from learnErrorsAdj
    names(dd_F) <- SampleNames
    names(drp_F) <- SampleNames
    names(dd_R) <- SampleNames
    names(drp_R) <- SampleNames
    
    # only for the ReadSummary later
    Uniques_F <- sapply(drp_F, function(x) length(x$uniques))
    Uniques_R <- sapply(drp_R, function(x) length(x$uniques))
    Denoised_F <- sapply(dd_F, function(x) length(x$denoised))
    Denoised_R <- sapply(dd_R, function(x) length(x$denoised))
    NoFilteredReads <- sapply(drp_F, function(x) sum(x$uniques)) # would be the same for drp_R, or sum dd_F$denoised
    
    mergers <- mergePairs(dd_F, drp_F, dd_R, drp_R, minOverlap = minOverlap, maxMismatch = maxMismatch, justConcatenate=justConcatenate) #!justConcatenate=TRUE
    
    rm(drp_F, drp_R, dd_F, dd_R)
    
    
  } else if (exists("errorsFW", inherits = FALSE) && exists("errorsRV", inherits = FALSE) && # do not do dereplication and denoising again if all samples had been used for determining the error matrixes
             !is.null(errorsFW$dds) && !is.null(errorsRV$dds)) {
    
    dd_F <- errorsFW$dds
    drp_F <- errorsFW$drps
    dd_R <- errorsRV$dds
    drp_R <- errorsRV$drps
    
    # just to not save the huge files:
    errorsFW$dds <- NULL
    errorsFW$drps <- NULL
    errorsRV$dds <- NULL
    errorsRV$drps <- NULL
    
    # NB: the following demands that the samples were denoised in the given order, I therefore removed tha randomize option from learnErrorsAdj
    names(dd_F) <- SampleNames
    names(drp_F) <- SampleNames
    names(dd_R) <- SampleNames
    names(drp_R) <- SampleNames
    
    # only for the ReadSummary later
    Uniques_F <- sapply(drp_F, function(x) length(x$uniques))
    Uniques_R <- sapply(drp_R, function(x) length(x$uniques))
    Denoised_F <- sapply(dd_F, function(x) length(x$denoised))
    Denoised_R <- sapply(dd_R, function(x) length(x$denoised))
    NoFilteredReads <- sapply(drp_F, function(x) sum(x$uniques)) # would be the same for drp_R, or sum dd_F$denoised
    
    mergers <- mergePairs(dd_F, drp_F, dd_R, drp_R, minOverlap = minOverlap, maxMismatch = maxMismatch)
    
    rm(drp_F, drp_R, dd_F, dd_R)
    
  } else {
    
    if (exists("errorsFW", inherits = FALSE)) {
      errorsFW$dds <- NULL
      errorsFW$drps <- NULL
    }
    if (exists("errorsRV", inherits = FALSE)) {
      errorsRV$dds <- NULL
      errorsRV$drps <- NULL
    }
    
    mergers <- vector("list", length(SampleNames))
    names(mergers) <- SampleNames
    NoFilteredReads <- vector("numeric", length(SampleNames))
    names(NoFilteredReads) <- SampleNames
    Uniques_F <- vector("numeric", length(SampleNames))
    names(Uniques_F) <- SampleNames
    Uniques_R <- vector("numeric", length(SampleNames))
    names(Uniques_R) <- SampleNames
    Denoised_F <- vector("numeric", length(SampleNames))
    names(Denoised_F) <- SampleNames
    Denoised_R <- vector("numeric", length(SampleNames))
    names(Denoised_R) <- SampleNames
    
    for(sam in SampleNames) {
      cat("Processing:", sam, "\n")
      cat(paste("\nDenoising sample:", sam), file = LogFile, append = TRUE)
      derepF <- derepFastq(filtFs[[sam]])
      NoFilteredReads[sam] <- sum(derepF$uniques)
      ddF <- dada(derepF, err=err_F, multithread=TRUE) 
      Uniques_F[sam] <- length(derepF$uniques)
      Denoised_F[sam] <- length(ddF$denoised)
      
      derepR <- derepFastq(filtRs[[sam]])
      ddR <- dada(derepR, err=err_R, multithread=TRUE)
      Uniques_R[sam] <- length(derepR$uniques)
      Denoised_R[sam] <- length(ddR$denoised)
      merger <- mergePairs(ddF, derepF, ddR, derepR, minOverlap = minOverlap, maxMismatch = maxMismatch)
      mergers[[sam]] <- merger
      
      rm(derepF, derepR, ddF, ddR, merger)
    }
    
  }
  
  
  
  
  if(!all((sapply(1:length(mergers), function(i) dim(mergers[[i]])[1]))!=0)){
    message("**Not in all samples merged amplicons could be found! maybe minOverlap too strict??**")
    cat("\n*** Not in all samples merged amplicons could be found! maybe minOverlap too strict?? ***", file = LogFile, append = TRUE)
  }
  
  if(all((sapply(1:length(mergers), function(i) dim(mergers[[i]])[1]))==0)){
    cat("\n*** ERROR: In none of the samples merged amplicons could be found! maybe minOverlap too strict?? ***", file = LogFile, append = TRUE)
    stop("**In none of the samples merged amplicons could be found! maybe minOverlap too strict??**")
  }
  
  message("*********************** all samples denoised mergerd amplicons generated ***********************
          ********************************************************************")
  cat("\n*** all samples denoised and mergerd amplicons generated ***", file = LogFile, append = TRUE)
  TimePassed <- proc.time()-ptm
  cat("\nTime after denoising: ", file = LogFile, append = TRUE)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
  cat("\n*** Start generating Sequence table and removing bimeras ***", file = LogFile, append = TRUE)
  
  
  ##############################
  ### Generating sequence table
  ##############################
  
  seqtab <- makeSequenceTable(mergers)
  
  
  ##############################
  ### Bimera removal
  ##############################
  
  # removes denoised sequence FW RV combi that was paired in merger for which the FW or the RV sequence was identified as bimera
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  # 1-(dim(seqtab.nochim)[2]/dim(seqtab)[2]) # The fraction chimeras made up of detected sequence variants
  # 1-(sum(seqtab.nochim)/sum(seqtab)) # The fraction of total reads that were chimeras (over all samples)
  
  message("*********************** Seq table generated and Bimeras removed ***********************
          ********************************************************************")
  cat("\n*** seqtable generated, bimeras removed ***", file = LogFile, append = TRUE)
  cat("\n*** Start saving the data ***", file = LogFile, append = TRUE)
  
  
  # generate also a read summary data frame
  ReadSummary <- data.frame(Sample = SampleNames, NoReads = sapply(F_QualityStats, function(df){df$NoReads[1]}))
  # NB, No total reads would be the same from R_QualityStats, could be used for sanity check
  ReadSummary$FilteredReads <- NoFilteredReads
  ReadSummary$MergedReads <- rowSums(seqtab)
  ReadSummary$MergedReadsWOBimera <- rowSums(seqtab.nochim)
  ReadSummary$UniqueSequences_F <- Uniques_F
  ReadSummary$DenoisedSequences_F <- Denoised_F
  ReadSummary$UniqueSequences_R <- Uniques_R
  ReadSummary$DenoisedSequences_R <- Denoised_R
  ReadSummary$Unique_Amplicons <- rowSums(seqtab != 0)
  ReadSummary$UniqueAmpliconsWOBimera <- rowSums(seqtab.nochim != 0)
  rownames(ReadSummary) <- NULL
  
  # I now also save errorsFW and RV to be able to recapitulate the error plots in later analyses
  if (!exists("errorsFW", inherits = FALSE)) {
    errorsFW <- NULL
  }
  if (!exists("errorsRV", inherits = FALSE)) {
    errorsRV <- NULL
  }
  
  save(errorsFW, errorsRV, seqtab.nochim, seqtab, mergers, ReadSummary, file = file.path(DataFolder, "DenoisedData.RData"))
  
  message("*********************** Sequence table generated, Data saved ***********************
          ********************************************************************")
  cat("\n*** Sequence table generated, Data saved ***", file = LogFile, append = TRUE)
  cat("\n*** Start generating Summary Plots ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Summary Plots
  ##############################
  
  # the width of the plots probably needs adjustment
  width = 5 + 0.5*(length(SampleNames)/10)
  
  # Plot the number of reads at the different steps for each sample (see also ReadSummary)
  pdf(file = file.path(PlotFolder, "NoReads_AllSamples.pdf"), width = width, height = 6)
  print(NoReads_StepsSimple(ReadSummary = ReadSummary, SampleNames = SampleNames, sort = TRUE))
  dev.off()
  
  # PLot the total number of amplicons against the number of unique amplicons for each sample
  FinalNumbers <- dplyr::select(ReadSummary, Sample, UniqueAmplicons = UniqueAmpliconsWOBimera, NoAmplicons = MergedReadsWOBimera)
  # FinalNumbers <- data.frame(Sample = rownames(seqtab.nochim), UniqueAmplicons = rowSums(seqtab.nochim != 0), NoAmplicons = rowSums(seqtab.nochim))
  
  Tr <- TotalandUniqueAmplicons(FinalNumbers = FinalNumbers, seqtab = seqtab.nochim)
  pdf(file = file.path(PlotFolder, "TotalvsUniqueAmplicons.pdf"), width = width, height = 6)
  print(Tr)
  dev.off()
  
  # Plot a linear regression line to illustrate the possible association between the total number of amplicons and the number of 
  # unique amplicons
  Tr2 <- ggplot(FinalNumbers, aes(x = NoAmplicons, y = UniqueAmplicons)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    xlab("Total Amplicons") +
    ylab("Unique Amplicons") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = "#999999", size = .15),
          panel.grid.major.x = element_line(color = "#999999", size = .15))
  
  
  pdf(file = file.path(PlotFolder, "AssociationTotaltoUniqueAmplicons.pdf"), width = 7, height = 6)
  print(Tr2)
  dev.off()
  
  # Plot a histogram illustrating in how many samples the amplicons are present
  if(length(SampleNames) < 10) {binwidth = 1}
  if(length(SampleNames) > 10 && length(SampleNames) < 100) {binwidth = 2}
  if(length(SampleNames) > 100) {binwidth = 3}
  
  FinalNumbers2 <- data.frame(InNoSamples = colSums(seqtab.nochim != 0))
  
  Trr <- ggplot(data = FinalNumbers2, aes(x = InNoSamples))  + 
    geom_histogram(binwidth = binwidth, col = "black", fill = "#E69F00")+
    geom_rug() +
    xlab("Present in No Samples") + 
    ylab("Count") +
    theme_bw() + 
    ggtitle(paste("Total No of unique amplicons:", dim(seqtab.nochim)[2], "Only in 1 sample:", sum(FinalNumbers2 ==1))) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(color = "#999999", size = .15)) +
    coord_cartesian(ylim = c(0,150))
  
  pdf(file = file.path(PlotFolder, "HistogramAmpliconsinNoSamples.pdf"), width = 7, height = 6)
  print(Trr)
  dev.off()
  cat("\n*** Summary Plots generated, Function done ***", file = LogFile, append = TRUE)
  message("*********************** Dada_wrap done ***********************
          ********************************************************************")
}





#######################################
### FUNCTION: Dada2_QualityCheck
#######################################

# runs the first part of the Dada2_wrap function, creating the quality plots and data. Based on these one can decide on the filtering
# parameters to subsequently use the Dada2_wrap function
## Input
# path: The path to the folder containing the sample folders with the fastq files: 
# F_pattern: a regular expression to find the fastq files with the forward reads in the sample folders
# R_pattern: a regular expression to find the fastq files with the reverse reads in the sample folders
# path2: default = NULL, the function creates the folders Dada_Plots and Dada_Data, it creates them by default in the path folder, when path2 is given 
# the folders will instead be generated in path2
## Output
# PLOTS:
# Quality Plots saved as pdfs in the generated Dada_Plots folder # NB: the plots are currently assuming 250 nt
# DATA:
# saved in the Dada_Data folder:
# QualityStats.RData contains PackageVersions, F_QualityStats, and R_QualityStats

Dada2_QualityCheck <- function(path, F_pattern, R_pattern, path2 = NULL) {
  
  ##############################
  ### call the required packages
  ##############################
  
  ## dada2:
  # source("https://bioconductor.org/biocLite.R")
  try(library(dada2), biocLite("dada2"))
  
  ## Short Read
  try(library(ShortRead), biocLite("ShortRead"))
  
  ## ggplot2
  try(library(ggplot2), install.packages("ggplot2"))
  
  ## dplyr
  try(library(dplyr), install.packages("dplyr"))
  
  ## dplyr
  try(library(tidyr), install.packages("tidyr"))
  
  
  message("*********************** All packages loaded ***********************
          ********************************************************************")
  
  ##############################
  ### save the package Versions
  ##############################
  # NB: outputs an error and stops function if a Package is not installed
  PackageVersions <- data.frame(Package = c("dada2", "ShortRead", "ggplot2", "dplyr", "tidyr"),
                                Version = c(packageVersion("dada2"),
                                            packageVersion("ShortRead"),
                                            packageVersion("ggplot2"),
                                            packageVersion("dplyr"),
                                            packageVersion("tidyr")))
  
  message(paste(PackageVersions$Package, ": ", PackageVersions$Version, "; ", sep= ""))
  
  ##############################
  ### Construct character vectors to the FW and RV fastq files
  ##############################
  if(is.null(path2)){
    
    path2 = path
  }
  
  
  folders <- list.dirs(path, recursive = FALSE, full.names = FALSE)
  
  if(length(folders) != 0) {
    
    SampleNames <- folders
    
    if(sum(grepl("^Dada", SampleNames)) != 0){
      warning("**There are folders starting with \"Dada\" in your path folder, maybe you have run the function on this path folder before\n.
              The folders starting with Dada will not be considered sample folders!!\n Files within Dada_Data Dada_FilteredFastqs and Dada_Plots will be overwritten**")
    }
    
    # exclude Folders starting with "Dada" from the folders considered as sample fodlers
    if(length(grep("^Dada", SampleNames))!=0) {
      SampleNames <- SampleNames[-grep("^Dada", SampleNames)]
    }
    
    F_fastq <- character(length = length(SampleNames))
    R_fastq <- character(length = length(SampleNames))
    
    for (i in 1:length(SampleNames)) {
      CurrentPath <- file.path(path, SampleNames[i])
      
      if(sum(grepl(F_pattern, list.files(CurrentPath))) == 0) {
        stop(paste("F_pattern fits no file in ", CurrentPath))
      }
      if(sum(grepl(F_pattern, list.files(CurrentPath))) > 1) {
        stop(paste("F_pattern fits several files in ", CurrentPath))
      }
      if(sum(grepl(R_pattern, list.files(CurrentPath))) == 0) {
        stop(paste("R_pattern fits no file in ", CurrentPath))
      }
      if(sum(grepl(R_pattern, list.files(CurrentPath))) > 1) {
        stop(paste("R_pattern fits several files in ", CurrentPath))
      }
      F_fastq[i] <- file.path(CurrentPath, list.files(CurrentPath)[grepl(F_pattern, list.files(CurrentPath))])
      R_fastq[i] <- file.path(CurrentPath, list.files(CurrentPath)[grepl(R_pattern, list.files(CurrentPath))])
      
    } 
    
    } else {
      
      stop("No sample folders were found in the given path! Currently the Dada2_wrap function can only handle the situation where the fastq files are in separate folders for each sample.
           These folders have to be in the \"path\" folder. Other situations have to be added.")
    }
  
  
  ##############################
  ### Determine the quality scores and save the stats in Data folder
  ##############################
  
  DataFolder <- file.path(path2, "Dada_Data")
  
  # Not sure if this is wanted, deleting all files that are already in the DataFolder folder
  if(file.exists(DataFolder)){
    file.remove(list.files(DataFolder, full.names = TRUE))
  }
  
  dir.create(DataFolder, showWarnings = TRUE)
  
  if(!file.exists(DataFolder)){
    stop("DataFolder could not be created! something wrong with path2, maybe permission denied.")
  }
  
  # Collect quality Score data of the fastQ files and store a data.frame for each sample in the lists FW_QualityStats and RV_QualityStats ####################
  
  F_QualityStats <- list()
  R_QualityStats <- list()
  
  for (i in seq_along(F_fastq)) {
    
    Current_FWfq <- F_fastq[i]
    Current_dfFW <- qa(Current_FWfq, n = 1e06)[["perCycle"]]$quality
    # df is a data frame containing for each cycle (nt) the distribution of Quality scores,
    # e.g. Cycle 1 had 7 different quality scores, then 7 rows of cycle one, for each score the count says how many reads had this score
    Current_dfFW <- dplyr::group_by(Current_dfFW, Cycle)
    
    Current_dfQStatsFW <- dplyr::summarise(
      Current_dfFW,
      NoReads = sum(Count),
      Mean_QS = sum(Count*Score)/sum(Count),
      SD_QS = sqrt(sum(Count*((Score-Mean_QS)^2))/(NoReads-1)),
      Median_QS = Score[which(cumsum(Count)/sum(Count) >= .5)][[1]],
      q25_QS = Score[which(cumsum(Count)/sum(Count) >= .25)][[1]],
      q75_QS = Score[which(cumsum(Count)/sum(Count) >= .75)][[1]])
    
    ## Check that all reads are of same length
    # x <- range(Current_dfQStatsFW$NoReads)
    # if (!all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5)) {
    #         stop("Not all reads of same length in file", F_fastq[i])
    # }
    
    F_QualityStats[[i]] <- as.data.frame(Current_dfQStatsFW)
    # as.data.frame to un-dplyr the data.frame
    
    ####### collect the same stats for the RV FastQ files
    Current_RVfq <- R_fastq[i]
    Current_dfRV <- qa(Current_RVfq, n = 1e06)[["perCycle"]]$quality
    
    Current_dfRV <- dplyr::group_by(Current_dfRV, Cycle)
    
    Current_dfQStatsRV <- dplyr::summarise(
      Current_dfRV,
      NoReads = sum(Count),
      Mean_QS = sum(Count*Score)/sum(Count),
      SD_QS = sqrt(sum(Count*((Score-Mean_QS)^2))/(NoReads-1)),
      Median_QS = Score[which(cumsum(Count)/sum(Count) >= .5)][[1]],
      q25_QS = Score[which(cumsum(Count)/sum(Count) >= .25)][[1]],
      q75_QS = Score[which(cumsum(Count)/sum(Count) >= .75)][[1]])
    
    ## Check that all reads are of same length
    # x <- range(Current_dfQStatsRV$NoReads)
    # if (!all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5)) {
    #         stop("Not all reads of same length in file", R_fastq[i])
    # }
    
    R_QualityStats[[i]] <- as.data.frame(Current_dfQStatsRV)
    
    rm(Current_FWfq, Current_dfFW, Current_dfQStatsFW,
       Current_RVfq, Current_dfRV, Current_dfQStatsRV, x)
    
  }
  
  # add the sample names as names to the lists:
  names(F_QualityStats) <- SampleNames
  names(R_QualityStats) <- SampleNames
  
  save(PackageVersions, F_QualityStats, R_QualityStats, file = file.path(DataFolder, "QualityStats.RData"))
  
  message("*********************** Quality Stats Collected ***********************
          ********************************************************************")
  
  ##############################
  ### Generate and save some quality plots
  ##############################
  
  PlotFolder <- file.path(path2, "Dada_Plots")
  
  ## Not sure if this is wanted, deleting all files that are already in the PlotFolder folder
  if(file.exists(PlotFolder)){
    file.remove(list.files(PlotFolder, full.names = TRUE))
  }
  
  dir.create(PlotFolder, showWarnings = FALSE)
  
  
  pdf(file = file.path(PlotFolder, "MedianQScore_F_AllNucleotides.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(F_QualityStats, SampleNames, Prefix = "FW"))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_F_150to240.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(F_QualityStats, SampleNames, xlim_low = 150, xlim_high = 240, Prefix = "FW"))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_R_AllNucleotides.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(R_QualityStats, SampleNames, Prefix = "RV"))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_R_150to240.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(R_QualityStats, SampleNames, xlim_low = 150, xlim_high = 240, Prefix = "RV"))
  dev.off()
  
  message("*********************** Quality plots generated ***********************
          ********************************************************************")
}


#######################################
### FUNCTION: learnErrorsAdj
#######################################
# is basically the dada2:::learnErrors function unchanged, only that the NREADS used are added to the output list,
# in addition: if dada command has been run on all samples, dds is also saved to not have to run it again in 
# DadaWrapper.
# NB: this makes this function a bit dangerous if fls were already dereplicated samples, so keep fls to be file names!


learnErrorsAdj <- function (fls, nreads = 1e+06, errorEstimationFunction = loessErrfun, 
                            multithread = FALSE) #, randomize = FALSE 
{
  NREADS <- 0
  drps <- vector("list", length(fls))
  # if (randomize) {
  #         fls <- sample(fls)
  # }
  for (i in seq_along(fls)) {
    if (dada2:::is.list.of(fls, "derep")) {
      drps[[i]] <- fls[[i]]
    }
    else {
      drps[[i]] <- derepFastq(fls[[i]])
    }
    NREADS <- NREADS + sum(drps[[i]]$uniques)
    if (NREADS > nreads) {
      break
    }
  }
  drps <- drps[1:i]
  dds <- dada(drps, err = NULL, selfConsist = TRUE, multithread = multithread)
  # cat("Total reads used: ", NREADS, "\n")
  ErrorList <- getErrors(dds, detailed = TRUE)
  ErrorList$NREADS <- NREADS
  if (length(drps) == length(fls)){
    ErrorList$dds <- dds 
    ErrorList$drps <- drps
  } else {
    ErrorList$dds <- NULL
    ErrorList$drps <- NULL
  }
  return(ErrorList)
}
































##################################### alternatives ########################################

#######################################
### FUNCTION: Dada2_wrap_BimFWRV
#######################################

# basically the same as Dada2_wrap, just that the bimeras are calculated individually for FW and RV reads.
# That was done in the original dada2 paper but has been replaced by removeBimeraDenovo on the final seqtab
# in the dada2 tutorials. 
# Removal of bimeras individually on FW and RV reads as in this function usually removes slightly more amplicons as
# bimeras than does Dada2_wrap

Dada2_wrap_BimFWRV <- function(path, F_pattern, R_pattern, path2 = NULL,
                               trimLeft = c(10,10), truncLen = c(220, 160), 
                               maxN = 0, maxEE = 2, truncQ = 2,
                               nreadsLearn = 1e+06,
                               err_F = NULL,
                               err_R = NULL,
                               minOverlap = 20,
                               maxMismatch = 0,
                               F_QualityStats = NULL,
                               R_QualityStats = NULL,
                               filtFs = NULL,
                               filtRs = NULL) {
  
  ##############################
  ### call the required packages
  ##############################
  
  ## dada2:
  # source("https://bioconductor.org/biocLite.R")
  try(library(dada2), biocLite("dada2"))
  
  ## Short Read
  try(library(ShortRead), biocLite("ShortRead"))
  
  ## ggplot2
  try(library(ggplot2), install.packages("ggplot2"))
  
  ## dplyr
  try(library(dplyr), install.packages("dplyr"))
  
  ## dplyr
  try(library(tidyr), install.packages("tidyr"))
  
  
  message("*********************** All packages loaded ***********************
          ********************************************************************")
  
  ##############################
  ### save the package Versions
  ##############################
  # NB: outputs an error and stops function if a Package is not installed
  PackageVersions <- data.frame(Package = c("dada2", "ShortRead", "ggplot2", "dplyr", "tidyr"),
                                Version = c(packageVersion("dada2"),
                                            packageVersion("ShortRead"),
                                            packageVersion("ggplot2"),
                                            packageVersion("dplyr"),
                                            packageVersion("tidyr")))
  
  message(paste(PackageVersions$Package, ": ", PackageVersions$Version, "; ", sep= ""))
  
  
  ##############################
  ### Construct character vectors to the FW and RV fastq files
  ##############################
  
  if(is.null(path2)){
    
    path2 = path
  }
  
  
  folders <- list.dirs(path, recursive = FALSE, full.names = FALSE)
  
  if(length(folders) != 0) {
    
    SampleNames <- folders
    
    if(sum(grepl("^Dada", SampleNames)) != 0){
      warning("**There are folders starting with \"Dada\" in your path folder, maybe you have run the function on this path folder before\n.
              The folders starting with Dada will not be considered sample folders!!\n Files within Dada_Data Dada_FilteredFastqs and Dada_Plots will be overwritten**")
    }
    
    # exclude Folders starting with "Dada" from the folders considered as sample fodlers
    if(length(grep("^Dada", SampleNames))!=0) {
      SampleNames <- SampleNames[-grep("^Dada", SampleNames)]
    }
    
    F_fastq <- character(length = length(SampleNames))
    R_fastq <- character(length = length(SampleNames))
    
    for (i in 1:length(SampleNames)) {
      CurrentPath <- file.path(path, SampleNames[i])
      
      if(sum(grepl(F_pattern, list.files(CurrentPath))) == 0) {
        stop(paste("F_pattern fits no file in ", CurrentPath))
      }
      if(sum(grepl(F_pattern, list.files(CurrentPath))) > 1) {
        stop(paste("F_pattern fits several files in ", CurrentPath))
      }
      if(sum(grepl(R_pattern, list.files(CurrentPath))) == 0) {
        stop(paste("R_pattern fits no file in ", CurrentPath))
      }
      if(sum(grepl(R_pattern, list.files(CurrentPath))) > 1) {
        stop(paste("R_pattern fits several files in ", CurrentPath))
      }
      F_fastq[i] <- file.path(CurrentPath, list.files(CurrentPath)[grepl(F_pattern, list.files(CurrentPath))])
      R_fastq[i] <- file.path(CurrentPath, list.files(CurrentPath)[grepl(R_pattern, list.files(CurrentPath))])
      
    }
    
    } else {
      
      stop("No sample folders were found in the given path! Currently the Dada2_wrap function can only handle the situation where the fastq files are in separate folders for each sample.
           These folders have to be in the \"path\" folder. Other situations have to be added.")
    }
  
  ##############################
  ### Start the log file
  ##############################
  
  DataFolder <- file.path(path2, "Dada_Data")
  
  # ATTENTION: FUNCTION DELETES ALL FILES ALREADY PRESRENT in DataFolder
  if(file.exists(DataFolder)){
    file.remove(list.files(DataFolder, full.names = TRUE))
  }
  
  
  dir.create(DataFolder, showWarnings = FALSE)
  
  if(!file.exists(DataFolder)){
    stop("DataFolder could not be created! something wrong with path2, maybe permission denied.")
  }
  
  ptm <- proc.time()
  
  LogFile <- file.path(DataFolder, "DadaWrapper.log")
  cat("Time after package installation: ", file = LogFile)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  cat("Package Versions: ", file = LogFile, append = TRUE)
  cat(paste(PackageVersions$Package, ": ", PackageVersions$Version, "; ", sep= ""), file = LogFile, append = TRUE)
  
  
  # Collect inputs in a data frame for saving it later
  Input <- list(path = path, F_pattern = F_pattern, R_pattern = R_pattern, path2 = path2,
                trimLeft = trimLeft, truncLen = truncLen, 
                maxN = maxN, maxEE = maxEE, truncQ = truncQ,
                nreadsLearn = nreadsLearn, err_F = err_F, err_R = err_R,
                minOverlap = minOverlap, maxMismatch = maxMismatch, F_QualityStats = F_QualityStats,
                R_QualityStats = R_QualityStats, filtFs = filtFs, filtRs = filtRs)
  
  
  ##############################
  ### Determine the quality scores and save the stats in Data folder
  ##############################
  
  # Collect quality Score data of the fastQ files and store a data.frame for each sample in the lists FW_QualityStats and RV_QualityStats ####################
  
  if (is.null(F_QualityStats) || is.null(R_QualityStats)) {
    
    F_QualityStats <- list()
    R_QualityStats <- list()
    
    for (i in seq_along(F_fastq)) {
      
      Current_FWfq <- F_fastq[i]
      Current_dfFW <- qa(Current_FWfq, n = 1e06)[["perCycle"]]$quality
      # df is a data frame containing for each cycle (nt) the distribution of Quality scores
      # e.g. Cycle 1 had 7 different quality scores, then 7 rows of cycle one, for each score the count says how many reads had this score
      Current_dfFW <- dplyr::group_by(Current_dfFW, Cycle)
      
      Current_dfQStatsFW <- dplyr::summarise(
        Current_dfFW,
        NoReads = sum(Count),
        Mean_QS = sum(Count*Score)/sum(Count),
        SD_QS = sqrt(sum(Count*((Score-Mean_QS)^2))/(NoReads-1)),
        Median_QS = Score[which(cumsum(Count)/sum(Count) >= .5)][[1]],
        q25_QS = Score[which(cumsum(Count)/sum(Count) >= .25)][[1]],
        q75_QS = Score[which(cumsum(Count)/sum(Count) >= .75)][[1]])
      
      # Check that all reads are of same length
      # x <- range(Current_dfQStatsFW$NoReads)
      # if (!all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5)) {
      #         stop("Not all reads of same length in file", F_fastq[i])
      # }
      # 
      F_QualityStats[[i]] <- as.data.frame(Current_dfQStatsFW)
      # as.data.frame to un-dplyr the data.frame
      
      # collect the same stats for the RV FastQ files
      Current_RVfq <- R_fastq[i]
      Current_dfRV <- qa(Current_RVfq, n = 1e06)[["perCycle"]]$quality
      
      Current_dfRV <- dplyr::group_by(Current_dfRV, Cycle)
      
      Current_dfQStatsRV <- dplyr::summarise(
        Current_dfRV,
        NoReads = sum(Count),
        Mean_QS = sum(Count*Score)/sum(Count),
        SD_QS = sqrt(sum(Count*((Score-Mean_QS)^2))/(NoReads-1)),
        Median_QS = Score[which(cumsum(Count)/sum(Count) >= .5)][[1]],
        q25_QS = Score[which(cumsum(Count)/sum(Count) >= .25)][[1]],
        q75_QS = Score[which(cumsum(Count)/sum(Count) >= .75)][[1]])
      
      ## Check that all reads are of same length
      # x <- range(Current_dfQStatsRV$NoReads)
      # if (!all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5)) {
      #         stop("Not all reads of same length in file", R_fastq[i])
      # }
      # 
      R_QualityStats[[i]] <- as.data.frame(Current_dfQStatsRV)
      
      rm(Current_FWfq, Current_dfFW, Current_dfQStatsFW,
         Current_RVfq, Current_dfRV, Current_dfQStatsRV, x)
      
    }
    
    # add the sample names as names to the lists:
    names(F_QualityStats) <- SampleNames
    names(R_QualityStats) <- SampleNames
    
    message("*********************** Quality Stats Collected ***********************
            ********************************************************************")
    cat("\n*** Quality Stats Collected ***", file = LogFile, append = TRUE)
    
  } else {
    
    if((names(F_QualityStats) != SampleNames) || (names(R_QualityStats) != SampleNames)) {
      
      stop("The given F_QualityStats or R_QualityStats do not fit to the SampleNames in path!")
    }
    
    message("*********************** Quality Stats given ***********************
            ********************************************************************")
    cat("\n*** Quality Stats given***", file = LogFile, append = TRUE)
    
  }
  
  save(PackageVersions, F_QualityStats, R_QualityStats, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  TimePassed <- proc.time()-ptm
  cat("\nTime after Quality Stats collection: ", file = LogFile, append = TRUE)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
  cat("\n*** Start generating Quality Plots ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Generate and save some quality plots
  ##############################
  
  PlotFolder <- file.path(path2, "Dada_Plots")
  
  # ATTENTION: FUNCTION DELETES ALL FILES ALREADY PRESRENT in PlotFolder
  if(file.exists(PlotFolder)){
    file.remove(list.files(PlotFolder, full.names = TRUE))
  }
  
  dir.create(PlotFolder, showWarnings = FALSE)
  
  
  pdf(file = file.path(PlotFolder, "MedianQScore_F_AllNucleotides.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(F_QualityStats, SampleNames))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_F_150to240.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(F_QualityStats, SampleNames, xlim_low = 150, xlim_high = 240))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_R_AllNucleotides.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(R_QualityStats, SampleNames))
  dev.off()
  
  pdf(file = file.path(PlotFolder, "MedianQScore_R_150to240.pdf"), width = 7, height = 6)
  print(QS_Median_OverviewPlot(R_QualityStats, SampleNames, xlim_low = 150, xlim_high = 240))
  dev.off()
  
  message("*********************** Plots generated start filtering ***********************
          ********************************************************************")
  
  cat("\n*** Quality Plots generated ***", file = LogFile, append = TRUE)
  cat("\n*** Start Filtering ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Filtering
  ##############################
  
  FilteredFolder <- file.path(path2, "Dada_FilteredFastqs")
  
  
  if(is.null(filtFs) || is.null(filtRs)){
    
    ## Not sure if this is wanted, deleting all files that are already in the DataFolder folder
    if(file.exists(FilteredFolder)){
      file.remove(list.files(FilteredFolder, full.names = TRUE))
    }
    
    dir.create(FilteredFolder, showWarnings = FALSE)
    
    filtFs <- file.path(FilteredFolder, paste0(SampleNames, "_F_Filtered.fastq.gz"))
    filtRs <- file.path(FilteredFolder, paste0(SampleNames, "_R_Filtered.fastq.gz"))
    names(filtFs) <- SampleNames
    names(filtRs) <- SampleNames
    
    
    for(i in seq_along(F_fastq)) {
      
      message(paste("**Filtering sample No:", i, "called:", SampleNames[i]), " **")
      cat(paste("\n**Filtering sample No:", i, "called:", SampleNames[i]), file = LogFile, append = TRUE)
      fastqPairedFilter(c(F_fastq[i], R_fastq[i]), c(filtFs[i], filtRs[i]), trimLeft = trimLeft, 
                        truncLen = truncLen, maxN = maxN, maxEE = maxEE, truncQ = truncQ, 
                        compress=TRUE, verbose=TRUE)
    }
    
    # check if files have been created
    if(!all(file.exists(filtFs)) || !all(file.exists(filtRs))) {
      cat("\n*** ERROR: Not all filtered files were created, maybe trimming impossible***", file = LogFile, append = TRUE)
      #stop("** Not all filtered files were created, maybe trimming impossible")
      
    }
    
    message("*********************** Filtering Done ***********************
            ********************************************************************")
    cat("\n*** Filtering Done ***", file = LogFile, append = TRUE)
    
  } else {
    
    if((names(filtFs) != SampleNames) || (names(filtRs) != SampleNames)) {
      
      cat("\n*** ERROR: The given filtFs or filtRs do not have sample names!***", file = LogFile, append = TRUE)
      stop("The given filtFs or filtRs do not have sample names!")
    }
    
    if(!all(file.exists(filtFs)) || !all(file.exists(filtRs))) {
      cat("\n*** ERROR: Not all files in the given filtFs or filtRs existed***", file = LogFile, append = TRUE)
      stop("** Not all files in the given filtFs or filtRs existed")
      
    }
    
    message("*********************** Filtered files were given ***********************
            ********************************************************************")
    cat("\n*** Filtered files were given ***", file = LogFile, append = TRUE)
    
    
  }
  
  
  save(PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
  
  cat("\nTime after filtering step: ", file = LogFile, append = TRUE)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  TimePassed <- proc.time()-ptm
  cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
  cat("\n*** Start estimating err_F if not given ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Estimate err_F matrix
  ##############################
  
  if(is.null(err_F)){
    
    # the new learnErrors from dada2 allows setting nreads (usually 1 million), it then reads in as many samples until
    # the number of unique sequences (drp$uniques) reaches nreads.
    # learnErrorsAdj is basically dada2::learnErrors, just adds the NREADS from the function to the output list and also adds
    # dreplicated (drps) and denoised (dds) datasets if all samples were used for error matrix estimation. If not drps and dds
    # remain NULL
    errorsFW <- learnErrorsAdj(fls = filtFs, nreads = nreadsLearn, multithread = TRUE)
    
    message(paste("**", errorsFW$NREADS, "reads have been used for err_F estimation **"))
    
    err_F <- errorsFW$err_out
    
    pdf(file = file.path(PlotFolder, "errorRates_F.pdf"), width = 10, height = 10)
    print(plotErrors(errorsFW, nominalQ=TRUE))
    dev.off()
    
    save(err_F, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
    
    message("*********************** err_F has been estimated ***********************
            ********************************************************************")
    cat("\n*** err_F has been estimated ***", file = LogFile, append = TRUE)
    cat(paste("\nReads used for err_F estimation: ", errorsFW$NREADS), file = LogFile, append = TRUE)
    TimePassed <- proc.time()-ptm
    cat("\nTime after err_F estimation: ", file = LogFile, append = TRUE)
    cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
    cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
    cat("\n*** Start estimating err_R if not given ***", file = LogFile, append = TRUE)
    
  }
  
  ##############################
  ### Estimate err_R matrix
  ##############################
  
  if(is.null(err_R)){
    
    
    errorsRV <- learnErrorsAdj(fls = filtRs, nreads = nreadsLearn, multithread = TRUE)
    
    message(paste("**", errorsRV$NREADS, "reads have been used for err_R estimation **"))
    
    err_R <- errorsRV$err_out
    
    pdf(file = file.path(PlotFolder, "errorRates_R.pdf"), width = 10, height = 10)
    print(plotErrors(errorsRV, nominalQ=TRUE))
    dev.off()
    
    save(err_F, err_R, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, Input, file = file.path(DataFolder, "QualityStats.RData"))
    
    message("*********************** err_R has been estimated ***********************
            ********************************************************************")
    cat("\n*** err_R has been estimated ***", file = LogFile, append = TRUE)
    cat(paste("\nReads used for err_R estimation: ", errorsRV$NREADS), file = LogFile, append = TRUE)
    TimePassed <- proc.time()-ptm
    cat("\nTime after err_R estimation: ", file = LogFile, append = TRUE)
    cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
    cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
    cat("\n*** Start denoising data, bimera detection, and merging of reads into amplicons***", file = LogFile, append = TRUE)
    
  }
  
  
  ##############################
  ### Denoising (dada command for all samples) and bimeara identification
  ##############################
  
  message("*********************** Start denoising and bimera detection ***********************
          ********************************************************************")
  
  # do not do dereplication and denoising again if all samples had been used for determining the error matrixes
  if (exists("errorsFW", inherits = FALSE) && exists("errorsRV", inherits = FALSE) &&
      !is.null(errorsFW$dds) && !is.null(errorsRV$dds)) {
    
    dd_F <- errorsFW$dds
    drp_F <- errorsFW$drps
    dd_R <- errorsRV$dds
    drp_R <- errorsRV$drps
    
    rm(errorsRV, errorsFW)
    
    # NB: the following demands that the samples were denoised in the given order, I therefore removed tha randomize option from learnErrorsAdj
    names(dd_F) <- SampleNames
    names(drp_F) <- SampleNames
    names(dd_R) <- SampleNames
    names(drp_R) <- SampleNames
    
    # only for the ReadSummary later
    Uniques_F <- sapply(drp_F, function(x) length(x$uniques))
    Uniques_R <- sapply(drp_R, function(x) length(x$uniques))
    Denoised_F <- sapply(dd_F, function(x) length(x$denoised))
    Denoised_R <- sapply(dd_R, function(x) length(x$denoised))
    NoFilteredReads <- sapply(drp_F, function(x) sum(x$uniques)) # would be the same for drp_R, or sum dd_F$denoised
    
    mergers <- mergePairs(dd_F, drp_F, dd_R, drp_R, minOverlap = minOverlap, maxMismatch = maxMismatch)
    
    bimFs <- vector("list", length(SampleNames))
    names(bimFs) <- SampleNames
    bimRs <- vector("list", length(SampleNames))
    names(bimRs) <- SampleNames
    
    for(sam in SampleNames) {
      cat("Processing:", sam, "\n")
      #cat(paste("\nDenoising sample:", sam), file = LogFile, append = TRUE)
      ddF <- dd_F[[sam]] 
      bimFs[[sam]] <- isBimeraDenovo(ddF, verbose=TRUE)
      ddR <- dd_R[[sam]]
      bimRs[[sam]] <- isBimeraDenovo(ddR, verbose=TRUE)
      
      rm(ddF, ddR)
    }
    
    rm(drp_F, drp_R, dd_F, dd_R)
    
  } else {
    
    if (exists("errorsFW", inherits = FALSE)) {
      rm(errorsFW)
    }
    if (exists("errorsRV", inherits = FALSE)) {
      rm(errorsRV)
    }
    
    mergers <- vector("list", length(SampleNames))
    names(mergers) <- SampleNames
    NoFilteredReads <- vector("numeric", length(SampleNames))
    names(NoFilteredReads) <- SampleNames
    bimFs <- vector("list", length(SampleNames))
    names(bimFs) <- SampleNames
    bimRs <- vector("list", length(SampleNames))
    names(bimRs) <- SampleNames
    Uniques_F <- vector("numeric", length(SampleNames))
    names(Uniques_F) <- SampleNames
    Uniques_R <- vector("numeric", length(SampleNames))
    names(Uniques_R) <- SampleNames
    Denoised_F <- vector("numeric", length(SampleNames))
    names(Denoised_F) <- SampleNames
    Denoised_R <- vector("numeric", length(SampleNames))
    names(Denoised_R) <- SampleNames
    
    for(sam in SampleNames) {
      cat("Processing:", sam, "\n")
      cat(paste("\nDenoising sample:", sam), file = LogFile, append = TRUE)
      derepF <- derepFastq(filtFs[[sam]])
      NoFilteredReads[sam] <- sum(derepF$uniques)
      ddF <- dada(derepF, err=err_F, multithread=TRUE) 
      bimFs[[sam]] <- isBimeraDenovo(ddF, verbose=TRUE)
      Uniques_F[sam] <- length(derepF$uniques)
      Denoised_F[sam] <- length(ddF$denoised)
      
      derepR <- derepFastq(filtRs[[sam]])
      ddR <- dada(derepR, err=err_R, multithread=TRUE)
      bimRs[[sam]] <- isBimeraDenovo(ddR, verbose=TRUE)
      Uniques_R[sam] <- length(derepR$uniques)
      Denoised_R[sam] <- length(ddR$denoised)
      merger <- mergePairs(ddF, derepF, ddR, derepR, minOverlap = minOverlap, maxMismatch = maxMismatch)
      mergers[[sam]] <- merger
      
      rm(derepF, derepR, ddF, ddR, merger)
    }
    
  }
  
  if(!all((sapply(1:length(mergers), function(i) dim(mergers[[i]])[1]))!=0)){
    message("**Not in all samples merged amplicons could be found! maybe minOverlap too strict??**")
    cat("\n*** Not in all samples merged amplicons could be found! maybe minOverlap too strict?? ***", file = LogFile, append = TRUE)
  }
  
  if(all((sapply(1:length(mergers), function(i) dim(mergers[[i]])[1]))==0)){
    cat("\n*** ERROR: In none of the samples merged amplicons could be found! maybe minOverlap too strict?? ***", file = LogFile, append = TRUE)
    stop("**In none of the samples merged amplicons could be found! maybe minOverlap too strict??**")
  }
  
  message("*********************** all samples denoised, bimeras identified, mergerd amplicons generated ***********************
          ********************************************************************")
  cat("\n*** all samples denoised, bimeras identified, mergerd amplicons generated ***", file = LogFile, append = TRUE)
  TimePassed <- proc.time()-ptm
  cat("\nTime after denoising: ", file = LogFile, append = TRUE)
  cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
  cat(paste("Time Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
  cat("\n*** Start removing bimera ***", file = LogFile, append = TRUE)
  
  
  ##############################
  ### Bimera removal
  ##############################
  
  # removes denoised sequence FW RV combi that was paired in merger for which the FW or the RV sequence was identified as bimera
  mergers.nochim <- mergers
  for (i in seq_along(mergers)) {
    mergers.nochim[[i]] <- mergers[[i]][!bimFs[[i]][mergers[[i]]$forward] & !bimRs[[i]][mergers[[i]]$reverse],]
  }
  
  message("*********************** Bimeras removed ***********************
          ********************************************************************")
  cat("\n*** bimeras removed ***", file = LogFile, append = TRUE)
  cat("\n*** Start generating Sequence table ***", file = LogFile, append = TRUE)
  
  
  ##############################
  ### Generating sequence table
  ##############################
  seqtab <- makeSequenceTable(mergers)
  seqtab.nochim <- makeSequenceTable(mergers.nochim)
  # 1-(dim(seqtab.nochimFWRW)[2]/dim(seqtab)[2]) # The fraction chimeras made up of detected sequence variants
  # 1-(sum(seqtab.nochimFWRW)/sum(seqtab)) # The fraction of total reads that were chimeras (over all samples)
  
  # generate also a read summary data frame
  ReadSummary <- data.frame(Sample = SampleNames, NoReads = sapply(F_QualityStats, function(df){df$NoReads[1]}))
  # NB, No total reads would be the same from R_QualityStats, could be used for sanity check
  ReadSummary$FilteredReads <- NoFilteredReads
  ReadSummary$MergedReads <- rowSums(seqtab)
  ReadSummary$MergedReadsWOBimera <- rowSums(seqtab.nochim)
  ReadSummary$UniqueSequences_F <- Uniques_F
  ReadSummary$DenoisedSequences_F <- Denoised_F
  ReadSummary$UniqueSequences_R <- Uniques_R
  ReadSummary$DenoisedSequences_R <- Denoised_R
  ReadSummary$Unique_Amplicons <- rowSums(seqtab != 0)
  ReadSummary$UniqueAmpliconsWOBimera <- rowSums(seqtab.nochim != 0)
  rownames(ReadSummary) <- NULL
  
  save(seqtab, seqtab.nochim, mergers, mergers.nochim, bimFs, bimRs, ReadSummary,
       file = file.path(DataFolder, "DenoisedData.RData"))
  
  message("*********************** Sequence table generated, Data saved ***********************
          ********************************************************************")
  cat("\n*** Sequence table generated, Data saved ***", file = LogFile, append = TRUE)
  cat("\n*** Start generating Summary Plots ***", file = LogFile, append = TRUE)
  
  ##############################
  ### Summary Plots
  ##############################
  
  # the width of the plots probably needs adjustment
  width = 5 + 0.5*(length(SampleNames)/10)
  
  # Plot the number of reads at the different steps for each sample (see also ReadSummary)
  pdf(file = file.path(PlotFolder, "NoReads_AllSamples.pdf"), width = width, height = 6)
  print(NoReads_StepsSimple(ReadSummary = ReadSummary, SampleNames = SampleNames, sort = TRUE))
  dev.off()
  
  # PLot the total number of amplicons against the number of unique amplicons for each sample
  FinalNumbers <- dplyr::select(ReadSummary, Sample, UniqueAmplicons = UniqueAmpliconsWOBimera, NoAmplicons = MergedReadsWOBimera)
  # FinalNumbers <- data.frame(Sample = rownames(seqtab.nochim), UniqueAmplicons = rowSums(seqtab.nochim != 0), NoAmplicons = rowSums(seqtab.nochim))
  
  Tr <- TotalandUniqueAmplicons(FinalNumbers = FinalNumbers, seqtab = seqtab.nochim)
  pdf(file = file.path(PlotFolder, "TotalvsUniqueAmplicons.pdf"), width = width, height = 6)
  print(Tr)
  dev.off()
  
  # Plot a linear regression line to illustrate the possible association between the total number of amplicons and the number of 
  # unique amplicons
  Tr2 <- ggplot(FinalNumbers, aes(x = NoAmplicons, y = UniqueAmplicons)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    xlab("Total Amplicons") +
    ylab("Unique Amplicons") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = "#999999", size = .15),
          panel.grid.major.x = element_line(color = "#999999", size = .15))
  
  
  pdf(file = file.path(PlotFolder, "AssociationTotaltoUniqueAmplicons.pdf"), width = 7, height = 6)
  print(Tr2)
  dev.off()
  
  # Plot a histogram illustrating in how many samples the amplicons are present
  if(length(SampleNames) < 10) {binwidth = 1}
  if(length(SampleNames) > 10 && length(SampleNames) < 100) {binwidth = 2}
  if(length(SampleNames) > 100) {binwidth = 3}
  
  FinalNumbers2 <- data.frame(InNoSamples = colSums(seqtab.nochim != 0))
  
  Trr <- ggplot(data = FinalNumbers2, aes(x = InNoSamples))  + 
    geom_histogram(binwidth = binwidth, col = "black", fill = "#E69F00")+
    geom_rug() +
    xlab("Present in No Samples") + 
    ylab("Count") +
    theme_bw() + 
    ggtitle(paste("Total No of unique amplicons:", dim(seqtab.nochim)[2], "Only in 1 sample:", sum(FinalNumbers2 ==1))) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(color = "#999999", size = .15)) +
    coord_cartesian(ylim = c(0,150))
  
  pdf(file = file.path(PlotFolder, "HistogramAmpliconsinNoSamples.pdf"), width = 7, height = 6)
  print(Trr)
  dev.off()
  cat("\n*** Summary Plots generated, Function done ***", file = LogFile, append = TRUE)
}
