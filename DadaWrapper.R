 # ---- Sourcing the plot and wrapper functions ----
 working_dir <- ""
 path_funct <- paste0(working_dir,"Clean_Data/")
 
# CHANGE pathToFunctions here:
pathToFunctions <- "Functions/"

source(file.path(pathToFunctions, "Dada_PlotFunctions.R"))
source(file.path(pathToFunctions, "Dada_WrapFunctions.R"))
# ----

# ---- Calling the wrap function (Adjust INPUTS) ----
Dada2_wrap(path = paste0(path_funct),
           F_pattern = "S*/*_1.fq.gz", 
           R_pattern = "S*/*_2.fq.gz",
           path2 = NULL, 
           trimLeft = 10, # how many nucleotides will be clipped from the beginning of the FW and RV reads, respectively
           truncLen = 0, # how many nucleotides will be clipped from the end of the FW and RV reads, respectively
           maxEE = 0.5, # After truncation, reads with higher than maxEE "expected errors" will be discarded
           maxN = 0, # After truncation, sequences with more than maxN Ns will be discarded, NB: Dada2 requires no Ns!
           truncQ = 2,# Truncate reads at the first instance of a quality score less than or equal to truncQ 
           # NB: see: https://github.com/benjjneb/dada2/issues/140 
           # ? I did not understand this because I thought it should clash with dada2 not allowing sequences of variable length, but this i now supported:
           # https://github.com/benjjneb/dada2/issues/55 (at the end)
           nreadsLearn = 1.2e+06, # the number of reads (distributed over far less unique reads) used to learn the error matrixes, i.e. nreads in dada2:::learnErrors
           err_F = NULL, # when error matrix given, the error matrix estimation is skipped
           err_R = NULL,
           minOverlap = 20, # minOverlap from the mergePairs command
           maxMismatch = 0, # maxMismatch from the mergePairs command
           F_QualityStats = NULL, # if given e.g. from Dada_QualityCheck the quality stats collection part is jumped over
           R_QualityStats = NULL,
           filtFs = NULL, # if given and FilteredFolder with files exist, Filtering can be jumped over
           filtRs = NULL,
           pool = FALSE) # look for a ASV in all the samples

# Then call on terminal: Rscript DadaWrapper.R
