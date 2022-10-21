# ---- Source the functions ----

# ATTENTION: change pathToFunctions here if necessary#
working_dir <- ""
path_funct <- paste0(working_dir,"Clean_Data/")
path2_funct <- paste0(working_dir,"QualityCheck/")

setwd("Cdiff_toxins_github")
pathToFunctions <- "Functions/"


source(file.path(pathToFunctions, "Dada_PlotFunctions.R"))
source(file.path(pathToFunctions, "Dada_WrapFunctions.R"))

# ----
# ########################## SAMPLES ##########################
dir.create(file.path(path2_funct, showWarnings = FALSE))

# ---- Call the wrap function (Adjust INPUTS) ----
Dada2_QualityCheck(path = path_funct, 
           F_pattern = "S*/*_1.fq.gz", 
           R_pattern = "S*/*_2.fq.gz",
           path2 = path2_funct)
# ----

# Then call on terminal Rscript Dada_QualityCheck.R
