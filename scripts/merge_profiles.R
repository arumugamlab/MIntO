#!/usr/bin/env Rscript

library(optparse)
parser = OptionParser()
parser = add_option(parser, c("-t", "--threads"), type="integer", default=4, help="number of threads [default: %default]")
parser = add_option(parser, c("-m", "--memory"), type="integer", default=10, help="maximum memory to be used")
parser = add_option(parser, c("-o", "--out"), type="character", default=NULL, help="output file")
parser = add_option(parser, c("-i", "--input"), type="character", default=NULL, help="files to be combined; as one comma-delimited string")
parser = add_option(parser, c("-k", "--keys"), type="character", default=NULL, help="columns shared across files and must be used for join; as one comma-delimited string")
parser = add_option(parser, c("-z", "--zeroes"), action="store_true", default=FALSE, help="flag whether there are zeroes in the profile output that should be removed")

opt = parse_args(parser)

# Check required arguments

if (is.null(opt$out)) {
    print_help(parser)
    stop("Missing required option '--out'")
}
if (is.null(opt$input)) {
    print_help(parser)
    stop("Missing required option '--input'")
}
if (is.null(opt$keys)) {
    print_help(parser)
    stop("Missing required option '--keys'")
}

# Store arguments

threads_n = opt$threads
memory_lim = opt$memory
out_file = opt$out
in_files = as.list(strsplit(opt$input, ",", fixed=TRUE))[[1]]
keys = as.vector(strsplit(opt$keys, ",", fixed=TRUE))[[1]]
zeroes = opt$zeroes

# Check input files

for (f in in_files) {
    if(!file.exists(f)) {
        stop(sprintf("File '%s' does not exist", f))
    }
    if(file.access(f, 4) == -1) {
        stop(sprintf("File '%s' cannot be read", f))
    }
}

# Check output file

if (file.exists(out_file) && file.access(out_file, 2) == -1) {
    stop(sprintf("Output file '%s' exists but cannot be overwritten", out_file))
}

# Command line check done

# Load libraries

library(R.utils)
library(data.table)
library(parallel)

setDTthreads(threads = threads_n)

# Log progress

logmsg <- function(...) {
    message("# ", format(Sys.time(), digits=0), " - ", paste(list(...)), collapse=" ")
}

# Read contents of a file into data.table, remove zeroes if necessary, index by given keys, and return it
get_dt_from_file <- function(filename, index_keys=c('ID'), remove_zeroes) {
    # Load libraries, for parallel
    library(R.utils)
    library(data.table)

    logmsg <- function(...) {
        message("# ", format(Sys.time(), digits=0), " - ", paste(list(...)), collapse=" ")
    }

    # Read file
    logmsg("FILE: ", filename)
    logmsg("  Reading")
    dt = fread(filename, skip='ID\t', header=T, data.table=TRUE, sep='\t')

    if (remove_zeroes == TRUE) {
        # Index by sample data
        logmsg("  Indexing table by data")
        sample_col = setdiff(colnames(dt), index_keys)
        setkeyv(dt, sample_col)

        # Remove 0-values
        logmsg("  Removing 0's ")
        logmsg("    before = ", nrow(dt))
        dt = dt[!J(0),]
        logmsg("     after = ", nrow(dt))
    }

    # Index by index_keys
    logmsg("  Indexing table by ", index_keys)
    setkeyv(dt, index_keys)

    return(dt)
}

# ID gets first key status.
# Others are sorted and sent as secondary keys.

secondary_keys = sort(setdiff(keys, c('ID')))
merged_dt = NULL

# Get list of DT from list of filenames
# Use get_dt_from_file() as the workhorse
logmsg("Reading ", length(in_files), " files using ", threads_n, " threads")
if (threads_n > 1) {
    cluster <- makeCluster(threads_n)
    dt_list = parLapply(cluster, in_files, get_dt_from_file, index_keys=c('ID', secondary_keys), remove_zeroes=zeroes)
    stopCluster(cluster)
} else {
    dt_list = lapply(in_files, get_dt_from_file, index_keys=c('ID', secondary_keys), remove_zeroes=zeroes)
}
logmsg("  done")

# Merge dt_list into a single dt
logmsg("Merging")
merged_dt = Reduce(function(...) merge(..., all=TRUE, sort=TRUE), dt_list)
logmsg("  done")

# Write output
logmsg("Writing merged output")
fwrite(merged_dt, file = out_file, row.names = F, col.names = T, sep = "\t", quote = F, na = 0)
logmsg("  done")

gc()
