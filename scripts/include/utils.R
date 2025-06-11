# '''
# Commonly used utility functions
#
# Authors: Mani Arumugam
#
# '''

logmsg <- function(..., decorate=FALSE) {
    msg_str  = paste(list(...), collapse="")
    time_str = format(Sys.time(), digits=0)
    if (decorate) {
        border_len = nchar(msg_str) + 10  # extra for padding and box
        top        = paste0("┌", strrep("─", border_len - 2), "┐")
        decorated  = paste0("│    ", msg_str, "    │")
        bottom     = paste0("└", strrep("─", border_len - 2), "┘")
        message("# ", time_str, " - ", top)
        message("# ", time_str, " - ", decorated)
        message("# ", time_str, " - ", bottom)
    } else {
        message("# ", time_str, " - ", msg_str)
    }
}

# Return the mem usage (in bytes) of a list of variables in a given environment (calling env by default)
# Combine into a data frame for better readability
report_mem_usage <- function(variables, env = parent.frame()) {
    library(pryr)

    # Is this variable of interest?
    var_of_interest <- function(x, env) {
        voi = get(x, envir = env)
        !(is.function(voi) || (is.atomic(voi) && length(voi)<=1))
    }

    # Get variables of interest: skip functions and scalars
    v = variables[sapply(variables, var_of_interest, env)]

    # Get a df with their sizes
    df = data.frame(VARIABLE = v, SIZE = sapply(v, function(x) pryr::object_size(get(x, env))))
    rownames(df) = NULL

    # Pad with hlines for readability
    return(rbind(data.frame(VARIABLE = "--------", SIZE = "--------"), df, data.frame(VARIABLE = c("--------", "TOTAL", "--------"), SIZE = c("--------", sum(df$SIZE), "--------"))))
}

# Return an md5hash of all columns in the data.table
md5hash <- function(dt) {
    library(openssl)
    column_hashes = dt[, lapply(.SD, function(col) openssl::md5(as.character(col)))]
    column_hashes = sapply(names(column_hashes), function(col) openssl::md5(paste0(column_hashes[[col]], collapse = "")))
    return(column_hashes)
}
