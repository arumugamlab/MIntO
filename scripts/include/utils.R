# '''
# Commonly used utility functions
#
# Authors: Mani Arumugam
#
# '''

logmsg <- function(...) {
    message("# ", format(Sys.time(), digits=0), " - ", paste(list(...)), collapse=" ")
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
