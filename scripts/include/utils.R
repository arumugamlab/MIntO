# '''
# Commonly used utility functions
#
# Authors: Mani Arumugam
#
# '''

logmsg <- function(...) {
    message("# ", format(Sys.time(), digits=0), " - ", paste(list(...)), collapse=" ")
}
