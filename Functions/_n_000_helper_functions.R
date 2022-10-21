# - define cbPalette and QuantColors15 schemes -
# R color blind palette is used for group comparisons (unless more than 8 groups)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # ref: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# QuantColors15 is used for phyla comparisons as long as there are < 15 phyla, each color is easy to distinguish from the two colors surrounding it
tol15rainbow=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA") # ref: https://tradeblotter.wordpress.com/2013/02/28/the-paul-tol-21-color-salute/
QuantColors15 <- tol15rainbow[c(1, 12, 4, 7, 13, 2, 11, 5, 8, 14, 3, 10, 6, 9, 15)]

pal <- function(col, border = "light gray", ...){
        n <- length(col)
        plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
             axes = FALSE, xlab = "", ylab = "", ...)
        rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}
# --


# --
###########################
### FUNCTION: areColors ###
############################
# Found here: #https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
# is used to test color_levels input, i.e. tests if each entry in color_levels is a character
# representation of a color

areColors <- function(x) {
        sapply(x, function(X) {
                tryCatch(is.matrix(grDevices::col2rgb(X)), 
                         error = function(e) FALSE)
        })
}
# --



# --
#######################################
### FUNCTION: assign_default_colors
#######################################
# Function that takes a df and a variable_name as input, and defines colors to the factor levels or unique entries based on cbPalette, QuantColors15 or viridis

assign_default_colors <- function(info_df, variable_name){
        
        if (!(variable_name %in% colnames(info_df))) {
                stop("variable_name must be a variable/column in info_df")   
        }
        
        if (is.factor(info_df[[variable_name]])) {
                uniqueEntries <- levels(info_df[[variable_name]])
        } else {
                uniqueEntries <- unique(info_df[[variable_name]])
        }
        
        
        if (length(uniqueEntries) < 9) {
                color_char <- cbPalette[1:length(uniqueEntries)]
        } else if (length(uniqueEntries) < 16) {
                color_char <- QuantColors15[1:length(uniqueEntries)]
        } else {
                color_char <- viridis(length(uniqueEntries))
        }
        
        names(color_char) <- uniqueEntries
        
        color_char
}
# --



# --
#######################################
### FUNCTION: make_color_vector
#######################################
# give a character vector or factor as input together with a color Palette
# outputs a named color vector, using the factor levels or the unique(character vector)
# The change to "NA" is to have NAs in the order I want them to be, if using na.value in ggplot it is always at the last position


make_color_vector <- function(in_vector, col_pal){
        
        if (is.factor(in_vector)) {
                vec_names <- levels(in_vector)
        } else {
                vec_names <- unique(as.character(in_vector))
        }
        
        vec_names[is.na(vec_names)] <- "NA"
        
        
        if (length(col_pal) < length(vec_names)) {
                stop("not enough colors in your color palette")   
        }
        
        col_vector <- col_pal[1:length(vec_names)]
        names(col_vector) <- vec_names
        col_vector
}
# --



# --
#######################################
### lmp (to get p_value from lm fit)
#######################################
# needed to get the p-value from-linear fit objects (from stackoverflow)

lmp <- function (modelobject) {
        if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
        f <- summary(modelobject)$fstatistic
        p <- pf(f[1],f[2],f[3],lower.tail=F)
        attributes(p) <- NULL
        return(p)
}
# --


# --
#######################################
### gm_own: calculate geometric mean
#######################################
# see commented below, this function comes from 
# <http://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in>

## Input:
# x numeric vector
# na.rm: if FALSE you get NA as soon as an NA is in your data, if TRUE the NA get basically treated as 0 (but NOTE these zeros always count also when zero.count = FALSE)
# zeros.count: This is IMPORTANT, if TRUE 0s count, so if x contains a 0 the GM will be lower than when zero.count = FALSE, in 
# which case the gm is calculated for the same vector in which all 0 have been removed.
## Output:
# the geometric mean of x
gm_own = function(x, na.rm=FALSE, zeros.count = TRUE){
        if(any(x < 0, na.rm = TRUE)){
                return(NaN)
        }
        if(zeros.count){
                exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
        } else {
                exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x[x > 0]))
        }
}
# --




# --
#######################################
### get_unique_facLevel_combinations
#######################################
# often needed to set the comparisons attribute in stat_compare_means from ggplot, it excludes comparison with each other and also removes 2 to 3 vs 3 to 2

get_unique_facLevel_combinations <- function(fac_levels){
        
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(jvec), function(x){
                        i <- jvec[x]
                })
        })
        j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        j <- ivec[x]
                })
        })
        i_s <- i_s[lower.tri(i_s)]
        j_s <- j_s[lower.tri(j_s)]
        
        comparisonList <- list()
        for (k in seq_along(i_s)){
                comparisonList[[k]] <- c(fac_levels[i_s[k]], fac_levels[j_s[k]])
        }
        
        comparisonList
        
}
# --




# --
#######################################
### get_unique_facLevel_combis
#######################################
# as get_unique_facLevel_combinations but outputting i_s and j_s

get_unique_facLevel_combis <- function(fac_levels){
        
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(jvec), function(x){
                        i <- jvec[x]
                })
        })
        j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        j <- ivec[x]
                })
        })
        i_s <- i_s[lower.tri(i_s)]
        j_s <- j_s[lower.tri(j_s)]
        
        list(i_s = i_s, j_s = j_s)
        
}
# --




