# '''
# Functions for plotting PCA
#
# Authors: Mani Arumugam
#
# '''

# Requirements:
#  'profile' should only have numeric columns - IDs should have been removed
#  'profile' should have removed all zerosum rows
#  'profile' should be sorted (desc) by rowSum
plot_PCA <- function(profile, label, color, metadata) {

    library(data.table)
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
    logmsg("  Replacement for 0 : ", min_value*1e-2)
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
