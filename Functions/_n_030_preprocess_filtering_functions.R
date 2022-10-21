# --
#######################################
#### visualize_filtering
#######################################


visualize_filtering <- function(physeq, prevalence, taxa_sums_quantile, phylum_colors = NULL){ 
        
        if (taxa_are_rows(physeq)) {
                physeq <- t(physeq)
        }
        
        df_ab_prev <- data.frame(Taxon_No = 1:ntaxa(physeq), 
                                 total_counts = taxa_sums(physeq),
                                 prevalence = colSums(as(otu_table(physeq), "matrix") != 0),
                                 sparsity = colSums(as(otu_table(physeq), "matrix") == 0),
                                 mean_count_nonzero = apply(as(otu_table(physeq), "matrix"), 2, function(x){mean(x[x > 0])}),
                                 median_count_nonzero = apply(as(otu_table(physeq), "matrix"), 2, function(x){median(x[x > 0])})
        )
        
        df_ab_prev <- cbind(df_ab_prev, tax_table(physeq))
        
        # - adjust color and order of the phyla in the following plots - 
        colori <- "Phylum"
        
        if (!is.null(phylum_colors)){
                df_ab_prev[[colori]] <- as.character(df_ab_prev[[colori]])
                df_ab_prev[[colori]][is.na(df_ab_prev[[colori]])] <- "NA" # NB: pools NAs
                if (!all(unique(df_ab_prev[[colori]]) %in% names(phylum_colors))){
                        stop("provided phylum_colors did not cover all Phyla in physeq")
                }
                df_ab_prev[[colori]] <- factor(df_ab_prev[[colori]], levels = names(phylum_colors), ordered = TRUE)
                custom_colors <- phylum_colors
        } else {
                CountOrder <- dplyr::group_by_(df_ab_prev, colori) %>% dplyr::summarise(total_count_sum = sum(total_counts)) %>% dplyr::arrange(desc(total_count_sum))
                
                CountOrder[[colori]] <- as.character(CountOrder[[colori]])
                CountOrder[[colori]][is.na(CountOrder[[colori]])] <- "NA"
                
                if (nrow(CountOrder) <= 15){
                        custom_colors <- make_color_vector(CountOrder[[colori]], QuantColors15)
                } else {
                        custom_colors <- make_color_vector(CountOrder[[colori]], viridis(nrow(CountOrder)))
                }
                
                df_ab_prev[[colori]] <- as.character(df_ab_prev[[colori]])
                df_ab_prev[[colori]][is.na(df_ab_prev[[colori]])] <- "NA"
                df_ab_prev[[colori]] <- factor(df_ab_prev[[colori]], levels = names(custom_colors), ordered = TRUE)
                
        }
        # --

        prev_thresh <- (prevalence/100)*nsamples(physeq)
        abund_thresh <- quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100)
        
        df_ab_prev_filt <- dplyr::filter(df_ab_prev, prevalence > prev_thresh | total_counts > abund_thresh)
        
        no_samples <- nsamples(physeq)
        shade_df <- data.frame(total_counts = 0, prevalence = 0)
        
        Tr_prev_vs_log10_ab <- ggplot(df_ab_prev, aes(x = total_counts, y = prevalence))
        Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab +
                geom_point(aes_string(col = colori), size = 2, alpha = 0.7) +
                scale_x_log10() +
                geom_rect(data = shade_df, xmin = -Inf, xmax = log10(abund_thresh), ymin = -Inf, ymax = prev_thresh, fill = "#660033", alpha = 0.4) +
                geom_hline(yintercept = (prevalence/100)*nsamples(physeq), col = cbPalette[1], lty = "dashed") +
                geom_vline(xintercept = quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100), col = cbPalette[1], lty = "dashed") +
                xlab("total counts (taxa_sums())") +
                scale_color_manual("", values = custom_colors) +
                theme_bw() +
                ggtitle(paste(nrow(df_ab_prev_filt), " of ", nrow(df_ab_prev), " taxa (", round(100*nrow(df_ab_prev_filt)/nrow(df_ab_prev), 1),
                              " %) and ", round(sum(df_ab_prev_filt$total_counts)), " of ", round(sum(df_ab_prev$total_counts)), " counts (",
                              round((sum(df_ab_prev_filt$total_counts)/sum(df_ab_prev$total_counts))*100, 1), " %) remain", sep = ""))
        
        
        
        Tr_prev_vs_log10_ab_wrap <- ggplot(df_ab_prev, aes(x = total_counts, y = prevalence))
        Tr_prev_vs_log10_ab_wrap <- Tr_prev_vs_log10_ab_wrap +
                geom_point(size = 2, alpha = 0.7) +
                scale_x_log10() +
                geom_rect(data = shade_df, xmin = -Inf, xmax = log10(abund_thresh), ymin = -Inf, ymax = prev_thresh, fill = "#660033", alpha = 0.4) +
                geom_hline(yintercept = (prevalence/100)*nsamples(physeq), col = cbPalette[1], lty = "dashed") +
                geom_vline(xintercept = quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100), col = cbPalette[1], lty = "dashed") +
                xlab("total counts (taxa_sums())") +
                facet_wrap(~ Phylum) +
                scale_color_manual("", values = custom_colors) +
                theme_bw() +
                ggtitle(paste(nrow(df_ab_prev_filt), " of ", nrow(df_ab_prev), " taxa (", round(100*nrow(df_ab_prev_filt)/nrow(df_ab_prev), 1),
                              " %) and ", round(sum(df_ab_prev_filt$total_counts)), " of ", round(sum(df_ab_prev$total_counts)), " counts (",
                              round((sum(df_ab_prev_filt$total_counts)/sum(df_ab_prev$total_counts))*100, 1), " %) remain", sep = ""))
        
        
        
        Tr_prev_vs_log10_ab_wrap <- Tr_prev_vs_log10_ab +
                facet_wrap(~ Phylum) +
                theme(legend.position = "none")
        
        out <- list(Tr_prev_vs_log10_ab = Tr_prev_vs_log10_ab,
                    Tr_prev_vs_log10_ab_wrap = Tr_prev_vs_log10_ab_wrap)
        
}
# --











