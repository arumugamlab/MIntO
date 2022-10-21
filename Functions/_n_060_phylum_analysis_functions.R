# --
#######################################
#### plot_sample_bars
#######################################
# based on plot_bar from phyloseq, difference, orders and colors the Samples based on group_var, and orders the fill based on abundance
# I guess inputs can be guessed on

plot_sample_bars <- function(physeq, x = "Sample", y = "Abundance", group_var, color_levels, fill = NULL,
                             color_sample_names = TRUE, col_vec = NULL, facet_grid = NULL, order_by_firmicutes = TRUE){
        
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the phyloseq object.")
        }
        
        if (!all(names(color_levels) %in% unique(sample_data(physeq)[[group_var]]))) {
                stop("Not all names in names(color_levels)are actually levels in the group_var column.")
        }
        
        
        # in case you do not want to see all samples
        if (!all(unique(sample_data(physeq)[[group_var]]) %in% names(color_levels))) {
                keepSamples <- sample_names(physeq)[sample_data(physeq)[[group_var]] %in% names(color_levels)]
                physeq <- prune_samples(samples = keepSamples, physeq)
        }
        sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], levels = names(color_levels), order = TRUE)
        
        # if (!is.factor(sample_data(physeq)[[group_var]])) {sample_data(physeq)[[group_var]] <- as.factor(sample_data(physeq)[[group_var]])}
        
        if (is.null(fill)) { fill = "Phylum"}
        
        mdf <- phyloseq::psmelt(physeq)
        
        
        # order fill levels according to abundance over all samples
        mdf[, fill] <- as.character(mdf[, fill])
        mdf[is.na(mdf[, fill]), fill] <- "NA" # NB: pools all AN
        sums <- group_by_(mdf, fill) %>% summarise(sum_abundance = sum(Abundance)) %>% arrange(sum_abundance)
        mdf[, fill] <- factor(mdf[, fill], levels = as.character(sums[[1]]), ordered = TRUE)
        
        # order samples according to levels or Firmicutes
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        
        if (order_by_firmicutes) {
                mdf_firmicutes <- dplyr::filter(mdf, Phylum == "Firmicutes") %>% arrange_(group_var, "Abundance")
                mdf$Sample <- factor(mdf$Sample, levels = mdf_firmicutes$Sample, ordered = TRUE)
        } else {
                mdf$Sample <- factor(mdf$Sample, levels = LookUpDF$Sample, ordered = TRUE)
        }
        
        
        
        # - define names of x axis using color_levels (which must be a named character vector) -        
        colxaxis <- color_levels[LookUpDF$Group]
        # --
        
        if (is.null(col_vec)){
                if (length(levels(mdf[, fill])) <= 15) {
                        fill_colors <- make_color_vector(mdf[, fill], rev(QuantColors15[1:length(levels(mdf[, fill]))]))
                } else {
                        fill_colors <- make_color_vector(mdf[, fill], viridis(length(levels(mdf[, fill]))))
                }
                
        } else {
                fill_colors <- col_vec
        }
        
        
        
        Tr <- ggplot(mdf, aes_string(x = x, y = y, fill = fill))
        Tr <- Tr + 
                geom_bar(stat = "identity", position = "stack") +
                theme_bw() +
                scale_fill_manual(values = fill_colors) +
                xlab("") +
                theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, colour = colxaxis))
        
        if (!is.null(facet_grid)) {
                formulation <- as.formula(paste("~ ", facet_grid, sep = ""))
                Tr <- Tr + facet_grid(formulation)
        }
        
        
        Tr
}
# --






# --
####################################
## plot_taxa_ratios_AllLevels 
###################################
# see plot_taxa_ratios_levelPairs: Here you directly calculate the count by count ratio matrix only for the tax_nom, you still facet by taxa_den
# (denominator) but you keep all levels in all plots. Makes extensive use of ggpubr, NB: ggpubr is so smart to adjust p-values when you use
# scale_y_log10

## Output: 
# - list(pVals = pVals, Tr = Tr, Tr1 = Tr1, pValsLog = pValsLog, Tr2 = Tr2, Tr3 = Tr3)


plot_taxa_ratios_AllLevels <- function(physeq, group_var, color_levels, tax_names = NULL,
                                   taxa_nom = "Firmicutes", taxa_den = NULL, test = "t.test", p_adjust_method = "fdr",
                                   tax_order = NULL,
                                   symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), hide.ns = FALSE) {
        
        
        if (taxa_are_rows(physeq)) { 
                physeq <- t(physeq)
        }
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the loaded phyloseq object.")
        }
        
        if (!all(names(color_levels) %in% unique(sample_data(physeq)[[group_var]]))) {
                stop("Not all names in names(color_levels)are actually levels in the group_var column.")
        }
        
        if (!(test %in% c("t.test", "wilcox.test"))) {
                stop("test should be t.test or wilcox.test")
        }
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        compare <- names(color_levels)
        
        if (!is.null(compare)) {
                group_var_levels <- compare
        } else {
                group_var_levels <- levels(group_fac)
        }
        
        # if (length(group_var_levels) != 2) {
        #         stop(paste0("compare (group_var_levels) must consist of two groups - you asked for ", 
        #                     paste(group_var_levels, collapse = ", ")))
        # }
        
        
        # - check that given tax_names fit to physeq and change taxa_names of physeq -
        if (is.null(tax_names)){
                tax_names <- paste("T", 1:ntaxa(physeq), sep = "_")
        } 
        
        if(!identical(ntaxa(physeq), length(tax_names))){stop("tax_names do not fit in length to physeq")}
        
        tax_names <- make.unique(tax_names)
        taxa_names(physeq) <- tax_names
        # --
        
        # - calculate the matrix taxa_nom/(all other taxa) -         
        CT <- t(as(otu_table(physeq), 'matrix')) # now taxa are rows and samples are columns
        # ONLY keep samples defined by group_var_levels (= names(color_levels))
        CT <- CT[, group_fac %in% group_var_levels]
        
        
        i <- which(rownames(CT) == taxa_nom)
        if (length(i) != 1) {stop("taxa_nom not found in tax_names or tax_names not unique!")}
        
        
        TbTmatrix <- apply(CT, 2, function(samp_cnts){samp_cnts[i]/samp_cnts})
        # produces for each taxon (= host taxon) a TbTMatrix
        # NB: there are possibly Inf, and NaN values in the matrix, specifically
        # 0/x = 0, x/0 = Inf; 0/0 = NaN!
        # --
        
        TbT_DF <- as.data.frame(TbTmatrix)
        TbT_DF$Taxon <- rownames(TbT_DF)
        # - use taxa_den (denominator) to restrict the taxa to which taxa_nom is compared to -
        if (is.null(taxa_den)) {taxa_den <- tax_names}
        TbT_DF <- TbT_DF[TbT_DF$Taxon %in% taxa_den, ]
        # --
        
        # - change to long DF -
        TbT_DF_l <- gather(TbT_DF, key = Sample, value = Ratio, -Taxon)
        # --
        
        # - add the group_var level information  -
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        TbT_DF_l$Group <- as.character(LookUpDF$Group[match(TbT_DF_l$Sample, LookUpDF$Sample)])
        TbT_DF_l$Group <- factor(TbT_DF_l$Group, levels = group_var_levels, ordered = T)
        # --
        
        # - change all ratios where either nominator taxon or denominator taxon had count = 0 to NA -
        # remember: 0/0 = NaN (not finite), 0/x = 0, x/0 = Inf
        TbT_DF_l$Ratio[!is.finite(TbT_DF_l$Ratio) | TbT_DF_l$Ratio == 0] <- NA
        # --
        
        # - find taxa that would throw an error in statistical test and remove those taxa from DF -
        # first find the taxa that would throw an error in t.test or wilcox.test
        var_plus_length_check <- group_by(TbT_DF_l, Taxon, Group) %>% summarise(Variance = var(Ratio, na.rm = T), NotNA = sum(!is.na(Ratio)))
        if (test == "t.test"){
                var_plus_length_check <- dplyr::filter(var_plus_length_check, !(Variance > 0) | NotNA < 2) # variance > 0 also to remove test where host_taxon == taxon
        } else if (test == "wilcox.test") {
                var_plus_length_check <- dplyr::filter(var_plus_length_check, !(Variance > 0) | NotNA < 1)
        }
        
        if (nrow(var_plus_length_check) != 0) { 
                TbT_DF_l <- filter(TbT_DF_l, !(Taxon %in% unique(var_plus_length_check$Taxon)))
        }
        # --
        
        # - use ggpubr::compare_menas to calculate all pValues of the Ratios for the different taxa_den between current group_levels -
        pVals <- ggpubr::compare_means(formula = Ratio ~ Group, data = TbT_DF_l, group.by = "Taxon", method = test, p.adjust.method = p_adjust_method, symnum.args = symnum.args)
        
        pVals <- dplyr::arrange(pVals, p)
        
        # NB: the pVals for t.test change when you log the ratios (scale_y_log10())
        TbT_DF_l$RatioLog10 <- log10(TbT_DF_l$Ratio)
        pValsLog <- ggpubr::compare_means(formula = RatioLog10 ~ Group, data = TbT_DF_l, group.by = "Taxon", method = test, p.adjust.method = p_adjust_method, symnum.args = symnum.args)
        pValsLog <- dplyr::arrange(pValsLog, p)
        # --
        
        # NB: I plot now for log and non log independently, even though ggpubr is so smart to change the p-values when you just use
        # Tr + scale_y_log10(). I plot independently because log might change the order in case of t.test!
        TbT_DF_l_log <- TbT_DF_l
        
        # - order taxa based on pVals result or based on tax_order -
        if (is.null(tax_order)){
                TbT_DF_l$Taxon <- factor(TbT_DF_l$Taxon, levels = unique(pVals$Taxon), ordered = TRUE)
                TbT_DF_l_log$Taxon <- factor(TbT_DF_l_log$Taxon, levels = unique(pValsLog$Taxon), ordered = TRUE)
                
        } else {
                if(!all(unique(pVals$Taxon) %in% tax_order)){
                        stop("given tax_order does not fit to tax_names")
                }
                TbT_DF_l$Taxon <- factor(TbT_DF_l$Taxon, levels = tax_order, ordered = TRUE)
                TbT_DF_l_log$Taxon <- factor(TbT_DF_l_log$Taxon, levels = tax_order, ordered = TRUE)
        }
        # --
        
        # - since you might have more than two levels in each plot you need to set the comparisons argument in stat_compare_means -
        comparisonList <- get_unique_facLevel_combinations(group_var_levels)
        # --
        
        # - plot: NB: also for non removed taxa some samples might have NA ratios that will be removed -
        
        Tr <- ggplot(TbT_DF_l, aes(x = Group, y = Ratio, col = Group))
        Tr <- Tr +
                geom_violin() +
                geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                scale_color_manual(values = color_levels) +
                facet_wrap(~ Taxon, scales = "free_y") +
                xlab("") +
                ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                theme_bw() +
                theme(legend.position = "none")
        
        Tr <- Tr + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = test, hide.ns = hide.ns)
        # Tr <- Tr + ggpubr::stat_compare_means(label = "p.signif", method = test, label.x = 1.5, hide.ns = hide.ns)
        
        Tr1 <- ggplot(TbT_DF_l, aes(x = Group, y = Ratio, col = Group))
        Tr1 <- Tr1 +
                geom_boxplot(outlier.color = NA) +
                geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                scale_color_manual(values = color_levels) +
                facet_wrap(~ Taxon, scales = "free_y") +
                xlab("") +
                ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                theme_bw() +
                theme(legend.position = "none")
        
        Tr1 <- Tr1 + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = test, hide.ns = hide.ns)
        #Tr1 <- Tr1 + ggpubr::stat_compare_means(label.x = 1.5, label = "p.signif", method = test, hide.ns = hide.ns) # p.format
        
        
        Tr2 <- ggplot(TbT_DF_l_log, aes(x = Group, y = Ratio, col = Group))
        Tr2 <- Tr2 +
                geom_violin() +
                geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                scale_color_manual(values = color_levels) +
                facet_wrap(~ Taxon, scales = "free_y") +
                xlab("") +
                ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                theme_bw() +
                theme(legend.position = "none")
        
        
        
        # Tr2 <- Tr2 + ggpubr::stat_compare_means(label.x = 1.5, label = "p.signif", method = test, hide.ns = hide.ns)
        Tr2 <- Tr2 + scale_y_log10()
        Tr2 <- Tr2 + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = test, hide.ns = hide.ns)
        
        
        Tr3 <- ggplot(TbT_DF_l_log, aes(x = Group, y = Ratio, col = Group))
        Tr3 <- Tr3 +
                geom_boxplot(outlier.color = NA) +
                geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                scale_color_manual(values = color_levels) +
                facet_wrap(~ Taxon, scales = "free_y") +
                xlab("") +
                ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                theme_bw() +
                theme(legend.position = "none")
        
        # Tr3 <- Tr3 + ggpubr::stat_compare_means(label.x = 1.5, label = "p.signif", method = test, hide.ns = hide.ns)
        Tr3 <- Tr3 + scale_y_log10()
        Tr3 <- Tr3 + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = test, hide.ns = hide.ns)
        
        
        list(pVals = pVals, Tr = Tr, Tr1 = Tr1, pValsLog = pValsLog, Tr2 = Tr2, Tr3 = Tr3)
        
}
#--





# --
####################################
## calculate_raw_TbTmatrixes:
###################################

calculate_raw_TbTmatrixes = function(physeq){
        
        if (taxa_are_rows(physeq)) {physeq <- t(physeq)}
        
       
        CT <- t(as(otu_table(physeq), 'matrix')) # now taxa are rows and samples are columns
        
        TbTmatrixes <- lapply(1:nrow(CT), function(i){apply(CT, 2, function(samp_cnts){samp_cnts[i]/samp_cnts})})
        # produces for each taxon (= host taxon) a TbTMatrix
        # NB: there are Inf, and NaN values in the matrixes, specifically
        # 0/x = 0, x/0 = Inf; 0/0 = NaN!
        
        names(TbTmatrixes) <- rownames(TbTmatrixes[[1]])
        
        
        TbTmatrixes
        
}
# --








# --
####################################
## create_raw_TbT_TilePlot: 
###################################
# NB: In this version TbTmatrixes should have been calculated on all samples in physeq. 
# names(color_levels) should only contain two levels in group_var!

create_raw_TbT_TilePlot <- function(TbTmatrixes, physeq, group_var, color_levels, tax_names = NULL, tax_order = NULL, 
                                           test = "wilcoxon", signi_level = 0.05, p_adjust_method = "none") {
        
        if(!identical(length(TbTmatrixes), ntaxa(physeq))){stop("TbTmatrixes don't fit to physeq")}
        
        if(ncol(TbTmatrixes[[1]]) != nsamples(physeq)){stop("TbTmatrixes don't fit to physeq. Not the same number of samples.")}
        
        
        if(test != "wilcoxon" & test != "t.test"){stop("test unknown, must be wilcoxon or t.test")}
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the loaded phyloseq object.")
        }
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        compare <- names(color_levels)
        
        if (!is.null(compare)) {
                group_var_levels <- compare
        } else {
                group_var_levels <- levels(group_fac)
        }
        
        if (length(group_var_levels) != 2) {
                stop(paste0("compare (names(color_levels)) must consist of two groups - you asked for ", 
                            paste(group_var_levels, collapse = ", ")))
        }
        
        if (!all(group_var_levels %in% levels(group_fac))) {
                stop("Not all given compare (group_var_levels) are actually levels in group_var column.")
        }
        
  
        
        # - check that given tax_names fit to physeq and change taxa_names of physeq -
        if (is.null(tax_names)){
                tax_names <- paste("T", 1:ntaxa(physeq), sep = "_")
        } 
        
        if(!identical(ntaxa(physeq), length(tax_names))){stop("tax_names do not fit in length to physeq")}
        
        tax_names <- make.unique(tax_names)
        # --
        
        
        
        names(TbTmatrixes) <- tax_names
        
        TbTmatrixes <- lapply(TbTmatrixes, function(mat){
                rownames(mat) <- tax_names
                mat
        })
        
        
        i <- group_var_levels[1]
        j <- group_var_levels[2]
        
        # ntaxa * ntaxa wilcoxon tests take time if you have a lot of taxa!
        pValMatrix <- sapply(TbTmatrixes, function(mat){
                apply(mat, 1, function(taxon_ratios){
                        x <- taxon_ratios[group_fac == i]
                        x <- x[is.finite(x) & x != 0] # removes all ratios in which one of the two taxa was not present!
                        y <- taxon_ratios[group_fac == j]
                        y <- y[is.finite(y) & y != 0] # removes 0/0 = NaN, 0/x = 0, x/0 = Inf
                        if (test == "wilcoxon"){
                                if (length(x) > 0 && length(y) > 0){
                                        pValue <- wilcox.test(x = x, y = y, alternative = "two", paired = F, exact = F)$p.value
                                        # NB: wilcox.test ignores in default setting (maybe see na.action) NA, NaN, Inf, -Inf
                                        # For plot: change sign of pValue to negative if taxon is more abundant in group 1. 
                                        Ranks <- rank(c(x[!is.na(x)], y[!is.na(y)]))
                                        n1 <- length(x[!is.na(x)])
                                        n2 <- length(y[!is.na(y)])
                                        Wx <- sum(Ranks[1:n1])-(n1*(n1+1)/2)
                                        Wy <- sum(Ranks[(n1+1):(n1+n2)])-(n2*(n2+1)/2)
                                        if(Wx > Wy){pValue <- -1*pValue}
                                        pValue
                                        
                                } else {
                                        pValue = 1
                                }
                                
                        } else if (test == "t.test") {
                                if (length(x) > 1 && length(y) > 1 && var(x) > 0 && var(y) > 0){
                                        pValue <- t.test(x = x, y = y, alternative = "two")$p.value
                                        if (mean(x, na.rm = T) > mean(y, na.rm = T)){pValue <- -1*pValue}
                                        pValue
                                } else {
                                        pValue <- 1
                                }
                                
                        }
                        
                })
        })
        
        # make sure diagonal is all NA (can be exceptions especially for t.test)
        diag(pValMatrix) <- NA
        
        # - adjust p-values if asked for -        
        signs <- pValMatrix < 0
        signs[is.na(signs)] <- FALSE
        
        pValMatrix <- abs(pValMatrix)
        for (e in 1:nrow(pValMatrix)){
                pValMatrix[e, ] <- p.adjust(pValMatrix[e, ], method = p_adjust_method)
        } # equal to t(apply(pValMatrix, 1, p.adjust, method = p_adjust))
        
        pValMatrix[signs] <- pValMatrix[signs]*(-1)
        # --
        
        # -- add a tile plot of the pValMatrix --
        
        DF <- as.data.frame(pValMatrix)
        DF[is.na(DF)] <- 2 # just to avoid missing values in plot and have a clear non-pValue value to mark self comparisons as black
        DF$HostTaxon <- rownames(pValMatrix)
        DF <- tidyr::gather(DF, key = Taxon , value = pValue, - HostTaxon)
        if (is.null(tax_order)) {
                DF$Taxon <- factor(DF$Taxon, levels = rownames(pValMatrix), ordered = TRUE)
                DF$HostTaxon <- factor(DF$HostTaxon, levels = rev(rownames(pValMatrix)), ordered = TRUE)
        } else {
                if(!all(rownames(pValMatrix) %in% tax_order)){
                        stop("given tax_order does not fit to tax_names")
                }
                DF$Taxon <- factor(DF$Taxon, levels = tax_order, ordered = TRUE)
                DF$HostTaxon <- factor(DF$HostTaxon, levels = rev(tax_order), ordered = TRUE)
                
        }
        
        fill_colors <- c(color_levels, ns = "gray98", " " = "black")
        
        DF$Fill <- "ns"
        DF$Fill[DF$pValue < signi_level & DF$pValue > 0] <- i
        DF$Fill[DF$pValue > -1*signi_level & DF$pValue < 0] <- j
        DF$Fill[DF$pValue == 2] <- " "
        DF$Fill <- factor(DF$Fill, levels = names(fill_colors), ordered = T)
        TileTr <- ggplot(DF, aes(x = Taxon, y = HostTaxon, fill = Fill))
        TileTr <- TileTr + 
                geom_raster() + 
                # ggtitle(paste(i, " vs ", j, sep = "")) +
                scale_fill_manual("", values = fill_colors) +
                scale_x_discrete(position = "top") +
                labs(x=NULL, y=NULL) +
                theme_bw() +
                #theme_tufte(base_family="Helvetica") +
                theme(panel.border = element_blank(),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
                      axis.ticks=element_blank())
        
        TileTr
        
        
}
# --



































########################### Functions currently not in use #######################################






# --
####################################
## plot_taxa_ratiosSingle 
###################################
# see plot_taxa_ratios_levelPairs: Here you directly calculate the count by count ratio matrix only for the tax_nom, you still facet by taxa_den
# (denominator) but you keep all levels in all plots. Makes extensive use of ggpubr, NB: ggpubr is so smart to adjust p-values when you use
# scale_y_log10

## Output: 
# - list(pVals = pVals, Tr = Tr, Tr1 = Tr1, pValsLog = pValsLog, Tr2 = Tr2, Tr3 = Tr3)

plot_taxa_ratiosSingle <- function(physeq, group_var, color_levels, tax_names = NULL,
                                   taxa_nom = "Firmicutes", taxa_den = NULL, test = "t.test", p_adjust_method = "fdr",
                                   tax_order = NULL,
                                   symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), hide.ns = FALSE) {
        
        
        if (taxa_are_rows(physeq)) { 
                physeq <- t(physeq)
        }
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the loaded phyloseq object.")
        }
        
        if (!all(names(color_levels) %in% unique(sample_data(physeq)[[group_var]]))) {
                stop("Not all names in names(color_levels)are actually levels in the group_var column.")
        }
        
        if (!(test %in% c("t.test", "wilcox.test"))) {
                stop("test should be t.test or wilcox.test")
        }
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        compare <- names(color_levels)
        
        if (!is.null(compare)) {
                group_var_levels <- compare
        } else {
                group_var_levels <- levels(group_fac)
        }
        
        # if (length(group_var_levels) != 2) {
        #         stop(paste0("compare (group_var_levels) must consist of two groups - you asked for ", 
        #                     paste(group_var_levels, collapse = ", ")))
        # }
        
        if (!all(group_var_levels %in% levels(group_fac))) {
                stop("Not all given compare (group_var_levels) are actually levels in group_var column.")
        }
        
        
        # - check that given tax_names fit to physeq and change taxa_names of physeq -
        if (is.null(tax_names)){
                tax_names <- paste("T", 1:ntaxa(physeq), sep = "_")
        } 
        
        if(!identical(ntaxa(physeq), length(tax_names))){stop("tax_names do not fit in length to physeq")}
        
        tax_names <- make.unique(tax_names)
        taxa_names(physeq) <- tax_names
        # --
        
        # - calculate the matrix taxa_nom/(all other taxa) -         
        CT <- t(as(otu_table(physeq), 'matrix')) # now taxa are rows and samples are columns
        i <- group_var_levels[1]
        j <- group_var_levels[2]
        
        CT <- CT[, group_fac %in% c(i, j)]
        
        
        
        i <- which(rownames(CT) == taxa_nom)
        if (length(i) != 1) {stop("taxa_nom not found in tax_names or tax_names not unique!")}
        
        
        TbTmatrix <- apply(CT, 2, function(samp_cnts){samp_cnts[i]/samp_cnts})
        # produces for each taxon (= host taxon) a TbTMatrix
        # NB: there are possibly Inf, and NaN values in the matrix, specifically
        # 0/x = 0, x/0 = Inf; 0/0 = NaN!
        # --
        
        TbT_DF <- as.data.frame(TbTmatrix)
        TbT_DF$Taxon <- rownames(TbT_DF)
        # - use taxa_den (denominator) to restrict the taxa to which taxa_nom is compared to -
        if (is.null(taxa_den)) {taxa_den <- tax_names}
        TbT_DF <- TbT_DF[TbT_DF$Taxon %in% taxa_den, ]
        # --
        
        # - change to long DF -
        TbT_DF_l <- gather(TbT_DF, key = Sample, value = Ratio, -Taxon)
        # --
        
        # - add the group_var level information  -
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        TbT_DF_l$Group <- as.character(LookUpDF$Group[match(TbT_DF_l$Sample, LookUpDF$Sample)])
        TbT_DF_l$Group <- factor(TbT_DF_l$Group, levels = group_var_levels, ordered = T)
        # --
        
        # - change all ratios where either nominator taxon or denominator taxon had count = 0 to NA -
        # remember: 0/0 = NaN (not finite), 0/x = 0, x/0 = Inf
        TbT_DF_l$Ratio[!is.finite(TbT_DF_l$Ratio) | TbT_DF_l$Ratio == 0] <- NA
        # --
        
        # - find taxa that would throw an error in statistical test and remove those taxa from DF -
        # first find the taxa that would throw an error in t.test or wilcox.test
        var_plus_length_check <- group_by(TbT_DF_l, Taxon, Group) %>% summarise(Variance = var(Ratio, na.rm = T), NotNA = sum(!is.na(Ratio)))
        if (test == "t.test"){
                var_plus_length_check <- dplyr::filter(var_plus_length_check, !(Variance > 0) | NotNA < 2) # variance > 0 also to remove test where host_taxon == taxon
        } else if (test == "wilcox.test") {
                var_plus_length_check <- dplyr::filter(var_plus_length_check, !(Variance > 0) | NotNA < 1)
        }
        
        if (nrow(var_plus_length_check) != 0) { 
                TbT_DF_l <- filter(TbT_DF_l, !(Taxon %in% unique(var_plus_length_check$Taxon)))
        }
        # --
        
        # - use ggpubr::compare_menas to calculate all pValues of the Ratios for the different taxa_den between current group_levels -
        pVals <- ggpubr::compare_means(formula = Ratio ~ Group, data = TbT_DF_l, group.by = "Taxon", method = test, p.adjust.method = p_adjust_method, symnum.args = symnum.args)
        
        pVals <- dplyr::arrange(pVals, p)
        
        # NB: the pVals change when you log the ratios (scale_y_log10())
        TbT_DF_l$RatioLog10 <- log10(TbT_DF_l$Ratio)
        pValsLog <- ggpubr::compare_means(formula = RatioLog10 ~ Group, data = TbT_DF_l, group.by = "Taxon", method = test, p.adjust.method = p_adjust_method, symnum.args = symnum.args)
        pValsLog <- dplyr::arrange(pValsLog, p)
        # --
        
        # NB: I plot now for log and non log independently, even though ggpubr is so smart to change the p-values when you just use
        # Tr + scale_y_log10(). I plot independently because log might change the order!
        TbT_DF_l_log <- TbT_DF_l
        
        # - order taxa based on pVals result or based on tax_order -
        if (is.null(tax_order)){
                TbT_DF_l$Taxon <- factor(TbT_DF_l$Taxon, levels = unique(pVals$Taxon), ordered = TRUE)
                TbT_DF_l_log$Taxon <- factor(TbT_DF_l_log$Taxon, levels = unique(pValsLog$Taxon), ordered = TRUE)
                
        } else {
                if(!all(unique(pVals$Taxon) %in% tax_order)){
                        stop("given tax_order does not fit to tax_names")
                }
                TbT_DF_l$Taxon <- factor(TbT_DF_l$Taxon, levels = tax_order, ordered = TRUE)
                TbT_DF_l_log$Taxon <- factor(TbT_DF_l_log$Taxon, levels = tax_order, ordered = TRUE)
        }
        # --
        
        # - since you might have more than two levels in each plot you need to set the comparisons argument in stat_compare_means -
        # comparisonList <- list(group_var_levels)
        # --
        
        # - plot: NB: also for non removed taxa some samples might have NA ratios that will be removed -
        
        Tr <- ggplot(TbT_DF_l, aes(x = Group, y = Ratio, col = Group))
        Tr <- Tr +
                geom_violin() +
                geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                scale_color_manual(values = color_levels) +
                facet_wrap(~ Taxon, scales = "free_y") +
                xlab("") +
                ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                theme_bw() +
                theme(legend.position = "none")
        
        # Tr <- Tr + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = "t.test", hide.ns = hide.ns)
        Tr <- Tr + ggpubr::stat_compare_means(label = "p.signif", method = test, label.x = 1.5, hide.ns = hide.ns)
        
        Tr1 <- ggplot(TbT_DF_l, aes(x = Group, y = Ratio, col = Group))
        Tr1 <- Tr1 +
                geom_boxplot(outlier.color = NA) +
                geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                scale_color_manual(values = color_levels) +
                facet_wrap(~ Taxon, scales = "free_y") +
                xlab("") +
                ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                theme_bw() +
                theme(legend.position = "none")
        
        
        Tr1 <- Tr1 + ggpubr::stat_compare_means(label.x = 1.5, label = "p.signif", method = test, hide.ns = hide.ns) # p.format
        
        
        Tr2 <- ggplot(TbT_DF_l_log, aes(x = Group, y = Ratio, col = Group))
        Tr2 <- Tr2 +
                geom_violin() +
                geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                scale_color_manual(values = color_levels) +
                facet_wrap(~ Taxon, scales = "free_y") +
                xlab("") +
                ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                theme_bw() +
                theme(legend.position = "none")
        
        
        Tr2 <- Tr2 + ggpubr::stat_compare_means(label.x = 1.5, label = "p.signif", method = test, hide.ns = hide.ns)
        Tr2 <- Tr2 + scale_y_log10()
        
        
        Tr3 <- ggplot(TbT_DF_l_log, aes(x = Group, y = Ratio, col = Group))
        Tr3 <- Tr3 +
                geom_boxplot(outlier.color = NA) +
                geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                scale_color_manual(values = color_levels) +
                facet_wrap(~ Taxon, scales = "free_y") +
                xlab("") +
                ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                theme_bw() +
                theme(legend.position = "none")
        
        Tr3 <- Tr3 + ggpubr::stat_compare_means(label.x = 1.5, label = "p.signif", method = test, hide.ns = hide.ns)
        Tr3 <- Tr3 + scale_y_log10()
        
        
        list(pVals = pVals, Tr = Tr, Tr1 = Tr1, pValsLog = pValsLog, Tr2 = Tr2, Tr3 = Tr3)
        
}
#--





# --
####################################
## plot_taxa_ratios_levelPairs 
###################################
## Input: 
# - TbTmatrixes_list: The list with the lists of TbTmatrixes for each level combi in group factor
# NB: here it should be raw_TbTmatrixes with 0 = 0/x, Inf = x/0, NaN = 0/0, all these values will be ignored in the ratio plots by being set to NA!!
# NB: NAs are also ignored in t.test or wilcoxon test
# - physeq: used for TbTmatrixes_list generation
# - group_var: the factor used to separate the samples
# - tax_names: the names of the taxa in physeq you want to use, e.g. taxa_names(physeq) if you like them, if NULL T1 to Tn will be used
# - taxa_nom: the nominator taxon of the abundance ratios, NB: only 1 allowed, must be included in tax_names
# - taxa_den: the denominator taxa of the abundance ratios, several allowed, all must be included in tax_names, if NULL all are used,
# i.e you get plots facet_wrapped around the taxa_den taxa
# NB: only taxa are plotted for which statistical test is possible!
# - test: either "t.test" or "wilcoxon"
# - p_adjust_method

## Output: 
# - list of pVals data frame (the pValues from t.test or wilcox.test after p.adjust) plus Violin plot plus Boxplot,
# so for each level combi in group_var a list of length 2.


plot_taxa_ratios_levelPairs <- function(TbTmatrixes_list, physeq, group_var, tax_names = NULL,
                                        taxa_nom = "Firmicutes", taxa_den = NULL, color_levels, test = "t.test", p_adjust_method = "fdr",
                                        symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) {
        
        if(!identical(length(TbTmatrixes_list[[1]]), ntaxa(physeq))){stop("TbTmatrixes can not fit to physeq")}
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # - get the level combis -
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_sAndj_s <- get_unique_facLevel_combis(fac_levels)
        i_s <- i_sAndj_s[["i_s"]]
        j_s <- i_sAndj_s[["j_s"]]
        # --
        
        result_list <- vector("list", length = length(i_s))
        
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        # order so the samples of a group are shown beside of each other (not really necessary for a violin plot)
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        
        
        for (k in seq_along(i_s)) {
                
                i <- i_s[k]
                j <- j_s[k]
                
                TbTmatrixes <- TbTmatrixes_list[[k]]
                
                if (is.null(tax_names)){
                        tax_names <- paste("T", 1:length(TbTmatrixes), sep = "_")
                } else {
                        if(!identical(length(TbTmatrixes), length(tax_names))){stop("tax_names do not fit in length to TbTmatrixes")}
                }
                
                names(TbTmatrixes) <- make.unique(tax_names)
                names(TbTmatrixes) <- tax_names #tax_names should be unique
                
                TbTmatrix <- TbTmatrixes[[taxa_nom]] # the only relevant matrix for this plot, it is where tax_nom counts have been divided by the counts of all other taxa
                
                if (is.null(TbTmatrix)) {stop("taxa_nom not found")}
                
                rownames(TbTmatrix) <- tax_names
                
                TbT_DF <- as.data.frame(TbTmatrix)
                TbT_DF$Taxon <- rownames(TbT_DF)
                # - use taxa_den (denominator) to restrict the taxa to which taxa_nom is compared to -
                if (is.null(taxa_den)) {taxa_den <- tax_names}
                TbT_DF <- TbT_DF[TbT_DF$Taxon %in% taxa_den, ]
                # --
                
                # - change to long DF -
                TbT_DF_l <- gather(TbT_DF, key = Sample, value = Ratio, -Taxon)
                # --
                # - add the group_var level information using the original level order in group_var -
                TbT_DF_l$Group <- as.character(LookUpDF$Group[match(TbT_DF_l$Sample, LookUpDF$Sample)])
                group_fac_current <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)]) # keeps the current levels in original order
                TbT_DF_l$Group <- factor(TbT_DF_l$Group, levels = levels(group_fac_current), ordered = T)
                # --
                
                # - change all ratios where either nominator taxon or denominator taxon had count = 0 to NA -
                # remember: 0/0 = NaN (not finite), 0/x = 0, x/0 = Inf
                TbT_DF_l$Ratio[!is.finite(TbT_DF_l$Ratio) | TbT_DF_l$Ratio == 0] <- NA
                # --
                
                # - find taxa that would throw an error in statistical test and remove those taxa from DF -
                # first find the taxa that would throw an error in t.test or wilcox.test
                var_plus_length_check <- group_by(TbT_DF_l, Taxon, Group) %>% summarise(Variance = var(Ratio, na.rm = T), NotNA = sum(!is.na(Ratio)))
                if (test == "t.test"){
                        var_plus_length_check <- filter(var_plus_length_check, !(Variance > 0) | NotNA < 2)
                } else if (test == "wilcox.test") {
                        var_plus_length_check <- filter(var_plus_length_check, NotNA < 1)
                }
                
                if (nrow(var_plus_length_check) != 0) { # should always remove at least the case where taxa_denominator = taxa_nom, since those ratios are all 1, hence no variance
                        TbT_DF_l <- filter(TbT_DF_l, !(Taxon %in% unique(var_plus_length_check$Taxon)))
                }
                # --
                
                # - use ggpubr::compare_menas to calculate all pValues of the Ratios for the different taxa_den between current group_levels -
                pVals <- ggpubr::compare_means(formula = Ratio ~ Group, data = TbT_DF_l, group.by = "Taxon", method = test, p.adjust.method = p_adjust_method, symnum.args = symnum.args)
                
                pVals <- dplyr::arrange(pVals, p)
                
                # NB: the pVals change when you log the ratios (scale_y_log10())
                TbT_DF_l$RatioLog10 <- log10(TbT_DF_l$Ratio)
                pValsLog <- ggpubr::compare_means(formula = RatioLog10 ~ Group, data = TbT_DF_l, group.by = "Taxon", method = test, p.adjust.method = p_adjust_method, symnum.args = symnum.args)
                pValsLog <- dplyr::arrange(pValsLog, p)
                # --
                
                # NB: I plot now for log and non log independently, even though ggpubr is so smart to change the p-values when you just use
                # Tr + scale_y_log10(). I plot independently because log might change the order!
                TbT_DF_l_log <- TbT_DF_l
                
                # - order taxa based on pVals result -
                TbT_DF_l$Taxon <- factor(TbT_DF_l$Taxon, levels = pVals$Taxon, ordered = TRUE)
                TbT_DF_l_log$Taxon <- factor(TbT_DF_l_log$Taxon, levels = pValsLog$Taxon, ordered = TRUE)
                # --
                
                
                # - plot: NB: also for non removed taxa some samples might have NA ratios that will be removed -
                Tr <- ggplot(TbT_DF_l, aes(x = Group, y = Ratio, col = Group))
                Tr <- Tr +
                        geom_violin() +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                        # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                        scale_color_manual(values = color_levels) +
                        facet_wrap(~ Taxon, scales = "free_y") +
                        xlab("") +
                        ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                        theme_bw() +
                        theme(legend.position = "none")
                
                
                Tr <- Tr + ggpubr::stat_compare_means(label = "p.signif", method = test, label.x = 1.5)
                
                
                Tr1 <- ggplot(TbT_DF_l, aes(x = Group, y = Ratio, col = Group))
                Tr1 <- Tr1 +
                        geom_boxplot(outlier.color = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                        # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                        scale_color_manual(values = color_levels) +
                        facet_wrap(~ Taxon, scales = "free_y") +
                        xlab("") +
                        ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                        theme_bw() +
                        theme(legend.position = "none")
                
                
                Tr1 <- Tr1 + ggpubr::stat_compare_means(label = "p.signif", method = test, label.x = 1.5) # p.format
                
                
                Tr2 <- ggplot(TbT_DF_l_log, aes(x = Group, y = Ratio, col = Group))
                Tr2 <- Tr2 +
                        geom_violin() +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                        # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                        scale_color_manual(values = color_levels) +
                        facet_wrap(~ Taxon, scales = "free_y") +
                        xlab("") +
                        ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                        theme_bw() +
                        theme(legend.position = "none")
                
                
                Tr2 <- Tr2 + ggpubr::stat_compare_means(label = "p.signif", method = test, label.x = 1.5)
                Tr2 <- Tr2 + scale_y_log10()
                
                
                Tr3 <- ggplot(TbT_DF_l_log, aes(x = Group, y = Ratio, col = Group))
                Tr3 <- Tr3 +
                        geom_boxplot(outlier.color = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                        # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                        scale_color_manual(values = color_levels) +
                        facet_wrap(~ Taxon, scales = "free_y") +
                        xlab("") +
                        ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                        theme_bw() +
                        theme(legend.position = "none")
                
                
                Tr3 <- Tr3 + ggpubr::stat_compare_means(label = "p.signif", method = test, label.x = 1.5)
                Tr3 <- Tr3 + scale_y_log10()
                
                
                # --
                
                result_list[[k]] <- list(pVals = pVals, Tr = Tr, Tr1 = Tr1, pValsLog = pValsLog, Tr2 = Tr2, Tr3 = Tr3)
                # rm(Tr, TbT_DF_l, TbT_DF, TbT_DF_l_test, pVals)
        }
        
        names(result_list) <- names(TbTmatrixes_list)
        result_list
}
# --





# --
####################################
## calculate_raw_TbTmatrixes_Pairs:
###################################
# see calculate_TbTmatrixes, this one just calculates the "raw" ratio matrixes i.e. without log and gm

## Input: 
# - physeq = a phyloseq object
# - group_var, name of the group_fac in sample_data(physeq)
## Output: 
# - list of TbTmatrixes, one list for each combi of levels in group_fac, list items are named by level_vs_level


calculate_raw_TbTmatrixes_Pairs = function(physeq, group_var){
        
        if (taxa_are_rows(physeq)) {physeq <- t(physeq)}
        
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # - get the level combis -
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_sAndj_s <- get_unique_facLevel_combis(fac_levels)
        i_s <- i_sAndj_s[["i_s"]]
        j_s <- i_sAndj_s[["j_s"]]
        # --
        
        TbTmatrixes_list <- vector("list", length = length(i_s))
        
        for (k in seq_along(i_s)) {
                
                i <- i_s[k]
                j <- j_s[k]
                
                
                CT <- t(as(otu_table(physeq), 'matrix')) # now taxa are rows and samples are columns
                group_fac_current <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])
                CT <- CT[, group_fac %in% group_fac_current]
                
                
                
                TbTmatrixes <- lapply(1:nrow(CT), function(i){apply(CT, 2, function(samp_cnts){samp_cnts[i]/samp_cnts})})
                # produces for each taxon (= host taxon) a TbTMatrix
                # NB: there are Inf, and NaN values in the matrixes, specifically
                # 0/x = 0, x/0 = Inf; 0/0 = NaN!
                
                names(TbTmatrixes) <- rownames(TbTmatrixes[[1]])
                TbTmatrixes_list[[k]] <- TbTmatrixes
                names(TbTmatrixes_list)[[k]] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        TbTmatrixes_list
        
}
# --



# --
####################################
## calculate_raw_TbTmatrixesSingle:
###################################

calculate_raw_TbTmatrixesSingle = function(physeq, group_var, compare){
        
        if (taxa_are_rows(physeq)) {physeq <- t(physeq)}
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the loaded phyloseq object.")
        }
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        if (!is.null(compare)) {
                group_var_levels <- compare
        } else {
                group_var_levels <- levels(group_fac)
        }
        
        if (length(group_var_levels) != 2) {
                stop(paste0("compare (group_var_levels) must consist of two groups - you asked for ", 
                            paste(group_var_levels, collapse = ", ")))
        }
        
        if (!all(group_var_levels %in% levels(group_fac))) {
                stop("Not all given compare (group_var_levels) are actually levels in group_var column.")
        }
        
        
        CT <- t(as(otu_table(physeq), 'matrix')) # now taxa are rows and samples are columns
        #group_fac_current <- droplevels(group_fac[group_fac %in% c(i, j)])
        CT <- CT[, group_fac %in% group_var_levels]
        
        
        
        TbTmatrixes <- lapply(1:nrow(CT), function(i){apply(CT, 2, function(samp_cnts){samp_cnts[i]/samp_cnts})})
        # produces for each taxon (= host taxon) a TbTMatrix
        # NB: there are Inf, and NaN values in the matrixes, specifically
        # 0/x = 0, x/0 = Inf; 0/0 = NaN!
        
        names(TbTmatrixes) <- rownames(TbTmatrixes[[1]])
        
        
        TbTmatrixes
        
}
# --







# --
####################################
## create_raw_TbT_TilePlots_Pairs: 
###################################
# wilcoxon.test is used
## Input: 
# - TbTmatrixes_list: The list with the lists of TbTmatrixes for each level combi in group factor
# NB: here it should be raw_TbTmatrixes with 0 = 0/x, Inf = x/0, NaN = 0/0, all these values will be ignored!!
# - physeq: used for TbTmatrixes_list generation
# - p_adjust: method for p.adjust, NB: if not bonferroni or none, the tile plots are not necessarily symmetric anymore
## Output: 
# - list of pValMatrixes plus TileTr for each level combination. 
# NB: I negated the p-values if host taxon was more abundant in grp2 compared to other taxon!!

create_raw_TbT_TilePlots_Pairs <- function(TbTmatrixes_list, physeq, group_var, tax_names = NULL,
                                     test = "wilcoxon", p_adjust = "none") {
        
        if(!identical(length(TbTmatrixes_list[[1]]), ntaxa(physeq))){stop("TbTmatrixes can not fit to physeq")}
        
        if(test != "wilcoxon" & test != "t.test"){stop("test unknown, must be wilcoxon or t.test")}
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # - get the level combis -
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_sAndj_s <- get_unique_facLevel_combis(fac_levels)
        i_s <- i_sAndj_s[["i_s"]]
        j_s <- i_sAndj_s[["j_s"]]
        # --
        
        result_list <- vector("list", length = length(i_s))
        
        for (k in seq_along(i_s)) {
                
                i <- i_s[k]
                j <- j_s[k]
                
                TbTmatrixes <- TbTmatrixes_list[[k]]
                
                if (is.null(tax_names)){
                        tax_names <- paste("T", 1:length(TbTmatrixes), sep = "_")
                } else {
                        if(!identical(length(TbTmatrixes), length(tax_names))){stop("tax_names do not fit in length to TbTmatrixes")}
                }
                
                tax_names <- make.unique(tax_names)
                names(TbTmatrixes) <- tax_names
                
                TbTmatrixes <- lapply(TbTmatrixes, function(mat){
                        rownames(mat) <- tax_names
                        mat
                })
                group_fac_current <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])
                
                # ntaxa * ntaxa wilcoxon tests take time!
                pValMatrix <- sapply(TbTmatrixes, function(mat){
                        apply(mat, 1, function(taxon_ratios){
                                x <- taxon_ratios[group_fac_current == fac_levels[i]]
                                x <- x[is.finite(x) & x != 0] # removes all ratios in which one of the two taxa was not present!
                                y <- taxon_ratios[group_fac_current == fac_levels[j]]
                                y <- y[is.finite(y) & y != 0] # removes 0/0 = NaN, 0/x = 0, x/0 = Inf
                                if (test == "wilcoxon"){
                                        if (length(x) > 0 && length(y) > 0){
                                                pValue <- wilcox.test(x = x, y = y, alternative = "two", paired = F, exact = F)$p.value
                                                # NB: wilcox.test ignores in default setting (maybe see na.action) NA, NaN, Inf, -Inf
                                                # For plot: change sign of pValue to negative if taxon is more abundant in group 1. 
                                                Ranks <- rank(c(x[!is.na(x)], y[!is.na(y)]))
                                                n1 <- length(x[!is.na(x)])
                                                n2 <- length(y[!is.na(y)])
                                                Wx <- sum(Ranks[1:n1])-(n1*(n1+1)/2)
                                                Wy <- sum(Ranks[(n1+1):(n1+n2)])-(n2*(n2+1)/2)
                                                if(Wx > Wy){pValue <- -1*pValue}
                                                pValue
                                                
                                        } else {
                                                pValue = 1
                                        }
                                        
                                } else if (test == "t.test") {
                                        if (length(x) > 1 && length(y) > 1 && var(x) > 0 && var(y) > 0){
                                                pValue <- t.test(x = x, y = y, alternative = "two")$p.value
                                                if (mean(x, na.rm = T) > mean(y, na.rm = T)){pValue <- -1*pValue}
                                                pValue
                                        } else {
                                                pValue <- 1
                                        }
                                        
                                }
                                
                        })
                })
                
                # make sure diagonal is all NA (can be exceptions especially for t.test)
                diag(pValMatrix) <- NA
                
                signs <- pValMatrix < 0
                signs[is.na(signs)] <- FALSE
                
                pValMatrix <- abs(pValMatrix)
                for (e in 1:nrow(pValMatrix)){
                        pValMatrix[e, ] <- p.adjust(pValMatrix[e, ], method = p_adjust)
                } # equal to t(apply(pValMatrix, 1, p.adjust, method = p_adjust))
                
                pValMatrix[signs] <- pValMatrix[signs]*(-1)
                
                # -- add a tile plot of the pValMatrix --
                
                DF <- as.data.frame(pValMatrix)
                DF[is.na(DF)] <- 2 # just to avoid missing values in plot and have a clear non-pValue value to mark self comparisons as black
                DF$HostTaxon <- rownames(pValMatrix)
                DF <- tidyr::gather(DF, key = Taxon , value = pValue, - HostTaxon)
                DF$Taxon <- factor(DF$Taxon, levels = rownames(pValMatrix), ordered = TRUE)
                DF$HostTaxon <- factor(DF$HostTaxon, levels = rev(rownames(pValMatrix)), ordered = TRUE)
                
                # # add color to the taxa names so you see up and down
                # colyaxis <- vector(mode = "character", length = nrow(DF))
                # colyaxis[] <- "black"
                # colyaxis[grepl("TP-U", levels(DF$HostTaxon))] <- "#E69F00"
                # colyaxis[grepl("TP-D", levels(DF$HostTaxon))] <- "#009E73"
                # colxaxis <- vector(mode = "character", length = nrow(DF))
                # colxaxis[] <- "black"
                # colxaxis[grepl("TP-U", levels(DF$Taxon))] <- "#E69F00"
                # colxaxis[grepl("TP-D", levels(DF$Taxon))] <- "#009E73"
                DF$Fill <- "not significant"
                DF$Fill[DF$pValue < 0.05 & DF$pValue > 0] <- "up (p < 0.05)"
                DF$Fill[DF$pValue > -0.05 & DF$pValue < 0] <- "down (p < 0.05)"
                DF$Fill[DF$pValue == 2] <- " "
                DF$Fill <- factor(DF$Fill, levels = c("up (p < 0.05)", "not significant", "down (p < 0.05)", " "), ordered = T)
                TileTr <- ggplot(DF, aes(x = Taxon, y = HostTaxon, fill = Fill))
                TileTr <- TileTr + 
                        geom_raster() + 
                        ggtitle(names(TbTmatrixes_list)[k]) +
                        scale_fill_manual("", values = c("not significant" = "gray98", "up (p < 0.05)" = "#AA4499", "down (p < 0.05)" = "#88CCEE", " " = "black")) +
                        scale_x_discrete(position = "top") +
                        labs(x=NULL, y=NULL) +
                        theme_tufte(base_family="Helvetica") +
                        theme(plot.title=element_text(hjust=0)) +
                        theme(axis.ticks=element_blank()) +
                        theme(axis.text=element_text(size=7)) +
                        theme(legend.title=element_blank()) +
                        theme(legend.text=element_text(size=6)) +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)) # colour = colxaxis
                result_list[[k]] <- list(pValMatrix = pValMatrix, TileTr = TileTr)
        }
        
        names(result_list) <- names(TbTmatrixes_list)
        result_list
}
# --







# --
####################################
## create_raw_TbT_TilePlotsSingle2: 
###################################
# NB: this version demands that TbTmatrixes only contain the samples defined by names(color_levels)

create_raw_TbT_TilePlotsSingle2 <- function(TbTmatrixes, physeq, group_var, color_levels, tax_names = NULL, tax_order = NULL, 
                                            test = "wilcoxon", signi_level = 0.05, p_adjust_method = "none") {
        
        if(!identical(length(TbTmatrixes), ntaxa(physeq))){stop("TbTmatrixes don't fit to physeq")}
        
        if(test != "wilcoxon" & test != "t.test"){stop("test unknown, must be wilcoxon or t.test")}
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the loaded phyloseq object.")
        }
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        compare <- names(color_levels)
        
        if (!is.null(compare)) {
                group_var_levels <- compare
        } else {
                group_var_levels <- levels(group_fac)
        }
        
        if (length(group_var_levels) != 2) {
                stop(paste0("compare (names(color_levels)) must consist of two groups - you asked for ", 
                            paste(group_var_levels, collapse = ", ")))
        }
        
        if (!all(group_var_levels %in% levels(group_fac))) {
                stop("Not all given compare (group_var_levels) are actually levels in group_var column.")
        }
        
        # special to this Single version
        if (ncol(TbTmatrixes[[1]]) != sum(group_fac %in% group_var_levels)){
                stop("TbTmatrixes do not fit to the group_fac in physeq restricted to samples defined by names(color_levels).")
        }
        
        
        # - check that given tax_names fit to physeq and change taxa_names of physeq -
        if (is.null(tax_names)){
                tax_names <- paste("T", 1:ntaxa(physeq), sep = "_")
        } 
        
        if(!identical(ntaxa(physeq), length(tax_names))){stop("tax_names do not fit in length to physeq")}
        
        tax_names <- make.unique(tax_names)
        # --
        
        
        
        names(TbTmatrixes) <- tax_names
        
        TbTmatrixes <- lapply(TbTmatrixes, function(mat){
                rownames(mat) <- tax_names
                mat
        })
        group_fac_current <- droplevels(group_fac[group_fac %in% group_var_levels])
        
        i <- group_var_levels[1]
        j <- group_var_levels[2]
        
        # ntaxa * ntaxa wilcoxon tests take time if you have a lot of taxa!
        pValMatrix <- sapply(TbTmatrixes, function(mat){
                apply(mat, 1, function(taxon_ratios){
                        x <- taxon_ratios[group_fac_current == i]
                        x <- x[is.finite(x) & x != 0] # removes all ratios in which one of the two taxa was not present!
                        y <- taxon_ratios[group_fac_current == j]
                        y <- y[is.finite(y) & y != 0] # removes 0/0 = NaN, 0/x = 0, x/0 = Inf
                        if (test == "wilcoxon"){
                                if (length(x) > 0 && length(y) > 0){
                                        pValue <- wilcox.test(x = x, y = y, alternative = "two", paired = F, exact = F)$p.value
                                        # NB: wilcox.test ignores in default setting (maybe see na.action) NA, NaN, Inf, -Inf
                                        # For plot: change sign of pValue to negative if taxon is more abundant in group 1. 
                                        Ranks <- rank(c(x[!is.na(x)], y[!is.na(y)]))
                                        n1 <- length(x[!is.na(x)])
                                        n2 <- length(y[!is.na(y)])
                                        Wx <- sum(Ranks[1:n1])-(n1*(n1+1)/2)
                                        Wy <- sum(Ranks[(n1+1):(n1+n2)])-(n2*(n2+1)/2)
                                        if(Wx > Wy){pValue <- -1*pValue}
                                        pValue
                                        
                                } else {
                                        pValue = 1
                                }
                                
                        } else if (test == "t.test") {
                                if (length(x) > 1 && length(y) > 1 && var(x) > 0 && var(y) > 0){
                                        pValue <- t.test(x = x, y = y, alternative = "two")$p.value
                                        if (mean(x, na.rm = T) > mean(y, na.rm = T)){pValue <- -1*pValue}
                                        pValue
                                } else {
                                        pValue <- 1
                                }
                                
                        }
                        
                })
        })
        
        # make sure diagonal is all NA (can be exceptions especially for t.test)
        diag(pValMatrix) <- NA
        
        # - adjust p-values if asked for -        
        signs <- pValMatrix < 0
        signs[is.na(signs)] <- FALSE
        
        pValMatrix <- abs(pValMatrix)
        for (e in 1:nrow(pValMatrix)){
                pValMatrix[e, ] <- p.adjust(pValMatrix[e, ], method = p_adjust_method)
        } # equal to t(apply(pValMatrix, 1, p.adjust, method = p_adjust))
        
        pValMatrix[signs] <- pValMatrix[signs]*(-1)
        # --
        
        # -- add a tile plot of the pValMatrix --
        
        DF <- as.data.frame(pValMatrix)
        DF[is.na(DF)] <- 2 # just to avoid missing values in plot and have a clear non-pValue value to mark self comparisons as black
        DF$HostTaxon <- rownames(pValMatrix)
        DF <- tidyr::gather(DF, key = Taxon , value = pValue, - HostTaxon)
        if (is.null(tax_order)) {
                DF$Taxon <- factor(DF$Taxon, levels = rownames(pValMatrix), ordered = TRUE)
                DF$HostTaxon <- factor(DF$HostTaxon, levels = rev(rownames(pValMatrix)), ordered = TRUE)
        } else {
                if(!all(rownames(pValMatrix) %in% tax_order)){
                        stop("given tax_order does not fit to tax_names")
                }
                DF$Taxon <- factor(DF$Taxon, levels = tax_order, ordered = TRUE)
                DF$HostTaxon <- factor(DF$HostTaxon, levels = rev(tax_order), ordered = TRUE)
                
        }
        
        fill_colors <- c(color_levels, ns = "gray98", " " = "black")
        
        DF$Fill <- "ns"
        DF$Fill[DF$pValue < signi_level & DF$pValue > 0] <- i
        DF$Fill[DF$pValue > -1*signi_level & DF$pValue < 0] <- j
        DF$Fill[DF$pValue == 2] <- " "
        DF$Fill <- factor(DF$Fill, levels = names(fill_colors), ordered = T)
        TileTr <- ggplot(DF, aes(x = Taxon, y = HostTaxon, fill = Fill))
        TileTr <- TileTr + 
                geom_raster() + 
                # ggtitle(paste(i, " vs ", j, sep = "")) +
                scale_fill_manual("", values = fill_colors) +
                scale_x_discrete(position = "top") +
                labs(x=NULL, y=NULL) +
                theme_bw() +
                #theme_tufte(base_family="Helvetica") +
                theme(panel.border = element_blank(),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
                      axis.ticks=element_blank())
        
        TileTr
        
        
}
# --