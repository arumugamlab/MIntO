# --
#######################################
### boxplot_sampleSums
#######################################
## INPUTs
# color_levels: not only defines the colors but also the samples that will be compared
# so the function allows more than two color_levels. 
# The order of the names in the color_levels defines the order of plotting.


# NB: makes use of ggpubr::stat_compare_means that allows you to adjust more settings, e.g. plotting symbols instead of numbers (change label to "p.signif")

boxplot_sampleSums <- function(physeq, group_var, color_levels, shape, test = "t.test", p_adjust_method = "fdr",
                               symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                               hide.ns = FALSE){
        
        DF <- cbind(sample_data(physeq), total_count = sample_sums(physeq))
        # keep only samples that belong to the names in color_levels
        DF <- DF[DF[[group_var]] %in% names(color_levels),]
        
        # make sure the order of the groups is as defined in color levels
        DF[[group_var]] <- factor(DF[[group_var]], levels = names(color_levels), ordered = TRUE)
        
        Tr <-  ggplot(DF, aes_string(x = group_var, y = "total_count", color = group_var)) 
        Tr <- Tr + 
                geom_boxplot(outlier.color = NA) +
                geom_jitter(aes_string(shape = shape), width = .2, height = 0, alpha = 0.65) +
                xlab("") +
                ylab("sample_sums()") +
                scale_color_manual("", values = color_levels) +
                theme_bw()
        if(is.null(shape)){
                Tr <- Tr + theme(legend.position = "none")
        }
        
        # - since you might have more than two levels in each plot you need to set the comparisons argument in stat_compare_means -
        group_fac <- factor(DF[[group_var]])
        fac_levels <- levels(group_fac)
        
        comparisonList <- get_unique_facLevel_combinations(fac_levels)
        
        Tr <- Tr + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = test, hide.ns = hide.ns)
        
        formulaa <- as.formula(paste("total_count ~", group_var, sep = " "))
        
        pVals <- compare_means(formula = formulaa, data = DF, method = test, p.adjust.method = p_adjust_method, symnum.args = symnum.args)
        
        list(pVals = pVals, Tr = Tr)
}
# --





# --
#######################################
### calc_alphadiv_plusLmResids
#######################################
## INPUTs:
# group_var and compare are used here to allow the analysis of only a subset of the samples in physeq.
# Output:
# outlist: 
#    [[1]]: DF of alpha_diversity values, plus residuals of the values to a linear fit against total counts (sample_sums)
#    [[2]]: list of the fit objects derived from lm of alpha diversity values against sample_sums()

# NB: is based on estimate_richness function from phyloseq
# You could easily calculate Shannon, Chao1, Observed self, see alphaDiversityMeasures.Rmd
# estimate_richness itself uses functions from the vegan package

calc_alphadiv_plusLmResids <- function(physeq, measures = c("Observed", "Shannon"), group_var, compare) {
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the phyloseq object.")
        }
        
        if (!all(compare %in% unique(sample_data(physeq)[[group_var]]))) {
                stop("Not all given compare levels are actually levels in the group_var column.")
        }
        
        
        DF_alpha <- suppressWarnings(phyloseq::estimate_richness(physeq, measures = measures))
        
        DF_alpha$Total <- sample_sums(physeq)
        
        rownames(DF_alpha) <- sample_names(physeq)
        
        # - restrict to the samples defined by group_var and compare -
        DF_alpha <- DF_alpha[sample_data(physeq)[[group_var]] %in% compare,]
        
        # because linear fits of alpha diversity measures to total counts are often highly significant, I add the residuals of these
        # linear fits. 

        fitlist <- list()
        ncol_df_alpha <- ncol(DF_alpha)
        
        for (i in 1:length(measures)){
                fit <- lm(DF_alpha[,measures[i]] ~ DF_alpha[,"Total"])
                fitlist[[i]] <- fit
                DF_alpha[, ncol_df_alpha + i] <- residuals(fit)
                colnames(DF_alpha)[ncol_df_alpha + i] <- paste(measures[i], "resids", sep = "_")
        }
        
        names(fitlist) <- measures
        
        
        DF_alpha <- data.frame(DF_alpha, sample_data(physeq)[sample_data(physeq)[[group_var]] %in% compare,])
        DF_alpha[[group_var]] <- factor(DF_alpha[[group_var]], levels = compare, ordered = T)
        outlist <- list(DF_alpha = DF_alpha, fitlist = fitlist)
}
# --




# --
#######################################
### calc_pVals_alphdiv
#######################################
# Output:
# a DF with the p-values of all pairwise comparisons between the levels in group_var for the different alpha diversity measures
# NB: uses ggpubr: compare_means

calc_pVals_alphdiv <- function(DF_alpha, measures, group_var, compare, test = "t.test", 
                               symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                               p.adjust.method = "BH", paired = FALSE){
        
        if(! group_var %in% colnames(DF_alpha)) {
                stop("The given group_var is not a variable in the DF_alpha.")
        }
        
        if (!all(compare %in% unique(DF_alpha[[group_var]]))) {
                stop("Not all given compare levels are actually levels in the group_var column.")
        }
        
        DF_alpha <- DF_alpha[DF_alpha[[group_var]] %in% compare,]
        DF_alpha[[group_var]] <- factor(DF_alpha[[group_var]], levels = compare, ordered = T)
        
        y_columns <- which(colnames(DF_alpha) %in% c(measures, paste0(measures, "_resids", sep = "")))
        ttestList <- list()
        for (i in 1:length(y_columns)){
                comparison <- as.formula(paste(colnames(DF_alpha)[y_columns[i]], " ~ ", group_var, sep = ""))
                ttestList[[i]] <- ggpubr::compare_means(formula = comparison, data = DF_alpha, method = test, p.adjust.method = p.adjust.method, symnum.args = symnum.args, paired = paired)
        }
        alpha_div_ttests <- do.call("rbind", ttestList) %>% arrange(.y.)
}
# --




# --
#######################################
### boxplots_alphdiv
#######################################
# Output:
# a list of boxplots for the alpha diversity measures given (including resid plots), the plots indicate the p-values between all levels in the group variable

# NB: makes use of ggpubr::stat_compare_means that allows you to adjust more settings, e.g. plotting symbols instead of numbers (change label to "p.signif")

boxplots_alphdiv <- function(DF_alpha, measures, group_var, shape, color_levels, test = "t.test", symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), hide.ns = FALSE, comparisonList ){
        
        if(! group_var %in% colnames(DF_alpha)) {
                stop("The given group_var is not a variable in the DF_alpha.")
        }
        
        if (!all(names(color_levels) %in% unique(DF_alpha[[group_var]]))) {
                stop("Not all given color levels are actually levels in the group_var column.")
        }
        
        DF_alpha <- DF_alpha[DF_alpha[[group_var]] %in% names(color_levels),]
        DF_alpha[[group_var]] <- factor(DF_alpha[[group_var]], levels = names(color_levels), ordered = T)
        
        boxplotList <- list()
        fac_levels <- levels(DF_alpha[[group_var]])
        if(is.null(comparisonList)){
        comparisonList <- get_unique_facLevel_combinations(fac_levels)}else{ comparisonList = comparisonList}
       
        y_columns <- which(colnames(DF_alpha) %in% c(measures, paste0(measures, "_resids", sep = "")))
        
        for (i in 1:length(y_columns)){
                aes_map <- aes_string(x = group_var, y = colnames(DF_alpha)[y_columns[i]], color = group_var)
                Tr <-  ggplot(DF_alpha, aes_map) + 
                        geom_boxplot(na.rm = TRUE, outlier.color = NA) +
                        geom_jitter(aes_string(shape = shape), width = .2, height = 0, alpha = 0.65) +
                        xlab("") +
                        scale_color_manual("", values = color_levels) +
                        theme_bw()
                Tr <- Tr + stat_compare_means(comparisons = comparisonList, method = test, label = "p.signif", symnum.args = symnum.args, hide.ns = hide.ns) # "p.signif"
                boxplotList[[i]] <- Tr
                names(boxplotList)[i] <- colnames(DF_alpha)[y_columns[i]]
        }
        boxplotList
}
# --




# --
#######################################
### lmPlots_alphdiv
#######################################
# Output:
# a list of linear fit plots for the alpha diversity measures given against Total = sample_sums()

# NB: requires the lmp function!

lmPlots_alphdiv <- function(DF_alpha, lm_fitlist, measures, group_var, shape, color_levels, test = "t.test"){
        
        if(! group_var %in% colnames(DF_alpha)) {
                stop("The given group_var is not a variable in DF_alpha.")
        }
        
        if(!all(measures %in% colnames(DF_alpha))) {
                stop("The given measures are not all a variable in DF_alpha.")
        }
        
        if (!all(names(color_levels) %in% unique(DF_alpha[[group_var]]))) {
                stop("Not all given color levels are actually levels in the group_var column.")
        }
        
        DF_alpha <- DF_alpha[DF_alpha[[group_var]] %in% names(color_levels),]
        DF_alpha[[group_var]] <- factor(DF_alpha[[group_var]], levels = names(color_levels), ordered = T)
        
        
        lm_plotList <- list()
        for (i in 1:length(measures)){
                aes_map <- aes_string(x = "Total", y = measures[i], color = group_var, shape = shape)
                Tr <-  ggplot(DF_alpha, aes_map) + 
                        geom_point(na.rm = TRUE) +
                        xlab("total counts (sample_sums())") +
                        scale_color_manual("", values = color_levels) +
                        theme_bw()
                fit <- lm_fitlist[[measures[i]]]
                adjR2 <- round(summary(fit)$adj.r.squared,3)
                if (adjR2 != 0){
                        pfit <- lmp(fit)
                } else {
                        pfit <- NA
                }
                Tr <- Tr + geom_smooth(aes_string(x = "Total", y = measures[i]), method = "lm", se = TRUE, inherit.aes = F)
                
                Tr <- Tr + ggtitle(paste("lm fit: p-value: ", format(pfit, digits = 4), ", adj.r2: ", adjR2, sep = ""))
                lm_plotList[[i]] <- Tr
                names(lm_plotList)[i] <- measures[i]
                
        }
        lm_plotList
}
# --




# --
#######################################
### raref_curve_richness
#######################################
# read rarefaction_curve_own, this fast version uses the vegan::rarefy function, so it only generates rarefaction
# for richness. NB: vegan::rarefy used here does averaging (read help)


raref_curve_richness <- function(physeq, group_var = NULL, max_total = NULL, step_size = 200, color_levels, seed = 123) {
        
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the phyloseq object.")
        }
        
        if (!all(names(color_levels) %in% unique(sample_data(physeq)[[group_var]]))) {
                stop("Not all given color levels are actually levels in the group_var column.")
        }
        
        keepSamples <- sample_names(physeq)[sample_data(physeq)[[group_var]] %in% names(color_levels)]
        physeq <- prune_samples(keepSamples, physeq)
        sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], levels = names(color_levels), ordered = T)

        if (taxa_are_rows(physeq)) {
                seqtab <- t(as(otu_table(physeq), "matrix"))
        } else {
                seqtab <- as(otu_table(physeq), "matrix")
        }
        
        
        if (is.null(max_total)) {
                max_total <- quantile(sample_sums(physeq), probs = .25)
        }
        
        steps <- seq(from = 0, to = max_total, by = step_size)
        NoSamples <- nsamples(physeq)
        NoSteps <- length(steps)
        richness_matrix <- matrix(data = NA, nrow = NoSamples, ncol = NoSteps)
        totalamplicons <- sample_sums(physeq)
        
        set.seed(seed)
        
        # ptm <- Sys.time()
        for (i in 1:length(steps)) {
                richness_matrix[,i] <- vegan::rarefy(seqtab, sample = steps[i], se = FALSE)
        }
        # Sys.time() - ptm
        
        # set richness for samples at which steps > totalamplicons to NA
        for (i in 1:NoSamples) {
                NaIndex <- which(totalamplicons[i] < steps)[1]
                if (!is.na(NaIndex)){
                        richness_matrix[i, NaIndex:ncol(richness_matrix)] <- NA
                }
        }
        
        rownames(richness_matrix) <- rownames(seqtab)
        colnames(richness_matrix) <- paste("step_", steps, sep = "")
        richness_df <- as.data.frame(richness_matrix)
        
        
        plot_div_df3 <- function (div_df, type = "richness") {
                
                div_df$Sample <- rownames(div_df)
                div_df$Total <- totalamplicons
                div_df <- tidyr::gather(div_df, key = step, value = div, -Sample, -Total)
                div_df <- tidyr::separate(div_df, col = step, into = c("term", "step"), sep = "_")
                div_df$step <- as.numeric(div_df$step)
                
                Tr <- ggplot(div_df, aes(x = step, y = div, group = Sample, col = Total))
                Tr <- Tr +
                        geom_line() +
                        scale_color_gradient2("total counts bef.", low = cbPalette[6], mid = cbPalette[1], high = cbPalette[2], midpoint = median(div_df$Total)) +
                        xlab("total counts in sample") +
                        ylab(type) +
                        theme_bw()
        }
        
        Tr_richness_grad <- plot_div_df3(richness_df, type = "richness")
        
        
        
        if (!is.null(group_var)){
                Group <- sample_data(physeq)[[group_var]] 
                if (!is.null(Group) && !is.factor(Group)) {
                        Group <- as.factor(Group)
                }
        } else { Group <- NULL}
        
        
        if (!is.null(Group)) {
                
                richness_df$Group <- Group
                
                # pairwise.tt_richness <- lapply(richness_df[, -ncol(richness_df)], function(step){ 
                #         ptt <- pairwise.t.test(x = step, g = richness_df$Group, alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
                #         ptt$p.value
                # })
                
                pairwise.tt_richness <- NA # NB: I currently do not even show these p-values and it is currently a source of error as soon as you do not
                # have enough samples within a group at a certain step for a t.test (i.e. only 1 sample in each group left.)
                
                plot_div_df2 <- function (div_df, type = "richness") {
                        
                        div_df$Sample <- rownames(div_df)
                        div_df <- tidyr::gather(div_df, key = step, value = div, -Sample, -Group)
                        div_df <- tidyr::separate(div_df, col = step, into = c("term", "step"), sep = "_")
                        div_df$step <- as.numeric(div_df$step)
                        
                        Tr <- ggplot(div_df, aes(x = step, y = div, group = Sample, col = Group))
                        Tr <- Tr +
                                geom_line() +
                                xlab("total counts in sample") +
                                scale_color_manual("", values = color_levels) +
                                ylab(type) +
                                theme_bw(12)
                }
                
                plot_div_df_group <- function (div_df, type = "alpha diversity") {
                        div_df <- tidyr::gather(div_df, key = step, value = div, - Group)
                        div_df <- tidyr::separate(div_df, col = step, into = c("term", "step"), sep = "_")
                        div_df$step <- as.numeric(div_df$step)
                        div_df <- group_by(div_df, Group, step)
                        div_df <- dplyr::summarise(div_df, Mean_div = mean(div, na.rm = T), SD_div = sd(div, na.rm = T), n = n(), SE_div = SD_div/sqrt(n))
                        
                        Tr <- ggplot(div_df, aes(x = step, y = Mean_div, col = Group))
                        Tr <- Tr + 
                                geom_line() +
                                geom_point(size = 1) +
                                geom_errorbar(aes(ymin = Mean_div-SE_div, ymax = Mean_div+SE_div)) +
                                scale_color_manual("", values = color_levels) +
                                xlab("total counts in sample") +
                                ylab(paste(type, "(group mean +- SE)")) +
                                theme_bw(12)
                }
                Tr_richness_col <- plot_div_df2(richness_df, type = "richness")
                
                Tr_richness_group <- plot_div_df_group(richness_df, type = "richness")
                
        } else {
                Tr_richness_col <- NA
                Tr_richness_group <- NA
                pairwise.tt_richness <- NA
        }
        
        outlist <- list(rarefaction_richness_df = richness_df, Tr_richness_grad = Tr_richness_grad, Tr_richness_col = Tr_richness_col,
                        Tr_richness_group = Tr_richness_group, pairwise.tt_richness = pairwise.tt_richness)
        
}
# --





# --
#######################################
### calc_breakaway_from_phyloseq
#######################################
## OUTPUT:
# breakawayResults: a data frame with the breakaway estimates and standard errors cbind to sample_data(physeq) with
# sample_data(physeq)[[group_var]] set to a factor with names(color_levels) as ordered levels.


calc_breakaway_from_phyloseq <- function(physeq, group_var, color_levels) {
        
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the phyloseq object.")
        }
        
        if (!all(names(color_levels) %in% unique(sample_data(physeq)[[group_var]]))) {
                stop("Not all names in color_levels are actually levels in the group_var column.")
        }
        
        
        keepSamples <- sample_names(physeq)[sample_data(physeq)[[group_var]] %in% names(color_levels)]
        physeq <- prune_samples(keepSamples, physeq)
        sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], levels = names(color_levels), ordered = T)
        
        # - get the count table with taxa as columns -
        if (taxa_are_rows(physeq)) {
                CT <- t(as(otu_table(physeq), "matrix"))
        } else {
                CT <- as(otu_table(physeq), "matrix")
        }
        # --
        
        # - get the frequency distributions for all samples -
        SampleFrequencies <- apply(CT, 1, make_frequency_count_table)
        
        # -- make sure there is at least one singleton in each sample, i.e. add one if not --
        # this is because breakaway does not work without singleton
        SampleFrequencies <- lapply(SampleFrequencies, function(freqDF){
                if (freqDF[1,1] != 1) {
                        freqDF <- rbind(data.frame(Var1 = 1, Freq = 1), freqDF)
                }
                freqDF
        })
        
        # ----
        # --
        
        # - run breakaway on each frequency distributions -
        estimate_richness_breakaway <- function(my_freq) {
                result <- try(breakaway(my_freq, output=F, plot = F, answers = T), silent = T)
                if (class(result) == "list")  {
                        df <- data.frame(name = result$name, observed = sum(my_freq$Freq), est = result$est, seest = result$seest)
                } else {
                        df <- data.frame(name = NA, observed = sum(my_freq$Freq), est = NA, seest = NA)
                }
        }
        
        resultList <- lapply(SampleFrequencies, estimate_richness_breakaway)
        breakawayResults <- do.call("rbind", resultList)
        breakawayResults <- cbind(breakawayResults, sample_data(physeq))
        # --
        
        breakawayResults
}
# --
        



# --
#######################################
### analyse_breakaway_results
#######################################
# generates a boxplot of the breakaway results with all levels defined by color_levels in, then uses
# betta function from breakaway to analyse the results pairwise for all levels defined by names(color_levels)
## OUTPUT
# list: bettaResultList (a list of the betta results of all pairwise comparisons), boxPlot: a simple boxplot of the alpha diversity estimates



analyse_breakaway_results <- function(breakawayResults, group_var, color_levels, test = "t.test", 
                                      symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), 
                                      hide.ns = FALSE) {
        
        if(! group_var %in% colnames(breakawayResults)) {
                stop("The given group_var is not a variable in breakawayResults.")
        }
        
        if (!all(names(color_levels) %in% unique(breakawayResults[[group_var]]))) {
                stop("Not all names in color_levels are actually levels in the group_var column.")
        }
        
        breakawayResults <- breakawayResults[breakawayResults[[group_var]] %in% names(color_levels),]
        
        # - make a simple boxplot ignoring seest -
        fac_levels <- levels(breakawayResults[[group_var]])
        comparisonList <- get_unique_facLevel_combinations(fac_levels)
        
        
        Tr <-  ggplot(breakawayResults, aes_string(x = group_var, y = "est", color = group_var)) + 
                geom_boxplot(na.rm = TRUE, outlier.colour = NA) +
                geom_jitter(aes_string(shape = shape), width = .2, height = 0, alpha = 0.65) +
                xlab("") +
                ylab("breakaway richness estimate") +
                scale_color_manual("", values = color_levels) +
                theme_bw()
        Tr <- Tr + stat_compare_means(comparisons = comparisonList, method = test, label = "p.format", symnum.args = symnum.args, hide.ns = hide.ns) # "p.signif"
        # -- 
        
        # - run betta test on all pairwise combinations in names(color_levels) -
        # see ?betta and _UnderstandBreakaway.Rmd You need a covariance matrix for the test
        bettaResultList <- list()
        
        for (i in 1:length(comparisonList)){
                current_fac_levels <- comparisonList[[i]]
                current_breakawayResults <- breakawayResults[breakawayResults[[group_var]] %in% current_fac_levels, ]
                current_breakawayResults[[group_var]] <- factor(current_breakawayResults[[group_var]], levels = current_fac_levels, ordered = T)
                current_fac_levels_num <- (setNames(seq_along(current_fac_levels), current_fac_levels))-1
                X <- cbind("Intercept" = 1, "group_var" = current_fac_levels_num[current_breakawayResults[[group_var]]])
                testResult <- breakaway::betta(chats = current_breakawayResults$est, ses = current_breakawayResults$seest, X = X)
                resultDF <- data.frame(pval = testResult$table[2,3], estimate = testResult$table[2,1], se = testResult$table[2,2],
                                       intercept = testResult$table[1,1])
                resultDF$direction <- current_fac_levels[1]
                if (resultDF$estimate > 0) {resultDF$direction <- current_fac_levels[2]}
                resultDF$comparison <- paste(current_fac_levels, collapse = " vs ")
                symnum.args$x <- resultDF$pval
                resultDF$signi <- do.call("symnum", symnum.args)
                resultDF <- dplyr::select(resultDF, comparison, pval, signi, direction, estimate, se, intercept)
                bettaResultList[[i]] <- resultDF
        }
        
        list("bettaResultList" = bettaResultList, "boxPlot" = Tr)
}
# --

