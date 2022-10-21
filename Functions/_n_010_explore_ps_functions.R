# --
#######################################
### calc_SFs
#######################################
# implementation of DESeq2 library size similar to estimateSizeFactorsForMatrix
# in difference to adjust_LS (obsolete) ignore.zero.ratios is always TRUE (so 0 ratios are always ignored)


calc_SFs <- function(physeq, zeros.count = FALSE, percentile = 50)  {
        
        # ---- Step 1: Calculate Geometric mean for each taxa over all samples ------
        # NB: these GM is basically the reference sample
        if(taxa_are_rows(physeq)){
                GM <- apply(otu_table(physeq), 1, gm_own, zeros.count = zeros.count)   
        } else {
                GM <- apply(otu_table(physeq), 2, gm_own, zeros.count = zeros.count) 
        }
        
        # ---- Step 2: Calculate Size factors --------
        
        # NB: x/y = exp(log(x) - log(y))
        
        if (taxa_are_rows(physeq)) {
                SFs <- apply(as(otu_table(physeq), "matrix"), 2, function(sample_cnts){exp(quantile((log(sample_cnts)-log(GM))[sample_cnts > 0], probs = percentile/100, na.rm = T))})
                SFs <- SFs/exp(mean(log(SFs)))
        } else {
                SFs <- apply(as(otu_table(physeq), "matrix"), 1, function(sample_cnts){exp(quantile((log(sample_cnts)-log(GM))[sample_cnts > 0], probs = percentile/100, na.rm = T))})
                SFs <- SFs/exp(mean(log(SFs)))
        } 
        
        
        if(min(SFs) == 0) { warning("in at least one sample the Size Factor was 0!") }
        
        SFs
}
# --

calc_SFs2 <- function(physeq, zeros.count = FALSE, percentile = 50)  {
        
        # ---- Step 1: Calculate Geometric mean for each taxa over all samples ------
        # NB: these GM is basically the reference sample
        ps_ra<-transform_sample_counts(physeq, function(x){x/sum(x)}) 
        RA<-as.data.frame(unclass(otu_table(ps_ra)))
        #---- Step 2: Calculate Size factors --------
        
        # NB: x/y = exp(log(x) - log(y))
        
        if (taxa_are_rows(physeq)) {
                SFs <- apply(as(otu_table(physeq), "matrix"), 2, function(sample_cnts){exp(quantile((log(sample_cnts)-log(RA))[sample_cnts > 0], probs = percentile/100, na.rm = T))})
                SFs <- SFs/exp(mean(log(SFs)))
        } else {
                SFs <- apply(as(otu_table(ps_ra), "matrix"), 1, function(sample_cnts){exp(quantile((log(sample_cnts)-log(RA))[sample_cnts > 0], probs = percentile/100, na.rm = T))})
                SFs <- SFs/exp(mean(log(SFs)))
        } 
        
        
        if(min(SFs) == 0) { warning("in at least one sample the Size Factor was 0!") }
        
        SFs
}



# --
#######################################
### simply_adjust_LS
#######################################
# implementation of DESeq2 library size similar to estimateSizeFactorsForMatrix
# in difference to adjust_LS (obsolete) ignore.zero.ratios is always TRUE (so 0 ratios are always ignored)

## Input:
# physeq
# zeros.count: zeros.count of gm_own, if TRUE 0 will be considered when calculating the geometric mean
# if FALSE not and thus the geometric means will be bigger (see gm_own)
# percentile: the percentile to select the size factor SF of a sample based on its count ratios to the reference sample. In
# DESeq percentile = 50, i.e. stats::median is used. 
# plots: if TRUE SFs and plots will be given in an output list, otherwise just the adjusted physeq is returned


simply_adjust_LS <- function(physeq, SFs = NULL, zeros.count = FALSE, percentile = 50)  {
        
        # ---- Step 1: Calculate Geometric mean for each taxa over all samples ------
        # NB: these GM is basically the reference sample
        if(taxa_are_rows(physeq)){
                GM <- apply(otu_table(physeq), 1, gm_own, zeros.count = zeros.count)   
        } else {
                GM <- apply(otu_table(physeq), 2, gm_own, zeros.count = zeros.count) 
        }
        
        # ---- Step 2: Calculate Size factors  unless given --------
        
        # NB: x/y = exp(log(x) - log(y))
        
        if (is.null(SFs)){
                
                if (taxa_are_rows(physeq)) {
                        SFs <- apply(as(otu_table(physeq), "matrix"), 2, function(sample_cnts){exp(quantile((log(sample_cnts)-log(GM))[sample_cnts > 0], probs = percentile/100, na.rm = T))})
                        SFs <- SFs/exp(mean(log(SFs)))
                } else {
                        SFs <- apply(as(otu_table(physeq), "matrix"), 1, function(sample_cnts){exp(quantile((log(sample_cnts)-log(GM))[sample_cnts > 0], probs = percentile/100, na.rm = T))})
                        SFs <- SFs/exp(mean(log(SFs)))
                }
                
        }
        
        
        if(min(SFs) == 0) { warning("in at least one sample the Size Factor was 0!") }
        
        
        # --- 3: calculate the new counts and put into a physeq object
        
        if(taxa_are_rows(physeq)){
                if (!identical(colnames(otu_table(physeq)), names(SFs))) {stop("names SFs do not fit to physeq")}
                Mat <- sweep(otu_table(physeq),2,SFs, "/")
                phynew <- physeq
                otu_table(phynew) <- otu_table(Mat, taxa_are_rows = TRUE)
        } else {
                if (!identical(rownames(otu_table(physeq)), names(SFs))) {stop("names SFs do not fit to physeq")}
                Mat <- sweep(otu_table(physeq),1,SFs, "/") 
                phynew <- physeq
                otu_table(phynew) <- otu_table(Mat, taxa_are_rows = FALSE)
        }
        
        
        
        list(physeq = phynew, SFs = SFs)
        
}
# --





# --
#######################################
### FUNCTION: check_phyla_distribution#
#######################################
# NB: throws error if "Phylum" is not in colnames(tax_table(physeq))
# outputs a data.frame, summarising the taxa distribution over the different phyla.
## INPUT:
# physeq: physeq object
## OUTPUT:
# data.frame, summarising the taxa distribution over the different phyla.
# the Phyla are ordered so the Phylum with most counts (PC_of_counts) is on top, no of taxa is used to break ties in the ordering
# columns: PC stands for percentage, 
# mean/median_taxa_sum are the mean/median of the taxa_sums (total counts over all samples) of the taxa in the respective phylum
# other columns should be clear

check_phyla_distribution <- function(physeq) {
        
        
        if (phyloseq::taxa_are_rows(physeq)) {
                physeq <- t(physeq)
        }
        
        df_ab_prev <- data.frame(Taxon_No = 1:ntaxa(physeq), 
                                 total_counts = taxa_sums(physeq),
                                 prevalence = colSums(as(unclass(otu_table(physeq)), "matrix") != 0))
        
        
        df_ab_prev <- cbind(df_ab_prev, unclass(tax_table(physeq)))
        
        PhylaDistribution <- dplyr::summarise(group_by(df_ab_prev, Phylum), 
                                              taxa = n(), 
                                              PC_of_taxa = round(100*taxa/ntaxa(physeq),1),
                                              PC_of_counts = round(100*sum(total_counts)/sum(otu_table(physeq)), 1),
                                              PC_of_prevalence = round(100*sum(prevalence)/sum(otu_table(physeq) != 0), 1),
                                              mean_taxa_sum = round(mean(total_counts)),
                                              median_taxa_sum = round(median(total_counts)),
                                              mean_prevalence_in_PC = round(100*mean(prevalence)/nsamples(ps), 1)) %>% 
                arrange(desc(PC_of_counts), desc(taxa), desc(PC_of_prevalence))
        
        PhylaDistribution
        
}
# --



# --
#######################################
### FUNCTION: check_phyla_distribution_NA#
#######################################
# NB: throws error if "Phylum" is not in colnames(tax_table(physeq))
# outputs a data.frame, summarising the taxa distribution over the different phyla.
## INPUT:
# physeq: physeq object
## OUTPUT:
# data.frame, summarising the taxa distribution over the different phyla.
# the Phyla are ordered so the Phylum with most counts (PC_of_counts) is on top, no of taxa is used to break ties in the ordering
# columns: PC stands for percentage, 
# mean/median_taxa_sum are the mean/median of the taxa_sums (total counts over all samples) of the taxa in the respective phylum
# other columns should be clear

check_phyla_distribution_NA <- function(physeq) {
        
        
        if (phyloseq::taxa_are_rows(physeq)) {
                physeq <- t(physeq)
        }
        
        df_ab_prev <- data.frame(Taxon_No = 1:ntaxa(physeq), 
                                 total_counts = taxa_sums(physeq),
                                 prevalence = colSums(as(otu_table(physeq), "matrix") != 0))
        
        
        df_ab_prev <- cbind(df_ab_prev, tax_table(physeq))
        df_ab_prev$Phylum <- as.character(df_ab_prev$Phylum)
        df_ab_prev$Phylum[is.na(df_ab_prev$Phylum)] <- paste(df_ab_prev$Kingdom[is.na(df_ab_prev$Phylum)], "_NA", sep = "")
        
        PhylaDistribution <- dplyr::summarise(group_by(df_ab_prev, Phylum), 
                                              taxa = n(), 
                                              PC_of_taxa = round(100*taxa/ntaxa(physeq),1),
                                              PC_of_counts = round(100*sum(total_counts)/sum(otu_table(physeq)), 1),
                                              PC_of_prevalence = round(100*sum(prevalence)/sum(otu_table(physeq) != 0), 1),
                                              mean_taxa_sum = round(mean(total_counts)),
                                              median_taxa_sum = round(median(total_counts)),
                                              mean_prevalence_in_PC = round(100*mean(prevalence)/nsamples(ps), 1)) %>% 
                arrange(desc(PC_of_counts), desc(taxa), desc(PC_of_prevalence))
        
        PhylaDistribution
        
}
# --






# --
#######################################
#### plot_sample_bars_compare
#######################################
# generates abundance barplots (see plot_bar_own) to compare ps to ps_tca, i.e. to see how SFs adjustment affects the abundances 

plot_sample_bars_compare <- function(physeq, physeq2, x = "Sample", y = "Abundance", group_var, color_levels, fill = NULL,
                                     color_sample_names = TRUE, col_vec = NULL, order_by_raw_counts = TRUE){
        
        # - prepare mdf of ps physeq -
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        
        # in case you do not want to see all samples
        if (!all(unique(sample_data(physeq)[[group_var]]) %in% names(color_levels))) {
                keepSamples <- sample_names(physeq)[sample_data(physeq)[[group_var]] %in% names(color_levels)]
                physeq <- prune_samples(samples = keepSamples, physeq)
        }
        sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], levels = names(color_levels), order = TRUE)
        
        # if (!is.factor(sample_data(physeq)[[group_var]])) {sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], levels = names(color_levels), order = TRUE)}
        
        if (is.null(fill)) { fill = "Phylum"}
        
        mdf <- phyloseq::psmelt(physeq)
        
        # order fill levels according to abundance over all samples
        mdf[, fill] <- as.character(mdf[, fill])
        mdf[is.na(mdf[, fill]), fill] <- "NA" # NB: pools all NA
        sums <- group_by_(mdf, fill) %>% summarise(sum_abundance = sum(Abundance)) %>% arrange(sum_abundance)
        mdf[, fill] <- factor(mdf[, fill], levels = as.character(sums[[1]]), ordered = TRUE)
        
        # --
        
        # - prepare mdf of ps_tca physeq -
        if(taxa_are_rows(physeq2)) { physeq2 <- t(physeq2) }
        
        if (!all(unique(sample_data(physeq2)[[group_var]]) %in% names(color_levels))) {
                keepSamples <- sample_names(physeq2)[sample_data(physeq2)[[group_var]] %in% names(color_levels)]
                physeq2 <- prune_samples(samples = keepSamples, physeq2)
                sample_data(physeq2)[[group_var]] <- factor(sample_data(physeq2)[[group_var]], levels = names(color_levels), order = TRUE)
        }
        
        if (!is.factor(sample_data(physeq2)[[group_var]])) {sample_data(physeq2)[[group_var]] <- factor(sample_data(physeq2)[[group_var]], levels = names(color_levels), order = TRUE)}
        
        
        mdf2 <- phyloseq::psmelt(physeq2)
        
        # order fill levels according to abundance over all samples
        mdf2[, fill] <- as.character(mdf2[, fill])
        mdf2[is.na(mdf2[, fill]), fill] <- "NA"
        # sums <- group_by_(mdf2, fill) %>% summarise(sum_abundance = sum(Abundance)) %>% arrange(sum_abundance)
        mdf2[, fill] <- factor(mdf2[, fill], levels = as.character(sums[[1]]), ordered = TRUE)
        # --
        
        mdf$Typer <- "before"
        mdf2$Typer <- "after SF adjustment"
        
        mdf <- rbind(mdf, mdf2)
        mdf$Typer <- factor(mdf$Typer, levels = c("before", "after SF adjustment"), ordered = TRUE)
        
        
        # order samples according to group_var levels and potentially by total counts
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        
        if (order_by_raw_counts) {
                rawSampleSizes <- data.frame("Sample" = sample_names(physeq), "Group" = sample_data(physeq)[[group_var]], "Total" = sample_sums(physeq))
                rawSampleSizes <- arrange(rawSampleSizes, Group, Total)
                mdf$Sample <- factor(mdf$Sample, levels = rawSampleSizes$Sample, ordered = TRUE)
                
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
                scale_fill_manual(values = fill_colors, na.value = "red") +
                xlab("") +
                facet_wrap(~ Typer, ncol = 1) +
                theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, colour = colxaxis))
        
        Tr
}
# --






# --
#######################################
### plot_sizeFactors
#######################################
# just to check if size factors differ between groups and how they are related with original sample sizes

plot_sizeFactors <- function(physeq, SFs, group_var, color_levels, shape, test = "t.test", p_adjust_method = "fdr",
                             symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                             hide.ns = FALSE){
        
        
        # in case you do not want to see all samples
        if (!all(unique(sample_data(physeq)[[group_var]]) %in% names(color_levels))) {
                keepSamples <- sample_names(physeq)[sample_data(physeq)[[group_var]] %in% names(color_levels)]
                physeq <- prune_samples(samples = keepSamples, physeq)
                sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], levels = names(color_levels), order = TRUE)
                SFs <- SFs[names(SFs) %in% keepSamples]
        }
        
        DF <- cbind(sample_data(physeq), SF = SFs)
        
        Tr <-  ggplot(DF, aes_string(x = group_var, y = "SF", color = group_var)) 
        Tr <- Tr + 
                geom_boxplot(outlier.color = NA) +
                geom_jitter(aes_string(shape = shape), width = .2, height = 0, alpha = 0.65) +
                xlab("") +
                ylab("size factor") +
                scale_color_manual("", values = color_levels) +
                theme_bw()
        if(is.null(shape)){
                Tr <- Tr + theme(legend.position = "none")
        }
        
        # - since you might have more than two levels in each plot you need to set the comparisons argument in stat_compare_means -
        group_fac <- factor(sample_data(physeq)[[group_var]])
        fac_levels <- levels(group_fac)
        
        comparisonList <- get_unique_facLevel_combinations(fac_levels)
        
        Tr <- Tr + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = test, hide.ns = hide.ns)
        
        formulaa <- as.formula(paste("SF ~", group_var, sep = " "))
        
        pVals <- compare_means(formula = formulaa, data = DF, method = test, p.adjust.method = p_adjust_method, symnum.args = symnum.args)
        
        # - Add total_counts vs SFs plot -
        DF$total_count <- sample_sums(physeq)
        Tr1 <-  ggplot(DF, aes_string(x = "total_count", y = "SF", color = group_var)) 
        Tr1 <- Tr1 + 
                geom_point(aes_string(shape = shape), alpha = 0.65) +
                xlab("sample_sums()") +
                ylab("size factor") +
                scale_color_manual("", values = color_levels) +
                theme_bw()
        if(is.null(shape)){
                Tr1 <- Tr1 + theme(legend.position = "none")
        }
        
        
        
        list(pVals = pVals, Tr = Tr, Tr1 = Tr1)
}
# --




# --
#######################################
### FUNCTION: get_assignemnt_distribution
#######################################
## INPUT:
# physeq
## OUTPUT:
# df: illustrating the percentage of taxa that have been assigned to the different taxonomic levels. Non assigned levels must be NA

get_assignemnt_distribution <- function(physeq){
        
        taxa <- tax_table(physeq)
        
        countNA <- apply(taxa, 2, function(x){sum(is.na(x))})
        countNonNA <- apply(taxa, 2, function(x){sum(!is.na(x))})
        # ambiguous <- apply(taxa, 2, function(x){sum(grepl("/", x))})
        total <- countNA + countNonNA
        assignment_distribution <- data.frame(assigned = countNonNA, 
                                              # assigned_unamb = countNonNA - ambiguous,
                                              total = total,
                                              PC_assigned = round(100*countNonNA/total, 1))
                                              # PC_assigned_unamb = round(100*(countNonNA - ambiguous)/total, 1))
        
}
# --





# --
#######################################
### FUNCTION: check_assignment_vs_abundance
#######################################
# checks if there is a trend for better assignment for more abundant SVs

check_assignment_vs_abundance <- function(physeq, abundanceQuantiles = seq(0, 90, by = 10)){
        
        if (phyloseq::taxa_are_rows(physeq)) {
                physeq <- t(physeq)
        }
        
        taxa = tax_table(physeq) 
        seqtab = as(otu_table(physeq), "matrix")
        
        
        total_counts <- colSums(seqtab)
        abQuantiles <- quantile(total_counts, probs = abundanceQuantiles/100)
        remaining_SVs <- lapply(abQuantiles, function(quant) {
                total_counts >= quant
        })
        No_ASVs <- sapply(remaining_SVs, sum)
        filtered.taxas <- lapply(remaining_SVs, function(indexes){
                taxa[indexes,]
        })
        
        assignment_distributions <- lapply(filtered.taxas, get_assignemnt_distribution)
        
        PC_assigned <- sapply(assignment_distributions, function(distri){distri[["PC_assigned"]]})
        
        rownames(PC_assigned) <- colnames(taxa)
        colnames(PC_assigned) <- paste("Ab_", abundanceQuantiles, "_", round(abQuantiles), "_", No_ASVs, sep = "")
        
        PC_assigned <- as.data.frame(PC_assigned)
        PC_assigned$Level <- rownames(PC_assigned)
        df_long <- gather(PC_assigned, key = Ab, value = PC, -Level)
        df_long <- separate(df_long, col = Ab, into = c("Type", "Quant", "abQuant", "No_ASVs"), sep = "_")
        df_long$Level <- factor(df_long$Level, levels = colnames(taxa), ordered = TRUE)
        df_long$Quant <- factor(df_long$Quant, levels = as.character(abundanceQuantiles), ordered = TRUE)
        
        tr <- ggplot(df_long, aes(x = Quant, y = PC, col = Level))
        tr <- tr + 
                geom_point(size = 2) +
                scale_colour_manual("", values = cbPalette[2:8]) +
                scale_x_discrete(breaks = abundanceQuantiles, labels = paste(abundanceQuantiles, " (", No_ASVs, ")", sep = "")) +
                ylab("percentage of taxa assigned") +
                xlab("total counts quantile (No of remaining taxa)") +
                theme_bw() +
                theme(axis.text.x = element_text(angle=90, vjust=0.5))
        
        list(PC_assigned, tr)
        
}
# --




# --
#######################################
### FUNCTION: check_assignment_vs_prevalence
#######################################
# checks if there is a trend for better assignment for more prealent SVs

check_assignment_vs_prevalence <- function(physeq, prevalences = seq(0, 90, by = 10)){
        
        if (phyloseq::taxa_are_rows(physeq)) {
                physeq <- t(physeq)
        }
        
        taxa = tax_table(physeq) 
        seqtab = as(otu_table(physeq), "matrix")
        
        prev <- colSums(seqtab != 0)
        remaining_SVs <- lapply(prevalences, function(preva) {
                prev > (preva/100)*nrow(seqtab)
        })
        No_ASVs <- sapply(remaining_SVs, sum)
        filtered.taxas <- lapply(remaining_SVs, function(indexes){
                taxa[indexes,]
        })
        
        assignment_distributions <- lapply(filtered.taxas, get_assignemnt_distribution)
        
        PC_assigned <- sapply(assignment_distributions, function(distri){distri[["PC_assigned"]]})
        
        rownames(PC_assigned) <- colnames(taxa)
        colnames(PC_assigned) <- paste("Prev_", prevalences, "_", No_ASVs, sep = "")
        
        PC_assigned <- as.data.frame(PC_assigned)
        PC_assigned$Level <- rownames(PC_assigned)
        df_long <- gather(PC_assigned, key = Prev, value = PC, -Level)
        df_long <- separate(df_long, col = Prev, into = c("Type", "Prevalence", "No_ASVs"), sep = "_")
        df_long$Level <- factor(df_long$Level, levels = colnames(taxa), ordered = TRUE)
        df_long$Prevalence <- factor(df_long$Prevalence, levels = as.character(prevalences), ordered = TRUE)
        
        
        tr <- ggplot(df_long, aes(x = Prevalence, y = PC, col = Level))
        tr <- tr + 
                geom_point(size = 2) +
                scale_colour_manual("", values = cbPalette[2:8]) +
                scale_x_discrete(breaks = prevalences, labels = paste(prevalences, " (", No_ASVs, ")", sep = "")) +
                ylab("percentage of taxa assigned") +
                xlab("prevalence (No of remaining taxa)") +
                theme_bw() +
                theme(axis.text.x = element_text(angle=90, vjust=0.5))
        
        list(PC_assigned, tr)
        
}
# --




# --
#######################################
#### plot_correlations_abundance_prev_sparsity
#######################################
# physeq
# col = here a string of a taxonomic level

plot_correlations_abundance_prev_sparsity <- function(physeq, col = NULL, col_vec = NULL){ 
        
        
        if (phyloseq::taxa_are_rows(physeq)) {
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
        
        nsamples <- nsamples(physeq)
        
        if (!is.null(col)) {
                
                if (!is.null(col_vec)){
                        df_ab_prev[[col]] <- as.character(df_ab_prev[[col]])
                        df_ab_prev[[col]][is.na(df_ab_prev[[col]])] <- "NA" # NB: pools all NA
                        if (length(col_vec) != length(unique(df_ab_prev[[col]]))){
                                stop("provided col_vec did not fit in length to the col variable")
                        }
                        df_ab_prev[[col]] <- factor(df_ab_prev[[col]], levels = names(col_vec), ordered = TRUE)
                        custom_colors <- col_vec
                } else {
                        CountOrder <- dplyr::group_by_(df_ab_prev, col) %>% dplyr::summarise(total_count_sum = sum(total_counts)) %>% dplyr::arrange(desc(total_count_sum))
                        # - make sure that also NAs are plotted, i.e. where coloring taxonomic rank is NA -
                        CountOrder[[col]] <- as.character(CountOrder[[col]])
                        CountOrder[[col]][is.na(CountOrder[[col]])] <- "NA"
                        
                        if (nrow(CountOrder) <= 15){
                                custom_colors <- make_color_vector(CountOrder[[col]], QuantColors15)
                        } else {
                                custom_colors <- make_color_vector(CountOrder[[col]], viridis(nrow(CountOrder)))
                        }
                        
                        df_ab_prev[[col]] <- as.character(df_ab_prev[[col]])
                        df_ab_prev[[col]][is.na(df_ab_prev[[col]])] <- "NA"
                        df_ab_prev[[col]] <- factor(df_ab_prev[[col]], levels = names(custom_colors), ordered = TRUE)
                        
                }
                
        }
        
        
        
        Tr_ab <- ggplot(df_ab_prev, aes(x = Taxon_No, y = total_counts))
        if (is.null(col)) {
                Tr_ab <- Tr_ab + geom_point(col = cbPalette[2], size = 2, alpha = 0.95) 
        } else {
                Tr_ab <- Tr_ab + geom_point(aes_string(col = col), size = 2, alpha = 0.95)
                Tr_ab <- Tr_ab + scale_color_manual("", values = custom_colors)
                
        }
        Tr_ab <- Tr_ab +
                ylab("total counts of taxon (= taxa_sums())") +
                theme_bw()
        
        
        Tr_prev <- ggplot(df_ab_prev, aes(x = Taxon_No, y = prevalence))
        if (is.null(col)) {
                Tr_prev <- Tr_prev + geom_point(col = cbPalette[2], size = 2, alpha = 0.95) 
        } else {
                Tr_prev <- Tr_prev + geom_point(aes_string(col = col), size = 2, alpha = 0.95)
                Tr_prev <- Tr_prev + scale_color_manual(values = custom_colors)
        }
        Tr_prev <- Tr_prev +
                theme_bw()
        
        
        
        fit_prev_log10 <- lm(formula = prevalence ~ log10(total_counts), data = df_ab_prev)
        pval_prev_log10 <- lmp(fit_prev_log10)
    
        fit_prev_mean <- lm(formula = prevalence ~ mean_count_nonzero, data = df_ab_prev)
        pval_prev_mean <- lmp(fit_prev_mean)
        fit_prev_mean_log10 <- lm(formula = prevalence ~ log10(mean_count_nonzero), data = df_ab_prev)
        pval_prev_mean_log10 <- lmp(fit_prev_mean_log10)
        fit_prev_median <- lm(formula = prevalence ~ median_count_nonzero, data = df_ab_prev)
        pval_prev_median <- lmp(fit_prev_median)
        fit_prev_median_log10 <- lm(formula = prevalence ~ log10(median_count_nonzero), data = df_ab_prev)
        pval_prev_median_log10 <- lmp(fit_prev_median_log10)
        
        
        Tr_prev_vs_log10_ab <- ggplot(df_ab_prev, aes(x = total_counts, y = prevalence))
        if (is.null(col)) {
                Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab + geom_point(col = cbPalette[2], size = 2, alpha = 0.95) 
        } else {
                Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab + geom_point(aes_string(col = col), size = 2, alpha = 0.95)
                Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab + scale_color_manual("", values = custom_colors)
        }
        Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab +
                scale_x_log10() +
                scale_y_continuous(limits = c(-1, max(df_ab_prev$prevalence) + 5)) +
                geom_smooth(method = "lm") +
                # annotate("text", x = max(df_ab_prev$abundance)/5, y = 5, label = paste("R.square of lm: ", as.character(round(summary(fit_prev)$r.squared,4), sep = ""))) +
                xlab("total counts of taxon (= taxa_sums())") +
                ggtitle(paste("lm fit: p.val: ", format(pval_prev_log10, digits = 4), "R.square: ", as.character(round(summary(fit_prev_log10)$r.squared,4), sep = ""))) +
                theme_bw(12)
        
        
        
        Tr_prev_vs_log10_meanab <- ggplot(df_ab_prev, aes(x = mean_count_nonzero, y = prevalence))
        if (is.null(col)) {
                Tr_prev_vs_log10_meanab <- Tr_prev_vs_log10_meanab + geom_point(col = cbPalette[2], size = 2, alpha = 0.95) 
        } else {
                Tr_prev_vs_log10_meanab <- Tr_prev_vs_log10_meanab + geom_point(aes_string(col = col), size = 2, alpha = 0.95)
                Tr_prev_vs_log10_meanab <- Tr_prev_vs_log10_meanab + scale_color_manual("", values = custom_colors)
        }
        Tr_prev_vs_log10_meanab <- Tr_prev_vs_log10_meanab +
                scale_x_log10() +
                geom_smooth(method = "lm") +
                scale_y_continuous(limits = c(-1, max(df_ab_prev$prevalence) + 5)) +
                xlab("mean count of taxon in non-zero samples") +
                ggtitle(paste("lm fit: p.val: ", format(pval_prev_mean_log10, digits = 4), "R.square: ", as.character(round(summary(fit_prev_mean_log10)$r.squared,4), sep = ""))) +
                theme_bw(12)
        
        Tr_prev_vs_log10_medianab <- ggplot(df_ab_prev, aes(x = median_count_nonzero, y = prevalence))
        if (is.null(col)) {
                Tr_prev_vs_log10_medianab <- Tr_prev_vs_log10_medianab + geom_point(col = cbPalette[2], size = 2, alpha = 0.95) 
        } else {
                Tr_prev_vs_log10_medianab <- Tr_prev_vs_log10_medianab + geom_point(aes_string(col = col), size = 2, alpha = 0.95)
                Tr_prev_vs_log10_medianab <- Tr_prev_vs_log10_medianab + scale_color_manual("", values = custom_colors)
        }
        Tr_prev_vs_log10_medianab <- Tr_prev_vs_log10_medianab +
                scale_x_log10() +
                geom_smooth(method = "lm") +
                scale_y_continuous(limits = c(-1, max(df_ab_prev$prevalence) + 5)) +
                xlab("median count of taxon in non-zero samples") +
                ggtitle(paste("lm fit: p.val: ", format(pval_prev_median_log10, digits = 4), "R.square: ", as.character(round(summary(fit_prev_median_log10)$r.squared,4), sep = ""))) +
                theme_bw(12)
        
        
        
        
        fitlist <- list(fit_prev_log10 = fit_prev_log10, fit_prev_mean = fit_prev_mean, 
                        fit_prev_mean_log10 = fit_prev_mean_log10, fit_prev_median = fit_prev_median,
                        fit_prev_median_log10 = fit_prev_median_log10)
        
        out <- list(Tr_ab = Tr_ab, 
                    Tr_prev = Tr_prev, 
                    Tr_prev_vs_log10_ab = Tr_prev_vs_log10_ab, 
                    Tr_prev_vs_log10_meanab = Tr_prev_vs_log10_meanab,
                    Tr_prev_vs_log10_medianab = Tr_prev_vs_log10_medianab,
                    fitlist = fitlist)
        
}
# --




# --
#######################################
#### plot_ab_pev_distributions
#######################################


plot_ab_pev_distributions <- function(physeq, prevalence = 5) {
        
        if (phyloseq::taxa_are_rows(physeq)) {
                physeq <- t(physeq)
        }
        
        seqtab <- as(otu_table(physeq), "matrix")
        
        FinalNumbersSeq <- data.frame(Sequence = colnames(seqtab), InNumberSamples = colSums(seqtab != 0), TotalCounts = colSums(seqtab))
        # a look at it shows you that seqtab is ordered after total counts
        FinalNumbersSeq <- group_by(FinalNumbersSeq, InNumberSamples)
        CountDistribution <- dplyr::summarise(FinalNumbersSeq, No_ASVs = n(), TotalCounts = sum(TotalCounts))
        CountDistribution$CumSumUnique <- rev(cumsum(rev(CountDistribution$No_ASVs)))
        CountDistribution$CumPerCUnique <- rev(cumsum(rev(CountDistribution$No_ASVs/ncol(seqtab))))
        CountDistribution$CumSumTotal <- rev(cumsum(rev(CountDistribution$TotalCounts)))
        CountDistribution$CumPerCTotal <- rev(cumsum(rev(CountDistribution$TotalCounts/sum(colSums(seqtab)))))
        
        PCValue <- ceiling((prevalence/100)*dim(seqtab)[1]) # tells you in how many samples a SV must be present to meet the prevalence 
        
        # Diff <- CountDistribution$InNumberSamples - PCValue
        # index <- which.max(Diff[Diff<0]) + which.min(Diff[Diff>=0])
        index <- which(CountDistribution$InNumberSamples >= PCValue)[1]
        PCKeptAtPCValue <- CountDistribution$CumPerCTotal[index]
        SVskeptAtPCValue <- CountDistribution$CumPerCUnique[index]
        
        # The number of samples the SVs are present in
        Tr <- ggplot(CountDistribution, aes(x = InNumberSamples, y = No_ASVs))
        Tr <- Tr + geom_point(col = "#E69F00", size = 3) +
                xlab("prevalence") +
                ylab("number of taxa") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        
        In1Index <- which(CountDistribution$InNumberSamples == 1)
        if (length(In1Index) != 0) {
                Tr <- Tr + ggtitle(paste(CountDistribution$No_ASVs[In1Index], " of ", CountDistribution$CumSumUnique[1], " taxa (", round(100*CountDistribution$No_ASVs[In1Index]/CountDistribution$CumSumUnique[1], 1), " %)", " were only found in 1 sample", sep = ""))
        } 
        
        
        # Cumulative Percentage of SVs
        Tr1 <- ggplot(CountDistribution, aes(x = InNumberSamples, y = CumPerCUnique))
        Tr1 <- Tr1 + geom_point(col = "#E69F00", size = 3) +
                xlab("prevalence") +
                ylab("cumulative percentage of taxa") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        Tr1 <- Tr1 + 
                geom_hline(yintercept = SVskeptAtPCValue, lty =  "dashed") +
                geom_vline(xintercept = (prevalence/100)*dim(seqtab)[1], lty = 'dashed') +
                ggtitle(paste("prevalence ", prevalence, " % = ", round((prevalence/100)*dim(seqtab)[1],1), "; ", CountDistribution$CumSumUnique[index], " of ", CountDistribution$CumSumUnique[1], 
                              " taxa (", round(100*CountDistribution$CumSumUnique[index]/CountDistribution$CumSumUnique[1], 1), " %) have higher prevalence", sep = ""))
        
        
        Tr2 <- ggplot(CountDistribution, aes(x = InNumberSamples, y = TotalCounts))
        Tr2 <- Tr2 + geom_point(col = "#E69F00", size = 3) +
                xlab("prevalence") +
                ylab("total counts of taxa with given prevalence") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        Tr2 <- Tr2 + ggtitle(paste(CountDistribution$CumSumTotal[1] - CountDistribution$CumSumTotal[index], " of ",
                                   CountDistribution$CumSumTotal[1], " (", round(100*(CountDistribution$CumSumTotal[1] - CountDistribution$CumSumTotal[index])/CountDistribution$CumSumTotal[1], 2),
                                   " %) counts are from taxa present in less than ", round((prevalence/100)*dim(seqtab)[1],1), " samples.", sep = ""))
        
        Tr3 <- ggplot(CountDistribution, aes(x = InNumberSamples, y = CumPerCTotal))
        Tr3 <- Tr3 + geom_point(col = "#E69F00", size = 3) +
                xlab("prevalence") +
                ylab("cumulative percentage of counts") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_line(color = "#999999", size = .15),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        Tr3 <- Tr3 + 
                geom_hline(yintercept = PCKeptAtPCValue, lty =  "dashed") +
                geom_vline(xintercept = (prevalence/100)*dim(seqtab)[1], lty = 'dashed') +
                ggtitle(paste("prevalence ", prevalence, " % = ", round((prevalence/100)*dim(seqtab)[1],1), "; ", CountDistribution$CumSumTotal[index], " of ", CountDistribution$CumSumTotal[1], 
                              " counts (", round(100*CountDistribution$CumSumTotal[index]/CountDistribution$CumSumTotal[1], 1), " %) would remain", sep = ""))
        
        
        TrList <- list(Tr, Tr1, Tr2, Tr3, CountDistribution)
        
        TrList
        
}
# --



#######################################
#### plot_sample_bars_compare
#######################################
# generates abundance barplots (see plot_bar_own) to compared different ps_object 

plot_sample_bars_compare_2 <- function(physeq, physeq2, physeq3, x = "Sample", y = "Abundance", group_var, color_levels, fill = NULL,
                                     color_sample_names = TRUE, col_vec = NULL,  taxa_Level1=taxa_Level1, taxa_Level2= taxa_Level2, taxa_Level3= taxa_Level3){
        
        # - prepare mdf of ps ASV -
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        
        # in case you do not want to see all samples
        if (!all(unique(sample_data(physeq)[[group_var]]) %in% names(color_levels))) {
                keepSamples <- sample_names(physeq)[sample_data(physeq)[[group_var]] %in% names(color_levels)]
                physeq <- prune_samples(samples = keepSamples, physeq)
        }
        sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], levels = names(color_levels), order = TRUE)
        
        # if (!is.factor(sample_data(physeq)[[group_var]])) {sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], levels = names(color_levels), order = TRUE)}
        
        if (is.null(fill)) { fill = "Phylum"}
        
        mdf <- phyloseq::psmelt(physeq)
        
        # order fill levels according to abundance over all samples
        mdf[, fill] <- as.character(mdf[, fill])
        mdf[is.na(mdf[, fill]), fill] <- "NA" # NB: pools all NA
        sums <- group_by_(mdf, fill) %>% summarise(sum_abundance = sum(Abundance)) %>% arrange(sum_abundance)
        mdf[, fill] <- factor(mdf[, fill], levels = as.character(sums[[1]]), ordered = TRUE)
        
        # --
        
        # - prepare mdf of ps Species -
        if(taxa_are_rows(physeq2)) { physeq2 <- t(physeq2) }
        
        if (!all(unique(sample_data(physeq2)[[group_var]]) %in% names(color_levels))) {
                keepSamples <- sample_names(physeq2)[sample_data(physeq2)[[group_var]] %in% names(color_levels)]
                physeq2 <- prune_samples(samples = keepSamples, physeq2)
                sample_data(physeq2)[[group_var]] <- factor(sample_data(physeq2)[[group_var]], levels = names(color_levels), order = TRUE)
        }
        
        if (!is.factor(sample_data(physeq2)[[group_var]])) {sample_data(physeq2)[[group_var]] <- factor(sample_data(physeq2)[[group_var]], levels = names(color_levels), order = TRUE)}
        
        
        mdf2 <- phyloseq::psmelt(physeq2)
        
        # order fill levels according to abundance over all samples
        mdf2[, fill] <- as.character(mdf2[, fill])
        mdf2[is.na(mdf2[, fill]), fill] <- "NA"
        # sums <- group_by_(mdf2, fill) %>% summarise(sum_abundance = sum(Abundance)) %>% arrange(sum_abundance)
        mdf2[, fill] <- factor(mdf2[, fill], levels = as.character(sums[[1]]), ordered = TRUE)
        # --
        
        
        # - prepare mdf of ps genera -
        if(taxa_are_rows(physeq3)) { physeq3 <- t(physeq3) }
        
        if (!all(unique(sample_data(physeq3)[[group_var]]) %in% names(color_levels))) {
                keepSamples <- sample_names(physeq3)[sample_data(physeq3)[[group_var]] %in% names(color_levels)]
                physeq3 <- prune_samples(samples = keepSamples, physeq3)
                sample_data(physeq3)[[group_var]] <- factor(sample_data(physeq3)[[group_var]], levels = names(color_levels), order = TRUE)
        }
        
        if (!is.factor(sample_data(physeq3)[[group_var]])) {sample_data(physeq3)[[group_var]] <- factor(sample_data(physeq3)[[group_var]], levels = names(color_levels), order = TRUE)}
        
        
        mdf3 <- phyloseq::psmelt(physeq3)
        mdf3$Species<-NA
        # order fill levels according to abundance over all samples
        mdf3[, fill] <- as.character(mdf3[, fill])
        mdf3[is.na(mdf3[, fill]), fill] <- "NA"
        # sums <- group_by_(mdf2, fill) %>% summarise(sum_abundance = sum(Abundance)) %>% arrange(sum_abundance)
        mdf3[, fill] <- factor(mdf3[, fill], levels = as.character(sums[[1]]), ordered = TRUE)
        # --
        
   
        mdf$Typer <- taxa_Level1
        mdf2$Typer <- taxa_Level2
        mdf3$Typer <- taxa_Level3
        
        mdf <- rbind(mdf, mdf2,mdf3)
        mdf$Typer <- factor(mdf$Typer, levels = c(taxa_Level1, taxa_Level2, taxa_Level3), ordered = TRUE)
        
        
        
        # order samples according to group_var levels and potentially by total counts
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        mdf$Sample <- factor(mdf$Sample, levels = sort(LookUpDF$Sample, decreasing=FALSE), ordered = TRUE)
       
        
        colNames = (sapply(strsplit(as.character(levels(mdf$Sample)), split = "_"),  `[`,1))
        
        colxaxis = color_levels[colNames]
        
        # - define names of x axis using color_levels (which must be a named character vector) -        
        #colxaxis <- color_levels[LookUpDF$Group]
        # --
        mdf2$Sample <- factor(mdf2$Sample, levels = sort(LookUpDF$Sample, decreasing=FALSE), ordered = TRUE)
        mdf3$Sample <- factor(mdf3$Sample, levels = sort(LookUpDF$Sample, decreasing=FALSE), ordered = TRUE)
        
        if (is.null(col_vec)){
                if (length(levels(mdf[, fill])) <= 15) {
                        fill_colors <- make_color_vector(mdf[, fill], rev(QuantColors15[1:length(levels(mdf[, fill]))]))
                } else {
                        fill_colors <- make_color_vector(mdf[, fill], viridis(length(levels(mdf[, fill]))))
                }
                
        } else {
                fill_colors <- col_vec
        }
        
        fill_colors2 <- make_color_vector(mdf2[, "Species"], rev(QuantColors15[1:length(levels(mdf2[, "Species"]))]))
        
        Tr <- ggplot(mdf, aes_string(x = x, y = y, fill = fill))
        Tr <- Tr + 
                geom_bar(stat = "identity", position = "stack") +
                theme_bw() +
                scale_fill_manual(values = fill_colors, na.value = "grey") +
                xlab("") +
                facet_wrap(~ Typer, ncol = 1) +
                theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, colour = colxaxis))
        
        #Tr
        
        Tr1<- ggplot(mdf3, aes_string(x = x, y = y, fill = fill))
        Tr1 <- Tr1 + 
                geom_bar(stat = "identity", position = "stack") +
                theme_bw() +
                scale_fill_manual(values = fill_colors, na.value = "grey") +
                xlab("") +
                facet_wrap(~ Method, ncol = 1) +
                theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, colour = colxaxis))
        
        #Tr1
        
        Tr2<- ggplot(mdf2, aes_string(x = x, y = y, fill = "Species"))
        Tr2 <- Tr2 + 
                geom_bar(stat = "identity", position = "stack") +
                theme_bw() +
                scale_fill_manual(values = fill_colors2, na.value = "grey") +
                xlab("") +
                facet_wrap(~ Method, ncol = 1) +
                theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, colour = colxaxis))
        
        out<-list(Tax_levels=Tr,tax_genus=Tr1, tax_Species=Tr2)
}
# --



#######################################
#### plot_sample_bars_compare
#######################################
# generates abundance barplots (see plot_bar_own) to compared different ps_object 

plot_sample_bars_compare_3 <- function(physeq, physeq2, physeq3, x = "Sample", y = "Abundance", group_var, color_levels, fill = NULL,
                                       color_sample_names = TRUE, col_vec = NULL,  taxa_Level1=taxa_Level1, taxa_Level2= taxa_Level2, taxa_Level3= taxa_Level3){
        
        
        # - prepare mdf of ps genera -
        if(taxa_are_rows(physeq3)) { physeq3 <- t(physeq3) }
        
        if (!all(unique(sample_data(physeq3)[[group_var]]) %in% names(color_levels))) {
                keepSamples <- sample_names(physeq3)[sample_data(physeq3)[[group_var]] %in% names(color_levels)]
                physeq3 <- prune_samples(samples = keepSamples, physeq3)
                sample_data(physeq3)[[group_var]] <- factor(sample_data(physeq3)[[group_var]], levels = names(color_levels), order = TRUE)
        }
        
        if (!is.factor(sample_data(physeq3)[[group_var]])) {sample_data(physeq3)[[group_var]] <- factor(sample_data(physeq3)[[group_var]], levels = names(color_levels), order = TRUE)}
        
        
        mdf3 <- phyloseq::psmelt(physeq3)
        
        # order fill levels according to abundance over all samples
        # mdf3$Sample <- factor(mdf3$Sample, levels = sort(mdf3$Sample, decreasing=FALSE), ordered = TRUE)
        # 
        # 
        colNames = (sapply(strsplit(as.character(levels(mdf$Sample)), split = "_"),  `[`,1))

        colxaxis = color_levels[colNames]
        
        sample_data(physeq3)$Sample <- factor(mdf3$Sample, levels = sort(mdf3$Sample, decreasing=FALSE), ordered = TRUE)

        
        
        Tr1<- ggplot(mdf3, aes_string(x = x, y = y, fill = fill))
        Tr1 <- Tr1 + 
                geom_bar(stat = "identity", position = "stack") +
                theme_bw() +
                scale_fill_manual(values = viridis(length(levels(mdf3[, fill]))), na.value = "grey") +
                xlab("") +
                facet_wrap(~ Method, ncol = 1) +
                theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, colour = colxaxis ))
       
        out<-list(Tr1)
}
# --
