cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#######################################
### FUNCTION: assignTaxonomyaddSpecies
#######################################

# Function that simply calls assignTaxonomy and then addSpecies from dada2
## Input
# seqtab: your sequence table from dada2
# minBoot: minBoot for assignTaxonomy
# allowMultiple: default 3, allowMultiple for addSpecies
# PathToRefs: path to the folder with bot the RefDataBase and SpeciesDB in
# RefDataBase: the reference database for assignTaxonomy
# SpeciesDB: the reference database for addSpecies
# PathToSave: the folder into which taxa and taxa.species will be saved, if NULL working directory is used
## Output:
# taxa, taxa.species are saved as Taxonomy.RData


assignTaxonomyaddSpecies <- function(seqtab, 
                                     minBoot = 50, 
                                     allowMultiple = TRUE,
                                     PathToRefs = NULL,
                                     RefDataBase = "silva_nr99_v138.1_wSpecies_train_set.fa.gz",#"silva_nr99_v138.1_train_set.fa.gz",
                                     SpeciesDB = "silva_species_assignment_v138.1.fa.gz",
                                     PathToSave = getwd(),
                                     tryRC = FALSE){
        
        ## Packages
        try(library(dada2), biocLite("dada2"))
        
        
        ## Reference Databases
        if(is.null(PathToRefs)){
                stop("PathToRefs must be given")
        }
        
        RefDB <- file.path(PathToRefs, RefDataBase)
        
        if(!file.exists(RefDB)){
                stop("could not find the reference database for assignTaxonomy")
        }
        
        SpecDB <- file.path(PathToRefs, SpeciesDB)
        
        # if(!file.exists(SpecDB)){
        #         stop("could not find the reference database for addSpecies")
        # }
        
        ## list the input
        InputSave <- list(seqtab = seqtab, 
                          minBoot = minBoot,
                          allowMultiple = allowMultiple,
                          tryRC = tryRC,
                          PathToRefs = PathToRefs, 
                          RefDataBase = RefDataBase,
                          SpeciesDB = SpeciesDB,
                          PathToSave = PathToSave)
        
        
        ## run assignTaxonomy and save
        ptm <- proc.time()
        taxa <- assignTaxonomy(seqtab, refFasta = RefDB, verbose = TRUE, minBoot = minBoot, tryRC = tryRC)
        TimeForassignTaxonomy <- (proc.time()-ptm)[3]
        
        InputSave[[length(InputSave) + 1]] <- TimeForassignTaxonomy
        
        save(taxa, InputSave, file = file.path(PathToSave, "Taxonomy.RData"))
        
        message("assignTaxonomy has been run")
        
        # run addSpecies and save
        if(!is.null(SpeciesDB)){
        ptm <- proc.time()
        taxa.species <- addSpecies(taxa, refFasta = SpecDB, verbose = TRUE, allowMultiple = allowMultiple)
        TimeForaddSpecies <- (proc.time()-ptm)[3]

        InputSave[[length(InputSave) + 1]] <- TimeForaddSpecies

        save(taxa, taxa.species, InputSave, file = file.path(PathToSave, "Taxonomy.RData"))

        message("Function Done")
        }

}


#######################################
### FUNCTION: construct_phylogenetic_tree
#######################################
## Background: 
# the entire function is basically a copy from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4955027/ (page 8)
# the function uses the packages DECIPHER (bioconductor) and phangorn (cran)
# I had once a clash of the DECIPHER package with dada2 when dada2 was loaded first, so maybe function needs to be run in a fresh session
# you should run the wrapper on the server because the last phangorn::optim.pml takes long time
## Inputs, 
# savepath = determines were the output data list will be saved
# all other inputs: refer to arguments for phangorn::optim.pml and I have honestly no clue about them at the moment
# so read code and help
## output: 
# a list with the GTR: generalized time-reversible with Gamma rate variation model fit containing the tree, and the Input data saved
# the list is saved as phylogenetic_tree.rds in savepath


construct_phylogenetic_tree <- function(seqtab.nochim, savepath,
                                        k = 4, inv = .2, 
                                        model = "GTR", rearrangement = "stochastic",
                                        trace = 0) {
        
        library(DECIPHER)
        library(phangorn); packageVersion("phangorn")
        
        RVer <- R.Version()
        RVer <- RVer$version.string
        PackageVersions <- data.frame(Package = c("R", "DECIPHER", "phangorn"),
                                      Version = c(RVer,
                                                  as.character(packageVersion("DECIPHER")),
                                                  as.character(packageVersion("phangorn"))))
        
        
        seqs <- colnames(seqtab.nochim) # simple character string
        names(seqs) <- seqs
        
        Inputs <- list(seqs = seqs,
                       savepath = savepath,
                       k = k,
                       inv = inv,
                       model = model,
                       rearrangement = rearrangement,
                       trace = trace)
        
        # ---- multiple alignment of SVs using "DECIPHER" (takes only 10 seconds) ----
        
        alignment <- DECIPHER::AlignSeqs(DNAStringSet(seqs), anchor = NA)
        
        # --------
        
        # ---- use alignment to construct phylogenetic tree (phangorn package) -----
        
        # -- First construct a neighbor-joining tree --
        
        phang.align <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA") # just transformation of DNA data into phyDat format
        dm <- phangorn::dist.ml(phang.align) # computes pairwise distances
        
        treeNJ <- phangorn::NJ(dm)
        
        # NB: They warn: tip order != sequence order
        # however, I got a TRUE here
        # all.equal(treeNJ$tip.label, seqs, check.attributes = FALSE)
        
        fit <- phangorn::pml(treeNJ, data = phang.align) # warns: negative edges length changed to 0!
        
        # ----
        
        # -- use the neighbor-joining tree as a start point to fit a GTR+G+I (generalized time-reversible with Gamma rate variation) maximum 
        # likelihood tree --
        
        fitGTR <- update(fit, k = k, inv = inv)
        message("starting the time consuming fit now")
        fitGTR <- phangorn::optim.pml(fitGTR, model = model, optInv = TRUE, optGamma = TRUE,
                                      rearrangement = rearrangement, control = pml.control(trace = trace))
        message("done with the time consuming fit:)")
        
        
        out <- list(fitGTR = fitGTR,
                    Inputs = Inputs,
                    PackageVersions = PackageVersions)
        
        saveRDS(out, file = file.path(savepath, "phylog_tree.rds"))
        # ----
        
        # --------
        
        return(out)
}



#######################################
### FUNCTION: get_assignemnt_distribution
#######################################

# Function to determine how many sequences could not be assigned using a given minBoot value
## Input
# taxa: the output of the assignTaxonomy (dada2) command. A matrix with nrow = number of sequences and ncol usually 6, the taxonomic levels. When 
# the minBoot criterion was not fulfilled at a taxonomic level, assignTaxonomy assigns an NA. 
## Output
# assignment_distribution: data frame that states for each taxonomic level in taxa how many SVs have been assigned 

get_assignemnt_distribution <- function(taxa){
        
        countNA <- apply(taxa, 2, function(x){sum(is.na(x))})
        countNonNA <- apply(taxa, 2, function(x){sum(!is.na(x))})
        ambiguous <- apply(taxa, 2, function(x){sum(grepl("/", x))})
        total <- countNA + countNonNA
        assignment_distribution <- data.frame(assigned = countNonNA, 
                                              assigned_unamb = countNonNA - ambiguous,
                                              total = total,
                                              PC_assigned = round(100*countNonNA/total, 1),
                                              PC_assigned_unamb = round(100*(countNonNA - ambiguous)/total, 1))
        
}


#######################################
### FUNCTION: check_assignment_vs_prevalence
#######################################
# checks if there is a trend for better assignment for more prealent SVs

check_assignment_vs_prevalence <- function(taxa, seqtab, prevalences = seq(0, 90, by = 10)){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
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
                scale_colour_manual("", values = cbPalette[2:11]) +
                scale_x_discrete(breaks = prevalences, labels = paste(prevalences, " (", No_ASVs, ")", sep = "")) +
                ylab("percentage of ASVs assigned") +
                xlab("prevalence (No ASVs remaining)") +
                theme_bw(12) +
                theme(axis.text.x = element_text(angle=90, vjust=0.5))
        
        list(PC_assigned, tr)
        
}


#######################################
### FUNCTION: check_assignment_vs_abundance
#######################################
# checks if there is a trend for better assignment for more abundant SVs

check_assignment_vs_abundance <- function(taxa, seqtab, abundanceQuantiles = seq(0, 90, by = 10)){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
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
                scale_colour_manual("", values = cbPalette[2:11]) +
                scale_x_discrete(breaks = abundanceQuantiles, labels = paste(abundanceQuantiles, " (", round(abQuantiles), ", ", No_ASVs, ")", sep = "")) +
                ylab("percentage of ASVs assigned") +
                xlab("total counts quantile (total counts filter, No ASVs remaining)") +
                theme_bw(12) +
                theme(axis.text.x = element_text(angle=90, vjust=0.5))
        
        list(PC_assigned, tr)
        
}


#######################################
### FUNCTION: plotTaxLevelvsAbundPrev
#######################################

# Function to determine how many sequences could be assigned to the different taxonomic levels compared to the abundance
# and prevalence of the SV
## Input
# taxa: the output of the assignTaxonomy (dada2) and maybe addSpecies command. nrow = number of sequences and ncol number taxonomic levels.
# seqtab: The abundance table output from Dada2_wrap
## Output
# list of plots ab vs assignment and list of plots prevalence vs assignment

plotTaxLevelvsAbundPrev <- function(taxa, seqtab, color = "#E69F00"){
        
        taxa_vs_ab_pref <- as.data.frame(cbind(Abundance = colSums(seqtab), Prevalence = colSums(seqtab != 0), !is.na(taxa)))
        
        TrList_ab <- list()
        for (i in 3:ncol(taxa_vs_ab_pref)){
                Tr <- ggplot(taxa_vs_ab_pref, aes_string(x = "Abundance", y = colnames(taxa_vs_ab_pref)[i]))
                Tr <- Tr + geom_jitter(col = color, size = 1.5, alpha = 0.5) +
                        scale_y_continuous(limits = c(-0.5, 1.5), breaks = c(0,1), labels = c('No', "Yes")) +
                        xlab("abundance of SV") +
                        ylab(paste("assigned to ", colnames(taxa_vs_ab_pref)[i], sep = "")) +
                        theme_bw() 
                TrList_ab[[i-2]] <- Tr
        }
        names(TrList_ab) <- colnames(taxa_vs_ab_pref)[3:ncol(taxa_vs_ab_pref)]
        
        TrList_prev <- list()
        for (i in 3:ncol(taxa_vs_ab_pref)){
                Tr <- ggplot(taxa_vs_ab_pref, aes_string(x = "Prevalence", y = colnames(taxa_vs_ab_pref)[i]))
                Tr <- Tr + geom_jitter(col = color, size = 1.5, alpha = 0.5) +
                        scale_y_continuous(limits = c(-0.55, 1.55), breaks = c(0,1), labels = c('No', "Yes")) +
                        xlab("SV is present in No of samples") +
                        ylab(paste("assigned to ", colnames(taxa_vs_ab_pref)[i], sep = "")) +
                        theme_bw() 
                TrList_prev[[i-2]] <- Tr
        }
        names(TrList_prev) <- colnames(taxa_vs_ab_pref)[3:ncol(taxa_vs_ab_pref)]
        
        list(TrList_ab = TrList_ab, TrList_prev = TrList_prev)
        
}


#######################################
### rarefy_sample
#######################################
# simple rarefy function based on sample(), see phyloseq:::rarefaction_subsample

rarefy_sample <- function (sample_cnts, size) {
        # attention, only works if nrow(ProbMat) is not NULL
        simsample <- integer(length(sample_cnts)) 
        suppressWarnings(draws <- sample(1:length(sample_cnts), size = size, replace = TRUE, prob = sample_cnts))
        drawtable <- table(draws)
        simsample[as(names(drawtable), "integer")] <- drawtable
        return(simsample)
}


#######################################
### rarefaction_curve_own
#######################################
# creates total amplicon steps from 0 to max_total. Then rarefies otu_table(physeq) to these steps and calculates each time
# richness and shannon. It plots these rarefaction curves. If group_var (factor in sample_data(physeq)) is given, it also
# plots grouped (Mean plus SE) and calculates the pairwise t.test p.values between the groups.
# type can be vegan or sample, rarefaction is then either based on vegan::rrarefy or rarefy_sample (sample())
# Outputs a named list with all results.
# NB: depending on step_size and size of otu_table(physeq) expect it to take 2-4 minutes!
# NB2: runs only one rarefaction per step, so depends on seed!

rarefaction_curve_own <- function(physeq, group_var = NULL, max_total = NULL, step_size = 200, type = "vegan", seed = 123) {
        
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
        if (taxa_are_rows(physeq)) {
                seqtab <- t(as(otu_table(physeq), "matrix"))
        } else {
                seqtab <- as(otu_table(physeq), "matrix")
        }
        
        Group <- sample_data(physeq)[[group_var]]
        
        if (!is.null(Group) && !is.factor(Group)) {
                Group <- as.factor(Group)
        }
        
        if (is.null(max_total)) {
                max_total <- quantile(rowSums(seqtab), probs = .25)
        }
        
        steps <- seq(from = 0, to = max_total, by = step_size)
        NoSamples <- nsamples(ps)
        NoSteps <- length(steps)
        richness_matrix <- matrix(data = NA, nrow = NoSamples, ncol = NoSteps)
        shannon_matrix <- matrix(data = NA, nrow = NoSamples, ncol = NoSteps)
        totalamplicons <- rowSums(seqtab)
        
        set.seed(seed)
        # ptm <- Sys.time()
        if (type == "vegan") {
                for (i in 1:length(steps)) {
                        rare_CT <- rrarefy(seqtab, sample = steps[i])
                        richness_matrix[,i] <- rowSums(rare_CT != 0)
                        shannon_matrix[,i] <- vegan::diversity(rare_CT, index = "shannon")
                }
        } else if (type == "sample") {
                #ptm <- Sys.time()
                for (i in 1:length(steps)) {
                        rare_CT <- t(apply(seqtab, 1, function(cnts){rarefy_sample(cnts, size = steps[i])}))
                        richness_matrix[,i] <- rowSums(rare_CT != 0)
                        shannon_matrix[,i] <- vegan::diversity(rare_CT, index = "shannon")
                }
                # Sys.time() - ptm
        } else {
                stop("type neither vegan nor sample")
        }
        
        # Sys.time() - ptm
        
        # set alpha-diversities for samples at which steps > totalamplicons to NA
        for (i in 1:NoSamples) {
                NaIndex <- which(totalamplicons[i] < steps)[1]
                if (!is.na(NaIndex)){
                        richness_matrix[i, NaIndex:ncol(richness_matrix)] <- NA
                        shannon_matrix[i, NaIndex:ncol(shannon_matrix)] <- NA
                }
        }
        
        rownames(richness_matrix) <- rownames(seqtab)
        colnames(richness_matrix) <- paste("step_", steps, sep = "")
        rownames(shannon_matrix) <- rownames(seqtab)
        colnames(shannon_matrix) <- paste("step_", steps, sep = "")
        
        richness_df <- as.data.frame(richness_matrix)
        shannon_df <- as.data.frame(shannon_matrix)
        
        plot_div_df <- function (div_df, type = "alpha diversity") {
                
                div_df$Sample <- rownames(div_df)
                div_df <- tidyr::gather(div_df, key = step, value = div, -Sample)
                div_df <- tidyr::separate(div_df, col = step, into = c("term", "step"), sep = "_")
                div_df$step <- as.numeric(div_df$step)
                
                Tr <- ggplot(div_df, aes(x = step, y = div, group = Sample))
                Tr <- Tr +
                        geom_line() +
                        xlab("total amplicons") +
                        ylab(type) +
                        theme_bw(12)
        }
        
        Tr_richness <- plot_div_df(richness_df, type = "richness")
        Tr_shannon <- plot_div_df(shannon_df, type = "shannon")
        
        
        if (!is.null(Group)) {
                
                richness_df$Group <- Group
                shannon_df$Group <- Group
                
                pairwise.tt_richness <- lapply(richness_df[, -ncol(richness_df)], function(step){ 
                        ptt <- pairwise.t.test(x = step, g = richness_df$Group, alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
                        ptt$p.value
                })
                
                pairwise.tt_shannon <- lapply(shannon_df[, -ncol(shannon_df)], function(step){ 
                        ptt <- pairwise.t.test(x = step, g = shannon_df$Group, alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
                        ptt$p.value
                })
                
                
                plot_div_df_group <- function (div_df, type = "alpha diversity") {
                        div_df <- tidyr::gather(div_df, key = step, value = div, -Group)
                        div_df <- tidyr::separate(div_df, col = step, into = c("term", "step"), sep = "_")
                        div_df$step <- as.numeric(div_df$step)
                        div_df <- group_by(div_df, Group, step)
                        div_df <- dplyr::summarise(div_df, Mean_div = mean(div, na.rm = T), SD_div = sd(div, na.rm = T), n = n(), SE_div = SD_div/sqrt(n))
                        
                        Tr <- ggplot(div_df, aes(x = step, y = Mean_div, col = Group))
                        Tr <- Tr + 
                                geom_line() +
                                geom_point(size = 1) +
                                geom_errorbar(aes(ymin = Mean_div-SE_div, ymax = Mean_div+SE_div)) +
                                scale_color_manual("", values = cbPalette[2:11]) +
                                xlab("total amplicons") +
                                ylab(paste(type, "(group mean +- SE)")) +
                                theme_bw(12)
                }
                
                Tr_richness_group <- plot_div_df_group(richness_df, type = "richness")
                Tr_shannon_group <- plot_div_df_group(shannon_df, type = "shannon")
                
        }
        
        outlist <- list(richness_df = richness_df, shannon_df = shannon_df,
                        Tr_richness = Tr_richness, Tr_shannon = Tr_shannon)
        
        if (!is.null(Group)) {
                outlist[[5]] <- pairwise.tt_richness
                outlist[[6]] <- pairwise.tt_shannon
                outlist[[7]] <- Tr_richness_group
                outlist[[8]] <- Tr_shannon_group
        }
        
        return(outlist)
        
}


#######################################
### rarefaction_curve_own_fast
#######################################
# read rarefaction_curve_own, this fast version uses the vegan::rarefy function, so it only generates rarefaction
# for richness. NB: rarefaction_curve_own only does one rarefaction for each step, so in principle you had to run it
# several times and average (since I do not it gives more bumpy curves). vegan::rarefy used here does averaging (read help)
# and is super fast. I still do not use the SE provided by rarefy, so that could be added in here

# # NB: in case you wanted to add sample labels to a Tr_richness_col
# df <- Tr_richness_col$data
# df <- filter(df, step == max(df$step))
# df$step <- df$step + 400
# Tr_richness_col <- Tr_richness_col + 
#         geom_label(data = df, aes(label = Sample)) +
#         scale_x_continuous(limits = c(-100, max(df$step) + 2500))

rarefaction_curve_own_fast <- function(physeq, group_var = NULL, max_total = NULL, step_size = 200, seed = 123) {
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
        
        if (taxa_are_rows(physeq)) {
                seqtab <- t(as(otu_table(physeq), "matrix"))
        } else {
                seqtab <- as(otu_table(physeq), "matrix")
        }
        
        if (!is.null(group_var)){
                Group <- sample_data(physeq)[[group_var]] 
                if (!is.null(Group) && !is.factor(Group)) {
                        Group <- as.factor(Group)
                }
        } else { Group <- NULL}
        
        
        if (is.null(max_total)) {
                max_total <- quantile(sample_sums(physeq), probs = .25)
        }
        
        steps <- seq(from = 0, to = max_total, by = step_size)
        NoSamples <- nsamples(physeq)
        NoSteps <- length(steps)
        richness_matrix <- matrix(data = NA, nrow = NoSamples, ncol = NoSteps)
        totalamplicons <- rowSums(seqtab)
        
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
        
        plot_div_df <- function (div_df, type = "alpha diversity") {
                
                div_df$Sample <- rownames(div_df)
                div_df <- tidyr::gather(div_df, key = step, value = div, -Sample)
                div_df <- tidyr::separate(div_df, col = step, into = c("term", "step"), sep = "_")
                div_df$step <- as.numeric(div_df$step)
                
                Tr <- ggplot(div_df, aes(x = step, y = div, group = Sample))
                Tr <- Tr +
                        geom_line() +
                        xlab("total counts in sample") +
                        ylab(type) +
                        theme_bw(12)
        }
        
        Tr_richness <- plot_div_df(richness_df, type = "richness")
        
        plot_div_df3 <- function (div_df, type = "richness") {
                
                div_df$Sample <- rownames(div_df)
                div_df$Total <- totalamplicons
                div_df <- tidyr::gather(div_df, key = step, value = div, -Sample, -Total)
                div_df <- tidyr::separate(div_df, col = step, into = c("term", "step"), sep = "_")
                div_df$step <- as.numeric(div_df$step)
                
                Tr <- ggplot(div_df, aes(x = step, y = div, group = Sample, col = Total))
                Tr <- Tr +
                        geom_line() +
                        scale_color_gradient2("total counts bef.", low = cbPalette[6], mid = cbPalette[1], high = cbPalette[2], midpoint = 
                                                      median(div_df$Total)) +
                        xlab("total counts in sample") +
                        ylab(type) +
                        theme_bw(12)
        }
        
        Tr_richness_grad <- plot_div_df3(richness_df, type = "richness")
        
        
        if (!is.null(Group)) {
                
                richness_df$Group <- Group
                
                # pairwise.tt_richness <- lapply(richness_df[, -ncol(richness_df)], function(step){ 
                #         ptt <- pairwise.t.test(x = step, g = richness_df$Group, alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
                #         ptt$p.value
                #})
                
                plot_div_df2 <- function (div_df, type = "richness") {
                        
                        div_df$Sample <- rownames(div_df)
                        div_df <- tidyr::gather(div_df, key = step, value = div, -Sample, -Group)
                        div_df <- tidyr::separate(div_df, col = step, into = c("term", "step"), sep = "_")
                        div_df$step <- as.numeric(div_df$step)
                        
                        Tr <- ggplot(div_df, aes(x = step, y = div, group = Sample, col = Group))
                        Tr <- Tr +
                                geom_line() +
                                xlab("total counts in sample") +
                                scale_color_manual("", values = cbPalette[2:11]) +
                                ylab(type) +
                                theme_bw(12)
                }
                
                plot_div_df_group <- function (div_df, type = "alpha diversity") {
                  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
                  
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
                                scale_color_manual("", values = cbPalette[2:11]) +
                                xlab("total counts in sample") +
                                ylab(paste(type, "(group mean +- SE)")) +
                                theme_bw(12)
                }
                Tr_richness_col <- plot_div_df2(richness_df, type = "richness")
                
                Tr_richness_group <- plot_div_df_group(richness_df, type = "richness")
                
        } 
        
        outlist <- list(rarefaction_richness_df = richness_df, Tr_richness = Tr_richness, Tr_richness_grad = Tr_richness_grad)
        outlist <- list(Tr_richness = Tr_richness, Tr_richness_grad = Tr_richness_grad)
        
        if (!is.null(Group)) {
                #outlist[[4]] <- pairwise.tt_richness
                outlist[[4]] <- Tr_richness_col
                outlist[[5]] <- Tr_richness_group
                names(outlist)[4:5] <- c("Tr_richness_col", "Tr_richness_group")
        }
        
        return(outlist)
        
}


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


#######################################
### calculate_alphadiversity
#######################################
# NB: I use the estimate_richness function here (that is also used in plot_richness from phyloseq).
# You could easily calculate Shannon, Chao1, Observed self, see alphaDiversityMeasures.Rmd
# estimate_richness itself uses functions from the vegan package

calculate_alphadiversity <- function(physeq, measures = c("Observed", "Shannon")) {
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
        DF_alpha <- suppressWarnings(phyloseq::estimate_richness(physeq, measures = measures))
        rownames(DF_alpha) <- sample_names(physeq)
        if ("Observed" %in% colnames(DF_alpha)){
                colnames(DF_alpha)[colnames(DF_alpha) == "Observed"] <- "Richness" 
                measures[measures == "Observed"] <- "Richness"
        }
        
        DF_alpha$Total <- sample_sums(physeq)
        
        ncol_df_alpha <- ncol(DF_alpha)
        # add residuals of lm fits to total_amplicons
        fitlist <- list()
        for (i in 1:length(measures)){
                fit <- lm(DF_alpha[,measures[i]] ~ DF_alpha[,"Total"])
                fitlist[[i]] <- fit
                DF_alpha[, ncol_df_alpha + i] <- residuals(fit)
                colnames(DF_alpha)[ncol_df_alpha + i] <- paste(measures[i], "resids", sep = "_")
        }
        
        names(fitlist) <- measures
        
        # if you have, add lm fits to FilteredReads
        fitlist_FilteredReads <- list()
        
        if (!is.null(sample_data(physeq)$FilteredReads)){
                DF_alpha$filtered_reads <- sample_data(physeq)$FilteredReads
                
                
                ncol_df_alpha <- ncol(DF_alpha)
                for (i in 1:length(measures)){
                        fit <- lm(DF_alpha[,measures[i]] ~ DF_alpha[,"filtered_reads"])
                        fitlist_FilteredReads[[i]] <- fit
                        DF_alpha[, ncol_df_alpha + i] <- residuals(fit)
                        colnames(DF_alpha)[ncol_df_alpha + i] <- paste(measures[i], "resids", "FilteredReads", sep = "_")
                }
        }
        
        names(fitlist_FilteredReads) <- measures
        
        DF_alpha <- data.frame(Sample = rownames(DF_alpha), DF_alpha, sample_data(physeq))
        
        outlist <- list(DF_alpha = DF_alpha, fitlist = fitlist, fitlist_FilteredReads = fitlist_FilteredReads)
}


calculate_alphadiversity1 <- function(physeq, measures = c("Observed", "Shannon")) {
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
  DF_alpha <- suppressWarnings(phyloseq::estimate_richness(physeq, measures = measures))
  rownames(DF_alpha) <- sample_names(physeq)
  if ("Observed" %in% colnames(DF_alpha)){
    colnames(DF_alpha)[colnames(DF_alpha) == "Observed"] <- "Richness" 
    measures[measures == "Observed"] <- "Richness"
  }
  
  DF_alpha$Total <- sample_sums(physeq)
  
  ncol_df_alpha <- ncol(DF_alpha)
  # add residuals of lm fits to total_amplicons
  fitlist <- list()
  for (i in 1:length(measures)){
    fit <- lm(DF_alpha[,measures[i]] ~ DF_alpha[,"Total"])
    fitlist[[i]] <- fit
    DF_alpha[, ncol_df_alpha + i] <- residuals(fit)
    colnames(DF_alpha)[ncol_df_alpha + i] <- paste(measures[i], "resids", sep = "_")
  }
  
 
  
  DF_alpha <- data.frame(Sample = rownames(DF_alpha), DF_alpha, sample_data(physeq))
  
  outlist <- list(DF_alpha = DF_alpha)
}

#######################################
### adj_LS
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


adj_LS <- function(physeq, zeros.count = FALSE, percentile = 50, plots = FALSE)  {
        
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
        
        if (plots){
                
          cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
          
                
                # compare calculated SFs to library sizes of the samples
                if(!identical(names(SFs), names(sample_sums(physeq)))){warning("Probably some Mix Up CHECK!")}
                comp <- data.frame(Sample = names(SFs), TotAmps = sample_sums(physeq), SFs = SFs)
                comp$SFsNormed <- comp$SFs/median(comp$SFs)
                # comp$TotAmpsNormed <- comp$TotAmps/mean(comp$TotAmps)
                comp$TotAmpsNormed <- comp$TotAmps/median(comp$TotAmps)
                comp <- dplyr::arrange(comp, desc(TotAmpsNormed))
                comp$Sample <- as.factor(comp$Sample)
                LevelsWant <- as.character(comp$Sample)
                for(i in 1:length(LevelsWant)){
                        comp$Sample <- relevel(comp$Sample, LevelsWant[i])
                }
                
                comp <- comp[c(1,4,5)]
                names(comp)[c(2,3)] <- c("SizeFactor_DESeq", "TotalAmplicons_relAb")
                comp <- tidyr::gather(comp, key = Corrector, value = NormedValue, -Sample)
                Tr <- ggplot(comp, aes(x = Sample, y = NormedValue, color = Corrector))
                Tr <- Tr + geom_point(size = 2) +
                        xlab("") +
                        ylab("correction value (normalized to median)") +
                        ggtitle(paste("Median SF: ", round(median(SFs),3), " Median TotAmp: ", round(median(sample_sums(physeq)),3), sep = "")) +
                        scale_color_manual(values = cbPalette[c(4,2)]) +
                        theme_bw() +
                        theme(panel.grid.minor = element_blank(),
                              panel.grid.major.y = element_blank(),
                              panel.grid.major.x = element_line(color = "#999999", size = .15),
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                              legend.title = element_blank())
                
                # -- 3c: save Histograms of the SFs and sample_sums
                histo2 <- function(x, xtitle, gtitle) {
                        x <- data.frame(x = x)
                        Tr <- ggplot(x, aes(x = x))
                        Tr <- Tr + geom_histogram(binwidth = diff(range(x))/60, col = "black", fill = "#E69F00") +
                                geom_rug() +
                                geom_vline(xintercept = median(x$x), col = "#009E73", size = 1) +
                                ylab("No Samples") + 
                                xlab(xtitle) +
                                ggtitle(gtitle) +
                                theme_bw() + 
                                theme(panel.grid.minor = element_blank(),
                                      panel.grid.major.y = element_blank(),
                                      panel.grid.major.x = element_line(color = "#999999", size = .15))
                }
                
                Tr2 <- histo2(SFs, xtitle = "Size Factors", gtitle = "Size Factors a la DESeq")
                Tr3 <- histo2(sample_sums(physeq), xtitle = "Total Amplicons", gtitle = "Size Factors a la relative abundance")
                
                List <- list(Physeq = phynew, SFs = SFs, SizeFactorCompare = Tr, SFHistos = list(Tr2, Tr3), RefSample = GM)
                
                
        } else {
                
                list(physeq = phynew, SFs = SFs)
                
        }
}

#######################################
### FUNCTION: alpha_diversity_wrapper
#######################################
# just wraps around other functions here to save writing time
# currently only works for alpha_div_measures = c("Observed", "Shannon")
alpha_diversity_wrapper <- function(physeq, alpha_div_measures = c("Observed", "Shannon"), group_var = group_var, shape = shape){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )

        DF_alpha_list <- calculate_alphadiversity(physeq = physeq, measures = alpha_div_measures)
        DF_alpha <- DF_alpha_list[[1]]

        # just I prefer Richness over Observed
        alpha_div_measures2 <- alpha_div_measures
        if ("Observed" %in% alpha_div_measures) {
                alpha_div_measures2[alpha_div_measures2 == "Observed"] <- "Richness"
        }

        pairwise.tt_rich <- pairwise.t.test(x = DF_alpha$Richness, g = DF_alpha[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
        # compare, to get the same results as individual t.test you need pool.sd = F
        # t.test(x = DF_alpha$Richness[DF_alpha$Group == "Old"], y = DF_alpha$Richness[DF_alpha$Group == "MiddleAged"], alternative = "two", var.equal = F)$p.value

        pairwise.tt_shannon <- pairwise.t.test(x = DF_alpha$Shannon, g = DF_alpha[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)

        if(var(DF_alpha$Total) > 0) {
                pairwise.tt_totalcounts <- pairwise.t.test(x = DF_alpha$Total, g = DF_alpha[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
        } else {
                pairwise.tt_totalcounts <- NA
        }

        #pairwise.tt_rich_resids <- pairwise.t.test(x = DF_alpha$Richness_resids, g = DF_alpha[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)

        pairwise.tt_shannon_resids <- pairwise.t.test(x = DF_alpha$Shannon_resids, g = DF_alpha[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)

        if(var(DF_alpha$filtered_reads) > 0) {
                pairwise.tt_filtered_reads <- pairwise.t.test(x = DF_alpha$filtered_reads, g = DF_alpha[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
        } else {
                pairwise.tt_filtered_reads <- NA
        }

       # pairwise.tt_rich_resids_filt <- pairwise.t.test(x = DF_alpha$Richness_resids_FilteredReads, g = DF_alpha[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)

        pairwise.tt_shannon_resids_filt <- pairwise.t.test(x = DF_alpha$Shannon_resids_FilteredReads, g = DF_alpha[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)

        # add the linear models correcting for total or filtered amplicons
        fitter_rich_totalcounts <- lm(DF_alpha[["Richness"]] ~ DF_alpha[["Total"]] + DF_alpha[[group_var]])
        fitter_shannon_totalcounts <- lm(DF_alpha[["Shannon"]] ~ DF_alpha[["Total"]] + DF_alpha[[group_var]])
        fitter_rich_filteredreads <- lm(DF_alpha[["Richness"]] ~ DF_alpha[["filtered_reads"]] + DF_alpha[[group_var]])
        fitter_shannon_filteredreads <- lm(DF_alpha[["Shannon"]] ~ DF_alpha[["filtered_reads"]] + DF_alpha[[group_var]])
        fitter.list <- list(fitter_rich_totalcounts = fitter_rich_totalcounts, fitter_shannon_totalcounts = fitter_shannon_totalcounts,
                            fitter_shannon_filteredreads = fitter_shannon_filteredreads,
                            fitter_shannon_filteredreads = fitter_shannon_filteredreads)


        TrListBP <- boxplot_alphaDiv_fromDF(DF = DF_alpha, color = group_var, group = group_var, shape = shape, measures = c(alpha_div_measures2, paste(alpha_div_measures2, "resids", sep = "_"), paste(alpha_div_measures2, "resids", "FilteredReads", sep = "_")))

        TrList_lm <- plot_alphaDivVstotalCounts_fromList(DF_List = DF_alpha_list, measures = alpha_div_measures2, color = group_var, shape = shape)

        TrList_lm_filteredReads <- plot_alphaDivVsfilteredReads_fromList(DF_List = DF_alpha_list, measures = alpha_div_measures2, color = group_var, shape = shape)

        # out <- list(DF_alpha_list = DF_alpha_list, TrListBP = TrListBP, TrList_lm = TrList_lm, TrList_lm_filteredReads = TrList_lm_filteredReads, pairwise.tt_rich = pairwise.tt_rich,
        #             fitter_list = fitter.list
        #             )

        out <- list(DF_alpha_list = DF_alpha_list, TrListBP = TrListBP, TrList_lm = TrList_lm, TrList_lm_filteredReads = TrList_lm_filteredReads, pairwise.tt_rich = pairwise.tt_rich,
                    fitter_list = fitter.list, pairwise.tt_shannon = pairwise.tt_shannon, pairwise.tt_totalcounts = pairwise.tt_totalcounts, pairwise.tt_rich_resids = pairwise.tt_rich_resids,
                    pairwise.tt_shannon_resids = pairwise.tt_shannon_resids, pairwise.tt_filtered_reads = pairwise.tt_filtered_reads,
                    pairwise.tt_rich_resids_filt = pairwise.tt_rich_resids_filt, pairwise.tt_shannon_resids_filt = pairwise.tt_shannon_resids_filt)

 }


alpha_diversity_wrapper <- function(physeq, alpha_div_measures = c("Observed", "Shannon"), group_var1 = group_var, shape1 = shape, color_by=group_var){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
  DF_alpha_list <- calculate_alphadiversity(physeq = physeq, measures = alpha_div_measures)
  DF_alpha <- DF_alpha_list[[1]]
  
  # just I prefer Richness over Observed
  alpha_div_measures2 <- alpha_div_measures
  if ("Observed" %in% alpha_div_measures) {
    alpha_div_measures2[alpha_div_measures2 == "Observed"] <- "Richness" 
  }
  
  pairwise.tt_rich <- pairwise.t.test(x = DF_alpha$Richness, g = DF_alpha[[group_var1]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
  # compare, to get the same results as individual t.test you need pool.sd = F
   #t.test(x = DF_alpha$Richness[DF_alpha$Group == "Old"], y = DF_alpha$Richness[DF_alpha$Group == "MiddleAged"], alternative = "two", var.equal = F)$p.value

  pairwise.tt_shannon <- pairwise.t.test(x = DF_alpha$Shannon, g = DF_alpha[[group_var1]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)

   if(var(DF_alpha$Total) > 0) {
   pairwise.tt_totalcounts <- pairwise.t.test(x = DF_alpha$Total, g = DF_alpha[[group_var1]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
   } else {
     pairwise.tt_totalcounts <- NA
   }
  
  pairwise.tt_rich_resids <- pairwise.t.test(x = DF_alpha$Richness_resids, g = DF_alpha[[group_var1]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)

  pairwise.tt_shannon_resids <- pairwise.t.test(x = DF_alpha$Shannon_resids, g = DF_alpha[[group_var1]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
  
   if(var(DF_alpha$filtered_reads) > 0) {
    pairwise.tt_filtered_reads <- pairwise.t.test(x = DF_alpha$filtered_reads, g = DF_alpha[[group_var1]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
   } else {
     pairwise.tt_filtered_reads <- NA
   }

 pairwise.tt_rich_resids_filt <- pairwise.t.test(x = DF_alpha$Richness_resids_FilteredReads, g = DF_alpha[[group_var1]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)

 pairwise.tt_shannon_resids_filt <- pairwise.t.test(x = DF_alpha$Shannon_resids_FilteredReads, g = DF_alpha[[group_var1]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)

  #add the linear models correcting for total or filtered amplicons
  fitter_rich_totalcounts <- lm(DF_alpha[["Richness"]] ~ DF_alpha[["Total"]] + DF_alpha[[group_var1]])
  fitter_shannon_totalcounts <- lm(DF_alpha[["Shannon"]] ~ DF_alpha[["Total"]] + DF_alpha[[group_var1]])
  fitter_rich_filteredreads <- lm(DF_alpha[["Richness"]] ~ DF_alpha[["filtered_reads"]] + DF_alpha[[group_var1]])
  fitter_shannon_filteredreads <- lm(DF_alpha[["Shannon"]] ~ DF_alpha[["filtered_reads"]] + DF_alpha[[group_var1]])
  fitter.list <- list(fitter_rich_totalcounts = fitter_rich_totalcounts, fitter_shannon_totalcounts = fitter_shannon_totalcounts,
                      fitter_shannon_filteredreads = fitter_shannon_filteredreads, 
                      fitter_shannon_filteredreads = fitter_shannon_filteredreads)
  
  
  TrListBP <- boxplot_alphaDiv_fromDF(DF = DF_alpha, color = color_by, group = group_var1, shape = shape1, measures = c(alpha_div_measures2, paste(alpha_div_measures2, "resids", sep = "_"), paste(alpha_div_measures2, "resids", "FilteredReads", sep = "_")))
  
  TrList_lm <- plot_alphaDivVstotalCounts_fromList(DF_List = DF_alpha_list, measures = alpha_div_measures2, color = color_by, shape = shape1)
  
  TrList_lm_filteredReads <- plot_alphaDivVsfilteredReads_fromList(DF_List = DF_alpha_list, measures = alpha_div_measures2, color = color_by, shape = shape1)
  
  # out <- list(DF_alpha_list = DF_alpha_list, TrListBP = TrListBP, TrList_lm = TrList_lm, TrList_lm_filteredReads = TrList_lm_filteredReads,
  #             fitter_list = fitter.list)
  # 
   out <- list(DF_alpha_list = DF_alpha_list, TrListBP = TrListBP, TrList_lm = TrList_lm, TrList_lm_filteredReads = TrList_lm_filteredReads, pairwise.tt_rich = pairwise.tt_rich,
               fitter_list = fitter.list, pairwise.tt_shannon = pairwise.tt_shannon, pairwise.tt_totalcounts = pairwise.tt_totalcounts, pairwise.tt_rich_resids = pairwise.tt_rich_resids,
               pairwise.tt_shannon_resids = pairwise.tt_shannon_resids, pairwise.tt_filtered_reads = pairwise.tt_filtered_reads,
               pairwise.tt_rich_resids_filt = pairwise.tt_rich_resids_filt, pairwise.tt_shannon_resids_filt = pairwise.tt_shannon_resids_filt)
  
}







#######################################
### FUNCTION: arrange_p_values
#######################################
# see Generalized_Phyloseq_Analysis to understand it, it is just about putting a lot of p-value matrixes into one wide data.frame

arrange_p_values <- function(pairwise_list){
        pvalue_types <- strsplit(names(pairwise_list), split = "t_")
        pvalue_types <- sapply(pvalue_types, `[`, 2)
        pairwise_list <- lapply(pairwise_list, function(pairw) {
                pmat <- pairw$p.value
                rowCol <- expand.grid(rownames(pmat), colnames(pmat))
                labs <- rowCol[as.vector(lower.tri(pmat,diag=T)),]
                cbind(labs, p_value = pmat[lower.tri(pmat, diag = T)])
        })
        pairwise_list <- lapply(1:length(pairwise_list), function(i){cbind(pairwise_list[[i]], Type = pvalue_types[i])})
        pairwise_long <- do.call("rbind", pairwise_list)
        pairwise_wide <- spread(pairwise_long, key = Type, value = p_value)
}



#######################################
### FUNCTION: calc_distances
#######################################

calc_distances <- function(physeq, dist_methods = c("bray")) {
        
        dist_list <- vector("list", length(dist_methods))
        names(dist_list) = dist_methods
        
        for (i in dist_methods) {
                iDist <- phyloseq::distance(physeq, method=i)
                dist_list[[i]] = iDist
        }
        
        return(dist_list)
        
}

#######################################
### FUNCTION: calc_ordination_from_distances
#######################################

calc_ordination_from_distances <- function(physeq, dist_list, ordination_type = "PCoA", group_var = NULL, shape = NULL, coord_cor = FALSE, label=NULL, size=size){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
        ordination_list <- vector("list", length(dist_list))
        DFList <- vector("list", length(dist_list))
        DF_taxa_List <- vector("list", length(dist_list))
        # TrList <- vector("list", length(dist_list))
        TrList_own <- vector("list", length(dist_list))
        TrList_taxa <- vector("list", length(dist_list))
        
        axes <- 1:2 # currently only allowed to plot first and second
        
        for (i in seq_along(dist_list)) {
                
                ordination <- phyloseq::ordinate(physeq, method = ordination_type, distance = dist_list[[i]])
                ordination_list[[i]] <- ordination
                DF <- phyloseq::plot_ordination(physeq, ordination_list[[i]], color = group_var, justDF = TRUE, shape=shape, label=label)
                DFList[[i]] <- DF # just the first two axes cbind to sample_data in physeq
                
                x = colnames(DF)[1]
                y = colnames(DF)[2]
                Tr <- ggplot(DF, aes_string(x = x, y = y, col = group_var, shape = shape, label=label)) 
                Tr <- Tr + geom_point(size=size) +
                        scale_color_manual("", values = cbPalette[2:11]) +
                        theme_bw(12) +
                        ggtitle(names(dist_list)[i])
                
                # for labelling axes
                if (ordination_type == "PCoA" || ordination_type == "NMDS") {
                        if (ordination_type == "PCoA") {
                                eigvec <- phyloseq:::extract_eigenvalue.pcoa(ordination)
                        } else {
                                eigvec <- phyloseq:::extract_eigenvalue.default(ordination)
                        }
                        
                        if (length(eigvec[axes]) > 0){
                                fracvar = eigvec[axes]/sum(eigvec)
                                percvar = round(100 * fracvar, 1)
                                strivar = as(c(Tr$label$x, Tr$label$y), "character")
                                strivar = paste0(strivar, " (", percvar, " %)")
                                Tr <- Tr + xlab(strivar[1]) + ylab(strivar[2]) 
                        }
                        
                        if (!is.null(eigvec) && coord_cor) {
                                Tr <- Tr + coord_fixed(sqrt(eigvec[2] / eigvec[1]))
                        }
                        
                }
                 
                TrList_own[[i]] <- Tr
                rm(Tr)
                        
                
                # TrList[[i]] <- phyloseq::plot_ordination(physeq, ordination_list[[i]], color = group_var) + ggtitle(names(dist_list)[i])
                
                DF_taxa <- phyloseq::plot_ordination(physeq, ordination_list[[i]], type = "taxa", color = "Phylum", justDF = TRUE)
                DF_taxa_List[[i]] <- DF_taxa
                
                x = colnames(DF_taxa)[1]
                y = colnames(DF_taxa)[2]
                Tr <- ggplot(DF_taxa, aes_string(x = x, y = y, col = "Phylum")) 
                Tr <- Tr + geom_point(size=size) +
                        # scale_color_manual("", values = cbPalette[2:8]) +
                        theme_bw(12) +
                        ggtitle(names(dist_list)[i])
                
                # for labelling axes
                if (ordination_type == "PCoA" || ordination_type == "NMDS") {
                        if (ordination_type == "PCoA") {
                                eigvec <- phyloseq:::extract_eigenvalue.pcoa(ordination)
                        } else {
                                eigvec <- phyloseq:::extract_eigenvalue.default(ordination)
                        }
                        
                        if (length(eigvec[axes]) > 0){
                                fracvar = eigvec[axes]/sum(eigvec)
                                percvar = round(100 * fracvar, 1)
                                strivar = as(c(Tr$label$x, Tr$label$y), "character")
                                strivar = paste0(strivar, " (", percvar, " %)")
                                Tr <- Tr + xlab(strivar[1]) + ylab(strivar[2]) 
                        }
                        
                        if (!is.null(eigvec) && coord_cor) {
                                Tr <- Tr + coord_fixed(sqrt(eigvec[2] / eigvec[1]))
                        }
                        
                }
                
                TrList_taxa[[i]] <- Tr
                rm(Tr)
                
        }
        
        names(ordination_list) <- names(TrList_taxa) <- names(DFList) <- names(DF_taxa_List) <- names(TrList_own) <- names(dist_list)
        out <- list(ordination_list = ordination_list, DFList = DFList, DF_taxa_List = DF_taxa_List, ordination_Tr_samples = TrList_own, ordination_Tr_taxa = TrList_taxa)
}







#######################################
### FUNCTION: KeepAmplicons
#######################################

# Function to find the amplicons that are present in at least "Percentage" of samples
## Input
# taxa: the output of the assignTaxonomy (dada2) and maybe addSpecies command. nrow = number of sequences and ncol number taxonomic levels.
# seqtab: The abundance table output from Dada2_wrap
# percentage: only amplicons present in at least this percentage of samples will be kept
## Output
# taxseq: list with the shortened taxa and seqtab

KeepAmplicons <- function(taxa, seqtab, Percentage = 10){
        
        if(!identical(colnames(seqtab), rownames(taxa))){
                stop("taxa and seqtab do not fit togetehr")
        }
        
        
        PerCSampleValue <- ceiling((Percentage/100)*dim(seqtab)[1])
        KeepAmplis <- colnames(seqtab)[colSums(seqtab != 0) >= PerCSampleValue]
        seqtab.keep <- seqtab[,KeepAmplis]
        taxa.keep <- taxa[KeepAmplis,]
        taxseq <- list(taxa.keep, seqtab.keep)
        return(taxseq)
}






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

# Explanation of original gm function:
# Note zero.propagate just makes sure you get a 0 if there is any zero in your vector
# if there is no zero exp(mean(log(x), na.rm = TRUE)) == exp(sum(log(x[x > 0]), na.rm=TRUE) / length(x))
# I think this zero.propagate is silly, I therefore decided to make my own function, gm_own below

gm = function(x, na.rm=TRUE, zero.propagate = FALSE){
        if(any(x < 0, na.rm = TRUE)){
                return(NaN)
        }
        if(zero.propagate){
                if(any(x == 0, na.rm = TRUE)){
                        return(0)
                }
                exp(mean(log(x), na.rm = na.rm))
        } else {
                exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
        }
}



#######################################
### filttaxa_by_prevalence
#######################################
# filters taxa by prevalence and provides plots that illustrate how the sample_sums
# were changed by the filtering

## Input:
# physeq
# prevalence: in %, only take stay that are present in more than prevalence % of the samples unless MaxCountCheck is used
# MaxCountCheck, MaxCount: if max count is used a taxa can remain even if not fulfilling the prevalence
# criterion if its max count in one of the samples is bigger than MaxCount 
# the idea is that some samples might depend very much on a rare taxa

## Output:
# list with the filtered physeq, the sparsity measures, and 3 different plots (3-5)


filttaxa_by_prevalence <- function(physeq, prevalence = 30, MaxCountCheck = FALSE, MaxCount = 0.2*(min(sample_sums(physeq)))) {
        
        if (MaxCountCheck) {
                physeq_f <- filter_taxa(physeq, function(x){(sum(x != 0) > (prevalence/100)*length(x)) || (max(x) > MaxCount)}, prune = TRUE)
        } else {
                physeq_f <- filter_taxa(physeq, function(x){(sum(x != 0) > (prevalence/100)*length(x))}, prune = TRUE)
        }
        Sparsity <- 100*(sum(otu_table(physeq) == 0)/(ntaxa(physeq)*nsamples(physeq)))
        Sparsity_f <- 100*(sum(otu_table(physeq_f) == 0)/(ntaxa(physeq_f)*nsamples(physeq_f)))
        PCCountsRemoved <- 100*(1-(sample_sums(physeq_f)/sample_sums(physeq)))
        DF <- data.frame(Sample = names(PCCountsRemoved), PCRemoved = PCCountsRemoved)
        DF$Sample <- as.factor(DF$Sample)
        DF <- dplyr::arrange(DF, PCRemoved)
        # relevel so samples are shown from min NoReads to max NoReads
        LevelsWant <- as.character(DF$Sample) 
        for (i in 1:length(LevelsWant)) {
                DF$Sample <- relevel(DF$Sample, ref = LevelsWant[i])
        }
        
        DF2 <- data.frame(before = sample_sums(physeq), after = sample_sums(physeq_f))
        DF2$Sample <- rownames(DF2)
        DF2$Sample <- as.factor(DF2$Sample)
        DF2 <- dplyr::arrange(DF2, before)
        LevelsWant <- as.character(DF2$Sample) 
        for (i in 1:length(LevelsWant)) {
                DF2$Sample <- relevel(DF2$Sample, ref = LevelsWant[i])
        }
        DF2 <- tidyr::gather(DF2, key = physeq, value = count, -Sample)
        
        Tr <- ggplot(DF, aes(x = Sample, y = PCRemoved))
        Tr <- Tr + geom_point(col = "#E69F00", size = 2) +
                ylab("% of amplicon reads removed") +
                xlab("") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
        
        Tr2 <- ggplot(data = DF, aes(x = PCRemoved))
        Tr2 <- Tr2 + geom_histogram(binwidth = 2, col = "black", fill = "#E69F00") +
                geom_rug() +
                ylab("No Samples") + 
                xlab("% of amplicon reads removed") +
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
        
        Tr3 <- ggplot(data = DF2, aes(x = Sample, y = count, color = physeq))
        Tr3 <- Tr3 + geom_point(size = 2) +
                ylab("sample_sums") +
                xlab("") +
                scale_color_manual(values = cbPalette[c(4,2)]) +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                      legend.title = element_blank())
        
        
        list(FilteredPhyseq = physeq_f, Sparsity = data.frame(Before = Sparsity, After = Sparsity_f), PlotReadsRemoved = Tr,
             HistReadsRemoved = Tr2, PlotSampleSumsBeforeAfter = Tr3)
        
}




############################################
## simulate_totalabVSrichness_rarefaction ##
############################################
# Input:
# physeq is needed to simulate the sample, i.e. the composition of the original sample
# specifically: the sample in physeq with the highest richness is used as base_sample (only the SVs > 0)
# no_low_extra: defines the Number of low abundance SVs that will be added to the base_sample, the "counts"
# of these extra SVs are random numbers between 0 and min_count in base_sample
# NB: the low_extra SVs make sure that the generated samples have more low abundance counts than the base_sample, and
# therefore plateau slower in rarefaction curves
# NB2: you could of course think about other ways to generate the original proportion
# the idea of the simulation is then: 1 million amplicons of the original proportion go to the sequencer (DNA)
# S1 is size1 amplicons from this original
# S2 is size2 amplicons from this original
# S3 is size2 amplicons rarefied from S1
# richnesses are compared and rarefaction curves are generated, all is saved in the output list

# further Inputs:
# type: "sample" or "vegan" decides whether the rarefaction steps are based on sample() or vegan::rrarefy command
# sd: is used to spread the low_extra_SVs "counts" in an rnorm command, the bigger the higher the spread, but note anyway bewtween 0 and min_count
# no_DNA_total_amplicons: "the number of amplicons generated and put in the seqeuncer", the higher the more secure that the richness of the
# DNA sample was 100%


simulate_totalabVSrichness_rarefaction <- function(physeq, size1 = 50000, size2 = 15000, nsims = 100, no_low_extra = 92, type = "sample", seed = 1576,
                                                   sd = 2, no_DNA_total_amplicons = 1e6){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
        set.seed(seed)
        seqtab <- as(otu_table(physeq), "matrix")
        base_sample <- seqtab[which.max(rowSums(seqtab != 0)),]
        min_count <- min(base_sample[base_sample > 0])
        total_SVs <- no_low_extra + sum(base_sample > 0)
        
        DNA_richness <- vector("numeric", length = nsims)
        S1_richness <- vector("numeric", length = nsims)
        S2_richness <- vector("numeric", length = nsims)
        S3_richness <- vector("numeric", length = nsims)
        S1s <- matrix(nrow = nsims, ncol = total_SVs)
        
        for (i in 1:nsims) {
                
                factors <- abs(rnorm(n = no_low_extra, sd = 2))
                factors <- factors/max(factors)
                extra_low_SVs <- factors*min_count
                
                sample_prop <- c(extra_low_SVs, base_sample[base_sample > 0])
                sample_prop <- sample_prop/sum(sample_prop)
                
                DNA_sample <- round(sample_prop * no_DNA_total_amplicons)
                DNA_richness[i] <- sum(DNA_sample > 0)
                
                if (type == "vegan") {
                        S1 <- rrarefy(matrix(DNA_sample, nrow = 1), sample = size1)
                        S2 <- rrarefy(matrix(DNA_sample, nrow = 1), sample = size2)
                        S3 <- rrarefy(S1, sample = size2)
                        S1 <- as.vector(S1)
                        S2 <- as.vector(S2)
                        S3 <- as.vector(S3)
                        S1_richness[i] <- sum(S1 > 0)
                        S2_richness[i] <- sum(S2 > 0)
                        S3_richness[i] <- sum(S3 > 0)
                        S1s[i,] <- S1
                        
                } else if (type == "sample") {
                        
                        S1 <- rarefy_sample(DNA_sample, size = size1)
                        S2 <- rarefy_sample(DNA_sample, size = size2)
                        S3 <- rarefy_sample(S1, size = size2)
                        S1_richness[i] <- sum(S1 > 0)
                        S2_richness[i] <- sum(S2 > 0)
                        S3_richness[i] <- sum(S3 > 0)
                        S1s[i,] <- S1
                        
                } else {
                        stop("wrong type")
                }
                
        }
        
        DF <- data.frame(DNA = DNA_richness, S1 = S1_richness, S2 = S2_richness, S3 = S3_richness)
        
        DFl <- gather(DF, key = sample, value = richness)
        DFl$Sample <- factor(DFl$sample, levels <- c("DNA", "S1", "S2", "S3"), ordered = TRUE)
        
        pair_tt <- pairwise.t.test(x = DFl$richness, g = DFl$Sample, alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
        
        Tr <- ggplot(DFl, aes(x = sample, y = richness, col = sample))
        
        Tr <- Tr +
                geom_boxplot() +
                geom_jitter(alpha = .7) +
                scale_color_manual("", values = cbPalette[2:11]) +
                theme_bw(12)
        
        
        # add rarefaction curves for 5 of the S1s plus the base_sampled supplied with 0s for the low_extra
        CT <- rbind(S1s[sample(nsims, 5),], c(rep(0, no_low_extra), base_sample[base_sample > 0]))
        samdf <- sample_data(ps)
        samdf <- samdf[1:nrow(CT),]
        rownames(CT) <- rownames(samdf)
        colnames(CT) <- paste("T_", 1:ncol(CT), sep = "")
        pss <- phyloseq(otu_table(CT, taxa_are_rows = FALSE), 
                        sample_data(samdf))
        
        Curves <- rarefaction_curve_own_fast(physeq = pss, group_var = NULL)
        
        # show distribution of SVs with abundance below 20 and above 0
        low_abund_SVs <- apply(CT, 1, function(x){x[x > 0 & x < 20]})
        df_plot <- data.frame(Sample = rep(names(low_abund_SVs), sapply(low_abund_SVs, length)),
                              Abund = unlist(low_abund_SVs))
        Trr <- ggplot(df_plot, aes(x = Abund, group = Sample, col = Sample))
        Trr <- Trr + geom_density() +
                scale_color_manual("", values = cbPalette[2:11], labels = c(paste("Sim", 1:5), "dada2_sam")) +
                xlab("Abundances below 20") +
                theme_bw(20)
        
        out <- list(S1s = S1s, DF_richness = DF, pair_tt_of_richness = pair_tt, Tr = Tr, CT = CT, Curves = Curves, Trr = Trr, sample_prop = sample_prop, base_sample = base_sample)
        
}



#######################################
### FUNCTION: distance_t_analyse
#######################################
# INPUT:
# dist_list: named list of dist objects (so for each distance one object)
# physeq: phyloseq object used to make the dist objects
# group_var: character identifying the grouping variable in the sample_data of the phyloseq
# OUTPUT:
# list of two named lists, one with the plots and one with the data_frames of the p_values from pairwise t.tests

distance_t_analyse <- function(dist_list, physeq, group_var) {
        
        TrList <- vector(mode = "list", length = length(dist_list))
        pValList <- vector(mode = "list", length = length(dist_list))
        
        for (i in 1:length(dist_list)){
                DistMat <- as(dist_list[[i]], "matrix")
                rowCol <- expand.grid(rownames(DistMat), colnames(DistMat))
                labs <- rowCol[as.vector(lower.tri(DistMat, diag=F)),]
                df <- cbind(labs, DistMat[lower.tri(DistMat, diag=F)])
                colnames(df) <- c("Row","Col","Distance")
                samdf <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
                df$Row_Group <- samdf$Group[match(df$Row, samdf$Sample)]
                df$Col_Group <- samdf$Group[match(df$Col, samdf$Sample)]
                
                # add GroupvsGroup using sort (order + match) to make sure that Old vs Young and Young vs Old both become Young vs Old
                df$GroupvsGroup <- apply(df[, c("Row_Group", "Col_Group")], 1, function(x){
                        paste(x[order(match(x, levels(samdf$Group)))], collapse = " vs ")
                })
                
                # to make sure group comparisons are in the "order" of levels "group_var"
                df$GroupvsGroupOrder <- apply(df[, c("Row_Group", "Col_Group")], 1, function(x){
                        x1 <- match(x[1], levels(samdf$Group))
                        x2 <- match(x[2], levels(samdf$Group))
                        if (x1 <= x2) {
                                as.numeric(paste(x1, x2, sep = ""))
                        } else {
                                as.numeric(paste(x2, x1, sep = ""))
                        }
                })
                
                df$Type <- "between"
                df$Type[df$Row_Group == df$Col_Group] <- "within"
                # NB: for within group distances, you have samples*(samples-1)/2 distances, for between group distances
                # you have samplesgrp1 * samplesgrp2 distances.
                
                GvsGLeveldf <- unique(df[, c("GroupvsGroup", "GroupvsGroupOrder", "Type")])
                GroupvsGroupLevels <- GvsGLeveldf[order(GvsGLeveldf$GroupvsGroupOrder), ]$GroupvsGroup
                GvsGLeveldfBetween <- GvsGLeveldf[GvsGLeveldf$Type == "between", ]
                GroupvsGroupLevelsBetween <- GvsGLeveldfBetween[order(GvsGLeveldfBetween$GroupvsGroupOrder), ]$GroupvsGroup
                df$GroupvsGroup <- factor(df$GroupvsGroup, levels = GroupvsGroupLevels, ordered = TRUE)
                
                pairwise_dist_ttest <- pairwise.t.test(x = df$Distance, g = df$GroupvsGroup, alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
                pairwise_dist_ttest <- pairwise_dist_ttest$p.value
                rowCol <- expand.grid(rownames(pairwise_dist_ttest), colnames(pairwise_dist_ttest))
                labs <- rowCol[as.vector(lower.tri(pairwise_dist_ttest,diag=T)),]
                df_p <- cbind(labs, pairwise_dist_ttest[lower.tri(pairwise_dist_ttest,diag=T)])
                colnames(df_p) <- c("Distances_2","Distances_1","p_value")
                df_p <- df_p[, c("Distances_1", "Distances_2", "p_value")]
                df_p$Significance <- ""
                df_p$Significance[df_p$p_value <= 0.05] <- "*"
                df_p$Significance[df_p$p_value <= 0.01] <- "**"
                df_p$Significance[df_p$p_value <= 0.001] <- "***"
                
                pValList[[i]] <- df_p
                
                
                # in case there are more than two groups in group_var I want a faceted plot, i.e. one facet for each GroupvsGroupLevelsBetween, for this I need to duplicate some data
                
                if (length(GroupvsGroupLevels) > 1) {
                        df_list <- lapply(GroupvsGroupLevelsBetween, function(level){
                                grps <- unlist(strsplit(level, " vs "))
                                df_current <- df[df$Row_Group %in% grps & df$Col_Group %in% grps, ]
                                df_current$Level <- level
                                df_current
                        })
                        df_plot <- do.call(rbind, df_list)
                        df_plot$Level <- factor(df_plot$Level, levels = GroupvsGroupLevelsBetween, ordered = TRUE)
                        
                        Tr <- ggplot(df_plot, aes(x = GroupvsGroup, y = Distance, col = GroupvsGroup))
                        Tr <- Tr + geom_boxplot(outlier.color = NA) + 
                                #geom_jitter(position = position_jitter(width = 0.3, height = 0), alpha = 0.65) +
                                facet_wrap(~ Level, ncol = 2, scales = "free_x") +
                                theme_bw() +
                                xlab("") +
                                ylab(paste(names(dist_list)[i], "distance", sep = " ")) +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                              legend.title = element_blank())
                        if (length(GroupvsGroupLevels) <= 7) {
                                Tr <- Tr +
                                        scale_color_manual(values = cbPalette[2:11])
                        }
                        
                } else {
                        Tr <- ggplot(df, aes(x = GroupvsGroup, y = Distance, col = GroupvsGroup))
                        Tr <- Tr + geom_boxplot(outlier.color = NA) + 
                               # geom_jitter(position = position_jitter(width = 0.3, height = 0), alpha = 0.65) +
                                theme_bw() +
                                xlab("") +
                                ylab(paste(names(dist_list)[i], "distance", sep = " ")) + 
                                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                                      legend.title = element_blank())
                        if (length(GroupvsGroupLevels) <= 7) {
                                Tr <- Tr +
                                        scale_color_manual(values = cbPalette[2:11])
                        }
                        
                }
                
                TrList[[i]] <- Tr
        }
        
        names(TrList) <- names(dist_list)
        names(pValList) <- names(dist_list)
        out <- list(DistanceBoxplots = TrList, DistancePValues = pValList)
        
}

distance_t_analyse_shape <- function(dist_list, physeq, group_var, shape) {
  
  TrList <- vector(mode = "list", length = length(dist_list))
  pValList <- vector(mode = "list", length = length(dist_list))
  
  for (i in 1:length(dist_list)){
    DistMat <- as(dist_list[[i]], "matrix")
    rowCol <- expand.grid(rownames(DistMat), colnames(DistMat))
    labs <- rowCol[as.vector(lower.tri(DistMat, diag=F)),]
    df <- cbind(labs, DistMat[lower.tri(DistMat, diag=F)])
    colnames(df) <- c("Row","Col","Distance")
    samdf <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]], Group2 = sample_data(physeq)[[shape]])
    df$Row_Group <- samdf$Group[match(df$Row, samdf$Sample)]
    df$Col_Group <- samdf$Group[match(df$Col, samdf$Sample)]
    df$Row_Group2 <- samdf$Group2[match(df$Row, samdf$Sample)]
    df$Col_Group2 <- samdf$Group2[match(df$Col, samdf$Sample)]
    
    # add GroupvsGroup using sort (order + match) to make sure that Old vs Young and Young vs Old both become Young vs Old
    df$GroupvsGroup <- apply(df[, c("Row_Group", "Col_Group")], 1, function(x){
      paste(x[order(match(x, levels(samdf$Group)))], collapse = " vs ")
    })
    
    # to make sure group comparisons are in the "order" of levels "group_var"
    df$GroupvsGroupOrder <- apply(df[, c("Row_Group", "Col_Group")], 1, function(x){
      x1 <- match(x[1], levels(samdf$Group))
      x2 <- match(x[2], levels(samdf$Group))
      if (x1 <= x2) {
        as.numeric(paste(x1, x2, sep = ""))
      } else {
        as.numeric(paste(x2, x1, sep = ""))
      }
    })
    
    
    
    
    df$Type <- "between"
    df$Type[df$Row_Group == df$Col_Group] <- "within"
    # NB: for within group distances, you have samples*(samples-1)/2 distances, for between group distances
    # you have samplesgrp1 * samplesgrp2 distances.
    
    GvsGLeveldf <- unique(df[, c("GroupvsGroup", "GroupvsGroupOrder", "Type")])
    GroupvsGroupLevels <- GvsGLeveldf[order(GvsGLeveldf$GroupvsGroupOrder), ]$GroupvsGroup
    GvsGLeveldfBetween <- GvsGLeveldf[GvsGLeveldf$Type == "between", ]
    GroupvsGroupLevelsBetween <- GvsGLeveldfBetween[order(GvsGLeveldfBetween$GroupvsGroupOrder), ]$GroupvsGroup
    df$GroupvsGroup <- factor(df$GroupvsGroup, levels = GroupvsGroupLevels, ordered = TRUE)
    
    pairwise_dist_ttest <- pairwise.t.test(x = df$Distance, g = df$GroupvsGroup, alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
    pairwise_dist_ttest <- pairwise_dist_ttest$p.value
    rowCol <- expand.grid(rownames(pairwise_dist_ttest), colnames(pairwise_dist_ttest))
    labs <- rowCol[as.vector(lower.tri(pairwise_dist_ttest,diag=T)),]
    df_p <- cbind(labs, pairwise_dist_ttest[lower.tri(pairwise_dist_ttest,diag=T)])
    colnames(df_p) <- c("Distances_2","Distances_1","p_value")
    df_p <- df_p[, c("Distances_1", "Distances_2", "p_value")]
    df_p$Significance <- ""
    df_p$Significance[df_p$p_value <= 0.05] <- "*"
    df_p$Significance[df_p$p_value <= 0.01] <- "**"
    df_p$Significance[df_p$p_value <= 0.001] <- "***"
    
    pValList[[i]] <- df_p
    
    
    # in case there are more than two groups in group_var I want a faceted plot, i.e. one facet for each GroupvsGroupLevelsBetween, for this I need to duplicate some data
    
    if (length(GroupvsGroupLevels) > 1) {
      df_list <- lapply(GroupvsGroupLevelsBetween, function(level){
        grps <- unlist(strsplit(level, " vs "))
        df_current <- df[df$Row_Group %in% grps & df$Col_Group %in% grps, ]
        df_current$Level <- level
        df_current
      })
      df_plot <- do.call(rbind, df_list)
      df_plot$Level <- factor(df_plot$Level, levels = GroupvsGroupLevelsBetween, ordered = TRUE)
      
      Tr <- ggplot(df_plot, aes(x = GroupvsGroup, y = Distance, col = GroupvsGroup, shape=Group2))
      Tr <- Tr + geom_boxplot(outlier.color = NA) + 
        geom_jitter(position = position_jitter(width = 0.3, height = 0), alpha = 0.65) +
        facet_wrap(~ Level, ncol = 2, scales = "free_x") +
        theme_bw() +
        xlab("") +
        ylab(paste(names(dist_list)[i], "distance", sep = " ")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              legend.title = element_blank())
      if (length(GroupvsGroupLevels) <= 7) {
        Tr <- Tr +
          scale_color_manual(values = cbPalette[2:11])
      }
      
    } else {
      Tr <- ggplot(df, aes(x = GroupvsGroup, y = Distance, col = GroupvsGroup,  shape=Group2))
      Tr <- Tr + geom_boxplot(outlier.color = NA) + 
        geom_jitter(position = position_jitter(width = 0.3, height = 0), alpha = 0.65) +
        theme_bw() +
        xlab("") +
        ylab(paste(names(dist_list)[i], "distance", sep = " ")) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              legend.title = element_blank())
      if (length(GroupvsGroupLevels) <= 7) {
        Tr <- Tr +
          scale_color_manual(values = cbPalette[2:11])
      }
      
    }
    
    TrList[[i]] <- Tr
  }
  
  names(TrList) <- names(dist_list)
  names(pValList) <- names(dist_list)
  out <- list(DistanceBoxplots = TrList, DistancePValues = pValList)
  
}


##############
distance_t_analyse_within <- function(dist_list, physeq, group_var) {
  
  TrList <- vector(mode = "list", length = length(dist_list))
  TrList_bi <- vector(mode = "list", length = length(dist_list))
  pValList <- vector(mode = "list", length = length(dist_list))
  
  for (i in 1:length(dist_list)){
    DistMat <- as(dist_list[[i]], "matrix")
    rowCol <- expand.grid(rownames(DistMat), colnames(DistMat))
    labs <- rowCol[as.vector(lower.tri(DistMat, diag=F)),]
    df <- cbind(labs, DistMat[lower.tri(DistMat, diag=F)])
    colnames(df) <- c("Row","Col","Distance")
    samdf <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
    df$Row_Group <- samdf$Group[match(df$Row, samdf$Sample)]
    df$Col_Group <- samdf$Group[match(df$Col, samdf$Sample)]
    
    # add GroupvsGroup using sort (order + match) to make sure that Old vs Young and Young vs Old both become Young vs Old
    df$GroupvsGroup <- apply(df[, c("Row_Group", "Col_Group")], 1, function(x){
      paste(x[order(match(x, levels(samdf$Group)))], collapse = " vs ")
    })
    
    # to make sure group comparisons are in the "order" of levels "group_var"
    df$GroupvsGroupOrder <- apply(df[, c("Row_Group", "Col_Group")], 1, function(x){
      x1 <- match(x[1], levels(samdf$Group))
      x2 <- match(x[2], levels(samdf$Group))
      if (x1 <= x2) {
        as.numeric(paste(x1, x2, sep = ""))
      } else {
        as.numeric(paste(x2, x1, sep = ""))
      }
    })
    
    df$Type <- "between"
    df$Type[df$Row_Group == df$Col_Group] <- "within"
    GvsGLeveldf <- unique(df[, c("GroupvsGroup", "GroupvsGroupOrder", "Type")])
    GroupvsGroupLevels <- GvsGLeveldf[order(GvsGLeveldf$GroupvsGroupOrder), ]$GroupvsGroup
    GvsGLeveldfBetween <- GvsGLeveldf[GvsGLeveldf$Type == "within", ]
    GroupvsGroupLevelsBetween <- GvsGLeveldfBetween[order(GvsGLeveldfBetween$GroupvsGroupOrder), ]$GroupvsGroup
    df$GroupvsGroup <- factor(df$GroupvsGroup, levels = GroupvsGroupLevels, ordered = TRUE)
    
   
  
      df_list <- lapply(GroupvsGroupLevelsBetween, function(level){
        grps <- unlist(strsplit(level, " vs "))
        df_current <- df[df$Row_Group %in% grps & df$Col_Group %in% grps, ]
        df_current$Level <- level
        df_current
      })
      
      df_plot <- do.call(rbind, df_list)
      df_plot$SamplevsSample <- apply(df_plot[, c("Row", "Col")], 1, function(x){
        paste(x[order(match(x, levels(df_plot$Row)))], collapse = " vs ")
      })
      
      df_plot$BrayIndex<-(100-100*(df_plot$Distance))
      df_plot$Level <- factor(df_plot$Level, levels = GroupvsGroupLevelsBetween, ordered = TRUE)
      
      Tr <- ggplot(df_plot, aes(x = Col_Group, y = Distance, col = Col_Group))
      #Tr <- Tr + geom_boxplot(outlier.color = NA) + 
      Tr <- Tr + geom_point(aes(colour = factor(Col_Group))) +
        #geom_jitter(position = position_jitter(width = 0.3, height = 0), alpha = 0.65) +
       # facet_wrap(~ Level, ncol = 2, scales = "free_x") +
        theme_bw() +
        xlab("") +
        ylab(paste(names(dist_list)[i], "distance", sep = " ")) +
        theme(axis.text.x = element_text(angle = 360, hjust = 0.5, vjust = 0.5),
              legend.position = "none")  + scale_color_manual(values = cbPalette[2:11])#+scale_color_brewer(palette = "Dark2")
      if (length(GroupvsGroupLevels) <= 7) {
        Tr <- Tr + scale_color_manual(values = cbPalette[2:11]) #+scale_color_brewer(palette = "Dark2")
          #scale_color_manual(values = cbPalette[2:11])
      }
      
    
   
    TrList[[i]] <- Tr
    
  
    Tr_bi <- ggplot(df_plot, aes(x = Col_Group, y = BrayIndex, col = Col_Group))
    #Tr_bi <- Tr_bi + geom_boxplot(outlier.color = NA) + 
    Tr_bi <- Tr_bi + geom_point(aes(colour = factor(Col_Group))) +
      #geom_jitter(position = position_jitter(width = 0.3, height = 0), alpha = 0.65) +
      # facet_wrap(~ Level, ncol = 2, scales = "free_x") +
      theme_bw() +
      xlab("") +
      ylab(paste(names(dist_list)[i], "Index [% similarity]", sep = " ")) +
      theme(axis.text.x = element_text(angle = 360, hjust = 0.5, vjust = 0.5),
            legend.position = "none")  + scale_color_manual(values = cbPalette[2:11])#+scale_color_brewer(palette = "Dark2")
    if (length(GroupvsGroupLevels) <= 7) {
      Tr_bi <- Tr_bi  + scale_color_manual(values = cbPalette[2:11])#+scale_color_brewer(palette = "Dark2")
        #scale_color_manual(values = cbPalette[2:11])
    }
  
  TrList_bi[[i]] <- Tr_bi
    
  }
  names(TrList) <- names(dist_list)
  names(pValList) <- names(dist_list)
  out <- list(DistanceBoxplots = TrList[[1]], DistanceBoxplots_brayIndex=TrList_bi[[1]],  brayIndex=df_plot)
  #DistancePValues = pValList, 
}







#######################################
### FUNCTION: pairwise.perm.manova.own
#######################################
# Function is very much based on pairwise.perm.manova {RVAideMemoire}
# but it also records R2 while looping through vegan::adonis, and generates a
# result data frame in which the results are shown in the order of the group_fac levels
# INPUT:
# dist_obj: dist object
# group_fac: the factor that groups the samples
# nperm: permutations in adonis
# p.adj.method: method to adjust p.values
# OUTPUT:
# data.frame showing p.values and R2 and adjusted p.values for the different between group comparisons


pairwise.perm.manova.own <- function(dist_obj, group_fac, nperm = 999, 
                                     p.adj.method = "none") {
        
        if (!("dist" %in% class(dist_obj))){
                stop("dist_obj must be of class dist")
        }
        
        group_fac <- factor(group_fac)
        
        fac_levels <- levels(group_fac)
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) # see pairwise.table
        
        i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        i <- ivec[x]
                })
        })
        j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        j <- jvec[x]
                })
        })
        i_s <- i_s[upper.tri(i_s)]
        j_s <- j_s[upper.tri(j_s)]
        
        p_vals <- vector(mode = "numeric", length = length(i_s))
        r2s <- vector(mode = "numeric", length = length(i_s))
        
        for (k in seq_along(i_s)){
                i <- i_s[k]
                j <- j_s[k]
                group_fac2 <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])
                dist_obj_mat <- as.matrix(dist_obj)
                rows <- which(group_fac %in% levels(group_fac2))
                dist_obj2 <- as.dist(dist_obj_mat[rows, rows])
                fit <- vegan::adonis(dist_obj2 ~ group_fac2, permutations = nperm)
                p_vals[k] <- fit$aov.tab[1, "Pr(>F)"]
                r2s[k] <- fit$aov.tab[1, "R2"]
        }
        
        
        
        result_df <- data.frame(Comparison = paste0(names(fac_levels_num[i_s]), "_vs_", names(fac_levels_num[j_s]), sep = ""),
                                addonis_p_value = p_vals, adonis_R2 = r2s, p_val_adj = p.adjust(p_vals, p.adj.method))
        
}


#######################################
### get_overview_of_physeq##
#################
# The function came to live while doing simulations: it became clear that the count variation (both across samples and taxa) is
# usually higher in real data than in simulated data. 
# Input: just physeq
# Output: a DF summarising features of the abundance table. 

get_overview_of_physeq <- function(physeq){
        
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        # I prefer taxa_are_rows = FALSE so rows (= observations = samples), and colums = variables = taxa
        
        CT <- as(otu_table(physeq), "matrix") # taxa are columns
        CT_NA <- CT
        CT_NA[CT_NA == 0] <- NA
        # also get a relative abundance Table
        CT_RA <- CT/rowSums(CT)
        CT_RA_NA <- CT_RA
        CT_RA_NA[CT_RA_NA == 0] <- NA
        
        
        Overview <- data.frame(NoSamples = nsamples(physeq),
                               NoTaxa = ntaxa(physeq),
                               MedianSampleSum = round(median(sample_sums(physeq))),
                               MedianTaxaSum = round(median(taxa_sums(physeq))),
                               Sparsity = round(100*(sum(otu_table(physeq) == 0)/(ntaxa(physeq)*nsamples(physeq))), 3),
                               MaxCount = max(CT),
                               MedianTaxaSD = round(median(apply(CT, 2, sd), na.rm = TRUE), 3),
                               MedTaxaSDNoZ = round(median(apply(CT_NA, 2, sd, na.rm = TRUE), na.rm = TRUE), 3),
                               
                               # NB: SD correlates with taxa_sums (= colSums(CT))
                               # plot(colSums(CT), apply(CT, 2, sd))
                               # plot(colSums(CT_NA, na.rm = TRUE), apply(CT_NA, 2, sd, na.rm = TRUE))
                               # in principle you could correct by taxa_sums
                               # MedianTaxaSD_Cor = round(median(apply(CT*(max(colSums(CT))/colSums(CT)), 2, sd), na.rm = TRUE), 3),
                               # just keep in mind it is confounded
                               
                               # get same taxa variation for relative abundance table
                               MedianTaxaSD_RA = round(median(apply(CT_RA, 2, sd), na.rm = TRUE), 6),
                               MedTaxaSDNoZ_RA = round(median(apply(CT_RA_NA, 2, sd, na.rm = TRUE), na.rm = TRUE), 6),
                               
                               MedianSampleSD = round(median(apply(CT, 1, sd), na.rm = TRUE), 3),
                               MedianSampleSDNoZ = round(median(apply(CT_NA, 1, sd, na.rm = TRUE), na.rm = TRUE), 3),
                               # makes in principle only sense on relative abundance
                               MedianSampleSD_RA = round(median(apply(CT_RA, 1, sd), na.rm = TRUE), 6),
                               MedianSampleSDNoZ_RA = round(median(apply(CT_RA_NA, 1, sd, na.rm = TRUE), na.rm = TRUE), 3)
        )
        
}


#######################################
### make_heat_map_physeq##
#################
## Input:
# physeq object
# group_var: the name of the group_fac column in sample_data used to order the samples in the heat map
# max_abundance_for_color: if null the 90% quantile of count/relative abundance in the data is used. all counts above this value will be
# shown yellow in the heat map
# tax_order: character vector of the original taxon names in the order you want them shown. if NULL >> tax_order = taxa_names(physeq)
# tax_names: the names that will be used for the taxons, if Null Taxon_1, Taxon_2 and so on will be used. NB: renaming of course after
# ordering. 
# color_sample_names: if TRUE and if you have less than 7 levels, the sample names will be colored using cbPalette
# gradient_steps: the steps the blue to green to yellow will be distributed in the viridis gradient: 4 numbers ranging from 1e-14 to 1,
# see default, you might wanna try c(0.25, 0.5, 0.75, 1) as well
## Output:
# the heat map trellis

make_heat_map_physeq <- function(physeq, group_var, max_abundance_for_color = NULL, tax_order = NULL,
                                 tax_names = NULL, color_sample_names = TRUE, gradient_steps = c(0.15, 0.3, 0.45, 1)){
        
        gradient_steps <- c(0, 1e-14, gradient_steps)
        
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        
        if (!is.factor(sample_data(physeq)[[group_var]])) {sample_data(physeq)[[group_var]] <- as.factor(sample_data(physeq)[[group_var]])}
        
        
        CT <- as(otu_table(physeq), "matrix") # taxa are columns
        DF_CT <- as.data.frame(t(CT))
        
        if (is.null(tax_order)){tax_order <- taxa_names(physeq)}
        
        if (length(tax_order) != nrow(DF_CT) || !all(tax_order %in% rownames(DF_CT))) {stop("the given tax_order must contain all taxa_names(physeq). If you want to subset,
                                                                                            do it before with prune_taxa")}
        DF_CT <- DF_CT[tax_order, ]
        
        if (is.null(tax_names)) {tax_names <- paste("Taxon_", 1:nrow(DF_CT), sep = "")}
        
        if (length(tax_names) != nrow(DF_CT) ) {stop("the given tax_names must be a character vector of length = ntaxa(physeq)")}
        
        
        rownames(DF_CT) <- make.unique(tax_names)
        DF_CT$Taxa <- rownames(DF_CT)
        DF_CT <- tidyr::gather(DF_CT, key = Sample , value = Count, -Taxa)
        DF_CT$Taxa <- factor(DF_CT$Taxa, levels = rev(unique(DF_CT$Taxa)), ordered = TRUE)
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        DF_CT$Sample <- factor(DF_CT$Sample, levels = LookUpDF$Sample, ordered = TRUE)
        
        if (is.null(max_abundance_for_color)) {
                max_abundance_for_color <- quantile(DF_CT$Count, .9)
        }
        
        # Color the sample names based on levels in group fac
        if (length(levels(LookUpDF$Group)) <= 10 && color_sample_names){
                color_lookup <- data.frame(level = levels(LookUpDF$Group), color = cbPalette[2:(length(levels(LookUpDF$Group)) + 1)])
                colxaxis <- as.character(color_lookup$color[match(LookUpDF$Group, color_lookup$level)])
        } else {
                colxaxis <- rep("black", nrow(LookUpDF))
        }
        
        
        hmTr <- ggplot(DF_CT, aes(x = Sample, y = Taxa, fill = Count))
        hmTr <- hmTr + 
                geom_raster() + 
                scale_fill_gradientn("", limits = c(0, max_abundance_for_color), colors = c("white", viridis(5)), values = gradient_steps, oob = squish) +
                scale_x_discrete(position = "top") +
                #coord_equal() +
                labs(x=NULL, y=NULL) +
                theme_tufte(base_family = "Helvetica") +
                theme(axis.ticks=element_blank(),
                      axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0,
                                                 colour = colxaxis))
        hmTr
}

#######################################
### make_heat_map_physeq_levels##
#################
# same heat map plots as from make_heat_map_physeq but for each combination of levels in your physeq
## Input:
# physeq object
# group_var: the name of the group_fac column in sample_data used to order the samples in the heat map
# max_abundance_for_color: if null the 90% percentile count/relative abundance in the data is used. all counts above this value will be
# shown yellow in the heat map
# tax_order: character vector of the original taxon names in the order you want them shown. if NULL >> tax_order = taxa_names(physeq)
# tax_names: the names that will be used for the taxons, if Null Taxon_1, Taxon_2 and so on will be used. NB: renaming of course after
# ordering. 
# color_sample_names: if TRUE and if you have less than 7 levels, the sample names will be colored using cbPalette
# gradient_steps: the steps the blue to green to yellow will be distributed in the viridis gradient: 4 numbers ranging from 1e-14 to 1,
# see default, you might wanna try c(0.25, 0.5, 0.75, 1) as well
# Output: list of heat maps for each level combination

make_heat_map_physeq_levels <- function(physeq, group_var, max_abundance_for_color = NULL, tax_order = NULL,
                                        tax_names = NULL, color_sample_names = TRUE, gradient_steps = c(0.15, 0.3, 0.45, 1)){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
        gradient_steps <- c(0, 1e-14, gradient_steps)
        
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        
        if (!is.factor(sample_data(physeq)[[group_var]])) {sample_data(physeq)[[group_var]] <- as.factor(sample_data(physeq)[[group_var]])}
        
        CT <- as(otu_table(physeq), "matrix") # taxa are columns
        DF_CT <- as.data.frame(t(CT))
        
        if (is.null(tax_order)){tax_order <- taxa_names(physeq)}
        
        if (length(tax_order) != nrow(DF_CT) || !all(tax_order %in% rownames(DF_CT))) {stop("the given tax_order must contain all taxa_names(physeq). If you want to subset,
                                                                                            do it before with prune_taxa")}
        DF_CT <- DF_CT[tax_order, ]
        
        if (is.null(tax_names)) {tax_names <- paste("Taxon_", 1:nrow(DF_CT), sep = "")}
        
        if (length(tax_names) != nrow(DF_CT) ) {stop("the given tax_names must be a character vector of length = ntaxa(physeq)")}
        
        tax_names[is.na(tax_names)] <- "NA"
        rownames(DF_CT) <- make.unique(tax_names)
        DF_CT$Taxa <- rownames(DF_CT)
        DF_CT <- tidyr::gather(DF_CT, key = Sample , value = Count, -Taxa)
        DF_CT$Taxa <- factor(DF_CT$Taxa, levels = rev(unique(DF_CT$Taxa)), ordered = TRUE)
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        DF_CT$Sample <- factor(DF_CT$Sample, levels = LookUpDF$Sample, ordered = TRUE)
        DF_CT$Group <- LookUpDF$Group[match(DF_CT$Sample, LookUpDF$Sample)]
        
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # -- get the level combis --
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        i <- ivec[x]
                })
        })
        j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        j <- jvec[x]
                })
        })
        i_s <- i_s[upper.tri(i_s)]
        j_s <- j_s[upper.tri(j_s)]
        # ----
        
        plot_list <- vector("list", length = length(i_s))
        
        
        
        for (k in seq_along(i_s)){
                i <- i_s[k]
                j <- j_s[k]
                
                group_fac_current <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])
                DF_CT_current <- filter(DF_CT, Group %in% group_fac_current)
                DF_CT_current$Group <- factor(DF_CT_current$Group, levels = c(fac_levels[i], fac_levels[j]), ordered = TRUE)
                LookUpDF_current <- LookUpDF[LookUpDF$Sample %in% DF_CT_current$Sample, ]
                DF_CT_current$Sample <- factor(DF_CT_current$Sample, levels = LookUpDF_current$Sample, ordered = TRUE)
                
                if (is.null(max_abundance_for_color)) {
                        max_abundance_for_color_current <- quantile(DF_CT_current$Count, .9)
                } else {
                        max_abundance_for_color_current <- max_abundance_for_color
                }
                
                # Color the sample names based on levels in group fac
                if (length(levels(LookUpDF$Group)) <= 10 && color_sample_names){
                        color_lookup <- data.frame(level = levels(LookUpDF$Group), color = cbPalette[2:(length(levels(LookUpDF$Group)) + 1)])
                        #LookUpDF$Group[match(levels(DF_CT_current$Sample), LookUpDF$Sample)]
                        colxaxis <- as.character(color_lookup$color[match(LookUpDF$Group[match(levels(DF_CT_current$Sample), LookUpDF$Sample)], color_lookup$level)])
                } else {
                        colxaxis <- rep("black", nrow(LookUpDF))
                }
                
                
                hmTr <- ggplot(DF_CT_current, aes(x = Sample, y = Taxa, fill = Count))
                hmTr <- hmTr + 
                        geom_raster() + 
                        scale_fill_gradientn("", limits = c(0, max_abundance_for_color_current), colors = c("red", viridis(5)), values = gradient_steps, oob = squish) +
                        scale_x_discrete(position = "top") +
                        #coord_equal() +
                        labs(x=NULL, y=NULL) +
                        theme_tufte(base_family = "Helvetica") +
                        theme(axis.ticks=element_blank(),
                              axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0,
                                                         colour = colxaxis), axis.text.y =   element_text(size = 13))
                
                plot_list[[k]] <- hmTr
                names(plot_list)[k] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        plot_list
}




#######################################
### plot_top_abundances_boxAndviolin##
#################
# generates box and violin plots of the taxons in physeq in the order given by tax_order and named by tax_names
# if tax_order and tax_names are NULL the order and names of physeq will be used
# facet_cols determines ncol in facet_wrap

# OUTPUT:
# generates for each level combination of the grouping factor (defined by group_var) seven plots, so for each combi a list of 7 plots

plot_top_abundances_boxAndviolin <- function(physeq, group_var, tax_order = NULL, tax_names = NULL, facet_cols = 5){
        
  cbPalette <- c("#037896","#a59c15", "#1fccef", "#77a515","#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        # I prefer taxa_are_rows = FALSE so rows (= observations = samples), and colums = variables = taxa
        
        if (!is.factor(sample_data(physeq)[[group_var]])) {sample_data(physeq)[[group_var]] <- as.factor(sample_data(physeq)[[group_var]])}
        
        CT <- as(otu_table(physeq), "matrix") # taxa are columns
        DF_CT <- as.data.frame(t(CT))
        
        if (is.null(tax_order)){tax_order <- taxa_names(physeq)}
        
        if (length(tax_order) != nrow(DF_CT) || !all(tax_order %in% rownames(DF_CT))) {stop("the given tax_order must contain all taxa_names(physeq). If you want to subset,
                                                                                            do it before with prune_taxa")}
        DF_CT <- DF_CT[tax_order, ]
        
        if (is.null(tax_names)) {tax_names <- paste("Taxon_", 1:nrow(DF_CT), sep = "")}
        
        if (length(tax_names) != nrow(DF_CT) ) {stop("the given tax_names must be a character vector of length = ntaxa(physeq)")}
        
        rownames(DF_CT) <- make.unique(tax_names)
        DF_CT$Taxa <- rownames(DF_CT)
        
        DF_CT <- tidyr::gather(DF_CT, key = Sample , value = Count, -Taxa)
        DF_CT$Taxa <- factor(DF_CT$Taxa, levels = unique(DF_CT$Taxa), ordered = TRUE)
        
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        DF_CT$Group <- LookUpDF$Group[match(DF_CT$Sample, LookUpDF$Sample)]
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # -- get the level combis --
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        i <- ivec[x]
                })
        })
        j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        j <- jvec[x]
                })
        })
        i_s <- i_s[upper.tri(i_s)]
        j_s <- j_s[upper.tri(j_s)]
        # ----
        
        if (length(levels(LookUpDF$Group)) <= 7){
                color_lookup <- data.frame(level = levels(LookUpDF$Group), color = cbPalette[2:(length(levels(LookUpDF$Group)) + 1)])
        } else {
                color_lookup <- data.frame(level = levels(LookUpDF$Group), color = viridis(length(levels(LookUpDF$Group))))
        }
        for (e in seq_along(color_lookup)) {
                color_lookup[, e] <- as.character(color_lookup[,e])
        }
        colors_to_use <- color_lookup$color
        names(colors_to_use) <- color_lookup$level
        
        
        plot_list <- vector("list", length = length(i_s))
        
        for (k in seq_along(i_s)){
                i <- i_s[k]
                j <- j_s[k]
                
                group_fac_current <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])
                DF_CT_current <- filter(DF_CT, Group %in% group_fac_current)
                DF_CT_current$Group <- factor(DF_CT_current$Group, levels = c(fac_levels[i], fac_levels[j]), ordered = TRUE)
                
                Tr <- ggplot(DF_CT_current, aes(x = Taxa, y = Count, col = Group))
                Tr <- Tr +
                        geom_boxplot(outlier.color = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        scale_color_manual("", values = colors_to_use) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0))
                
                Tr1 <- ggplot(DF_CT_current, aes(x = Group, y = Count, col = Group))
                Tr1 <- Tr1 +
                        geom_boxplot(outlier.color = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        facet_wrap(~ Taxa, ncol = facet_cols, scales = "free_y") +
                        scale_color_manual("", values = colors_to_use) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw()
                
                Tr2 <- ggplot(DF_CT_current, aes(x = Taxa, y = Count, col = Group))
                Tr2 <- Tr2 +
                        geom_violin(fill = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        scale_color_manual("", values = colors_to_use) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0))
                
                Tr3 <- ggplot(DF_CT_current, aes(x = Group, y = Count, col = Group))
                Tr3 <- Tr3 +
                        geom_violin(fill = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        facet_wrap(~ Taxa, ncol = facet_cols, scales = "free_y") +
                        scale_color_manual("", values = colors_to_use) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        theme(legend.position = "none")
                
                # same for log10
                # if (min(DF_CT_current$Count[DF_CT_current$Count > 0]) > 1e-6) {
                #         DF_CT_current$Count[DF_CT_current$Count == 0] <- 1e-6
                # } else {
                DF_CT_current$Count[DF_CT_current$Count == 0] <- min(DF_CT_current$Count[DF_CT_current$Count > 0])
                # }
                
                Tr4 <- ggplot(DF_CT_current, aes(x = Taxa, y = Count, col = Group))
                Tr4 <- Tr4 +
                        geom_boxplot(outlier.color = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        scale_color_manual("", values = colors_to_use) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0)) +
                        scale_y_log10()
                
                Tr5 <- ggplot(DF_CT_current, aes(x = Group, y = Count, col = Group))
                Tr5 <- Tr5 +
                        geom_boxplot(outlier.color = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        facet_wrap(~ Taxa, ncol = facet_cols, scales = "free_y") +
                        scale_color_manual("", values = colors_to_use) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        scale_y_log10() +
                        theme(legend.position = "none")
                
                Tr6 <- ggplot(DF_CT_current, aes(x = Taxa, y = Count, col = Group))
                Tr6 <- Tr6 +
                        geom_violin(fill = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        scale_color_manual("", values = colors_to_use) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0)) +
                        scale_y_log10()
                
                Tr7 <- ggplot(DF_CT_current, aes(x = Group, y = Count, col = Group))
                Tr7 <- Tr7 +
                        geom_violin(fill = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        facet_wrap(~ Taxa, ncol = facet_cols, scales = "free_y") +
                        scale_color_manual("", values = colors_to_use) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        scale_y_log10() +
                        theme(legend.position = "none")
                
                
                plot_list[[k]] <- list(Tr, Tr1, Tr2, Tr3, Tr4, Tr5, Tr6, Tr7)
                names(plot_list)[k] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        plot_list
}


#######################################
### FUNCTION: test_differential_prevalence
#######################################
# Function performs fisher exact test on prevalence (absence presence) between the levels
# in a grouping factor
# INPUT:
# physeq: phyloseq
# group_var: name of the column in sample_data(physeq) that defines the groups
# p.adj.method, used in p.adjust
# minCount: present are taxa in species with more counts than minCount
# OUTPUT:
# list of data.frames, one data frame for each combi of levels in your grouping factor
# The data frames are ordered by p_value, and the tax_table has been cbound:)

test_differential_prevalence <- function(physeq, group_var, p.adj.method = "fdr", minCount = 0L) {
        
        
        ### The Fisher test of pooled successes
        # MatrixFisher <- matrix(c(200, 80, 200, 320), nrow = 2)
        # rownames(MatrixFisher) <- c("Strain A", "Strain B")
        # colnames(MatrixFisher) <- c("Success", "Fail")
        # fisher.test(MatrixFisher, alternative = "two")
        # PValueFisher <- fisher.test(MatrixFisher, alternative = "two")$p.value
        
        if (taxa_are_rows(physeq)) { physeq <- t(physeq)}
        
        CT <- as(otu_table(physeq), "matrix")
        CT <- CT > minCount
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        
        prev_list <- lapply(fac_levels, function(level){
                matrix(c(Present = colSums(CT[group_fac == level, ]),
                           Absent = colSums(!CT[group_fac == level, ])), ncol = 2)
        })
        
        Count_list <- lapply(fac_levels, function(level){
          ps1<-prune_samples(sample_data(physeq)[[group_var]]==level, physeq)
          total_counts_of_ASV <-as.matrix(taxa_sums(ps1))
          total_counts_of_ASV
        })
        
       
        
        # -- get the level combis --
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        i <- ivec[x]
                })
        })
        j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        j <- jvec[x]
                })
        })
        i_s <- i_s[upper.tri(i_s)]
        j_s <- j_s[upper.tri(j_s)]
        # ----
        
        p_val_list <- vector("list", length = length(i_s))
        for (k in seq_along(i_s)){
                i <- i_s[k]
                j <- j_s[k]
                prev_PC_1 <- round(100*(prev_list[[i]][, 1]/(prev_list[[i]][1,1] + prev_list[[i]][1,2])), 1)
                prev_PC_2 <- round(100*(prev_list[[j]][, 1]/(prev_list[[j]][1,1] + prev_list[[j]][1,2])), 1)
                total_prev_PC_1<-prev_list[[i]][, 1]
                total_prev_PC_2<-prev_list[[j]][, 1] 
                count_PC_1<-Count_list[[i]][, 1]
                count_PC_2<-Count_list[[j]][, 1]
                direction <- rep("down", length(prev_PC_1))
                direction[prev_PC_1 > prev_PC_2] <- "up"
                rowwise_compare_matrix <- cbind(prev_list[[i]], prev_list[[j]])
                p_vals <- sapply(1:nrow(rowwise_compare_matrix), function(e){
                        mat_fisher <- matrix(c(rowwise_compare_matrix[e, 1],
                                        rowwise_compare_matrix[e, 3],
                                        rowwise_compare_matrix[e, 2],
                                        rowwise_compare_matrix[e, 4]), ncol = 2)
                        fisher.test(mat_fisher)$p.value
                })
                p_vals_adj <- p.adjust(p_vals, p.adj.method)
                significance <- rep("", length(prev_PC_1))
                significance[p_vals_adj <= 0.05] <- "*"
                significance[p_vals_adj <= 0.01] <- "**"
                significance[p_vals_adj <= 0.001] <- "***"
                p_val_df <- data.frame(p_vals, p_vals_adj, significance, prev_PC_1, prev_PC_2, direction, total_prev_PC_1, total_prev_PC_2, count_PC_1, count_PC_2)
                colnames(p_val_df)[1:2] <- paste(c("p_val_", "p_val_adj_"), fac_levels[i], "_vs_", fac_levels[j], sep = "")
                colnames(p_val_df)[4:5] <- paste(c("prev_PC_"), c(fac_levels[i], fac_levels[j]), sep = "")
                colnames(p_val_df)[7:8] <- paste(c("Prev_"), c(fac_levels[i], fac_levels[j]), sep = "")
                colnames(p_val_df)[9:10] <- paste(c("Counts_"), c(fac_levels[i], fac_levels[j]), sep = "")
                tax.table<-as.data.frame(unclass(tax_table(physeq)))
                 p_val_df <- cbind(as.data.frame(p_val_df),tax.table )
                p_val_df$Taxon <- colnames(CT)
                p_val_df <- arrange(p_val_df, p_val_df[,1])
                p_val_df <- select(p_val_df, Taxon, 1:(ncol(p_val_df)-1))
                p_val_list[[k]] <- p_val_df
                names(p_val_list)[[k]] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        p_val_list
        
}









#######################################
### wilcoxTestApply_physeq
#######################################
# Accepts several groups in group_var and results in a list of DF for the pairwise comparisons
# directly related to wilcoxTestApply from SimulationFunctions.R just not looping through a list of phyloseq objects
# but doing the job on a single phyloseq object
# does wilcoxon test and adds further info
# in addition it performs a fisher exact test on the sparsity proportions
# NB: Tested: 
# with excludeZeros you can decide on whether 0 counts should be excluded for the wilcox.test, the fisher sparsity test is
# of course not affected. 
# NB: in case in one of the two groups all counts are 0 and excludeZeros = T, then NA is given for all wilcoxon 
# statistics!
# The teststatistic is based on the standardized teststatistic, equation provided by multtest::mt.minP (compare with mtApply)
# (see equation for standStat2 in the code, results in exactly the same as standStat when uncomment)
## Input
# physeq object
# group_var: refers to a factor in sample_data(phyloseq) that defines the groups
# excludeZeros: decides on whether 0s should be considered when comparing the groups in a wilcox.test
# p.adjust.method, used as method in p.adjust
## Output
# list of dataframes with the results for each pairwise group comparison in group_var associated group_fac

wilcoxTestApply_physeq <- function(physeq, group_var, excludeZeros = FALSE, p.adjust.method = "fdr") {
        
        
        if (taxa_are_rows(physeq)) { physeq <- t(physeq)}
        
        CT <- as(otu_table(physeq), "matrix")
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # -- get the level combis --
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        i <- ivec[x]
                })
        })
        j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        j <- jvec[x]
                })
        })
        i_s <- i_s[upper.tri(i_s)]
        j_s <- j_s[upper.tri(j_s)]
        # ----
        
        result_list <- vector("list", length = length(i_s))
        
        for (k in seq_along(i_s)){
                i <- i_s[k]
                j <- j_s[k]
                
                res_mat <- apply(CT, 2, function(taxon_counts){
                        x <- taxon_counts[as.numeric(group_fac) == i]
                        Zeros_grp1 <- sum(x == 0)
                        Sparsity_grp1 <- 100*(Zeros_grp1/length(x))
                        Present_grp1 <- length(x)-Zeros_grp1
                        if(excludeZeros){
                                x <- x[x != 0]
                        }
                        Median_grp1 <- median(x, na.rm = T) # NA in case all 0
                        Mean_grp1 <- mean(x, na.rm = T) # NaN in case all 0
                        if (is.na(Mean_grp1)){ Mean_grp1 = NA }
                        y <- taxon_counts[as.numeric(group_fac) == j]
                        Zeros_grp2 <- sum(y == 0)
                        Present_grp2 <- length(y)-Zeros_grp2
                        Sparsity_grp2 <- 100*(Zeros_grp2/length(y))
                        if(excludeZeros){
                                #if(all(y == 0)){y[1] <- ceiling(mean(taxon_counts))+1}
                                y <- y[y != 0]
                        }
                        Median_grp2 <- median(y, na.rm = T)
                        Mean_grp2 <- mean(y, na.rm = T)
                        if (is.na(Mean_grp2)){ Mean_grp2 = NA }
                        
                        if (length(x) != 0 && length(y) != 0){
                                wilcTest <- wilcox.test(x = x, y = y, alternative = "two", paired = F, exact = F)
                                pValue <- wilcTest$p.value
                                W <- wilcTest$statistic
                                # calculate standardized rank sum Wilcoxon statistics as in multtest::mt.minP
                                Ranks <- rank(c(x, y))
                                n1 <- length(x)
                                n2 <- length(y)
                                # Wx <- sum(Ranks[1:n1])-(n1*(n1+1)/2) # would be the same as W
                                # how about the other W?
                                # Wy <- sum(Ranks[(n1+1):n2]) - (n2*(n2+1)/2)
                                standStat <- -1*((sum(Ranks[1:n1]) - n1*(n1+n2+1)/2)/sqrt(n1*n2*(n1+n2+1)/12))
                                
                                # # if you want to check that multtest::mt.minP would give the same statistic
                                # mati <- matrix(c(x,y), nrow = 1)
                                # grFac <- c(rep(fac_levels[i], n1), rep(fac_levels[j], n2))
                                # grFac <- factor(grFac, levels = c(fac_levels[i], fac_levels[j]))
                                # standStat2 <- multtest::mt.minP(mati, grFac, test = "wilcoxon")$teststat
                                # # identical(standStat, standStat2) # TRUE
                                # uncomment all with standStat2 to test all the way
                                
                        } else {
                                pValue = NA
                                W <- NA
                                standStat = NA
                                n1 <- length(x)
                                n2 <- length(y)
                                # standStat2 = NA
                        }
                        
                        
                        # -- add fisher exact test of presence differences (should be none in simulation) --
                        fisherMat <- matrix(c(Present_grp1, Zeros_grp1, Present_grp2, Zeros_grp2), 
                                            ncol = 2, dimnames = list(c("Present", "Absent"), c("grp1", "grp2")) )
                        Test <- fisher.test(fisherMat)
                        
                        c(teststat = standStat, p_val = pValue, Median_grp1 = Median_grp1, Median_grp2 = Median_grp2, 
                          Mean_grp1 = Mean_grp1, Mean_grp2 = Mean_grp2, n1 = n1, n2 = n2, Present_grp1 = Present_grp1, 
                          Present_grp2 = Present_grp2, Zeros_grp1 = Zeros_grp1, Zeros_grp2 = Zeros_grp2, 
                          Sparsity_grp1 = Sparsity_grp1, Sparsity_grp2 = Sparsity_grp2, W, 
                          p_val_Fisher = Test$p.value, Test$estimate) #, teststat2 = standStat2
                })
                
                res_mat <- t(res_mat)
                colnames(res_mat) <- c("teststat", "p_val", "Median_grp1", "Median_grp2", "Mean_grp1", "Mean_grp2", "n1",
                                       "n2", "Present_grp1", "Present_grp2", "Zeros_grp1", "Zeros_grp2", "Sparsity_grp1", "Sparsity_grp2", "W",
                                       "p_val_Fisher", "oddsRatioFisher") #, , "teststat2"
                wilc.adjusted <- p.adjust(res_mat[,"p_val"], method = p.adjust.method)
                fisher.adjusted <- p.adjust(res_mat[,"p_val_Fisher"], method = p.adjust.method)
                significance <- rep("", length(wilc.adjusted))
                significance[wilc.adjusted <= 0.05] <- "*"
                significance[wilc.adjusted <= 0.01] <- "**"
                significance[wilc.adjusted <= 0.001] <- "***"
                direction <- rep("down", length(wilc.adjusted))
                direction[res_mat[, "Median_grp1"] > res_mat[, "Median_grp2"]] <- "up"
                significance_fisher <- rep("", length(fisher.adjusted))
                significance_fisher[fisher.adjusted <= 0.05] <- "*"
                significance_fisher[fisher.adjusted <= 0.01] <- "**"
                significance_fisher[fisher.adjusted <= 0.001] <- "***"
                direction_fisher <- rep("down", length(fisher.adjusted))
                direction_fisher[res_mat[, "Sparsity_grp1"] < res_mat[, "Sparsity_grp2"]] <- "up"
                
                
                DF <- data.frame(Taxon = rownames(res_mat), res_mat, p_val_adj = wilc.adjusted, significance = significance,
                                 direction = direction, p_val_Fisher_adj = fisher.adjusted,
                                 significance_fisher = significance_fisher,
                                 direction_fisher = direction_fisher)
                
                
                DF <- dplyr::select(DF, 1:3, 19:21, 17, 22:24, 4:7, 10:15, 8:9, 16, 18)
                # teststat/standStat2 version:
                # DF <- dplyr::select(DF, 1:2, 19, 3, 20:22, 17, 23:25, 4:7, 10:15, 8:9, 16, 18)

                DF <- cbind(DF, tax_table(physeq))
                # DF <- dplyr::arrange(DF, desc(abs(teststat)))
                DF <- dplyr::arrange(DF, p_val)
                
                result_list[[k]] <- DF
                names(result_list)[[k]] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        result_list
}



#######################################
### DESeq2Apply_physeq
#######################################
## Inputs
# physeq: phyloseq object
# group_var: name of column that defines group fac in sample_data
# SFs: often you might want to give the SizeFactors already because you wanted to calculate them on non-filtered data,
# when SFs are not NULL, type is ignored
# type: type in estimateSizeFactors, ignored when Size factors given
## OUTPUT:
# list of two lists: the first: List of DESeq2 results plus of fisher.exact test, for each level combination in group factor one data_frame = list entry
# in the second list is just the size factor adjusted physeq object


DESeq2Apply_physeq <- function(physeq, group_var, adjust, SFs = NULL, type = "ratio", p.adjust.method = "fdr"){
        
        if (taxa_are_rows(physeq)) {physeq <- t(physeq)}
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        fac_levels <- levels(group_fac)
        
        # NB: DESeq2 can not deal with ordered factors, sees them somehow as just one level, therefore
        sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], levels = levels(sample_data(physeq)[[group_var]]), ordered = FALSE)
        #sample_data(physeq)[[adjust]] <- factor(sample_data(physeq)[[adjust]], levels = levels(sample_data(physeq)[[adjust]]), ordered = FALSE)
        
        #DES = phyloseq_to_deseq2(physeq, formula(paste("~", adjust +group_var)))
        DES = phyloseq_to_deseq2(physeq, formula(paste("~", "+", group_var)))
        
        if (is.null(SFs)){
                if(type == "ratio"){
                        GM <- apply(otu_table(physeq), 2, gm_own, zeros.count = FALSE)
                }
                
                dds <- estimateSizeFactors(DES, type = type, geoMeans = GM)
                # NB: geoMeans is ignored when type = "iterate"
                # NB2: "iterate" takes much longer than "ratio", and when using GM via gm_own with "ratio" the size factors 
                # correlate with above 99% with size factors from "iterate"
                # SFs2 <- sizeFactors(dds)
                
        } else {
                dds <- DES
                sizeFactors(dds) = SFs
                # identical(sizeFactors(dds), SFs) 
        }
        
        
        dds <- estimateDispersions(dds, quiet = TRUE) 
        dds <- nbinomWaldTest(dds)
        
        # to get the size factor adjusted physeq object
        physeq_out <- physeq
        otu_table(physeq_out) <- otu_table(t(counts(dds, normalized = TRUE)), taxa_are_rows = FALSE)
        
        # -- get the level combis, NB: dds contains result infos on all combinations! --
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        i <- ivec[x]
                })
        })
        j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        j <- jvec[x]
                })
        })
        i_s <- i_s[upper.tri(i_s)]
        j_s <- j_s[upper.tri(j_s)]
        # ----
        
        result_list <- vector("list", length = length(i_s))
        
        for (k in seq_along(i_s)) {
                i <- i_s[k]
                j <- j_s[k]
                
                res <- as.data.frame(results(dds, contrast = c(group_var, fac_levels[i], fac_levels[j])))
                
                res$p_val_adj <- p.adjust(res$pvalue, method = p.adjust.method) # NB: in case of "fdr" same as default DESeq2
                
                CT <- counts(dds, normalized = TRUE)
                n1 <- sum(group_fac == fac_levels[i])
                n2 <- sum(group_fac == fac_levels[j])
                res$Median_grp1 <- apply(CT[, group_fac == fac_levels[i]], 1, median)
                res$Median_grp2 <- apply(CT[, group_fac == fac_levels[j]], 1, median)
                res$Mean_grp1 <- apply(CT[, group_fac == fac_levels[i]], 1, mean)
                res$Mean_grp2 <- apply(CT[, group_fac == fac_levels[j]], 1, mean)
                # res$baseMeanSelf <- apply(CT, 1, mean) # exactly the same as baseMean!
                res$Zeros_grp1 <- apply(CT[, group_fac == fac_levels[i]], 1, function(cnts){sum(cnts == 0)})
                res$Zeros_grp2 <- apply(CT[, group_fac == fac_levels[j]], 1, function(cnts){sum(cnts == 0)})
                res$Present_grp1 <- apply(CT[, group_fac == fac_levels[i]], 1, function(cnts){sum(cnts != 0)})
                res$Present_grp2 <- apply(CT[, group_fac == fac_levels[j]], 1, function(cnts){sum(cnts != 0)})
                res$prev_grp1 <- round(100*(res$Present_grp1/n1),1)
                res$prev_grp2 <- round(100*(res$Present_grp2/n2), 1)
                res$Sparsity_grp1 <- 100*(res$Zeros_grp1/n1)
                res$Sparsity_grp2 <- 100*(res$Zeros_grp2/n2)
                res$n1 <- n1
                res$n2 <- n2
                
                # -- add fisher exact test of sparsity/prevalence again --
                Fisher <- t(sapply(1:nrow(res), FUN = function(i){
                        fisherMat <- matrix(c(res$Present_grp1[i], res$Zeros_grp1[i], res$Present_grp2[i],
                                              res$Zeros_grp2[i]), ncol = 2, dimnames = list(c("Present", "Absent"),  c("grp1", "grp2")))
                        Test <- fisher.test(fisherMat)
                        cbind(Test$p.value, Test$estimate)
                }))
                
                p_val_adj_Fisher <- p.adjust(Fisher[,1], method = p.adjust.method)
                
                res$p_val_Fisher <- Fisher[,1]
                res$p_val_Fisher_adj <- p_val_adj_Fisher
                res$oddsRatioFisher <- Fisher[,2]
                
                significance <- rep("", nrow(res))
                significance[res$padj <= 0.05] <- "*"
                significance[res$padj <= 0.01] <- "**"
                
                significance[res$padj <= 0.001] <- "***"
                
                direction <- rep( fac_levels[j], nrow(res))
                direction[res$Mean_grp1 > res$Mean_grp2] <- fac_levels[i]
                significance_fisher <- rep("", nrow(res))
                significance_fisher[p_val_adj_Fisher <= 0.05] <- "*"
                significance_fisher[p_val_adj_Fisher <= 0.01] <- "**"
                significance_fisher[p_val_adj_Fisher <= 0.001] <- "***"
                direction_fisher <- rep("down", nrow(res))
                direction_fisher[res$Sparsity_grp1 < res$Sparsity_grp2] <- "up"
                
                res <- dplyr::mutate(res, Taxon = rownames(res), significance = significance, direction = direction,
                                     significance_fisher = significance_fisher, direction_fisher = direction_fisher)
                
                
                res <- dplyr::select(res, Taxon, teststat = stat, p_val = pvalue, p_val_adj,
                                     significance, direction, p_val_Fisher,
                                     p_val_Fisher_adj, significance_fisher,
                                     direction_fisher, Median_grp1, Median_grp2, Mean_grp1,
                                     Mean_grp2, Present_grp1, Present_grp2, Zeros_grp1, Zeros_grp2, prev_grp1, prev_grp2, 
                                     n1, n2, baseMean, log2FoldChange, lfcSE, oddsRatioFisher)
                # NB: I dropped here padj from DESeq since same as p_val_adj in case of p.adjust.method = "fdr"
                res <- cbind(res, unclass(tax_table(physeq)))
                res <- dplyr::arrange(res, (p_val_adj))
                # res <- res[order(res$p_val),]
                result_list[[k]] <- res
                names(result_list)[[k]] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        out <- list(result_list, physeq_out)
        
}
####################
#Adjust_multiple 
#####################
DESeq2Apply_physeq_adj_2 <- function(physeq, SFs = NULL, type = "ratio", p.adjust.method = "fdr",  
                                     adj1= adj1, adj2= adj2, group_var=group_var){
  
  if (taxa_are_rows(physeq)) {physeq <- t(physeq)}
  
  group_fac <- factor(sample_data(physeq)[[adj1]])
  fac_levels <- levels(group_fac)
  
  condition_fac <- factor(sample_data(physeq)[[group_var]])
  condition_fac_levels <- levels(condition_fac)
  
  # NB: DESeq2 can not deal with ordered factors, sees them somehow as just one level, therefore
  sample_data(physeq)[[adj1]] <- factor(sample_data(physeq)[[adj1]], levels = levels(sample_data(physeq)[[adj1]]), ordered = FALSE)
  sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], levels = levels(sample_data(physeq)[[group_var]]), ordered = FALSE)
  sample_data(physeq)[[adj2]] <- factor(sample_data(physeq)[[adj2]], levels = levels(sample_data(physeq)[[adj2]]), ordered = FALSE)
  
   # sample_data(physeq)[[condition2]] <- factor(sample_data(physeq)[[condition2]], levels = levels(sample_data(physeq)[[condition2]]), ordered = FALSE)
  
  
  #DES = phyloseq_to_deseq2(physeq, formula(paste("~", adjust +group_var)))
  DES = phyloseq_to_deseq2(physeq, formula(paste("~", adj1, "+",adj2, "+", group_var, "+",adj1, ":", group_var )))
  
  if (is.null(SFs)){
    if(type == "ratio"){
      GM <- apply(otu_table(physeq), 2, gm_own, zeros.count = FALSE)
    }
    
    dds <- estimateSizeFactors(DES, type = type, geoMeans = GM)
    # NB: geoMeans is ignored when type = "iterate"
    # NB2: "iterate" takes much longer than "ratio", and when using GM via gm_own with "ratio" the size factors 
    # correlate with above 99% with size factors from "iterate"
    # SFs2 <- sizeFactors(dds)
    
  } else {
    dds <- DES
    sizeFactors(dds) = SFs
    # identical(sizeFactors(dds), SFs) 
  }
  

  dds <- estimateDispersions(dds, quiet = TRUE) 
  dds <- nbinomWaldTest(dds)
  
  #plotSparsity(dds)
  
  
  
  # to get the size factor adjusted physeq object
  physeq_out <- physeq
  otu_table(physeq_out) <- otu_table(t(counts(dds, normalized = TRUE)), taxa_are_rows = FALSE)
  
  # -- get the level combis, NB: dds contains result infos on all combinations! --
  
  # ----
 
  results_names<-(resultsNames(dds))
  result_list <- vector("list", length = length(results_names))
  res1 <- as.data.frame(results(dds, contrast = c(condition1, condition_fac_levels[1], condition_fac_levels[2]))) # the condition effect for group I
  res1$p_val_adj <- p.adjust(res1$pvalue, method = p.adjust.method)
  res1$significance <- rep("", nrow(res1))
  res1$significance[res1$padj <= 0.05] <- "*"
  res1$significance[res1$padj <= 0.01] <- "**"
  res1$significance[res1$padj <= 0.001] <- "***"
  
  res2<- as.data.frame(results(dds, list(c(results_names[3], results_names[4])))) # the condition effect for group II
  res2$p_val_adj <- p.adjust(res2$pvalue, method = p.adjust.method)
  res2$significance <- rep("", nrow(res2))
  res2$significance[res2$padj <= 0.05] <- "*"
  res2$significance[res2$padj <= 0.01] <- "**"
  res2$significance[res2$padj <= 0.001] <- "***"
  
  res3<-as.data.frame(results(dds, name=results_names[4]))# is the condition effect different across groups?
  res3$p_val_adj <- p.adjust(res3$pvalue, method = p.adjust.method)
  res3$significance <- rep("", nrow(res3))
  res3$significance[res3$padj <= 0.05] <- "*"
  res3$significance[res3$padj <= 0.01] <- "**"
  res3$significance[res3$padj <= 0.001] <- "***"

  res4<- as.data.frame(results(dds, contrast = c(adj1, fac_levels[1], fac_levels[2]))) #effect of group
  res4$p_val_adj <- p.adjust(res4$pvalue, method = p.adjust.method)
  res4$significance <- rep("", nrow(res4))
  res4$significance[res4$padj <= 0.05] <- "*"
  res4$significance[res4$padj <= 0.01] <- "**"
  res4$significance[res4$padj <= 0.001] <- "***"
 
     # NB: in case of "fdr" same as default DESeq2
    res<-as.data.frame(rownames(res1))  
    CT <- counts(dds, normalized = TRUE)
    n1 <- sum(condition_fac == condition_fac_levels[1])
    n2 <- sum(condition_fac == condition_fac_levels[2])
    res$Median_grp1 <- apply(CT[, condition_fac == condition_fac_levels[1]], 1, median)
    res$Median_grp2 <- apply(CT[, condition_fac == condition_fac_levels[2]], 1, median)
    res$Mean_grp1 <- apply(CT[, condition_fac == condition_fac_levels[1]], 1, mean)
    res$Mean_grp2 <- apply(CT[, condition_fac == condition_fac_levels[2]], 1, mean)
    direction <- rep( condition_fac_levels[2], nrow(res))
    direction[res$Mean_grp1 > res$Mean_grp2] <- condition_fac_levels[1]
    res$Zeros_grp1 <- apply(CT[, condition_fac == condition_fac_levels[1]], 1, function(cnts){sum(cnts == 0)})
    res$Zeros_grp2 <- apply(CT[, condition_fac == condition_fac_levels[2]], 1, function(cnts){sum(cnts == 0)})
    res$Present_grp1 <- apply(CT[, condition_fac == condition_fac_levels[1]], 1, function(cnts){sum(cnts != 0)})
    res$Present_grp2 <- apply(CT[, condition_fac == condition_fac_levels[2]], 1, function(cnts){sum(cnts != 0)})
    res$prev_grp1 <- round(100*(res$Present_grp1/n1),1)
    res$prev_grp2 <- round(100*(res$Present_grp2/n2), 1)
    res$Sparsity_grp1 <- 100*(res$Zeros_grp1/n1)
    res$Sparsity_grp2 <- 100*(res$Zeros_grp2/n2)
    res$n1 <- n1
    res$n2 <- n2
    
    
    rownames(res)<-res$`rownames(res1)`
    res<-res[-1]
    res <- dplyr::mutate(res, Taxon = rownames(res), direction = direction)
    res <- dplyr::select(res, Taxon,  direction, Present_grp1, Present_grp2,  prev_grp1, prev_grp2 )
    # NB: 1 dropped here padj from DESeq since same as p_val_adj in case of p.adjust.method = "fdr"
    res <- cbind(res, tax_table(physeq))
    res1<-cbind(res1, res)
    res2<-cbind(res2, res)
    res3<-cbind(res3, res)
    
   ####Group effect
    res_g<-as.data.frame(rownames(res4)) 
    CT <- counts(dds, normalized = TRUE)
    n1 <- sum(group_fac == fac_levels[1])
    n2 <- sum(group_fac == fac_levels[2])
    res_g$Median_grp1 <- apply(CT[, group_fac == fac_levels[1]], 1, median)
    res_g$Median_grp2 <- apply(CT[, group_fac == fac_levels[2]], 1, median)
    res_g$Mean_grp1 <- apply(CT[, group_fac == fac_levels[1]], 1, mean)
    res_g$Mean_grp2 <- apply(CT[, group_fac == fac_levels[2]], 1, mean)
    direction_g <- rep(fac_levels[2], nrow(res))
    direction_g[res_g$Mean_grp1 > res_g$Mean_grp2] <- fac_levels[1]
    res_g$Zeros_grp1 <- apply(CT[, group_fac == fac_levels[1]], 1, function(cnts){sum(cnts == 0)})
    res_g$Zeros_grp2 <- apply(CT[, group_fac == fac_levels[2]], 1, function(cnts){sum(cnts == 0)})
    res_g$Present_grp1 <- apply(CT[, group_fac == fac_levels[1]], 1, function(cnts){sum(cnts != 0)})
    res_g$Present_grp2 <- apply(CT[, group_fac == fac_levels[2]], 1, function(cnts){sum(cnts != 0)})
    res_g$prev_grp1 <- round(100*(res$Present_grp1/n1),1)
    res_g$prev_grp2 <- round(100*(res$Present_grp2/n2), 1)
    
    res_g$n1 <- n1
    res_g$n2 <- n2
    rownames(res_g)<-res_g$`rownames(res4)`
    res_g<-res_g[-1]
    res_g <- dplyr::mutate(res_g, Taxon = rownames(res_g), direction = direction_g)
    res_g <- dplyr::select(res_g, Taxon,  direction, Present_grp1, Present_grp2,  prev_grp1, prev_grp2 )
    res_g <- cbind(res_g, tax_table(physeq))
    res4<-cbind(res4, res_g)
    
    #####
    
    # res <- res[order(res$p_val),]
    result_list[[1]] <- res1
    result_list[[2]] <- res2
    result_list[[3]] <- res3
    result_list[[4]] <- res4
    names(result_list) <-c(paste(condition_fac_levels[1], "_vs_", condition_fac_levels[2], "_", fac_levels[1], sep = ""),
                           paste(condition_fac_levels[1], "_vs_", condition_fac_levels[2], "_", fac_levels[2], sep = ""),
                           paste(condition_fac_levels[1], "_vs_", condition_fac_levels[2], "_", "Difference", sep = ""), 
                           paste(fac_levels[1], "_vs_", fac_levels[2], sep = "") 
                           ) 
 
  
  out <- list(result_list, physeq_out, dds)
  
}



##############
#Adjust
#############

DESeq2Apply_physeq_adj <- function(physeq, group_var,  SFs = NULL, type = "ratio", p.adjust.method = "fdr", adj1=adj1, adj2=adj2){
  
  if (taxa_are_rows(physeq)) {physeq <- t(physeq)}
  
  group_fac <- factor(sample_data(physeq)[[group_var]])
  fac_levels <- levels(group_fac)
  # NB: DESeq2 can not deal with ordered factors, sees them somehow as just one level, therefore
  sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], levels = levels(sample_data(physeq)[[group_var]]), ordered = FALSE)
  
  if(is.null(adj2)){
  sample_data(physeq)[[adj1]] <- factor(sample_data(physeq)[[adj1]], levels = levels(sample_data(physeq)[[adj1]]), ordered = FALSE)
  DES = phyloseq_to_deseq2(physeq, formula(paste("~", adj1, "+",group_var)))
  }else{
  sample_data(physeq)[[adj1]] <- factor(sample_data(physeq)[[adj1]], levels = levels(sample_data(physeq)[[adj1]]), ordered = FALSE)
  sample_data(physeq)[[adj2]] <- factor(sample_data(physeq)[[adj2]], levels = levels(sample_data(physeq)[[adj2]]), ordered = FALSE)
  DES = phyloseq_to_deseq2(physeq, formula(paste("~", adj1, "+", adj2, "+",group_var)))
  }
  
  
  
  
  #sample_data(physeq)[[adjust]] <- factor(sample_data(physeq)[[adjust]], levels = levels(sample_data(physeq)[[adjust]]), ordered = FALSE)
  
  #DES = phyloseq_to_deseq2(physeq, formula(paste("~", adjust +group_var)))
  
  #DES = phyloseq_to_deseq2(physeq, formula(paste("~", adj1, "+",group_var)))
  
  if (is.null(SFs)){
    if(type == "ratio"){
      GM <- apply(otu_table(physeq), 2, gm_own, zeros.count = FALSE)
    }
    
    dds <- estimateSizeFactors(DES, type = type, geoMeans = GM)
    # NB: geoMeans is ignored when type = "iterate"
    # NB2: "iterate" takes much longer than "ratio", and when using GM via gm_own with "ratio" the size factors 
    # correlate with above 99% with size factors from "iterate"
    # SFs2 <- sizeFactors(dds)
    
  } else {
    dds <- DES
    sizeFactors(dds) = SFs
    # identical(sizeFactors(dds), SFs) 
  }
  
  
  dds <- estimateDispersions(dds, quiet = TRUE) 
  dds <- nbinomWaldTest(dds)
  
  resultsNames(dds)
  
  # to get the size factor adjusted physeq object
  physeq_out <- physeq
  otu_table(physeq_out) <- otu_table(t(counts(dds, normalized = TRUE)), taxa_are_rows = FALSE)
  
  # -- get the level combis, NB: dds contains result infos on all combinations! --
  fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
  i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
    sapply(seq_along(ivec), function(x){
      i <- ivec[x]
    })
  })
  j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
    sapply(seq_along(ivec), function(x){
      j <- jvec[x]
    })
  })
  i_s <- i_s[upper.tri(i_s)]
  j_s <- j_s[upper.tri(j_s)]
  # ----
  
  result_list <- vector("list", length = length(i_s))
  
  for (k in seq_along(i_s)) {
    i <- i_s[k]
    j <- j_s[k]
    
    res <- as.data.frame(results(dds, contrast = c(group_var, fac_levels[i], fac_levels[j])))
    
    res$p_val_adj <- p.adjust(res$pvalue, method = p.adjust.method) # NB: in case of "fdr" same as default DESeq2
    
    CT <- counts(dds, normalized = TRUE)
    n1 <- sum(group_fac == fac_levels[i])
    n2 <- sum(group_fac == fac_levels[j])
    res$Median_grp1 <- apply(CT[, group_fac == fac_levels[i]], 1, median)
    res$Median_grp2 <- apply(CT[, group_fac == fac_levels[j]], 1, median)
    res$Mean_grp1 <- apply(CT[, group_fac == fac_levels[i]], 1, mean)
    res$Mean_grp2 <- apply(CT[, group_fac == fac_levels[j]], 1, mean)
    
    res$Zeros_grp1 <- apply(CT[, group_fac == fac_levels[i]], 1, function(cnts){sum(cnts == 0)})
    res$Zeros_grp2 <- apply(CT[, group_fac == fac_levels[j]], 1, function(cnts){sum(cnts == 0)})
    res$Present_grp1 <- apply(CT[, group_fac == fac_levels[i]], 1, function(cnts){sum(cnts != 0)})
    res$Present_grp2 <- apply(CT[, group_fac == fac_levels[j]], 1, function(cnts){sum(cnts != 0)})
    res$prev_grp1 <- round(100*(res$Present_grp1/n1),1)
    res$prev_grp2 <- round(100*(res$Present_grp2/n2), 1)
    
    res$n1 <- n1
    res$n2 <- n2
    
    # -- add fisher exact test of sparsity/prevalence again --
    Fisher <- t(sapply(1:nrow(res), FUN = function(i){
      fisherMat <- matrix(c(res$Present_grp1[i], res$Zeros_grp1[i], res$Present_grp2[i],
                            res$Zeros_grp2[i]), ncol = 2, dimnames = list(c("Present", "Absent"), c("grp1", "grp2")))
      Test <- fisher.test(fisherMat)
      cbind(Test$p.value, Test$estimate)
    }))
    
    p_val_adj_Fisher <- p.adjust(Fisher[,1], method = p.adjust.method)
    
    res$p_val_Fisher <- Fisher[,1]
    res$p_val_Fisher_adj <- p_val_adj_Fisher
    res$oddsRatioFisher <- Fisher[,2]
    
    significance <- rep("", nrow(res))
    significance[res$padj <= 0.05] <- "*"
    significance[res$padj <= 0.01] <- "**"
    significance[res$padj <= 0.001] <- "***"
    direction <- rep( fac_levels[j], nrow(res))
    direction[res$Mean_grp1 > res$Mean_grp2] <- fac_levels[i]
    significance_fisher <- rep("", nrow(res))
    significance_fisher[p_val_adj_Fisher <= 0.05] <- "*"
    significance_fisher[p_val_adj_Fisher <= 0.01] <- "**"
    significance_fisher[p_val_adj_Fisher <= 0.001] <- "***"
    direction_fisher <- rep("down", nrow(res))
    direction_fisher[res$Sparsity_grp1 < res$Sparsity_grp2] <- "up"
    
    res <- dplyr::mutate(res, Taxon = rownames(res), significance = significance, direction = direction,
                         significance_fisher = significance_fisher, direction_fisher = direction_fisher)
    
    res <- dplyr::select(res, Taxon, teststat = stat, p_val = pvalue, p_val_adj,
                         significance, direction, p_val_Fisher,
                         p_val_Fisher_adj, significance_fisher,
                         direction_fisher, Median_grp1, Median_grp2, Mean_grp1,
                         Mean_grp2, Present_grp1, Present_grp2, Zeros_grp1, Zeros_grp2, prev_grp1, prev_grp2, 
                         n1, n2, baseMean, log2FoldChange, lfcSE, oddsRatioFisher)
    # NB: I dropped here padj from DESeq since same as p_val_adj in case of p.adjust.method = "fdr"
    res <- cbind(res, unclass(tax_table(physeq)))
    res <- dplyr::arrange(res, (p_val_adj))
    # res <- res[order(res$p_val),]
    result_list[[k]] <- res
    names(result_list)[[k]] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
  }
  
  out <- list(result_list, physeq_out, dds)
  
}

####################################
## calculate_TbTmatrixes:
###################################

# The method is based on direct comparisons between Taxa (one to one ratio comparisons), thus avoiding compositionality effects.
# For each "host taxon" a matrix is created, where the host taxons counts are directly compared to the
# counts of each other taxon.
# Specifically, in steps, example for count table of 100 taxons and 20 samples
# Step 1: calculate taxa by taxa ratios (one matrix per host taxon) to get rid of compositionality. 
# The ratios show for example count taxon 1/ count taxon 2 over all samples, illustrating the non-compositionality affected
# ratios of taxon 1 to taxon 2 in the different samples
# Step 2: divide the ratios by geometric mean over all samples, to make ratios of taxon1/taxon2 comparable to taxon1/taxon3 ratios and so on
# Step 3: log the ratios, so the sum of ratios/values over all samples = 0. 
# (see also implementation steps below to understand the values in the final TbT matrixes even better)
# NB: Zero Counts are fully ingored here: Reasoning on this treatment of zero counts: 
# Values/Ratios where one of the two taxa is missing (count = 0) should not be used to determine
# differential abundance (instead sparsity should be tested separately)
# so if the host taxon is 0 in a sample, then the entire sample will basically be ignored for that taxon.
# Values where the host taxon or the other taxon is 0 should be set to NA, so that all 0 values come
# from cases log(x/y) where x = y but both x and y != 0
# NB: further explanations within the code

## Input: 
# - physeq = a phyloseq object
# - group_var, name of the group_fac in sample_data(physeq)
## Output: 
# - list of TbTmatrixes, one list for each combi of levels in group_fac, list items are named by level_vs_level

# Calculation/implementation Steps (assume a count table of 100 taxons and 20 samples):
# Step 1: For each host taxon (e.g. taxon 1) calculate the within-sample taxon 1/taxon ratios and log these ratios,
# i.e. log(count taxon 1/count taxon 1:100) which is simplyd log(count taxon 1) - log(count taxon 1:100).
# Step 2: ratios/values where on of the counts were 0 are set to NA, and the log(geometric mean) over all samples
# is subtracted. Done

calculate_TbTmatrixes = function(physeq, group_var){
        
        if (taxa_are_rows(physeq)) {physeq <- t(physeq)}
        
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # -- get the level combis --
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        i <- ivec[x]
                })
        })
        j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        j <- jvec[x]
                })
        })
        i_s <- i_s[upper.tri(i_s)]
        j_s <- j_s[upper.tri(j_s)]
        # ----
        
        TbTmatrixes_list <- vector("list", length = length(i_s))
        
        for (k in seq_along(i_s)) {
                
                i <- i_s[k]
                j <- j_s[k]
                
                
                CT <- t(as(otu_table(physeq), 'matrix')) # now taxa are rows and samples are columns
                group_fac_current <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])
                CT <- CT[, group_fac %in% group_fac_current]
                
                
                # == Take in-sample taxa by taxa ratios and take the log ==
                # NB: log(x/y) = log(x) - log(y)
                TbTmatrixes <- lapply(1:nrow(CT), function(i){apply(CT, 2, function(samp_cnts){log(samp_cnts[i]) - log(samp_cnts)})})
                # produces for each taxon (= host taxon) a TbTMatrix
                # NB: there are -Inf, Inf, and NaN values in the matrixes, specifically
                # 0/x = log(0) - log(x) = -Inf, x/0 = log(x) - log(0) = Inf; 0/0 = log(0) - log(0) = NaN!
                # Inf, -Inf, and NaN will be ignored for calculating the rowMeans (geometric means) in the
                # next step
                # ====
                
                # == dividing each row of a TbTMatrix by its geometric mean, and taking log ==
                # NB: because the ratios are already "logged" it is just subtracting the rowMeans, see NB2
                # NB2: remember (see keep Note: #statistics #work geometric mean): geometric_mean(x) with x = x1, ... xn = (x1*...*xn)^1/n = exp(mean(log(x))), i.e. exp((1/n)*(log(x1)+ ... + log(xn))). 
                # Therefore, log(geometric_mean(x)) = (1/n)*(log(x1)+ ... + log(xn)) (= rowMeans when row = log(x1), ... log(xn))
                # Therefore: the rowMeans of the "logged" ratios are log(geometric_mean), and thus:
                # log(Ratio/geometric_mean) = log(Ratio) - rowMean! Remember log(Ratio) is currently in the TbTMatrixes
                # NB3: all Inf, -Inf, and NaN are set to NA, thus all 0 in the matrixes are from ratios x/y where x = y and x AND y != 0.
                TbTmatrixes <- lapply(TbTmatrixes, function(mat) {
                        mat[!is.finite(mat)] <- NA # puts all Inf, -Inf, and also NaN to NA!
                        mat-rowMeans(mat, na.rm = TRUE)
                        #newM[is.na(newM)] <- 0
                })
                # ====
                names(TbTmatrixes) <- rownames(TbTmatrixes[[1]])
                TbTmatrixes_list[[k]] <- TbTmatrixes
                names(TbTmatrixes_list)[[k]] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        TbTmatrixes_list
        
}


####################################
## evaluate_TbTmatrixes: 
###################################
# NB: you need to give same physeq as in calculate_TbTmatrixes, of course two functions could be combined, would
# save some double calculations
# the function calculates the TbT_groupSums as "teststatistic" in the TbTmatrixes 
# it adds count infos from physeq, NB: Median_grp1 and so on are calculated only from non-zero
# values because TbT method ignores zeros
## Input: 
# - TbTmatrixes_list: The list with the list of TbTmatrixes for each level combination in group_var defined group factor
# - physeq: same physeq object as used in calculate_TbTmatrixes
# - group_var: defining the group factor in physeq, same again as in calculate_TbTmatrixes
# - p.adjust.method: here only for p.adjust on fisher p_values
## Output: 
# - list of result DF for each combination of levels in group_var defined group fac

evaluate_TbTmatrixes <- function(TbTmatrixes_list, physeq, group_var, p.adjust.method = "fdr") {

        if(!identical(length(TbTmatrixes_list[[1]]), ntaxa(physeq))){stop("TbTmatrixes can not fit to physeq")}

        CT <- t(as(otu_table(physeq), "matrix"))
        CT[CT == 0] <- NA # because the method ignores 0 counts, I want to ignore them when calculating Mean_grp1 and so on
        
        
        group_fac <- factor(sample_data(physeq)[[group_var]])

        fac_levels <- levels(group_fac)

        # -- get the level combis --
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels)
        i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        i <- ivec[x]
                })
        })
        j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        j <- jvec[x]
                })
        })
        i_s <- i_s[upper.tri(i_s)]
        j_s <- j_s[upper.tri(j_s)]
        # ----

        result_list <- vector("list", length = length(i_s))

        for (k in seq_along(i_s)) {
                
                i <- i_s[k]
                j <- j_s[k]

                TbTmatrixes <- TbTmatrixes_list[[k]]
                group_fac_current <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])
                
                sumGrp1 <- sapply(TbTmatrixes, function(mat){sum(mat[, group_fac_current == fac_levels[i]], na.rm = T)})
                sumGrp2 <- sapply(TbTmatrixes, function(mat){sum(mat[, group_fac_current == fac_levels[j]], na.rm = T)})
                sumAll <- sapply(TbTmatrixes, function(mat){sum(mat, na.rm = T)})
                DF <- data.frame(Taxon = names(sumAll), sumGrp1 = sumGrp1, sumGrp2 = sumGrp2, sumAll = sumAll, TbT_groupSum = pmax(sumGrp1, sumGrp2))
                
                CT_current <- CT[, group_fac %in% group_fac_current]
                
                DF$Median_grp1 <- apply(CT_current[, group_fac_current == fac_levels[i]], 1, median, na.rm = T)
                DF$Median_grp2 <- apply(CT_current[, group_fac_current == fac_levels[j]], 1, median, na.rm = T)
                DF$Mean_grp1 <- apply(CT_current[, group_fac_current == fac_levels[i]], 1, mean, na.rm = T)
                DF$Mean_grp2 <- apply(CT_current[, group_fac_current == fac_levels[j]], 1, mean, na.rm = T)
                DF$Present_grp1 <- apply(CT_current[, group_fac_current == fac_levels[i]], 1, function(cnts){sum(!is.na(cnts))})
                DF$Present_grp2 <- apply(CT_current[, group_fac_current == fac_levels[j]], 1, function(cnts){sum(!is.na(cnts))})
                DF$Zeros_grp1 <- apply(CT_current[, group_fac_current == fac_levels[i]], 1, function(cnts){sum(is.na(cnts))})
                DF$Zeros_grp2 <- apply(CT_current[, group_fac_current == fac_levels[j]], 1, function(cnts){sum(is.na(cnts))})
                n1 <- sum(group_fac_current == fac_levels[i])
                n2 <- sum(group_fac_current == fac_levels[j])
                DF$Sparsity_grp1 <- 100*(DF$Zeros_grp1/n1)
                DF$Sparsity_grp2 <- 100*(DF$Zeros_grp2/n2)
                DF$n1 <- n1
                DF$n2 <- n2
                # -- add fisher exact test of presence differences (should be none in simulation) --
                Fisher <- t(sapply(1:nrow(DF), FUN = function(i){
                        fisherMat <- matrix(c(DF$Present_grp1[i], DF$Zeros_grp1[i], DF$Present_grp2[i],
                                              DF$Zeros_grp2[i]), ncol = 2, dimnames = list(c("Present", "Absent"), c("grp1", "grp2")))
                        Test <- fisher.test(fisherMat)
                        cbind(Test$p.value, Test$estimate)
                }))
                
                p_val_adj_Fisher <- p.adjust(Fisher[,1], method = p.adjust.method)
                
                DF$p_val_Fisher <- Fisher[,1]
                DF$p_val_Fisher_adj <- p_val_adj_Fisher
                DF$oddsRatioFisher <- Fisher[,2]
                
                direction <- rep("down", nrow(DF))
                direction[DF$Median_grp1 > DF$Median_grp2] <- "up"
                significance_fisher <- rep("", nrow(DF))
                significance_fisher[p_val_adj_Fisher <= 0.05] <- "*"
                significance_fisher[p_val_adj_Fisher <= 0.01] <- "**"
                significance_fisher[p_val_adj_Fisher <= 0.001] <- "***"
                direction_fisher <- rep("down", nrow(DF))
                direction_fisher[DF$Sparsity_grp1 < DF$Sparsity_grp2] <- "up"
                
                DF <- dplyr::mutate(DF, direction = direction, significance_fisher = significance_fisher, direction_fisher = direction_fisher)
                
                DF <- dplyr::select(DF, Taxon, TbT_groupSum, direction, p_val_Fisher,
                                     p_val_Fisher_adj, significance_fisher,
                                     direction_fisher, Median_grp1, Median_grp2, Mean_grp1,
                                     Mean_grp2, Present_grp1, Present_grp2, Zeros_grp1, Zeros_grp2, Sparsity_grp1, Sparsity_grp2, 
                                     n1, n2, sumGrp1, sumGrp2, sumAll, oddsRatioFisher)
                
                
                DF <- cbind(DF, tax_table(physeq))
                DF <- dplyr::arrange(DF, desc(abs(TbT_groupSum)))
                result_list[[k]] <- DF
                names(result_list)[[k]] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        result_list
}


####################################
## evaluate_TbTmatrixes_wilcoxTest: 
###################################
# An attempt to analyse TbTMatrixes with wilcox.test, similar to wilcoxTestApply. Basically simply applying wilcox.test on the colSums or colMeans
# of the TbTMatrixes, i.e. the sum or the mean over all taxa for each sample
## Input: 
# - TbTMatrixes_list: The list with the list of TbTMatrixes for each simulation
# - simlist: The list of physeq objects that have been used to calculate the TBTMatrixes_list
# - classlabel: the name of the column in the simlist physeq objects with the factor that defines the two sample groups
#  usually in your simulations = sample_data(simlist[[1]])$postfix
# - colSums: if TRUE colSums are used, if FALSE colMeans
# - excludeZeros: as in wilcoxTestApply, probably here should always be TRUE, samples where the host taxon was not present are excluded from the wilcox.test.
## Output: 
# - list of result DF for each simulation ordered after the wilcox.test standardized test statistic
# - also here a fisher exact test is added to compare sparsity levels.


evaluate_TbTmatrixes_wilcoxTest <- function(TbTmatrixes_list, physeq, group_var, p.adjust.method = "fdr", type = "median") {

        if(!identical(length(TbTmatrixes_list[[1]]), ntaxa(physeq))){stop("TbTmatrixes can not fit to physeq")}

        CT <- t(as(otu_table(physeq), "matrix"))
        CT[CT == 0] <- NA # because the method ignores 0 counts, I want to ignore them when calculating Mean_grp1 and so on


        group_fac <- factor(sample_data(physeq)[[group_var]])

        fac_levels <- levels(group_fac)

        # -- get the level combis --
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels)
        i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        i <- ivec[x]
                })
        })
        j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        j <- jvec[x]
                })
        })
        i_s <- i_s[upper.tri(i_s)]
        j_s <- j_s[upper.tri(j_s)]
        # ----

        result_list <- vector("list", length = length(i_s))

        for (k in seq_along(i_s)) {

                i <- i_s[k]
                j <- j_s[k]

                TbTmatrixes <- TbTmatrixes_list[[k]]
                group_fac_current <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])

                if(type == "sum"){
                        measureMatrix <- sapply(TbTmatrixes, colSums, na.rm = T)
                        # NB: samples in which host taxon was 0 are all NA in the TbTMatrix, colSums of all NA samples = 0!
                        measureMatrix[measureMatrix == 0] <- NA # assuming that it is practically impossible to get 0 if not the host taxon was absent in the sample
                } else if (type == "median") {
                        measureMatrix <- sapply(1:length(TbTmatrixes), function(e){
                                mat <- TbTmatrixes[[e]]
                                mat <- mat[-e, ] # removes the all 0 row where host taxon was compared to itself
                                apply(mat, 2, median, na.rm = T)

                        })
                        colnames(measureMatrix) <- taxa_names(physeq)

                        # measureMatrix[is.na(measureMatrix)] <- 0
                } else {
                        stop("incorrect type, must be sum or median")
                }

                res_mat <- apply(measureMatrix, 2, function(taxon_measures){
                        x <- taxon_measures[group_fac_current == fac_levels[i]]
                        Zeros_grp1 <- sum(is.na(x))
                        Sparsity_grp1 <- 100*(Zeros_grp1/length(x))
                        Present_grp1 <- length(x)-Zeros_grp1
                        x <- x[!is.na(x)]
                        Median_grp1 <- median(x, na.rm = T)
                        Mean_grp1 <- mean(x, na.rm = T)
                        if (is.na(Mean_grp1)){ Mean_grp1 = NA }

                        y <- taxon_measures[group_fac_current == fac_levels[j]]
                        Zeros_grp2 <- sum(is.na(y))
                        Sparsity_grp2 <- 100*(Zeros_grp2/length(y))
                        Present_grp2 <- length(y)-Zeros_grp2
                        y <- y[!is.na(y)]
                        Median_grp2 <- median(y, na.rm = T)
                        Mean_grp2 <- mean(y, na.rm = T)
                        if (is.na(Mean_grp2)){ Mean_grp2 = NA }

                        # I think you might just as well do a t.test here at least for median version, because after log transform data should almost be normal
                        # you could add a normality test
                        if (length(x) != 0 && length(y) != 0){
                                wilcTest <- wilcox.test(x = x, y = y, alternative = "two", paired = F, exact = F)
                                pValue <- wilcTest$p.value
                                W <- wilcTest$statistic
                                # calculate standardized rank sum Wilcoxon statistics as in multtest::mt.minP
                                Ranks <- rank(c(x, y))
                                n1 <- length(x) # should be always equal to Present_grp1 here
                                n2 <- length(y)
                                # Wx <- sum(Ranks[1:n1])-(n1*(n1+1)/2) # would be the same as W
                                # how about the other W?
                                # Wy <- sum(Ranks[(n1+1):n2]) - (n2*(n2+1)/2)
                                standStat <- -1*((sum(Ranks[1:n1]) - n1*(n1+n2+1)/2)/sqrt(n1*n2*(n1+n2+1)/12))

                                # # if you want to check that multtest::mt.minP would give the same statistic
                                # mati <- matrix(c(x,y), nrow = 1)
                                # grFac <- c(rep(fac_levels[i], n1), rep(fac_levels[j], n2))
                                # grFac <- factor(grFac, levels = c(fac_levels[i], fac_levels[j]))
                                # standStat2 <- multtest::mt.minP(mati, grFac, test = "wilcoxon")$teststat
                                # # identical(standStat, standStat2) # TRUE
                                # uncomment all with standStat2 to test all the way

                        } else {
                                pValue = NA
                                W <- NA
                                standStat = NA
                                n1 <- length(x)
                                n2 <- length(y)
                                # standStat2 = NA
                        }
                        
                        fisherMat <- matrix(c(Present_grp1, Zeros_grp1, Present_grp2, Zeros_grp2), 
                                            ncol = 2, dimnames = list(c("Present", "Absent"), c("grp1", "grp2")) )
                        Test <- fisher.test(fisherMat)
                        
                        c(teststat = standStat, p_val = pValue, Median_grp1 = Median_grp1, Median_grp2 = Median_grp2, 
                          Mean_grp1 = Mean_grp1, Mean_grp2 = Mean_grp2, n1 = n1, n2 = n2, Present_grp1 = Present_grp1, 
                          Present_grp2 = Present_grp2, Zeros_grp1 = Zeros_grp1, Zeros_grp2 = Zeros_grp2, 
                          Sparsity_grp1 = Sparsity_grp1, Sparsity_grp2 = Sparsity_grp2, W, 
                          p_val_Fisher = Test$p.value, Test$estimate) #, teststat2 = standStat2
                        
                })
                
                res_mat <- t(res_mat)
                colnames(res_mat) <- c("teststat", "p_val", "Median_grp1", "Median_grp2", "Mean_grp1", "Mean_grp2", "n1",
                                       "n2", "Present_grp1", "Present_grp2", "Zeros_grp1", "Zeros_grp2", "Sparsity_grp1", "Sparsity_grp2", "W",
                                       "p_val_Fisher", "oddsRatioFisher") #, , "teststat2"
                wilc.adjusted <- p.adjust(res_mat[,"p_val"], method = p.adjust.method)
                fisher.adjusted <- p.adjust(res_mat[,"p_val_Fisher"], method = p.adjust.method)
                significance <- rep("", length(wilc.adjusted))
                significance[wilc.adjusted <= 0.05] <- "*"
                significance[wilc.adjusted <= 0.01] <- "**"
                significance[wilc.adjusted <= 0.001] <- "***"
                direction <- rep("down", length(wilc.adjusted))
                direction[res_mat[, "Median_grp1"] > res_mat[, "Median_grp2"]] <- "up"
                significance_fisher <- rep("", length(fisher.adjusted))
                significance_fisher[fisher.adjusted <= 0.05] <- "*"
                significance_fisher[fisher.adjusted <= 0.01] <- "**"
                significance_fisher[fisher.adjusted <= 0.001] <- "***"
                direction_fisher <- rep("down", length(fisher.adjusted))
                direction_fisher[res_mat[, "Sparsity_grp1"] < res_mat[, "Sparsity_grp2"]] <- "up"
                
                
                DF <- data.frame(Taxon = rownames(res_mat), res_mat, p_val_adj = wilc.adjusted, significance = significance,
                                 direction = direction, p_val_Fisher_adj = fisher.adjusted,
                                 significance_fisher = significance_fisher,
                                 direction_fisher = direction_fisher)
                
                
                DF <- dplyr::select(DF, 1:3, 19:21, 17, 22:24, 4:7, 10:15, 8:9, 16, 18)
                # teststat/standStat2 version:
                # DF <- dplyr::select(DF, 1:2, 19, 3, 20:22, 17, 23:25, 4:7, 10:15, 8:9, 16, 18)
                
                DF <- cbind(DF, tax_table(physeq))
                # DF <- dplyr::arrange(DF, desc(abs(teststat)))
                DF <- dplyr::arrange(DF, p_val)
                
                result_list[[k]] <- list(DF, measureMatrix)
                names(result_list)[[k]] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        result_list
}



####################################
## calculate_raw_TbTmatrixes:
###################################
# see calculate_TbTmatrixes, this one just calculates the "raw" ratio matrixes

## Input: 
# - physeq = a phyloseq object
# - group_var, name of the group_fac in sample_data(physeq)
## Output: 
# - list of TbTmatrixes, one list for each combi of levels in group_fac, list items are named by level_vs_level


calculate_raw_TbTmatrixes = function(physeq, group_var){
        
        if (taxa_are_rows(physeq)) {physeq <- t(physeq)}
        
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # -- get the level combis --
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        i <- ivec[x]
                })
        })
        j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        j <- jvec[x]
                })
        })
        i_s <- i_s[upper.tri(i_s)]
        j_s <- j_s[upper.tri(j_s)]
        # ----
        
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


####################################
## create_TbT_TilePlots: 
###################################
# NB for a lot of taxa this function takes a lot of time, see 100 x 100 taxa is already 10000 wilcoxon tests
# the function calculates ntaxa x ntaxa matrixes of wilcoxon test p-values for TbTmatrixes list in TbTmatrixes_list
# i.e. for each level comparison within the group factor defined by group_var in physeq:)
# Based on the TbTMatrixes (probably ratios better than gm log ones) the p-value matrix
# indicates which taxa is enriched compared to which other taxa.
# in addition it adds a tile plot for each pValueMatrix
# wilcoxon.test is used
## Input: 
# - TbTmatrixes_list: The list with the lists of TbTmatrixes for level combi in group factor
# - physeq: used for TbTmatrixes_list generation
# - NB: so far no pvalue adjustment
## Output: 
# - list of pValMatrixes plus TileTr for each level combination. 
# NB: I negated the p-values if host taxon was more abundant in grp2 compared to other taxon!!

create_TbT_TilePlots <- function(TbTmatrixes_list, physeq, group_var, tax_names = NULL,
                                 test = "wilcoxon") {
        
        if(!identical(length(TbTmatrixes_list[[1]]), ntaxa(physeq))){stop("TbTmatrixes can not fit to physeq")}
        
        if(test != "wilcoxon" & test != "t.test"){stop("test unknown, must be wilcoxon or t.test")}
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # -- get the level combis --
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels)
        i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        i <- ivec[x]
                })
        })
        j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        j <- jvec[x]
                })
        })
        i_s <- i_s[upper.tri(i_s)]
        j_s <- j_s[upper.tri(j_s)]
        # ----
        
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
                                # NB: it is possible that x is completely NA (if taxon was only present when host taxon was not!), therefore
                                if(sum(!is.na(x)) == 0){x[1] <- 0}
                                y <- taxon_ratios[group_fac_current == fac_levels[j]]
                                if(sum(!is.na(y)) == 0){y[1] <- 0}
                                if (test == "wilcoxon"){
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
                                } else if (test == "t.test") {
                                        if (sum(!is.na(x)) < 2 || sum(!is.na(y)) < 2) {
                                                pValue <- 1
                                        } else {
                                                pValue <- t.test(x = x, y = y, alternative = "two")$p.value
                                        }
                                        if (mean(x, na.rm = T) > mean(y, na.rm = T)){pValue <- -1*pValue}
                                        pValue
                                }
                                
                        })
                })
                
                # make sure diagonal is all NA (can be exceptions especially for t.test)
                diag(pValMatrix) <- NA
                
                
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

####################################
## create_raw_TbT_TilePlots: 
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

create_raw_TbT_TilePlots <- function(TbTmatrixes_list, physeq, group_var, tax_names = NULL,
                                 test = "wilcoxon", p_adjust = "none") {
        
        if(!identical(length(TbTmatrixes_list[[1]]), ntaxa(physeq))){stop("TbTmatrixes can not fit to physeq")}
        
        if(test != "wilcoxon" & test != "t.test"){stop("test unknown, must be wilcoxon or t.test")}
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # -- get the level combis --
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels)
        i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        i <- ivec[x]
                })
        })
        j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        j <- jvec[x]
                })
        })
        i_s <- i_s[upper.tri(i_s)]
        j_s <- j_s[upper.tri(j_s)]
        # ----
        
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


####################################
## create_taxa_ratio_violin_plots: 
###################################
# wilcoxon.test is used
## Input: 
# - TbTmatrixes_list: The list with the lists of TbTmatrixes for each level combi in group factor
# NB: here it should be raw_TbTmatrixes with 0 = 0/x, Inf = x/0, NaN = 0/0, all these values will be ignored in the ratio plots by being set to NA!!
# - physeq: used for TbTmatrixes_list generation
# - group_var: the factor used to separate the samples
# - tax_names: the names of the taxa in physeq you want to use, e.g. taxa_names(physeq) if you like them, if NULL T1 to Tn will be used
# - taxa_nom: the nominator taxon of the abundance ratios, NB: only 1 allowed, must be included in tax_names
# - taxa_den: the denominator taxa of the abundance ratios, several allowed, all must be included in tax_names, if NULL all are used,
# i.e you get plots facet_wrapped around the taxa_den taxa
# - test: either "t.test" or "wilcoxon"
# - p_adjust: default "BH", method in p.adjust

## Output: 
# - list of pVals data frame (the pValues from t.test or wilcox.test after p.adjust) plus Violin plot,
# so for each level combi in group_var a list of length 2


create_taxa_ratio_violin_plots <- function(TbTmatrixes_list, physeq, group_var, tax_names = NULL,
                                     taxa_nom = "Firmicutes", taxa_den = NULL, test = "t.test", p_adjust = "BH") {
        
        if(!identical(length(TbTmatrixes_list[[1]]), ntaxa(physeq))){stop("TbTmatrixes can not fit to physeq")}
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # -- get the level combis --
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels)
        i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        i <- ivec[x]
                })
        })
        j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        j <- jvec[x]
                })
        })
        i_s <- i_s[upper.tri(i_s)]
        j_s <- j_s[upper.tri(j_s)]
        # ----
        
        result_list <- vector("list", length = length(i_s))
        
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        
        # Color the sample names based on levels in group fac
        if (length(levels(LookUpDF$Group)) <= 7){
                color_lookup <- data.frame(level = levels(LookUpDF$Group), color = cbPalette[2:(length(levels(LookUpDF$Group)) + 1)])
        } else {
                color_lookup <- data.frame(level = levels(LookUpDF$Group), color = viridis(length(levels(LookUpDF$Group))))
        }
        for (e in seq_along(color_lookup)) {
                color_lookup[, e] <- as.character(color_lookup[,e])
        }
        colors_to_use <- color_lookup$color
        names(colors_to_use) <- color_lookup$level
        
        for (k in seq_along(i_s)) {
                
                i <- i_s[k]
                j <- j_s[k]
                
                TbTmatrixes <- TbTmatrixes_list[[k]]
                
                if (is.null(tax_names)){
                        tax_names <- paste("T", 1:length(TbTmatrixes), sep = "_")
                } else {
                        if(!identical(length(TbTmatrixes), length(tax_names))){stop("tax_names do not fit in length to TbTmatrixes")}
                }
                
                names(TbTmatrixes) <- tax_names
                
                TbTmatrix <- TbTmatrixes[[taxa_nom]]
                
                if (is.null(TbTmatrix)) {stop("taxa_nom not found")}
                
                rownames(TbTmatrix) <- tax_names
                
                TbT_DF <- as.data.frame(TbTmatrix)
                TbT_DF$Taxon <- rownames(TbT_DF)
                if (is.null(taxa_den)) {taxa_den <- tax_names}
                TbT_DF <- TbT_DF[TbT_DF$Taxon %in% taxa_den, ]
                TbT_DF_l <- gather(TbT_DF, key = Sample, value = Ratio, -Taxon)
                TbT_DF_l$Group <- as.character(LookUpDF$Group[match(TbT_DF_l$Sample, LookUpDF$Sample)])
                group_fac_current <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])
                TbT_DF_l$Group <- factor(TbT_DF_l$Group, levels = levels(group_fac_current), ordered = T)
                TbT_DF_l$Ratio[!is.finite(TbT_DF_l$Ratio) | TbT_DF_l$Ratio == 0] <- NA
                var_plus_length_check <- group_by(TbT_DF_l, Taxon, Group) %>% summarise(Variance = var(Ratio, na.rm = T), NotNA = sum(!is.na(Ratio)))
                if (test == "t.test"){
                        var_plus_length_check <- filter(var_plus_length_check, !(Variance > 0) | NotNA < 2)
                } else if (test == "wilcoxon") {
                        var_plus_length_check <- filter(var_plus_length_check, NotNA < 1)
                }
                
                if (nrow(var_plus_length_check) != 0) {
                        TbT_DF_l_test <- filter(TbT_DF_l, !(Taxon %in% unique(var_plus_length_check$Taxon)))
                } else {
                        TbT_DF_l_test <- TbT_DF_l
                }
                
                # needs probably some sanity checks here on TbT_DF_l_test
                
                TbT_DF_l_test <- group_by(TbT_DF_l_test, Taxon)
                
                if (test == "t.test"){
                        pVals <- summarise(TbT_DF_l_test, pVal = t.test(Ratio ~ Group, alternative = "two")$p.value)
                } else if (test == "wilcoxon"){
                        pVals <- summarise(TbT_DF_l_test, pVal = wilcox.test(Ratio ~ Group, alternative = "two", paired = F, exact = F)$p.value)
                }
                
                pVals$padj <- p.adjust(pVals$pVal, method = p_adjust)
                pVals$significance <- ""
                pVals$significance[pVals$padj <= 0.05] <- "*"
                pVals$significance[pVals$padj <= 0.01] <- "**"
                pVals$significance[pVals$padj <= 0.001] <- "***"
                pVals$Name <- paste(pVals$Taxon, pVals$significance)
                
                TbT_DF_l$Taxon[TbT_DF_l$Taxon %in% pVals$Taxon] <- pVals$Name[match(TbT_DF_l$Taxon[TbT_DF_l$Taxon %in% pVals$Taxon], pVals$Taxon)]
                # order levels of TbT_DF_l$Taxon based on pValues
                pVals <- arrange(pVals, pVal)
                
                TbT_DF_l$Taxon <- factor(TbT_DF_l$Taxon, levels = c(pVals$Name, unique(TbT_DF_l$Taxon)[!(unique(TbT_DF_l$Taxon) %in% pVals$Name)]), ordered = T)
                
                Tr <- ggplot(TbT_DF_l, aes(x = Group, y = Ratio, col = Group))
                Tr <- Tr +
                        geom_violin() +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                        # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                        scale_color_manual(values = colors_to_use) +
                        facet_wrap(~ Taxon, scales = "free_y") +
                        xlab("") +
                        ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                        theme_bw() +
                        theme(legend.position = "none")
  
                
                result_list[[k]] <- list(pVals, Tr)
                rm(Tr, TbT_DF_l, TbT_DF, TbT_DF_l_test, pVals)
        }
        
        names(result_list) <- names(TbTmatrixes_list)
        result_list
}



####################################
## rank_evaluate_TbTmatrixes: 
###################################
# NB: you need to give same physeq as in calculate_TbTmatrixes, of course two functions could be combined, would
# save some double calculations
# the function calculates the TbT_groupSums as "teststatistic" in the TbTmatrixes 
# it adds count infos from physeq, NB: Median_grp1 and so on are calculated only from non-zero
# values because TbT method ignores zeros
## Input: 
# - TbTmatrixes_list: The list with the list of TbTmatrixes for each level combination in group_var defined group factor
# - physeq: same physeq object as used in calculate_TbTmatrixes
# - group_var: defining the group factor in physeq, same again as in calculate_TbTmatrixes
# - p.adjust.method: here only for p.adjust on fisher p_values
## Output: 
# - list of result DF for each combination of levels in group_var defined group fac

rank_evaluate_TbTmatrixes <- function(TbTmatrixes, group_fac, p.adjust.method = "fdr") {
        
        for (i in seq_along(TbTmatrixes)) {
                
                mat <- TbTmatrixes[[i]]
                
                # calculate ranks, and transform NaN (i.e. where both host taxon and denominator taxon were absent to NA)
                
                mat_rank <- t(apply(mat, 1, rank, na.last = "keep"))
                # NB: in case no NA then rowSums(mat_rank) is equal to nsamples*(nsamples + 1)/2
                
                # next steps:
                # - divide by geometric mean and log so rowSums are 0 in each row
                # - set all NA to 0
                # - get colSums to get evaluation matrix1
                # - run wilcoxon test
                
               
        }
        
        
}



#######################################
### gm_own_tbt: calculate geometric mean
#######################################
# see above for gm_own, just a simpler version and based on the fact that NA might be in the data


gm_own_tbt = function(x){
        
        exp(sum(log(x), na.rm = TRUE) / length(x[!is.na(x)]))

}

get_assignment_distribution <- function(taxa){
  
  countNA <- apply(taxa, 2, function(x){sum(is.na(x))})
  countNonNA <- apply(taxa, 2, function(x){sum(!is.na(x))})
  total <- countNA + countNonNA
  assignment_distribution <- data.frame(assigned = countNonNA, 
                                        total = total,
                                        PC_assigned = round(100*countNonNA/total, 1)
  )
  return(assignment_distribution)
}


get_pretty_taxon_name = function(x) {
  name <- "Unknown";
  if (!is.na(x['Kingdom'])) {
    name <- x['Kingdom']
  }
  if(!is.na(x['Phylum'])) {
    name <- x['Phylum']
  }
  if(!is.na(x['Class'])) {
    name <- x['Class']
  }
  if(!is.na(x['Order'])) {
    name <- x['Order']
  }
  if(!is.na(x['Family'])) {
    name <- x['Family']
  }
  if(!is.na(x['Genus'])) {
    name <- x['Genus']
  }
  
  if (!is.na(x['Species'])) {
    if(x['Species']=="bacterium"){name<-x['Genus']}
  else{
    name <-paste(x['Genus'], x['Species'])
    name<- strsplit(name, split = "/")
    if(length(name[[1]])>1){
    name <- paste(x['Genus'], "(Multi-Affiliation)")
    }else{name <-paste(x['Genus'], x['Species'])}
    
  }
  }
id<-paste("[", x["Taxon"], "]",sep = '')
name1<-paste(id,name, x["sig"])
 
  return(name1)
}


get_pretty_taxon_name1 = function(x) {
  name <- "Unknown";
  if (!is.na(x['Kingdom'])) {
    name <- x['Kingdom']
  }
  if(!is.na(x['Phylum'])) {
    name <- x['Phylum']
  }
  if(!is.na(x['Class'])) {
    name <- x['Class']
  }
  if(!is.na(x['Order'])) {
    name <- x['Order']
  }
  if(!is.na(x['Family'])) {
    name <- x['Family']
  }
  if(!is.na(x['Genus'])) {
    name <- x['Genus']
  }
  if (!is.na(x['Species'])) {
    if(x['Species']=="bacterium"){name<-x['Genus']}
    else{
      name <-paste(x['Genus'], x['Species'])
      name<- strsplit(name, split = "/")
      if(length(name[[1]])>1){
        name <- paste(x['Genus'], "(Multi-Affiliation)")
      }else{name <-paste(x['Genus'], x['Species'])}
      
    }
  }
  #id<-paste("[", x["Taxon"], "]",sep = '')
  #id<-x["Taxon"]
  #name1<-paste(name,id,sep = " ")
  name1<-paste(name, x["sig"])
  return(name1)
}


get_pretty_taxon_name2 = function(x) {
  name <- "Unknown";
  if (!is.na(x['Kingdom'])) {
    name <- x['Kingdom']
  }
  if(!is.na(x['Phylum'])) {
    name <- x['Phylum']
  }
  if(!is.na(x['Class'])) {
    name <- x['Class']
  }
  if(!is.na(x['Order'])) {
    name <- x['Order']
  }
  if(!is.na(x['Family'])) {
    name <- x['Family']
  }
  if(!is.na(x['Genus'])) {
    name <- x['Genus']
  }
  if (!is.na(x['Species'])) {
    if(x['Species']=="bacterium"){name<-x['Genus']}
    else{
      name <-paste(x['Genus'], x['Species'])
      name<- strsplit(name, split = "/")
      if(length(name[[1]])>1){
        name <- paste(x['Genus'], "(Multi-Affiliation)")
      }else{name <-paste(x['Genus'], x['Species'])}
      
    }
  }
  id<-paste("[", x["Taxon"], "]",sep = '')
  name1<-paste(id,name)
  
  return(name1)
}


get_pretty_taxon_name_fun = function(x) {
  name <-x['Module']
  id<-paste("[", x["Taxon"], "]",sep = '')
  name1<-paste(id,name, x["sig"])
  
  return(name1)
}




get_pretty_taxon_name_path = function(x) {
  name <-x["Pathway"]
  id<-paste("[", x["Taxon"], "]",sep = '')
  name1<-paste(id,name, x["sig"])
  
  return(name1)
}

get_pretty_taxon_name_fun_2 = function(x) {
  name <-paste(x['Level_2'],"/",x['Module'], sep = '')
  id<-paste("[", x["Taxon"], "]",sep = '')
  name1<-paste(id,x['Module'])
  
  return(name1)
}

get_pretty_taxon_name_path_2 = function(x) {
  name <-x["Pathway"]
  id<-paste("[", x["Taxon"], "]",sep = '')
  name1<-paste(id,x['Pathway'])
  
  return(name1)
}


####pairwise.adonis #https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R

pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='fdr')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0.0001 *** 0.001 ** 0.01 * 0.05 .")
  return(pairw.res)
  
} 


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


loop_vegan_adonis <- function(dist_obj, group_fac, nperm = 999, 
                              p.adj.method = "none", symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) {
  
  if (!("dist" %in% class(dist_obj))){
    stop("dist_obj must be of class dist")
  }
  
  group_fac <- factor(group_fac)
  
  # - get all level combinations within group_fac - # changed this so that I think order is more logical, check regularly if correct
  fac_levels <- levels(group_fac)
  fac_levels_num <- setNames(seq_along(fac_levels), fac_levels)
  i_sAndj_s <- get_unique_facLevel_combinations(fac_levels)
  i_s <- i_sAndj_s[["i_s"]]
  j_s <- i_sAndj_s[["j_s"]]
  # -- 
  
  p_vals <- vector(mode = "numeric", length = length(i_s))
  r2s <- vector(mode = "numeric", length = length(i_s))
  
  # - loop vegan::adonis for all level combinations and generate df -
  for (k in seq_along(i_s)){
    i <- i_s[k]
    j <- j_s[k]
    group_fac2 <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])
    dist_obj_mat <- as.matrix(dist_obj)
    rows <- which(group_fac %in% levels(group_fac2))
    dist_obj2 <- as.dist(dist_obj_mat[rows, rows])
    fit <- vegan::adonis(dist_obj2 ~ group_fac2, permutations = nperm)
    p_vals[k] <- fit$aov.tab[1, "Pr(>F)"]
    r2s[k] <- fit$aov.tab[1, "R2"]
  }
  
  
  
  result_df <- data.frame(Comparison = paste0(names(fac_levels_num[i_s]), "_to_", names(fac_levels_num[j_s]), sep = ""),
                          adonis_p_value = p_vals, adonis_R2 = r2s, p_val_adj = p.adjust(p_vals, p.adj.method))
  # --
  
  # - add vegan::adonis including all group levels -
  fit <- vegan::adonis(dist_obj ~ group_fac, permutations = nperm)
  df <- data.frame(Comparison = "Overall",
                   adonis_p_value = fit$aov.tab[1, "Pr(>F)"], 
                   adonis_R2 = fit$aov.tab[1, "R2"], 
                   p_val_adj = fit$aov.tab[1, "Pr(>F)"])
  
  result_df <- rbind(df, result_df)
  # -- 
  
  # - add significance levels -
  symnum.args$x <- result_df$adonis_p_value
  result_df$significance <- do.call(stats::symnum, symnum.args) %>% as.character()
  result_df <- mutate(result_df, adonis_R2_PC = round(100*adonis_R2, 2))
  # --
  
}


Deseq2SummaryResults<- function(DESeq2_outPut, p_val_adj, ps, Functional=FALSE, group_var, pathway=FALSE){
  
  library(DESeq2)
  DESeq2_result_list <- DESeq2_outPut[[1]]
  physeq_DESadjust <- DESeq2_outPut[[2]]
  
  suppressWarnings(head_values <- sapply(DESeq2_result_list, function(df){
    max(which(df$p_val_adj <= p_val_adj))
  }))
  
  
  to_keep<-which(is.finite(head_values))
  
  if(length(to_keep)>0){
  DESeq2_result_list<-DESeq2_result_list[to_keep]
  
  suppressWarnings(head_values <- sapply(DESeq2_result_list, function(df){
    max(which(df$p_val_adj <= p_val_adj))
  }))
  
  
  original_head_values <- data.frame(Comparison = c(names(head_values), "Total"), 
                                     NoSignificant = c(head_values, ntaxa(ps)), 
                                     PC_Significant = 100*c(head_values,ntaxa(ps))/ntaxa(ps))
  
  if(Functional){
    if(pathway){DES_plot_list <- lapply(1:length(DESeq2_result_list), function(i){
      DEseq2_Control<-Deseq_plot_fun_path(DESeq2_result_list[[i]], p_val_adj,label=FALSE)+ theme(axis.text.x = element_text(angle = 360, hjust = 0, vjust=0.5))
    })}else{
    DES_plot_list <- lapply(1:length(DESeq2_result_list), function(i){
      DEseq2_Control<-Deseq_plot_fun(DESeq2_result_list[[i]], p_val_adj,label=FALSE)+ theme(axis.text.x = element_text(angle = 360, hjust = 0, vjust=0.5))
    })}}else{
      DES_plot_list <- lapply(1:length(DESeq2_result_list), function(i){
      DEseq2_Control<-Deseq_plot(DESeq2_result_list[[i]], p_val_adj,label=FALSE)+ theme(axis.text.x = element_text(angle = 360, hjust = 0, vjust=0.5))})
    }
  names(DES_plot_list) <- names(DESeq2_result_list)
  
  
  
  if(Functional){
    
      if(pathway){DES_table_list <- lapply(1:length(DESeq2_result_list), function(i){
        df <- dplyr::filter(DESeq2_result_list[[i]],  p_val_adj<=0.05)
        df <- dplyr::select(df, 1,  Pathway,  4:6,24)
        colnames(df)[3:6] <- c("padj", "sig", "dir", "log2FC")
        df<-dplyr::arrange(df, (log2FC))
        df
      })}else{
            DES_table_list <- lapply(1:length(DESeq2_result_list), function(i){
              df <- dplyr::filter(DESeq2_result_list[[i]],  p_val_adj<=0.05)
              df <- dplyr::select(df, 1, Level_1, Level_2, Module,  4:6,24)
              colnames(df)[5:8] <- c("padj", "sig", "dir", "log2FC")
              df<-dplyr::arrange(df, (log2FC))
              df
            })}} else{
      
      DES_table_list<- lapply(1:length(DESeq2_result_list), function(i){
        df <- dplyr::filter(DESeq2_result_list[[i]],  p_val_adj<=0.05)
        df <- dplyr::select(df, 1, Phylum, Family, Genus, Species, 2, 4:6,24)
        colnames(df)[7:10] <- c("padj", "sig", "dir", "log2FC")
        df
      })
      
    }
  names(DES_table_list) <- names(DESeq2_result_list)
  
  DES_table_list <- lapply(1:length(DES_table_list), function(i){
    df <- DES_table_list[[i]]
  })
  
  names(DES_table_list) <- names(DESeq2_result_list)
  tax_orders <- lapply(1:length(DES_table_list), function(i){
    df <- DES_table_list[[i]]
    df$Taxon
  })
  names(tax_orders) <- names(DES_table_list)
  
  
  row_names_for_heat_maps_new <- lapply(1:length(DES_table_list), function(i){
    row_names <- apply(DES_table_list[[i]] , 1, if(Functional){if(pathway){get_pretty_taxon_name_path}else{
      get_pretty_taxon_name_fun}}else{get_pretty_taxon_name})
    
  })
  names(row_names_for_heat_maps_new) <- names(DES_table_list)
  
  
  row_names_df<- lapply(1:length(DES_table_list), function(i){
    row_names <- apply(DES_table_list[[i]] , 1, if(Functional){if(pathway){get_pretty_taxon_name_path_2}else{
      get_pretty_taxon_name_fun_2}}else{get_pretty_taxon_name2})
  })
  names(row_names_df) <- names(DES_table_list)
  
  DES_table_list_print <- lapply(1:length(DES_table_list), function(i){
    df<-cbind(row_names_df[[i]],DES_table_list[[i]][-c(1:5)])
    colnames(df)[1]<-"Taxon"
    df
  })
  names(DES_table_list_print) <- names(DES_table_list)
  
  
  pruned_physeqs <- lapply(1:length(DES_table_list), function(i){
    
    prune_taxa(tax_orders[[i]], ps)
  })
  
  names(pruned_physeqs) <- names(DES_table_list)
  
  
  heat_maps_DES <- lapply(1:length(pruned_physeqs), function(i){
    maps <- make_heat_map_physeq_levels(pruned_physeqs[[i]], group_var = group_var, max_abundance_for_color = NULL, tax_order = tax_orders[[i]], tax_names = row_names_for_heat_maps_new[[i]], color_sample_names = TRUE, gradient_steps = c(0.15, 0.3, 0.45, 1)) 
    maps[[i]]
  })
  names(heat_maps_DES) <- names(DESeq2_result_list)
  
  violin_plots_DES <- lapply(1:length(tax_orders), function(i){
    plotlist <- plot_top_abundances_boxAndviolin(physeq = pruned_physeqs[[i]], group_var = group_var, tax_order = tax_orders[[i]], tax_names = row_names_for_heat_maps_new[[i]])
    plotlist[[i]]
  })
  
  # there are 8 plots per list, only pick violin plot faceted, and logged abundance:
  violin_plots_DES <- lapply(violin_plots_DES, `[[`, 8)

  
  
  
  
  out <- list(original_head_values=original_head_values,DES_table_list=DES_table_list, DES_table_list_print= DES_table_list_print, DES_plot_list=DES_plot_list, heat_maps_DES=heat_maps_DES, violin_plots_DES=violin_plots_DES)}else{out<- NULL
                                                                                                                                                                                                                                     print("No significant results were found")}

   #return(out)
  }
#overlaping results
overlap_results2=function(df1,df2,compare=compare){
  A<-df1[,compare]
  B<-df2[,compare]
  C= A[(A %in% B)]
  if(length(C)!=0){
  A<-df1[C,c(2,7)]
  B<-df2[C,c(2,7)]
  D<-merge(A,B,by=0)
  colnames(D)<-c(compare, "HR.1st","p.value.1st" ,"HR.2nd","p.value.2nd")
  D}else(print("No overlaping")) 
}


overlap_results=function(df1,df2){
  A<-rownames(df1)
  B<-rownames(df2)
  C= A[(A %in% B)]
  if(length(C)!=0){
     D}else(print("No overlaping")) 
}

overlap_results3=function(df1,df2,compare=compare){
  A<-df1[,compare]
  B<-df2[,compare]
  C= A[(A %in% B)]
  row.names(df1)<-A
  row.names(df2)<-B
  if(length(C)!=0){
    A<-df1[C,c(2,7,8)]
    B<-df2[C,c(2,7,8)]
    D<-merge(A,B,by=0)
    row.names(D)<-D[,4]
    D<-D[,-c(4,7)]
    colnames(D)<-c(compare, "HR.1st","p.value.1st" ,"HR.2nd","p.value.2nd")
    D}else(print("No overlaping")) 
}

#tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")


writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

###################################################################################
#' Convert phyloseq OTU count data into DGEList for edgeR package
#' https://joey711.github.io/phyloseq-extensions/edgeR.html
#' Further details.
#' 
#' @param physeq (Required).  A \code{\link{phyloseq-class}} or
#'  an \code{\link{otu_table-class}} object. 
#'  The latter is only appropriate if \code{group} argument is also a 
#'  vector or factor with length equal to \code{nsamples(physeq)}.
#'  
#' @param group (Required). A character vector or factor giving the experimental
#'  group/condition for each sample/library. Alternatively, you may provide
#'  the name of a sample variable. This name should be among the output of
#'  \code{sample_variables(physeq)}, in which case
#'  \code{get_variable(physeq, group)} would return either a character vector or factor.
#'  This is passed on to \code{\link[edgeR]{DGEList}},
#'  and you may find further details or examples in its documentation.
#'  
#' @param method (Optional). The label of the edgeR-implemented normalization to use.
#'  See \code{\link[edgeR]{calcNormFactors}} for supported options and details. 
#'  The default option is \code{"RLE"}, which is a scaling factor method 
#'  proposed by Anders and Huber (2010).
#'  At time of writing, the \link[edgeR]{edgeR} package supported 
#'  the following options to the \code{method} argument:
#'  
#'  \code{c("TMM", "RLE", "upperquartile", "none")}.
#'
#' @param ... Additional arguments passed on to \code{\link[edgeR]{DGEList}}
#' #' 
#' @examples
#'

phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
  }
################################################################################

















# ----------------------------- obsolete functions ----------------------------------------------------

# #######################################
# ### adjust_LS
# #######################################
# # own implementation of DESeq2 library size adjustment with little tweaks
# 
# ## Input:
# # physeq
# # zeros.count: zeros.count of gm_own, if TRUE 0 will be considered when calculating the geometric mean
# # if FALSE not and thus the geometric means will be bigger (see gm_own)
# # percentile: the percentile to select the size factor SF of a sample based on its count ratios to the reference sample. In
# # DESeq percentile = 50, i.e. the median is used. NOTE when ignore.zero.ratios = FALSE, you might get SF = 0 in case that more 
# # than 50% of the taxa of in a sample are 0. Function will through a warning then. s
# # ignore.zero.ratios: if TRUE, the SF of each sample (i.e. the percentile) will be calculated based on the pool of non-zero ratios of the sample
# # plots: if TRUE SFs and plots will be given in an output list, otherwise just the adjusted physeq is returned
# 
# 
# adjust_LS <- function(physeq, zeros.count = FALSE, percentile = 50, ignore.zero.ratios = TRUE, plots = FALSE)  {
#         
#         # ---- Step 1: Use the geometric means for each taxa over all samples to calculate a "representative" Reference Sample ------
#         if(taxa_are_rows(physeq)){
#                 GM <- apply(otu_table(physeq), 1, gm_own, zeros.count = zeros.count)   
#         } else {
#                 GM <- apply(otu_table(physeq), 2, gm_own, zeros.count = zeros.count) 
#         }
#         
#         # ---- Step 2: Calculate Size factors --------
#         # -- 2a: Calculate Ratio Matrix by dividing counts of each sample by Reference Sample
#         
#         if(taxa_are_rows(physeq)){
#                 RatioMatrix <- sweep(otu_table(physeq), 1, GM, "/")
#         } else {
#                 RatioMatrix <- sweep(otu_table(physeq), 2, GM, "/") 
#         }
#         
#         # -- 2b: Insert: Calculate 0 Percentage for each sample
#         
#         if(taxa_are_rows(physeq)){
#                 zeroPC <- 100*(colSums(otu_table(physeq) == 0)/ntaxa(physeq))
#         } else {
#                 zeroPC <- 100*(rowSums(otu_table(physeq) == 0)/ntaxa(physeq))
#         }
#         
#         
#         # -- 2c: calculate SFs for each sample depending on given percentile and whether to ignore.zero.ratios
#         # if ignore.zero.ratios = TRUE, the percentile refers only to the non zero ratios in each sample
#         # NB: DESEQ uses ignore.zero.ratios = TRUE
#         if(ignore.zero.ratios){
#                 if (taxa_are_rows(physeq)) {
#                         SFs <- apply(RatioMatrix, 2, function(x){x <- x[x>0]; quantile(x, probs = percentile/100, na.rm = T)})
#                 } else {
#                         SFs <- apply(RatioMatrix, 1, function(x){x <- x[x>0]; quantile(x, probs = percentile/100, na.rm = T)})
#                 }  
#         } else {
#                 if (taxa_are_rows(physeq)) {
#                         SFs <- apply(RatioMatrix, 2, quantile, probs = percentile/100, na.rm = T)
#                 } else {
#                         SFs <- apply(RatioMatrix, 1, quantile, probs = percentile/100, na.rm = T) 
#                 }   
#         }
#         
#         if(min(SFs) == 0) { warning("in at least one sample the Size Factor was 0!") }
#         
#         
#         # --- 3: calculate the new counts and put into a physeq object
#         
#         if(taxa_are_rows(physeq)){
#                 if (!identical(colnames(otu_table(physeq)), names(SFs))) {stop("names SFs do not fit to physeq")}
#                 Mat <- sweep(otu_table(physeq),2,SFs, "/")
#                 phynew <- phyloseq(otu_table(Mat, taxa_are_rows = TRUE), sample_data(physeq), tax_table(physeq))
#         } else {
#                 if (!identical(rownames(otu_table(physeq)), names(SFs))) {stop("names SFs do not fit to physeq")}
#                 Mat <- sweep(otu_table(physeq),1,SFs, "/") 
#                 phynew <- phyloseq(otu_table(Mat, taxa_are_rows = FALSE), sample_data(physeq), tax_table(physeq))
#         }
#         
#         if (plots){
#                 
#                 # ---- step 4: Generate Plots that illustrate the process
#                 # -- 4a: calculate for each sample a histogram of the ratios used to determine the SF
#                 
#                 cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#                 
#                 # Define the histos function
#                 histos <- function(x, SF, zPC, SampleName) {
#                         x <- data.frame(x = x)
#                         Tr <- ggplot(x, aes(x = x))
#                         Tr <- Tr + geom_histogram(binwidth = 0.3, col = "black", fill = "#E69F00") +
#                                 geom_rug() +
#                                 geom_vline(xintercept = SF, col = "#009E73") +
#                                 ylab("Frequency") + 
#                                 xlab("Count/(Count of Ref Sample)") +
#                                 ggtitle(paste("Sample: ", SampleName, "; zeroPC: ", round(zPC, 2), "; size factor: ", round(SF,2), sep = "")) +
#                                 theme_bw() + 
#                                 theme(panel.grid.minor = element_blank(),
#                                       panel.grid.major.y = element_blank(),
#                                       panel.grid.major.x = element_line(color = "#999999", size = .15))
#                 }
#                 
#                 
#                 if(taxa_are_rows(physeq) && (!all.equal(colnames(RatioMatrix), names(SFs), names(zeroPC)))) {warning("names messed up!")}
#                 if(!taxa_are_rows(physeq) && (!all.equal(rownames(RatioMatrix), names(SFs), names(zeroPC)))) {warning("names messed up!")}
#                 
#                 if(ignore.zero.ratios){
#                         if(taxa_are_rows(physeq)){
#                                 HistoList <- lapply(1:nsamples(physeq), function(i){histos(x = RatioMatrix[,i][RatioMatrix[,i] > 0], SF = SFs[i], zPC = zeroPC[i], SampleName = names(SFs)[i])})
#                         } else {
#                                 HistoList <- lapply(1:nsamples(physeq), function(i){histos(x = RatioMatrix[i,][RatioMatrix[i,] > 0], SF = SFs[i], zPC = zeroPC[i], SampleName = names(SFs)[i])})
#                         }  
#                 } else {
#                         if(taxa_are_rows(physeq)){
#                                 HistoList <- lapply(1:nsamples(physeq), function(i){histos(x = RatioMatrix[,i], SF = SFs[i], zPC = zeroPC[i], SampleName = names(SFs)[i])})
#                         } else {
#                                 HistoList <- lapply(1:nsamples(physeq), function(i){histos(x = RatioMatrix[i,], SF = SFs[i], zPC = zeroPC[i], SampleName = names(SFs)[i])}) 
#                         }   
#                 }
#                 
#                 # -- 3b: compare calculated SFs to library sizes of the samples
#                 if(!identical(names(SFs), names(sample_sums(physeq)))){warning("Probably some Mix Up CHECK!")}
#                 comp <- data.frame(Sample = names(SFs), TotAmps = sample_sums(physeq), SFs = SFs)
#                 comp$SFsNormed <- comp$SFs/median(comp$SFs)
#                 comp$TotAmpsNormed <- comp$TotAmps/mean(comp$TotAmps)
#                 comp <- dplyr::arrange(comp, desc(TotAmpsNormed))
#                 comp$Sample <- as.factor(comp$Sample)
#                 LevelsWant <- as.character(comp$Sample)
#                 for(i in 1:length(LevelsWant)){
#                         comp$Sample <- relevel(comp$Sample, LevelsWant[i])
#                 }
#                 
#                 comp <- comp[c(1,4,5)]
#                 names(comp)[c(2,3)] <- c("SizeFactor_DESeq", "TotalAmplicons_relAb")
#                 comp <- tidyr::gather(comp, key = Corrector, value = NormedValue, -Sample)
#                 Tr <- ggplot(comp, aes(x = Sample, y = NormedValue, color = Corrector))
#                 Tr <- Tr + geom_point(size = 2) +
#                         xlab("") +
#                         ylab("correction value (normalized to median)") +
#                         ggtitle(paste("Median SF: ", round(median(SFs),3), " Median TotAmp: ", round(median(sample_sums(physeq)),3), sep = "")) +
#                         scale_color_manual(values = cbPalette[c(4,2)]) +
#                         theme_bw() +
#                         theme(panel.grid.minor = element_blank(),
#                               panel.grid.major.y = element_blank(),
#                               panel.grid.major.x = element_line(color = "#999999", size = .15),
#                               axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
#                               legend.title = element_blank())
#                 
#                 # -- 3c: save Histograms of the SFs and sample_sums
#                 histo2 <- function(x, xtitle, gtitle) {
#                         x <- data.frame(x = x)
#                         Tr <- ggplot(x, aes(x = x))
#                         Tr <- Tr + geom_histogram(binwidth = diff(range(x))/60, col = "black", fill = "#E69F00") +
#                                 geom_rug() +
#                                 geom_vline(xintercept = median(x$x), col = "#009E73", size = 1) +
#                                 ylab("No Samples") + 
#                                 xlab(xtitle) +
#                                 ggtitle(gtitle) +
#                                 theme_bw() + 
#                                 theme(panel.grid.minor = element_blank(),
#                                       panel.grid.major.y = element_blank(),
#                                       panel.grid.major.x = element_line(color = "#999999", size = .15))
#                 }
#                 
#                 Tr2 <- histo2(SFs, xtitle = "Size Factors", gtitle = "Size Factors a la DESeq")
#                 Tr3 <- histo2(sample_sums(physeq), xtitle = "Total Amplicons", gtitle = "Size Factors a la relative abundance")
#                 
#                 List <- list(Physeq = phynew, SFs = SFs, SizeFactorCompare = Tr, SFHistos = list(Tr2, Tr3), RefSample = GM, RatioMatrix = RatioMatrix, zeroPC = zeroPC, HistoList = HistoList)
#                 
#                 
#         } else {
#                 
#                 phynew
#                 
#         }
#         
#         
# }


# #######################################
# ### FUNCTION: plotTaxLevelvsprevalence
# #######################################
# 
# # Function to determine how many sequences could be assigned to the different taxonomic levels compared to the number of smaples the amplicons are present
# ## Input
# # taxa: the output of the assignTaxonomy (dada2) and maybe addSpecies command. nrow = number of sequences and ncol number taxonomic levels.
# # seqtab: The abundance table output from Dada2_wrap
# ## Output
# # TLvsprevalence: list of 2 data frames that state for each prevalence level the assigned taxonomic levels, first data frame as numbers, second as percentage 
# ## requires: 
# # dplyr!
# 
# plotTaxLevelvsprevalence <- function(taxa, seqtab){
#         
#         TLvsprevalence <- as.data.frame(cbind(colSums(seqtab != 0), apply(taxa, 2, function(x){!is.na(x)})))
#         colnames(TLvsprevalence)[1] <- "InNoSamples"
#         TLvsprevalence <- dplyr::group_by(TLvsprevalence, InNoSamples)
#         
#         if("Species" %in% colnames(taxa)){
#                 
#                 TLvsprevalence <- as.data.frame(dplyr::summarise(TLvsprevalence, NoAmplicons = n(), Kingdom = sum(Kingdom), Phylum = sum(Phylum), Class = sum(Class),
#                                                                  Order = sum(Order), Family = sum(Family), Genus = sum(Genus), Species = sum(Species)))
#                 TLPC <- dplyr::mutate(TLvsprevalence, Kingdom = Kingdom/NoAmplicons, Phylum = 100*(Phylum/NoAmplicons), Class = 100*(Class/NoAmplicons), Order =
#                                               100*(Order/NoAmplicons), Family = 100*(Family/NoAmplicons), Genus = 100*(Genus/NoAmplicons), Species = 100*(Species/NoAmplicons))
#                 
#         } else {
#                 
#                 TLvsprevalence <- as.data.frame(dplyr::summarise(TLvsprevalence, NoAmplicons = n(), Kingdom = sum(Kingdom), Phylum = sum(Phylum), Class = sum(Class),
#                                                                  Order = sum(Order), Family = sum(Family), Genus = sum(Genus)))
#                 TLPC <- dplyr::mutate(TLvsprevalence, Kingdom = Kingdom/NoAmplicons, Phylum = 100*(Phylum/NoAmplicons), Class = 100*(Class/NoAmplicons), Order =
#                                               100*(Order/NoAmplicons), Family = 100*(Family/NoAmplicons), Genus = 100*(Genus/NoAmplicons))
#                 
#         }
#         
#         TLvsprevalence <- list(Counts = TLvsprevalence, PerC = TLPC)
#         
#         return(TLvsprevalence)
#         
# }