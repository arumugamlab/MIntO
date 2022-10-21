tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
tol15rainbow=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA")
tol18rainbow=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
# ...and finally, the Paul Tol 21-color salute
tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

#######################################
#### IndivQualPlot
#######################################

IndivQualPlot <- function(QStatsList, SampleName, xlim_low = -5, xlim_high = 305) {
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        df <- QStatsList[[SampleName]]
        ggplot(data = df, aes(x = Cycle, y = Mean_QS))  + 
                geom_hline(yintercept = 25, color = 'darkred', linetype = "dashed", size = .35) +
                geom_line(color = "#009E73") + 
                geom_line(aes(y = Median_QS), color = "#E69F00") +
                #geom_line(aes(y = q25_QS), color = "#FC8D62", linetype = "dashed") + 
                #geom_line(aes(y = q75_QS), color = "#FC8D62", linetype = "dashed") +
                ylab("Quality Score") + 
                xlab("Cycle") +
                coord_cartesian(xlim=c(xlim_low, xlim_high),
                                ylim = c(8,42), expand = FALSE) +
                ## I added this
                ggtitle(paste("Sample:", SampleName, " No Reads: ", df$NoReads[1])) +
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.x = element_blank(),
                      panel.grid.major.y = element_line(color = "#999999", size = .15))
        
}

#######################################
#### QS_Median_OverviewPlot
#######################################

QS_Median_OverviewPlot <- function(QStatsList, SampleNames, Prefix = "FW", xlim_low = -5, xlim_high = 305, ylim_low=8, ylim_high=42) {
        # Input:
        # QStatsList: list of QStats data frames, such as FW_QualityStats
        # SampleNames: List of Sample names that must be the names of the QStatsList
        
        if (!all(SampleNames %in% names(QStatsList))) {
                stop("not all SampleNames were found in the QStatsList")
        }
        
        ReadLengths <- sapply(1:length(SampleNames), function(x) dim(QStatsList[[SampleNames[x]]])[1])
        x <- range(ReadLengths)/mean(ReadLengths)
        # if (!isTRUE(all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5))) {
        #         
        #         stop("Samples are of different read lengths") 
        # }
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        QList <- QStatsList[SampleNames]
        # add sample names to the data frames
        for (i in seq_along(QList)) {
                QList[[i]]$Sample <- names(QList[i])
        }
        # combine all df in the list to one big df
        df <- do.call("rbind", QList)
        
        ggplot(data = df, aes(x = Cycle, y = Median_QS, colour = NoReads, group = Sample))  + 
                geom_hline(yintercept = 25, color = 'darkred', linetype = "dashed", size = .35) +
                geom_line() +
                #geom_point() +
                ylab("Quality Score") + 
                xlab("Cycle") +
                coord_cartesian(xlim=c(xlim_low, xlim_high),
                                ylim = c(ylim_low,ylim_high), expand = FALSE) +
                scale_color_gradient2(name = "No reads", limits = c(10000, 60000), midpoint = 35000, low = "#009E73", high = "#E69F00", mid = "#999999") +
                #scale_color_gradient2(name = "No reads", limits = c(min(df$NoReads), max(df$NoReads)), midpoint = min(df$NoReads) + (range(df$NoReads)[2]- range(df$NoReads)[1])/2, low = "#009E73", high = "#E69F00", mid = "#999999")
                ggtitle(paste(Prefix, "reads, Median_QS, NoSamples:", length(SampleNames))) +
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.x = element_blank(),
                      panel.grid.major.y = element_line(color = "#999999", size = .15))
}


#####################################
QS_Median_OverviewPlot_ITS2 <- function(QStatsList, SampleNames, Prefix = "FW", xlim_low = -5, xlim_high = 305) {
  # Input:
  # QStatsList: list of QStats data frames, such as FW_QualityStats
  # SampleNames: List of Sample names that must be the names of the QStatsList
  
  if (!all(SampleNames %in% names(QStatsList))) {
    stop("not all SampleNames were found in the QStatsList")
  }
  
  ReadLengths <- sapply(1:length(SampleNames), function(x) dim(QStatsList[[SampleNames[x]]])[1])
  x <- range(ReadLengths)/mean(ReadLengths)
  # if (!isTRUE(all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5))) {
  #   
  #   stop("Samples are of different read lengths") 
  # }
  
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  QList <- QStatsList[SampleNames]
  # add sample names to the data frames
  for (i in seq_along(QList)) {
    QList[[i]]$Sample <- names(QList[i])
  }
  # combine all df in the list to one big df
  df <- do.call("rbind", QList)
  
  ggplot(data = df, aes(x = Cycle, y = Median_QS, colour = NoReads, group = Sample))  + 
    geom_hline(yintercept = 25, color = 'darkred', linetype = "dashed", size = .35) +
    geom_line() +
    #geom_point() +
    ylab("Quality Score") + 
    xlab("Cycle") +
    coord_cartesian(xlim=c(xlim_low, xlim_high),
                    ylim = c(8,42), expand = FALSE) +
    scale_color_gradient2(name = "No reads", limits = c(10000, 60000), midpoint = 35000, low = "#009E73", high = "#E69F00", mid = "#999999") +
    #scale_color_gradient2(name = "No reads", limits = c(min(df$NoReads), max(df$NoReads)), midpoint = min(df$NoReads) + (range(df$NoReads)[2]- range(df$NoReads)[1])/2, low = "#009E73", high = "#E69F00", mid = "#999999")
    ggtitle(paste(Prefix, "reads, Median_QS, NoSamples:", length(SampleNames))) +
    theme_bw() + 
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = "#999999", size = .15))
}
###
#######################################
#### QS_Mean_OverviewPlot
#######################################

QS_Mean_OverviewPlot <- function(QStatsList, SampleNames, Prefix = "FW", xlim_low = -5, xlim_high = 255) {
        # Input:
        # QStatsList: list of QStats data frames, such as FW_QualityStats
        # SampleNames: List of Sample names that must be the names of the QStatsList
        
        if (!all(SampleNames %in% names(QStatsList))) {
                stop("not all SampleNames were found in the QStatsList")
        }
        
        ReadLengths <- sapply(1:length(SampleNames), function(x) dim(QStatsList[[SampleNames[x]]])[1])
        x <- range(ReadLengths)/mean(ReadLengths)
        if (!isTRUE(all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5))) {
                
                stop("Samples are of different read lengths") 
        }
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        QList <- QStatsList[SampleNames]
        
        # add sample names to the data frames
        for (i in seq_along(QList)) {
                QList[[i]]$Sample <- names(QList[i])
        }
        
        df <- do.call("rbind", QList)
        ggplot(data = df, aes(x = Cycle, y = Mean_QS, colour = NoReads, group = Sample))  + 
                geom_hline(yintercept = 25, color = 'darkred', linetype = "dashed", size = .35) +
                geom_line() + 
                #geom_point() +
                ylab("Quality Score") + 
                xlab("Cycle") +
                coord_cartesian(xlim=c(xlim_low, xlim_high),
                                ylim = c(8,42), expand = FALSE)  +
                scale_color_gradient2(name = "No reads", limits = c(10000, 60000), midpoint = 35000, low = "#009E73", high = "#E69F00", mid = "#999999") +
                #scale_color_gradient2(name = "No reads", limits = c(min(df$NoReads), max(df$NoReads)), midpoint = min(df$NoReads) + (range(df$NoReads)[2]- range(df$NoReads)[1])/2, low = "#009E73", high = "#E69F00", mid = "#999999")
                ggtitle(paste(Prefix, " reads, Mean_QS, NoSamples", length(SampleNames))) +
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.x = element_blank(),
                      panel.grid.major.y = element_line(color = "#999999", size = .15))
        
}


#######################################
#### NoReads_DotPlot
#######################################

NoReads_DotPlot <- function(QStatsList, SampleNames) {
        # Input:
        # QStatsList: list of QStats data frames, such as FW_QualityStats
        # SampleNames: List of Sample names that must be the names of the QStatsList
        
        if (!all(SampleNames %in% names(QStatsList))) {
                stop("not all SampleNames were found in the QStatsList")
        }
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        QList <- QStatsList[SampleNames]
        # add sample names to the data frames
        for (i in seq_along(QList)) {
                QList[[i]]$Sample <- names(QList[i])
        }
        # combine all df in the list to one big df
        df <- do.call("rbind", QList)
        df <- df[!duplicated(df$Sample), c("Sample", "NoReads")]
        df <- dplyr::arrange(df, desc(NoReads))
        df$Sample <- factor(df$Sample)
        
        # relevel so samples are shown from min NoReads to max NoReads
        LevelsWant <- as.character(df$Sample) 
        for (i in 1:length(LevelsWant)) {
                df$Sample <- relevel(factor(df$Sample), ref = LevelsWant[i])
        }
        
        ggplot(data = df, aes(x = Sample, y = NoReads))  + 
                geom_hline(yintercept = 10000, color = 'darkred', linetype = "dashed", size = .35) +
                geom_point(color = "#E69F00", size = 1) +
                ylab("No Reads") + 
                xlab("") +
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
}


#######################################
#### NoReads_Histogram
#######################################

NoReads_Histogram <- function(QStatsList, SampleNames) {
        # Input:
        # QStatsList: list of QStats data frames, such as FW_QualityStats
        # SampleNames: List of Sample names that must be the names of the QStatsList
        
        if (!all(SampleNames %in% names(QStatsList))) {
                stop("not all SampleNames were found in the QStatsList")
        }
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        QList <- QStatsList[SampleNames]
        # add sample names to the data frames
        for (i in seq_along(QList)) {
                QList[[i]]$Sample <- names(QList[i])
        }
        # combine all df in the list to one big df
        df <- do.call("rbind", QList)
        df <- df[!duplicated(df$Sample), c("Sample", "NoReads")]
        
        ggplot(data = df, aes(x = NoReads))  + 
                geom_histogram(binwidth = 1000, col = "black", fill = "#E69F00")+
                geom_vline(xintercept = 10000, color = 'darkred', linetype = "dashed", size = .35) +
                geom_rug() +
                ylab("No Samples") + 
                xlab("No Reads") +
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
}

#######################################
#### NoReads_DotPlot_Type
#######################################

NoReads_DotPlot_Type <- function(QStatsList = NULL, derepFs = NULL, mergers = NULL, mergers.nochim = NULL, SampleNames, sort = TRUE) {
        
        # Input:
        # QStatsList: list of QStats data frames, such as FW_QualityStats
        # SampleNames: List of Sample names that must be names of the QStatsList, derepFs, mergers, and mergers.nochim
        # sort: if the samples should be sorted in the plot based on NoReads
        
        if (all(c(is.null(QStatsList), is.null(derepFs), is.null(mergers), is.null(mergers.nochim)))) {
                stop("QStatsList, derepFs, mergers, mergers.nochim can not all be NULL")
        }
        
        if (!all(SampleNames %in% c(names(QStatsList), names(derepFs), names(mergers), names(mergers.nochim)))) {
                stop("not all SampleNames were found in the given data")
        }
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        ## Construct the plotting data frame
        
        df.all <- NULL
        
        if (!is.null(QStatsList) & any(SampleNames %in% names(QStatsList))) {
                
                QStatsList <- QStatsList[SampleNames[which(SampleNames %in% names(QStatsList))]]
                
                for (i in seq_along(QStatsList)) {
                        QStatsList[[i]]$Sample <- names(QStatsList[i])
                }
                
                df.all <- do.call("rbind",QStatsList)
                df.all <- df.all[!duplicated(df.all$Sample), c("Sample", "NoReads")]
                df.all$Type <- "all"
                
        }
        
        df.filtered <- NULL
        
        if (!is.null(derepFs) & any(SampleNames %in% names(derepFs))) {
                
                derepFs <- derepFs[SampleNames[which(SampleNames %in% names(derepFs))]]
                
                df.filtered <- data.frame(Sample = names(derepFs), NoReads = 0, Type = "filtered")
                for (i in seq_along(derepFs)) {
                        df.filtered$NoReads[i] <- sum(derepFs[[i]]$uniques) 
                }
                
        }
        
        df.merged <- NULL
        
        if (!is.null(mergers) & any(SampleNames %in% names(mergers))) {
                
                mergers <- mergers[SampleNames[which(SampleNames %in% names(mergers))]]
                
                df.merged <- data.frame(Sample = names(mergers), NoReads = 0, Type = "merged")
                for(i in seq_along(mergers)) {
                        df.merged$NoReads[i] <- sum(mergers[[i]]$abundance)
                }
                
        } 
        
        df.nochim <- NULL
        
        if (!is.null(mergers.nochim) & any(SampleNames %in% names(mergers.nochim))) {
                
                mergers.nochim <- mergers.nochim[SampleNames[which(SampleNames %in% names(mergers.nochim))]]
                
                df.nochim <- data.frame(Sample = names(mergers.nochim), NoReads = 0, Type = "nochim")
                for(i in seq_along(mergers.nochim)) {
                        df.nochim$NoReads[i] <- sum(mergers.nochim[[i]]$abundance)
                }
                
        }
        
        df.plot <- rbind(df.all, df.filtered, df.merged, df.nochim)
        
        if (sort) {
                # sort always by highest level: all, filtered, merged, nochim
                if (!is.null(df.all)) {
                        df.all <- dplyr::arrange(df.all, desc(NoReads))
                        LevelsWant <- as.character(df.all$Sample)
                } else if (!is.null(df.filtered)) {
                        df.filtered <- dplyr::arrange(df.filtered, desc(NoReads))
                        LevelsWant <- as.character(df.filtered$Sample)
                } else if (!is.null(df.merged)) {
                        df.merged <- dplyr::arrange(df.merged, desc(NoReads))
                        LevelsWant <- as.character(df.merged$Sample)
                } else if (!is.null(df.nochim)) {
                        df.nochim <- dplyr::arrange(df.nochim, desc(NoReads))
                        LevelsWant <- as.character(df.nochim$Sample)
                }
                
                
                df.plot$Sample <- factor(df.plot$Sample)
                for (i in seq_along(LevelsWant)) {
                        df.plot$Sample <- relevel(factor(df.plot$Sample), ref = LevelsWant[i])
                }
                
        }
        
        
        ggplot(data = df.plot, aes(x = Sample, y = NoReads, color = Type))  + 
                geom_hline(yintercept = 10000, color = 'darkred', linetype = "dashed", size = .35) +
                geom_point(size = 2) +
                scale_color_manual(values = c(cbPalette[2:4],cbPalette[7])) +
                ylab("No Reads") + 
                xlab("") +
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                      legend.title = element_blank())
        
}



#######################################
#### NoReads_StepsSimple
#######################################

NoReads_StepsSimple <- function(ReadSummary, SampleNames, sort = TRUE) {
        
        
        if (!all(SampleNames %in% ReadSummary$Sample)) {
                stop("not all SampleNames were found in the given ReadSummary")
        }
        
        ReadSummary <- ReadSummary[ReadSummary$Sample %in% SampleNames,]
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
        
        ## Construct the plotting data frame
        
        ReadSummary <- subset(ReadSummary, select = c("Sample", "NoReads", "FilteredReads", "MergedReads", "MergedReadsWOBimera"))
        
        
        if (sort) {
                
                ReadSummary <- dplyr::arrange(ReadSummary, desc(NoReads))
                LevelsWant <- as.character(ReadSummary$Sample)
                for (i in seq_along(LevelsWant)) {
                        ReadSummary$Sample <- relevel(factor(ReadSummary$Sample), ref = LevelsWant[i])
                }
                
        }
        
        colnames(ReadSummary)[2:5] <- c("all", "filtered", "merged", "nochim")
        ReadSummary <- tidyr::gather(ReadSummary, key = Type, value = NoReads, -Sample)
        
        Tr <- ggplot(data = ReadSummary, aes(x = Sample, y = NoReads, color = Type))  
        Tr <- Tr + 
                geom_hline(yintercept = 10000, color = 'darkred', linetype = "dashed", size = .35) +
                geom_point(size = 2) +
                scale_color_manual(values = c(cbPalette[2:4],cbPalette[7])) +
                ylab("No Reads") + 
                xlab("") +
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                      legend.title = element_blank())
        Tr
        
}



#######################################
#### TotalandUniqueAmplicons
#######################################
# Input:
# FinalNumers: data frame with Sample, UniqueAmplicons and NoAmplicons columns
# sort: default TRUE: sort samples after number total amplicons
# Output:
# the Tr structure
# NB: requires tidyr


TotalandUniqueAmplicons <- function(FinalNumbers, seqtab, sort =TRUE) {
        
        ## show the final amplicon numbers in a plot
        FinalNumbers$TotalAmplicons <- FinalNumbers$NoAmplicons/100 
        
        FinalNumbersL <- tidyr::gather(FinalNumbers, key = Type, value = Number, -Sample, -NoAmplicons)
        
        if(sort) {
                
                FinalNumbersL <- dplyr::arrange(FinalNumbersL, desc(NoAmplicons))
                LevelsWant <- as.character(FinalNumbersL$Sample)
                for (i in seq_along(LevelsWant)) {
                        FinalNumbersL$Sample <- relevel(factor(FinalNumbersL$Sample), ref = LevelsWant[i])
                }
                
        }
        
        Tr <- ggplot(FinalNumbersL, aes(x = Sample, y = Number, col = Type)) +
                geom_point() +
                scale_color_manual(values = c("#E69F00", "#009E73"), labels = c("TotalAmplicons/100", "UniqueAmplicons"))+
                ylab("No Amplicons") + 
                xlab("") +
                #scale_y_continuous(breaks = c(100, 200, 300, 400), labels = c("100\n(10000)", "200\n(20000)", "300\n(30000)", "400\n(40000)")) +
                ggtitle(paste("No unique amplicons in all samples",dim(seqtab)[2]))+
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                      legend.title = element_blank())
        
        return(Tr)
        
}


#######################################
#### plotSVdistributions
#######################################
## SORRY FOR THE TERMINOLOGY CONFUSION< from now on counts is used instead of amplicons in the plots
## REQUIRES dplyr
## Input:
# Seqtab: seqtab from dada2 wrapper
# prevalence: a percentage of samples (bw 0 and 100)
## Output: 
# TrList: List of four Trelis objects, plus distribution data.frame


plotSVdistributions <- function(seqtab, prevalence = 10) {
        
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
                ylab("number of ASVs") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        
        In1Index <- which(CountDistribution$InNumberSamples == 1)
        if (length(In1Index) != 0) {
                Tr <- Tr + ggtitle(paste(CountDistribution$No_ASVs[In1Index], " of ", CountDistribution$CumSumUnique[1], " ASVs (", round(100*CountDistribution$No_ASVs[In1Index]/CountDistribution$CumSumUnique[1], 1), " %)", " were only found in 1 sample", sep = ""))
        } 
        
        
        # Cumulative Percentage of SVs
        Tr1 <- ggplot(CountDistribution, aes(x = InNumberSamples, y = CumPerCUnique))
        Tr1 <- Tr1 + geom_point(col = "#E69F00", size = 3) +
                xlab("prevalence") +
                ylab("cumulative percentage of ASVs") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        Tr1 <- Tr1 + 
                geom_hline(yintercept = SVskeptAtPCValue, lty =  "dashed") +
                geom_vline(xintercept = (prevalence/100)*dim(seqtab)[1], lty = 'dashed') +
                ggtitle(paste("prevalence ", prevalence, " % = ", round((prevalence/100)*dim(seqtab)[1],1), "; ", CountDistribution$CumSumUnique[index], " of ", CountDistribution$CumSumUnique[1], 
                              " ASVs (", round(100*CountDistribution$CumSumUnique[index]/CountDistribution$CumSumUnique[1], 1), " %) have higher prevalence", sep = ""))
        
        
        Tr2 <- ggplot(CountDistribution, aes(x = InNumberSamples, y = TotalCounts))
        Tr2 <- Tr2 + geom_point(col = "#E69F00", size = 3) +
                xlab("prevalence") +
                ylab("total counts of ASVs with given prevalence") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        Tr2 <- Tr2 + ggtitle(paste(CountDistribution$CumSumTotal[1] - CountDistribution$CumSumTotal[index], " of ",
                                   CountDistribution$CumSumTotal[1], " (", round(100*(CountDistribution$CumSumTotal[1] - CountDistribution$CumSumTotal[index])/CountDistribution$CumSumTotal[1], 2),
                                    " %) counts are from ASVs present in less than ", round((prevalence/100)*dim(seqtab)[1],1), " samples.", sep = ""))
        
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
        
        return(TrList)
        
}

#######################################
#### plotAlphaDiversity
#######################################
# plotAlphaDiversity is related to phyloseq::plot_richness but also offers boxplots (when group != NULL) and
# always renders each plot individually in a list instead of using facet_wraps

plotAlphaDiversity <- function(physeq, measures = NULL, x = "samples", color = NULL, shape = NULL,
                               group = NULL, defCol = "#E69F00"){
        
        
        AlphaDiv <- suppressWarnings(estimate_richness(physeq, measures = measures))
        measures = colnames(AlphaDiv)
        ses = colnames(AlphaDiv)[grep("^se\\.", colnames(AlphaDiv))]
        measures = measures[!measures %in% ses]
        if("Observed" %in% measures){
                
                measures[measures == "Observed"] <- "Richness"
                colnames(AlphaDiv)[colnames(AlphaDiv) == "Observed"] <- "Richness"
        }
        
        if (!is.null(sample_data(physeq, errorIfNULL = FALSE))) {
                DF <- data.frame(AlphaDiv, sample_data(physeq))
        } else {
                DF <- data.frame(AlphaDiv)
        }
        
        DF$samples <- sample_names(physeq)
        
        if (!is.null(x)) {
                if (x %in% c("sample", "samples", "sample_names", "sample.names", "Sample", "Samples")) {
                        x <- "samples"
                }
        } else {
                x <- "samples"
        }
        
        if(!is.null(group)){

                
                TrList <- list()
                
                for (i in 1:length(measures)) {
                        
                        aes_map = aes_string(x = group, y = measures[i], shape = shape, color = color, group = group)
                        
                        if(is.null(color)){
                                Tr = ggplot(DF, aes_map) + geom_boxplot(na.rm = TRUE, col = defCol)
                        } else {
                                
                                Tr = ggplot(DF, aes_map) + geom_boxplot(na.rm = TRUE)
                        }
                        
                        Tr = Tr + geom_jitter(width = .2, height = 0, alpha = 0.65)
                        
                        Tr <- Tr + scale_colour_manual(values = cbPalette[2:11])
                        
                        Tr <- Tr + theme_bw() + xlab("")
                        
                        Tr = Tr + theme(panel.grid.minor = element_blank(),
                                        panel.grid.major.y = element_blank(),
                                        panel.grid.major.x = element_line(color = "#999999", size = .15))
                        
                        TrList[[i]] <- Tr
                        
                }
                
                

        } else {
                
                TrList <- list()
                
                for (i in 1:length(measures)) {
                        
                        aes_map = aes_string(x = x, y = measures[i], shape = shape, color = color)
                        
                        if(is.null(color)){
                                Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE, col = defCol)
                        } else {
                                
                                Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE)
                        }
                        
                        
                        if (measures[i] == "Chao1") {
                                Tr = Tr + geom_errorbar(aes(ymax = Chao1 + se.chao1, ymin = Chao1 -
                                                                    se.chao1), width = 0.1)
                        }
                        if (measures[i] == "ACE") {
                                Tr = Tr + geom_errorbar(aes(ymax = ACE + se.ACE, ymin = ACE -
                                                                    se.ACE), width = 0.1)
                        }
                        
                        Tr <- Tr + scale_colour_manual(values = cbPalette[2:11])
                        
                        Tr <- Tr + theme_bw()
                        
                        Tr = Tr + theme(panel.grid.minor = element_blank(),
                                        panel.grid.major.y = element_blank(),
                                        panel.grid.major.x = element_line(color = "#999999", size = .15),
                                        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
                        
                        TrList[[i]] <- Tr
                        
                }
                
                
                
        }
        
        names(TrList) <- measures
        return(TrList)
        
}


#######################################
#### boxplot_alphaDiv_fromDF
#######################################
# see plotAlphaDiversity

boxplot_alphaDiv_fromDF <- function(DF, measures, color = NULL, shape = NULL,
                                    group = NULL, defCol = "#E69F00") {
        
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
        TrList <- list()
        
        if (!is.null(DF$Total) && !is.null(DF$filtered_reads)) {
                measures <- c(measures, "Total", "filtered_reads")
        }
        
        for (i in 1:length(measures)) {
                
                aes_map = aes_string(x = group, y = measures[i], color = color)
                
                if(is.null(color)){
                        Tr = ggplot(DF, aes_map) + geom_boxplot(na.rm = TRUE, col = defCol)
                } else {
                        
                        Tr = ggplot(DF, aes_map) + geom_boxplot(na.rm = TRUE)
                }
                
                Tr = Tr + geom_jitter(aes_string(shape = shape), width = .2, height = 0, alpha = 0.65)
                
                Tr <- Tr + scale_colour_manual("", values = cbPalette[2:11])
                
                Tr <- Tr + theme_bw() + xlab("") +theme(legend.position = "none")
                
                Tr = Tr + theme(panel.grid.minor = element_blank(),
                                panel.grid.major.y = element_blank(),
                                panel.grid.major.x = element_line(color = "#999999", size = .15))
                
                TrList[[i]] <- Tr
                
        }
        
        names(TrList) <- measures
        return(TrList)
}



#######################################
#### plot_alphaDivVstotalCounts
#######################################
# Function plots the alpha diversity measures from estimate_richness against the total number of reads/amplicons per sample
# and adds a linear fit.
# if color is given the fit is individual on the different colors, but the p-value in the title is still to an overall fit!!

plot_alphaDivVstotalCounts <- function(physeq, measures = NULL, color = NULL, shape = NULL, defCol = "#E69F00"){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
        # needed to get the p-value from-linear fit objects (from stackoverflow)
        lmp <- function (modelobject) {
                if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
                f <- summary(modelobject)$fstatistic
                p <- pf(f[1],f[2],f[3],lower.tail=F)
                attributes(p) <- NULL
                return(p)
        }

        AlphaDiv <- suppressWarnings(estimate_richness(physeq, measures = measures))
        measures = colnames(AlphaDiv)
        ses = colnames(AlphaDiv)[grep("^se\\.", colnames(AlphaDiv))]
        measures = measures[!measures %in% ses]
        
        if("Observed" %in% measures){
                
                measures[measures == "Observed"] <- "Richness"
                colnames(AlphaDiv)[colnames(AlphaDiv) == "Observed"] <- "Richness"
        }
        
        if (!is.null(sample_data(physeq, errorIfNULL = FALSE))) {
                DF <- data.frame(AlphaDiv, sample_data(physeq))
        } else {
                DF <- data.frame(AlphaDiv)
        }
        
        DF$samples <- sample_names(physeq)
        
        if(taxa_are_rows(physeq)){
                DF$TotalReads <- colSums(otu_table(physeq))
                
        } else {
                DF$TotalReads <- rowSums(otu_table(physeq))
        }
        
        
        TrList <- list()
        
        for (i in 1:length(measures)) {
                
                aes_map = aes_string(x = "TotalReads", y = measures[i], shape = shape, color = color)
                
                if(is.null(color)){
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE, col = defCol)
                } else {
                        
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE)
                }
                
                
                if (measures[i] == "Chao1") {
                        Tr = Tr + geom_errorbar(aes(ymax = Chao1 + se.chao1, ymin = Chao1 -
                                                            se.chao1), width = 0.1)
                }
                if (measures[i] == "ACE") {
                        Tr = Tr + geom_errorbar(aes(ymax = ACE + se.ACE, ymin = ACE -
                                                            se.ACE), width = 0.1)
                }
                
                Tr <- Tr + scale_colour_manual(values = cbPalette[2:11])
                
                # add the regression line and the p-values
                fit <- lm(DF[,measures[i]] ~ DF[,"TotalReads"])
                adjR2 <- round(summary(fit)$adj.r.squared,3)
                pfit <- lmp(fit)
                
                Tr <- Tr + geom_smooth(aes_string(x = "TotalReads", y = measures[i]), method = "lm", se = TRUE, inherit.aes = F)
                
                Tr <- Tr + ggtitle(paste("lm fit: p-value: ", format(pfit, digits = 4), ", adj.r.squared: ", adjR2, sep = ""))
                
                Tr <- Tr + theme_bw() + xlab("total counts in sample (sample_sums())")
                
                Tr = Tr + theme(panel.grid.minor = element_blank())
                
                TrList[[i]] <- Tr
                
        }
  
        names(TrList) <- measures
        return(TrList)
        
}    


#######################################
#### plot_alphaDivVstotalCounts
#######################################
# Function plots the alpha diversity measures from estimate_richness against the total number of reads/amplicons per sample
# and adds a linear fit.
# if color is given the fit is individual on the different colors, but the p-value in the title is still to an overall fit!!

plot_alphaDivVstotalCounts_fromList <- function(DF_List, measures = NULL, color = NULL, shape = NULL, defCol = "#E69F00"){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
        TrList <- list()
        
        DF <- DF_List[[1]]
        fitlist <- DF_List[["fitlist"]]
        
        for (i in 1:length(measures)) {
                
                aes_map = aes_string(x = "Total", y = measures[i], shape = shape, color = color)
                
                if(is.null(color)){
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE, col = defCol)
                } else {
                        
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE)
                }
                
                
                if (measures[i] == "Chao1") {
                        Tr = Tr + geom_errorbar(aes(ymax = Chao1 + se.chao1, ymin = Chao1 -
                                                            se.chao1), width = 0.1)
                }
                if (measures[i] == "ACE") {
                        Tr = Tr + geom_errorbar(aes(ymax = ACE + se.ACE, ymin = ACE -
                                                            se.ACE), width = 0.1)
                }
                
                Tr <- Tr + scale_colour_manual("", values = cbPalette[2:11])
                
                # add the regression line and the p-values
                fit <- fitlist[[measures[i]]]
                adjR2 <- round(summary(fit)$adj.r.squared,3)
                if (adjR2 != 0){
                        pfit <- lmp(fit)
                } else {
                        pfit <- NA
                }
                
                Tr <- Tr + geom_smooth(aes_string(x = "Total", y = measures[i]), method = "lm", se = TRUE, inherit.aes = F)
                
                Tr <- Tr + ggtitle(paste("lm fit: p-value: ", format(pfit, digits = 4), ", adj.r.squared: ", adjR2, sep = ""))
                
                Tr <- Tr + theme_bw() + xlab("total counts in sample (sample_sums())")
                
                Tr = Tr + theme(panel.grid.minor = element_blank())
                
                TrList[[i]] <- Tr
                
        }
        
        names(TrList) <- measures
        return(TrList)
        
}

#######################################
#### plot_alphaDivVsfilteredReads_fromList
#######################################
# Function plots the alpha diversity measures from estimate_richness against the given number of filteredReads per sample
# and adds a linear fit.
# if color is given the fit is individual on the different colors, but the p-value in the title is still to an overall fit!!

plot_alphaDivVsfilteredReads_fromList <- function(DF_List, measures = NULL, color = NULL, shape = NULL, defCol = "#E69F00"){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
        TrList <- list()
        
        DF <- DF_List[[1]]
        fitlist <- DF_List[["fitlist_FilteredReads"]]
        
        for (i in 1:length(measures)) {
                
                aes_map = aes_string(x = "filtered_reads", y = measures[i], shape = shape, color = color)
                
                if(is.null(color)){
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE, col = defCol)
                } else {
                        
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE)
                }
                
                
                if (measures[i] == "Chao1") {
                        Tr = Tr + geom_errorbar(aes(ymax = Chao1 + se.chao1, ymin = Chao1 -
                                                            se.chao1), width = 0.1)
                }
                if (measures[i] == "ACE") {
                        Tr = Tr + geom_errorbar(aes(ymax = ACE + se.ACE, ymin = ACE -
                                                            se.ACE), width = 0.1)
                }
                
                Tr <- Tr + scale_colour_manual("", values = cbPalette[2:11])
                
                # add the regression line and the p-values
                fit <- fitlist[[measures[i]]]
                adjR2 <- round(summary(fit)$adj.r.squared,3)
                if (adjR2 != 0){
                        pfit <- lmp(fit)
                } else {
                        pfit <- NA
                }
                
                Tr <- Tr + geom_smooth(aes_string(x = "filtered_reads", y = measures[i]), method = "lm", se = TRUE, inherit.aes = F)
                
                Tr <- Tr + ggtitle(paste("lm fit: p-value: ", format(pfit, digits = 4), ", adj.r.squared: ", adjR2, sep = ""))
                
                Tr <- Tr + theme_bw() + xlab("number of filtered reads that entered dada command")
                
                Tr = Tr + theme(panel.grid.minor = element_blank())
                
                TrList[[i]] <- Tr
                
        }
        
        names(TrList) <- measures
        return(TrList)
        
}

#######################################
#### plot_alphaDivVSoriginalTotalSampleCounts
#######################################
# wanted to add this plot after rarefaction

plot_alphaDivVSoriginalTotalSampleCounts <- function(DF_alpha_rare, DF_alpha_no_rare, measures = NULL, color = NULL, shape = NULL, defCol = "#E69F00"){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
        measures2 <- measures
        if ("Observed" %in% measures2) {
                measures2[measures2 == "Observed"] <- "Richness" 
        }
        
        DF_alpha_rare$originalTotal <- DF_alpha_no_rare$Total
        DF <- DF_alpha_rare
        
        TrList <- list()
        
        for (i in 1:length(measures2)) {
                
                aes_map = aes_string(x = "originalTotal", y = measures2[i], shape = shape, color = color)
                
                if(is.null(color)){
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE, col = defCol)
                } else {
                        
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE)
                }
                
                
                if (measures[i] == "Chao1") {
                        Tr = Tr + geom_errorbar(aes(ymax = Chao1 + se.chao1, ymin = Chao1 -
                                                            se.chao1), width = 0.1)
                }
                if (measures[i] == "ACE") {
                        Tr = Tr + geom_errorbar(aes(ymax = ACE + se.ACE, ymin = ACE -
                                                            se.ACE), width = 0.1)
                }
                
                Tr <- Tr + scale_colour_manual("", values = cbPalette[2:11])
                
                # add the regression line and the p-values
                fit <- lm(DF[, measures2[i]] ~ DF[,"originalTotal"])
                adjR2 <- round(summary(fit)$adj.r.squared,3)
                if (adjR2 != 0){
                        pfit <- lmp(fit)
                } else {
                        pfit <- NA
                }
                
                Tr <- Tr + geom_smooth(aes_string(x = "originalTotal", y = measures2[i]), method = "lm", se = TRUE, inherit.aes = F)
                
                Tr <- Tr + ggtitle(paste("lm fit: p-value: ", format(pfit, digits = 4), ", adj.r.squared: ", adjR2, sep = ""))
                
                Tr <- Tr + theme_bw() + xlab("total counts in sample (sample_sums()) before rarefying")
                
                Tr = Tr + theme(panel.grid.minor = element_blank())
                
                TrList[[i]] <- Tr
                
        }
        
        names(TrList) <- measures2
        return(TrList)
        
}



#######################################
#### plot_alphaDivVSoriginalTotalSampleCounts2
#######################################
# wanted to add this plot after rarefaction

plot_alphaDivVSoriginalTotalSampleCounts2 <- function(DF_alpha_rare, DF_alpha_no_rare, measures = NULL, color = NULL, defCol = "#E69F00"){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
        measures2 <- measures
        if ("Observed" %in% measures2) {
                measures2[measures2 == "Observed"] <- "Richness" 
        }
        
        DF_alpha_no_rare <- DF_alpha_no_rare[, c("Sample", measures2, "Total", "filtered_reads", color)]
        DF_alpha_no_rare$Type <- "Before rarefaction"
        DF_alpha_rare <- DF_alpha_rare[, c("Sample", measures2, "Total", "filtered_reads", color)]
        to_level <- DF_alpha_rare$Total[1]
        DF_alpha_rare$Type <- paste("After rarefaction to ", to_level, " counts", sep = "")
        DF_alpha_rare$Total <- DF_alpha_no_rare$Total
        
        DF <- rbind(DF_alpha_no_rare, DF_alpha_rare)
        DF$Type <- factor(DF$Type, levels = c("Before rarefaction", paste("After rarefaction to ", to_level, " counts", sep = "")),
                          ordered = TRUE)
        
        TrList <- list()
        
        for (i in 1:length(measures2)) {
                
                aes_map = aes_string(x = "Total", y = measures2[i], color = color)
                
                if(is.null(color)){
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE, col = defCol)
                } else {
                        
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE)
                }
                
                
                Tr <- Tr + scale_colour_manual("", values = cbPalette[2:11])
                
                # for the regressin lines and the p-values
                fit <- lm(DF_alpha_no_rare[, measures2[i]] ~ DF_alpha_no_rare[,"Total"])
                adjR2 <- round(summary(fit)$adj.r.squared,3)
                if (adjR2 != 0){
                        pfit <- lmp(fit)
                } else {
                        pfit <- NA
                }
                
                fit_rare <- lm(DF_alpha_rare[, measures2[i]] ~ DF_alpha_rare[,"Total"])
                adjR2_rare <- round(summary(fit_rare)$adj.r.squared,3)
                if (adjR2_rare != 0){
                        pfit_rare <- lmp(fit)
                } else {
                        pfit_rare <- NA
                }
                
                
                
                Tr <- Tr + geom_smooth(aes_string(x = "Total", y = measures2[i]), method = "lm", se = TRUE, inherit.aes = F)
                
                Tr <- Tr + facet_grid(~ Type, scales = "free_y")
                
                Tr <- Tr + ggtitle(paste("fit before: p: ", format(pfit, digits = 4), ", R^2: ", adjR2, "; after: p: ",
                                         format(pfit_rare, digits = 4), ", R^2: ", adjR2_rare, sep = ""))
                
                Tr <- Tr + theme_bw() + xlab("total counts in sample (sample_sums()) before rarefying")
                
                Tr = Tr + theme(panel.grid.minor = element_blank())
                
                TrList[[i]] <- Tr
                
        }
        
        names(TrList) <- measures2
        
        TrList2 <- list()
        
        for (i in 1:length(measures2)) {
                
                aes_map = aes_string(x = "filtered_reads", y = measures2[i], color = color)
                
                if(is.null(color)){
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE, col = defCol)
                } else {
                        
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE)
                }
                
                
                Tr <- Tr + scale_colour_manual("", values = cbPalette[2:11])
                
                # for the regressin lines and the p-values
                fit <- lm(DF_alpha_no_rare[, measures2[i]] ~ DF_alpha_no_rare[,"filtered_reads"])
                adjR2 <- round(summary(fit)$adj.r.squared,3)
                if (adjR2 != 0){
                        pfit <- lmp(fit)
                } else {
                        pfit <- NA
                }
                
                fit_rare <- lm(DF_alpha_rare[, measures2[i]] ~ DF_alpha_rare[,"filtered_reads"])
                adjR2_rare <- round(summary(fit_rare)$adj.r.squared,3)
                if (adjR2_rare != 0){
                        pfit_rare <- lmp(fit)
                } else {
                        pfit_rare <- NA
                }
                
                
                
                Tr <- Tr + geom_smooth(aes_string(x = "filtered_reads", y = measures2[i]), method = "lm", se = TRUE, inherit.aes = F)
                
                Tr <- Tr + facet_grid(~ Type, scales = "free_y")
                
                Tr <- Tr + ggtitle(paste("fit before: p: ", format(pfit, digits = 4), ", R^2: ", adjR2, "; after: p: ",
                                         format(pfit_rare, digits = 4), ", R^2: ", adjR2_rare, sep = ""))
                
                Tr <- Tr + theme_bw() + xlab("No of filtered reads (that entered dada() command)")
                
                Tr = Tr + theme(panel.grid.minor = element_blank())
                
                TrList2[[i]] <- Tr
                
        }
        
        names(TrList2) <- measures2
        
        TrList3 <- list()
        
        for (i in 1:length(measures2)) {
                
                aes_map = aes_string(x = "Total", y = measures2[i], color = color, shape = "Type")
                
                if(is.null(color)){
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE, col = defCol)
                } else {
                        
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE)
                }
                
                
                Tr <- Tr + scale_colour_manual("", values = cbPalette[2:11])
                
                
                # for the regressin lines and the p-values
                fit <- lm(DF_alpha_no_rare[, measures2[i]] ~ DF_alpha_no_rare[,"Total"])
                adjR2 <- round(summary(fit)$adj.r.squared,3)
                if (adjR2 != 0){
                        pfit <- lmp(fit)
                } else {
                        pfit <- NA
                }
                
                fit_rare <- lm(DF_alpha_rare[, measures2[i]] ~ DF_alpha_rare[,"Total"])
                adjR2_rare <- round(summary(fit_rare)$adj.r.squared,3)
                if (adjR2_rare != 0){
                        pfit_rare <- lmp(fit)
                } else {
                        pfit_rare <- NA
                }
                
                
                
                Tr <- Tr + geom_smooth(aes_string(x = "Total", y = measures2[i], fill = "Type"), method = "lm", se = FALSE, inherit.aes = F)
                
                
                Tr <- Tr + ggtitle(paste("fit before: p: ", format(pfit, digits = 4), ", R^2: ", adjR2, "; after: p: ",
                                         format(pfit_rare, digits = 4), ", R^2: ", adjR2_rare, sep = ""))
                
                Tr <- Tr + theme_bw() + xlab("total counts in sample (sample_sums()) before rarefying")
                
                Tr = Tr + theme(panel.grid.minor = element_blank(),
                                legend.title = element_blank(),
                                legend.position = "top")
                
                TrList3[[i]] <- Tr
                
        }
        
        names(TrList3) <- measures2
        
        TrList4 <- list()
        
        for (i in 1:length(measures2)) {
                
                aes_map = aes_string(x = "filtered_reads", y = measures2[i], color = color, shape = "Type")
                
                if(is.null(color)){
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE, col = defCol)
                } else {
                        
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE)
                }
                
                
                Tr <- Tr + scale_colour_manual("", values = cbPalette[2:11])
                
                
                # for the regressin lines and the p-values
                fit <- lm(DF_alpha_no_rare[, measures2[i]] ~ DF_alpha_no_rare[,"filtered_reads"])
                adjR2 <- round(summary(fit)$adj.r.squared,3)
                if (adjR2 != 0){
                        pfit <- lmp(fit)
                } else {
                        pfit <- NA
                }
                
                fit_rare <- lm(DF_alpha_rare[, measures2[i]] ~ DF_alpha_rare[,"filtered_reads"])
                adjR2_rare <- round(summary(fit_rare)$adj.r.squared,3)
                if (adjR2_rare != 0){
                        pfit_rare <- lmp(fit)
                } else {
                        pfit_rare <- NA
                }
                
                
                
                Tr <- Tr + geom_smooth(aes_string(x = "filtered_reads", y = measures2[i], fill = "Type"), method = "lm", se = FALSE, inherit.aes = F)
                
                
                Tr <- Tr + ggtitle(paste("fit before: p: ", format(pfit, digits = 4), ", R^2: ", adjR2, "; after: p: ",
                                         format(pfit_rare, digits = 4), ", R^2: ", adjR2_rare, sep = ""))
                
                Tr <- Tr + theme_bw() + xlab("No of filtered reads (that entered dada() command)")
                
                Tr = Tr + theme(panel.grid.minor = element_blank(),
                                legend.title = element_blank(),
                                legend.position = "top")
                
                TrList4[[i]] <- Tr
                
        }
        
        names(TrList4) <- measures2
        
        list(TrList_total = TrList, TrList_filtered_reads = TrList2, 
             TrList_total_one = TrList3, TrList_filtered_reads_one = TrList4)
}


#######################################
#### plot_correlations_abundance_prev_sparsity
#######################################
# df_ab_prev: data frame with "ASV_ID", "total_counts_of_ASV", "prevalence", "sparsity", "mean_count_nonzero",
# "median_count_nonzero"
# outputs a list with different plots and fits

plot_correlations_abundance_prev_sparsity <- function(df_ab_prev, col = NULL){ # NB: you could color by Phylum for example
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
        nsamples <- df_ab_prev$prevalence[1] + df_ab_prev$sparsity[1]
        
        Tr_ab <- ggplot(df_ab_prev, aes(x = ASV_ID, y = total_counts_of_ASV))
        if (is.null(col)) {
                Tr_ab <- Tr_ab + geom_point(col = cbPalette[2], size = 2, alpha = 0.7) 
        } else {
                Tr_ab <- Tr_ab + geom_point(aes_string(col = col), size = 2, alpha = 0.7)
        }
        Tr_ab <- Tr_ab +
                ylab("total counts of ASV (= taxa_sums())") +
                theme_bw(12)
        
        
        Tr_prev <- ggplot(df_ab_prev, aes(x = ASV_ID, y = prevalence))
        if (is.null(col)) {
                Tr_prev <- Tr_prev + geom_point(col = cbPalette[2], size = 2, alpha = 0.7) 
        } else {
                Tr_prev <- Tr_prev + geom_point(aes_string(col = col), size = 2, alpha = 0.7)
        }
        Tr_prev <- Tr_prev +
                theme_bw(12)
        
        
        # - associations of total abundance and log10(total_counts_of_ASV) to sparsity/prevalence -
        # NB: turned out: log10(total_counts_of_ASV) is far better correlated with prevalence/sparsity, so I stick with the log fits
        # Also NB: since prevalence = constant - sparsity, a fit prevalence ~ abundance is equal to fit sparsity ~ abundance just with
        # opposite coefficient signs
        
        # fit_spar <- lm(formula = sparsity ~ total_counts_of_ASV, data = df_ab_prev)
        # pval_spar <- lmp(fit_spar)
        # fit_prev <- lm(formula = prevalence ~ total_counts_of_ASV, data = df_ab_prev)
        # pval_prev <- lmp(fit_prev)
        fit_prev_log10 <- lm(formula = prevalence ~ log10(total_counts_of_ASV), data = df_ab_prev)
        pval_prev_log10 <- lmp(fit_prev_log10)
        # fit_spar_log10 <- lm(formula = sparsity ~ log10(total_counts_of_ASV), data = df_ab_prev)
        # pval_spar_log10 <- lmp(fit_spar_log10)
        # identical(pval_prev_log10, pval_spar_log10) # TRUE
        
        # it comes natural that total abundance and prevalence/sparsity are correlated, but how about mean abundance in non zero samples
        # here I prepare all combinations, but stick to prevalence for now (comment sparsity out) since correlations are the same
        
        
        # fit_spar_mean <- lm(formula = sparsity ~ mean_count_nonzero, data = df_ab_prev)
        # pval_spar_mean <- lmp(fit_spar_mean)
        # fit_spar_mean_log10 <- lm(formula = sparsity ~ log10(mean_count_nonzero), data = df_ab_prev)
        # pval_spar_mean_log10 <- lmp(fit_spar_mean_log10)
        # fit_spar_median <- lm(formula = sparsity ~ median_count_nonzero, data = df_ab_prev)
        # pval_spar_median <- lmp(fit_spar_median)
        # fit_spar_median_log10 <- lm(formula = sparsity ~ log10(median_count_nonzero), data = df_ab_prev)
        # pval_spar_median_log10 <- lmp(fit_spar_median_log10)
        
        fit_prev_mean <- lm(formula = prevalence ~ mean_count_nonzero, data = df_ab_prev)
        pval_prev_mean <- lmp(fit_prev_mean)
        fit_prev_mean_log10 <- lm(formula = prevalence ~ log10(mean_count_nonzero), data = df_ab_prev)
        pval_prev_mean_log10 <- lmp(fit_prev_mean_log10)
        fit_prev_median <- lm(formula = prevalence ~ median_count_nonzero, data = df_ab_prev)
        pval_prev_median <- lmp(fit_prev_median)
        fit_prev_median_log10 <- lm(formula = prevalence ~ log10(median_count_nonzero), data = df_ab_prev)
        pval_prev_median_log10 <- lmp(fit_prev_median_log10)
        
        
        # Tr_ab_vs_prev <- ggplot(df_ab_prev, aes(x = total_counts_of_ASV, y = prevalence))
        # Tr_ab_vs_prev <- Tr_ab_vs_prev +
        #         geom_point(col = cbPalette[2], size = 2, alpha = 0.7) +
        #         geom_smooth(method = "lm") +
        #         scale_y_continuous(limits = c(-1, max(df_ab_prev$prevalence) + 5)) +
        #         # annotate("text", x = max(df_ab_prev$abundance)/5, y = 5, label = paste("R.square of lm: ", as.character(round(summary(fit_prev)$r.squared,4), sep = ""))) +
        #         xlab("total abundance (taxa_sums())") +
        #         ggtitle(paste("lm fit: p.val: ", format(pval_prev, digits = 4), "R.square: ", as.character(round(summary(fit_prev)$r.squared,4), sep = ""))) +
        # theme_bw(12)
        # Tr_ab_vs_prev_75Q <- Tr_ab_vs_prev + coord_cartesian(xlim = c(-5, quantile(df_ab_prev$abundance, .75)))
        
        Tr_prev_vs_log10_ab <- ggplot(df_ab_prev, aes(x = total_counts_of_ASV, y = prevalence))
        if (is.null(col)) {
                Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab + geom_point(col = cbPalette[2], size = 2, alpha = 0.7) 
        } else {
                Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab + geom_point(aes_string(col = col), size = 2, alpha = 0.7)
        }
        Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab +
                scale_x_log10() +
                scale_y_continuous(limits = c(-1, max(df_ab_prev$prevalence) + 5)) +
                geom_smooth(method = "lm") +
                # annotate("text", x = max(df_ab_prev$abundance)/5, y = 5, label = paste("R.square of lm: ", as.character(round(summary(fit_prev)$r.squared,4), sep = ""))) +
                xlab("total counts of ASV (= taxa_sums())") +
                ggtitle(paste("lm fit: p.val: ", format(pval_prev_log10, digits = 4), "R.square: ", as.character(round(summary(fit_prev_log10)$r.squared,4), sep = ""))) +
                theme_bw(12)
        
        # NB: you could color or facet by phylum
        
        # Tr_prev_vs_log10_ab <- ggplot(df_ab_prev, aes(x = total_counts_of_ASV, y = prevalence))
        # Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab +
        #         geom_point(aes(col = Phylum), size = 2, alpha = 0.7) +
        #         scale_x_log10() +
        #         geom_smooth(method = "lm") +
        #         # annotate("text", x = max(df_ab_prev$abundance)/5, y = 5, label = paste("R.square of lm: ", as.character(round(summary(fit_prev)$r.squared,4), sep = ""))) +
        #         xlab("total abundance (taxa_sums())") +
        #         ggtitle(paste("lm fit: p.val: ", format(pval_prev_log10, digits = 4), "R.square: ", as.character(round(summary(fit_prev_log10)$r.squared,4), sep = ""))) +
        #         theme_bw(12)
        
        
        # Tr_spar_vs_meanab <- ggplot(df_ab_prev, aes(x = mean_count_nonzero, y = sparsity))
        # Tr_spar_vs_meanab <- Tr_spar_vs_meanab +
        #         geom_point(col = cbPalette[2], size = 2, alpha = 0.7) +
        #         geom_smooth(method = "lm") +
        #         scale_y_continuous(limits = c(-1, nsamples + 5)) +
        #         xlab("SVs mean abundance in non-zero samples") +
        #         ggtitle(paste("lm fit: p.val: ", format(pval_spar_mean, digits = 4), "R.square: ", as.character(round(summary(fit_spar_mean)$r.squared,4), sep = ""))) +
        #         theme_bw(12)
        # Tr_spar_vs_meanab_75Q <- Tr_spar_vs_meanab + coord_cartesian(xlim = c(-5, quantile(df_ab_prev$mean_count_nonzero, .75)))
        # 
        
        # Tr_spar_vs_log10_meanab <- ggplot(df_ab_prev, aes(x = mean_count_nonzero, y = sparsity))
        # Tr_spar_vs_log10_meanab <- Tr_spar_vs_log10_meanab +
        #         geom_point(col = cbPalette[2], size = 2, alpha = 0.7) +
        #         scale_x_log10() +
        #         geom_smooth(method = "lm") +
        #         xlab("SVs mean abundance in non-zero samples") +
        #         ggtitle(paste("lm fit: p.val: ", format(pval_spar_mean_log10, digits = 4), "R.square: ", as.character(round(summary(fit_spar_mean_log10)$r.squared,4), sep = ""))) +
        #         theme_bw(12)
        
        Tr_prev_vs_log10_meanab <- ggplot(df_ab_prev, aes(x = mean_count_nonzero, y = prevalence))
        if (is.null(col)) {
                Tr_prev_vs_log10_meanab <- Tr_prev_vs_log10_meanab + geom_point(col = cbPalette[2], size = 2, alpha = 0.7) 
        } else {
                Tr_prev_vs_log10_meanab <- Tr_prev_vs_log10_meanab + geom_point(aes_string(col = col), size = 2, alpha = 0.7)
        }
        Tr_prev_vs_log10_meanab <- Tr_prev_vs_log10_meanab +
                scale_x_log10() +
                geom_smooth(method = "lm") +
                scale_y_continuous(limits = c(-1, max(df_ab_prev$prevalence) + 5)) +
                xlab("mean count of ASV in non-zero samples") +
                ggtitle(paste("lm fit: p.val: ", format(pval_prev_mean_log10, digits = 4), "R.square: ", as.character(round(summary(fit_prev_mean_log10)$r.squared,4), sep = ""))) +
                theme_bw(12)
        
        Tr_prev_vs_log10_medianab <- ggplot(df_ab_prev, aes(x = median_count_nonzero, y = prevalence))
        if (is.null(col)) {
                Tr_prev_vs_log10_medianab <- Tr_prev_vs_log10_medianab + geom_point(col = cbPalette[2], size = 2, alpha = 0.7) 
        } else {
                Tr_prev_vs_log10_medianab <- Tr_prev_vs_log10_medianab + geom_point(aes_string(col = col), size = 2, alpha = 0.7)
        }
        Tr_prev_vs_log10_medianab <- Tr_prev_vs_log10_medianab +
                scale_x_log10() +
                geom_smooth(method = "lm") +
                scale_y_continuous(limits = c(-1, max(df_ab_prev$prevalence) + 5)) +
                xlab("median count of ASV in non-zero samples") +
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



#######################################
#### plot_abundance_prev_filter
#######################################
# df_ab_prev: data frame with "SV_ID", "total_counts_of_ASV", "prevalence", "sparsity", "mean_count_nonzero",
# "median_count_nonzero" plus tax_table

plot_abundance_prev_filter <- function(physeq, prevalence, taxa_sums_quantile){ 
        
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  df_ab_prev <- data.frame(SV_ID = 1:ntaxa(physeq), 
                                 total_counts_of_ASV = taxa_sums(physeq),
                                 prevalence = colSums(as(otu_table(physeq), "matrix") != 0),
                                 sparsity = colSums(as(otu_table(physeq), "matrix") == 0), 
                                 mean_count_nonzero = apply(as(otu_table(physeq), "matrix"), 2, function(x){mean(x[x > 0])}),
                                 median_count_nonzero = apply(as(otu_table(physeq), "matrix"), 2, function(x){median(x[x > 0])}))
        
        df_ab_prev <- cbind(df_ab_prev, tax_table(physeq))
        
        prev_thresh <- (prevalence/100)*nsamples(physeq)
        abund_thresh <- quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100)
        
        df_ab_prev_filt <- dplyr::filter(df_ab_prev, prevalence > prev_thresh | total_counts_of_ASV > abund_thresh)
        
        no_samples <- df_ab_prev$prevalence[1] + df_ab_prev$sparsity[1]
        shade_df <- data.frame(total_counts_of_ASV = 0, prevalence = 0)
        
        
        Tr_prev_vs_log10_ab <- ggplot(df_ab_prev, aes(x = total_counts_of_ASV, y = prevalence))
        Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab +
                geom_point(col = cbPalette[2], size = 2, alpha = 0.7) +
                scale_x_log10() +
                geom_rect(data = shade_df, xmin = -Inf, xmax = log10(abund_thresh), ymin = -Inf, ymax = prev_thresh, fill = "#660033", alpha = 0.4) +
                geom_hline(yintercept = (prevalence/100)*nsamples(physeq), col = cbPalette[1], lty = "dashed") +
                geom_vline(xintercept = quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100), col = cbPalette[1], lty = "dashed") +
                xlab("total counts of ASVs (taxa_sums())") + 
                theme_bw(12) +
                ggtitle(paste(nrow(df_ab_prev_filt), " of ", nrow(df_ab_prev), " ASVs (", round(100*nrow(df_ab_prev_filt)/nrow(df_ab_prev), 1),
                              " %) and ", round(sum(df_ab_prev_filt$total_counts_of_ASV)), " of ", round(sum(df_ab_prev$total_counts_of_ASV)), " counts (",
                              round((sum(df_ab_prev_filt$total_counts_of_ASV)/sum(df_ab_prev$total_counts_of_ASV))*100, 1), " %) remain", sep = ""))
        
        
        Tr_prev_vs_log10_ab_col <- ggplot(df_ab_prev, aes(x = total_counts_of_ASV, y = prevalence))
        Tr_prev_vs_log10_ab_col <- Tr_prev_vs_log10_ab_col +
                geom_point(aes(col = Phylum), size = 2, alpha = 0.7) +
                scale_x_log10() +
                geom_rect(data = shade_df, xmin = -Inf, xmax = log10(abund_thresh), ymin = -Inf, ymax = prev_thresh, fill = "#660033", alpha = 0.4) +
                geom_hline(yintercept = (prevalence/100)*nsamples(physeq), col = cbPalette[1], lty = "dashed") +
                geom_vline(xintercept = quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100), col = cbPalette[1], lty = "dashed") +
                xlab("total counts of ASVs (taxa_sums())") +
                facet_wrap(~Phylum) +
                theme_bw(12) +
                theme(legend.position = "none")
        
        
        # phylum_df <- df_ab_prev[, c("Phylum", "total_counts_of_ASV", "prevalence")]
        # phylum_df <- group_by(phylum_df, Phylum)
        # phylum_df <- dplyr::summarise(phylum_df, SVs = n(), abundance = round(sum(total_counts_of_ASV)))
        # phylum_df_filt <- df_ab_prev_filt[, c("Phylum", "total_counts_of_ASV", "prevalence")]
        # phylum_df_filt <- group_by(phylum_df_filt, Phylum)
        # phylum_df_filt <- dplyr::summarise(phylum_df_filt, SVs = n(), abundance = round(sum(total_counts_of_ASV)))
        # phylum_df_summary <- merge(phylum_df, phylum_df_filt, by = "Phylum")
        # colnames(phylum_df_summary) <- c("Phylum", "SVs_before", "abundance_before", "SVs_after", 'abundance_after')
        # phylum_df_summary <- mutate(phylum_df_summary, SV_r_PC = round(100*SVs_after/SVs_before, 1), abundance_r_PC = round(100*abundance_after/abundance_before, 1),
        #                             SV_PC = round(100*SVs_after/sum(SVs_after), 1), abundance_PC = round(100*abundance_after/sum(abundance_after), 1))
        
        Before <- dplyr::summarise(group_by(df_ab_prev, Phylum), ASVs_bef = n(), PC_ASV_bef = round(100*n()/nrow(df_ab_prev),1), PC_total_pre_bef = round(100*sum(prevalence)/sum(df_ab_prev$prevalence), 1),
                           PC_total_ab_bef = round(100*sum(total_counts_of_ASV)/sum(df_ab_prev$total_counts_of_ASV), 1), 
                           mean_pre_bef_inPC = round(100*mean(prevalence)/no_samples, 1),
                           mean_tot_ab_bef = round(mean(total_counts_of_ASV)), 
                           mean_mean_ab_nonzero_bef = round(mean(mean_count_nonzero)),
                           med_med_ab_nonzero_bef = round(median(median_count_nonzero)),
                           total_ab_bef = sum(total_counts_of_ASV))
        
        After <- dplyr::summarise(group_by(df_ab_prev_filt, Phylum), ASVs_aft = n(), PC_ASV_aft = round(100*n()/nrow(df_ab_prev_filt),1), PC_total_pre_aft = round(100*sum(prevalence)/sum(df_ab_prev_filt$prevalence), 1),
                           PC_total_ab_aft = round(100*sum(total_counts_of_ASV)/sum(df_ab_prev_filt$total_counts_of_ASV), 1), 
                           mean_pre_aft_inPC = round(100*mean(prevalence)/no_samples, 1),
                           mean_tot_ab_aft = round(mean(total_counts_of_ASV)),
                           mean_mean_ab_nonzero_aft = round(mean(mean_count_nonzero)),
                           med_med_ab_nonzero_aft = round(median(median_count_nonzero)),
                           total_ab_aft = sum(total_counts_of_ASV))
        
        Before_total <- dplyr::summarise(df_ab_prev, ASVs_bef = n(), PC_ASV_bef = round(100*n()/nrow(df_ab_prev),1), PC_total_pre_bef = round(100*sum(prevalence)/sum(df_ab_prev$prevalence), 1),
                            PC_total_ab_bef = round(100*sum(total_counts_of_ASV)/sum(df_ab_prev$total_counts_of_ASV), 1), 
                            mean_pre_bef_inPC = round(100*mean(prevalence)/no_samples, 1),
                            mean_tot_ab_bef = round(mean(total_counts_of_ASV)), 
                            mean_mean_ab_nonzero_bef = round(mean(mean_count_nonzero)),
                            med_med_ab_nonzero_bef = round(median(median_count_nonzero)),
                            total_ab_bef = sum(total_counts_of_ASV))
        
        Before_total <- data.frame(Phylum = "Total", Before_total)
        
        After_total <- dplyr::summarise(df_ab_prev_filt, ASVs_aft = n(), PC_ASV_aft = round(100*n()/nrow(df_ab_prev_filt),1), PC_total_pre_aft = round(100*sum(prevalence)/sum(df_ab_prev_filt$prevalence), 1),
                           PC_total_ab_aft = round(100*sum(total_counts_of_ASV)/sum(df_ab_prev_filt$total_counts_of_ASV), 1), 
                           mean_pre_aft_inPC = round(100*mean(prevalence)/no_samples, 1),
                           mean_tot_ab_aft = round(mean(total_counts_of_ASV)),
                           mean_mean_ab_nonzero_aft = round(mean(mean_count_nonzero)),
                           med_med_ab_nonzero_aft = round(median(median_count_nonzero)),
                           total_ab_aft = sum(total_counts_of_ASV))
        
        After_total <- data.frame(Phylum = "Total", After_total)
        
        Before <- rbind(Before, Before_total)
        
        After <- rbind(After, After_total)
        
        
        Merged <- merge(Before, After, by = "Phylum", all = TRUE, sort = FALSE)
        
        Merged <- dplyr::mutate(Merged, PC_ASV_rem = round(100*ASVs_aft/ASVs_bef, 1), PC_ab_rem = round(100*total_ab_aft/total_ab_bef, 1))
        
        Merged <- dplyr::select(Merged, 1, 20, 21, 2, 11, 3, 12, 4, 13, 5, 14, 6, 15, 7, 16, 8, 17, 9, 18, 10, 19)
        
        # in case some phyla have been kicked out: you need to reorder like it was "Before" to have total at the last place
        Merged[is.na(Merged)] <- 0 # Phyla that have been removed are NA
        colnames(Merged)[20:21] <- c("tot_ab_bef", "tot_ab_aft")
        Merged[, 20] <- round(Merged[, 20])
        Merged[, 21] <- round(Merged[, 21])
        
        #Merged$tot_ab_bef <- round(Before$total_ab_bef)
        #Merged$tot_ab_aft <- round(After$total_ab_aft)
        
        out <- list(Tr_prev_vs_log10_ab = Tr_prev_vs_log10_ab,
                    Tr_prev_vs_log10_ab_col = Tr_prev_vs_log10_ab_col,
                    phylum_df_summary = Merged)
        
}


####
#ASV with non-zero abundance in at least 2 protocols across all samples 
####
plot_abundance_protocols_filter <- function(physeq, prevalence, prevalence_methods, taxa_sums_quantile){ 
  
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  df_ab_prev <- data.frame(seqs=taxa_names(physeq),
                           SV_ID = 1:ntaxa(physeq), 
                           total_counts_of_ASV = taxa_sums(physeq),
                           prevalence = colSums(as(otu_table(physeq), "matrix") != 0),
                           sparsity = colSums(as(otu_table(physeq), "matrix") == 0), 
                           mean_count_nonzero = apply(as(otu_table(physeq), "matrix"), 2, function(x){mean(x[x > 0])}),
                           median_count_nonzero = apply(as(otu_table(physeq), "matrix"), 2, function(x){median(x[x > 0])}))
  
  prevalence_methods2<-prevalence_methods[,c(1,7)]
  df_ab_prev<-merge(prevalence_methods2, df_ab_prev, by="seqs")
  row.names(df_ab_prev)<-df_ab_prev$seqs
  df_ab_prev<-df_ab_prev[,-1]
  df_ab_prev <- cbind(df_ab_prev, tax_table(physeq))
  
  prev_thresh <- (prevalence/100)*nsamples(physeq)
  abund_thresh <- quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100)
  Meth_thresh<-2
  
  
  df_ab_prev_filt <- dplyr::filter(df_ab_prev,  Number_Methods >=Meth_thresh)
  df_ab_prev_filt1 <- dplyr::filter(df_ab_prev, (prevalence > prev_thresh | total_counts_of_ASV > abund_thresh) & Number_Methods >=Meth_thresh)
  
  no_samples <- df_ab_prev$prevalence[1] + df_ab_prev$sparsity[1]
  shade_df <- data.frame(total_counts_of_ASV = 0, prevalence = 0)
  
  
  Tr_prev_vs_log10_ab <- ggplot(df_ab_prev_filt, aes(x = total_counts_of_ASV, y = prevalence))
  Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab +
    #geom_point(col = cbPalette[2], size = 2, alpha = 0.7) +
    geom_point(aes(col = factor(Number_Methods)), size = 2, alpha = 0.7) +
    scale_x_log10() +
    geom_rect(data = shade_df, xmin = -Inf, xmax = log10(abund_thresh), ymin = -Inf, ymax = prev_thresh, fill = "#660033", alpha = 0.4) +
    geom_hline(yintercept = (prevalence/100)*nsamples(physeq), col = cbPalette[1], lty = "dashed") +
    geom_vline(xintercept = quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100), col = cbPalette[1], lty = "dashed") +
    xlab("total counts of ASVs (taxa_sums())") + 
    theme_bw(12) +
    ggtitle(paste(nrow(df_ab_prev_filt), " of ", nrow(df_ab_prev), " ASVs (", round(100*nrow(df_ab_prev_filt)/nrow(df_ab_prev), 1),
                  " %) and ", round(sum(df_ab_prev_filt$total_counts_of_ASV)), " of ", round(sum(df_ab_prev$total_counts_of_ASV)), " counts (",
                  round((sum(df_ab_prev_filt$total_counts_of_ASV)/sum(df_ab_prev$total_counts_of_ASV))*100, 1), " %) remain", sep = ""))+labs(subtitle=paste("All filters: ", nrow(df_ab_prev_filt1), " of ", nrow(df_ab_prev), " ASVs (", round(100*nrow(df_ab_prev_filt1)/nrow(df_ab_prev), 1),
                       " %) and ", round(sum(df_ab_prev_filt1$total_counts_of_ASV)), " of ", round(sum(df_ab_prev$total_counts_of_ASV)), " counts (",
                       round((sum(df_ab_prev_filt1$total_counts_of_ASV)/sum(df_ab_prev$total_counts_of_ASV))*100, 1), " %) remain", sep = ""))+
    theme(legend.position = "bottom")+guides(color=guide_legend(title="Number of methods"))
 
  
  
  
  Tr_prev_vs_log10_ab_col <- ggplot(df_ab_prev_filt, aes(x = total_counts_of_ASV, y = prevalence))
  Tr_prev_vs_log10_ab_col <- Tr_prev_vs_log10_ab_col +
    geom_point(aes(col = factor(Number_Methods)), size = 2, alpha = 0.7) +
    scale_x_log10() +
    geom_rect(data = shade_df, xmin = -Inf, xmax = log10(abund_thresh), ymin = -Inf, ymax = prev_thresh, fill = "#660033", alpha = 0.4) +
    geom_hline(yintercept = (prevalence/100)*nsamples(physeq), col = cbPalette[1], lty = "dashed") +
    geom_vline(xintercept = quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100), col = cbPalette[1], lty = "dashed") +
    xlab("total counts of ASVs (taxa_sums())") +
    facet_wrap(~Phylum) +
    theme_bw(12) +theme(legend.position = "bottom")+guides(color=guide_legend(title="Number of methods"))
  
  # phylum_df <- df_ab_prev[, c("Phylum", "total_counts_of_ASV", "prevalence")]
  # phylum_df <- group_by(phylum_df, Phylum)
  # phylum_df <- dplyr::summarise(phylum_df, SVs = n(), abundance = round(sum(total_counts_of_ASV)))
  # phylum_df_filt <- df_ab_prev_filt[, c("Phylum", "total_counts_of_ASV", "prevalence")]
  # phylum_df_filt <- group_by(phylum_df_filt, Phylum)
  # phylum_df_filt <- dplyr::summarise(phylum_df_filt, SVs = n(), abundance = round(sum(total_counts_of_ASV)))
  # phylum_df_summary <- merge(phylum_df, phylum_df_filt, by = "Phylum")
  # colnames(phylum_df_summary) <- c("Phylum", "SVs_before", "abundance_before", "SVs_after", 'abundance_after')
  # phylum_df_summary <- mutate(phylum_df_summary, SV_r_PC = round(100*SVs_after/SVs_before, 1), abundance_r_PC = round(100*abundance_after/abundance_before, 1),
  #                             SV_PC = round(100*SVs_after/sum(SVs_after), 1), abundance_PC = round(100*abundance_after/sum(abundance_after), 1))
  
  Before <- dplyr::summarise(group_by(df_ab_prev, Phylum), ASVs_bef = n(), PC_ASV_bef = round(100*n()/nrow(df_ab_prev),1), PC_total_pre_bef = round(100*sum(prevalence)/sum(df_ab_prev$prevalence), 1),
                             PC_total_ab_bef = round(100*sum(total_counts_of_ASV)/sum(df_ab_prev$total_counts_of_ASV), 1), 
                             mean_pre_bef_inPC = round(100*mean(prevalence)/no_samples, 1),
                             mean_tot_ab_bef = round(mean(total_counts_of_ASV)), 
                             mean_mean_ab_nonzero_bef = round(mean(mean_count_nonzero)),
                             med_med_ab_nonzero_bef = round(median(median_count_nonzero)),
                             total_ab_bef = sum(total_counts_of_ASV))
  
  After <- dplyr::summarise(group_by(df_ab_prev_filt, Phylum), ASVs_aft = n(), PC_ASV_aft = round(100*n()/nrow(df_ab_prev_filt),1), PC_total_pre_aft = round(100*sum(prevalence)/sum(df_ab_prev_filt$prevalence), 1),
                            PC_total_ab_aft = round(100*sum(total_counts_of_ASV)/sum(df_ab_prev_filt$total_counts_of_ASV), 1), 
                            mean_pre_aft_inPC = round(100*mean(prevalence)/no_samples, 1),
                            mean_tot_ab_aft = round(mean(total_counts_of_ASV)),
                            mean_mean_ab_nonzero_aft = round(mean(mean_count_nonzero)),
                            med_med_ab_nonzero_aft = round(median(median_count_nonzero)),
                            total_ab_aft = sum(total_counts_of_ASV))
  
  Before_total <- dplyr::summarise(df_ab_prev, ASVs_bef = n(), PC_ASV_bef = round(100*n()/nrow(df_ab_prev),1), PC_total_pre_bef = round(100*sum(prevalence)/sum(df_ab_prev$prevalence), 1),
                                   PC_total_ab_bef = round(100*sum(total_counts_of_ASV)/sum(df_ab_prev$total_counts_of_ASV), 1), 
                                   mean_pre_bef_inPC = round(100*mean(prevalence)/no_samples, 1),
                                   mean_tot_ab_bef = round(mean(total_counts_of_ASV)), 
                                   mean_mean_ab_nonzero_bef = round(mean(mean_count_nonzero)),
                                   med_med_ab_nonzero_bef = round(median(median_count_nonzero)),
                                   total_ab_bef = sum(total_counts_of_ASV))
  
  Before_total <- data.frame(Phylum = "Total", Before_total)
  
  After_total <- dplyr::summarise(df_ab_prev_filt, ASVs_aft = n(), PC_ASV_aft = round(100*n()/nrow(df_ab_prev_filt),1), PC_total_pre_aft = round(100*sum(prevalence)/sum(df_ab_prev_filt$prevalence), 1),
                                  PC_total_ab_aft = round(100*sum(total_counts_of_ASV)/sum(df_ab_prev_filt$total_counts_of_ASV), 1), 
                                  mean_pre_aft_inPC = round(100*mean(prevalence)/no_samples, 1),
                                  mean_tot_ab_aft = round(mean(total_counts_of_ASV)),
                                  mean_mean_ab_nonzero_aft = round(mean(mean_count_nonzero)),
                                  med_med_ab_nonzero_aft = round(median(median_count_nonzero)),
                                  total_ab_aft = sum(total_counts_of_ASV))
  
  After_total <- data.frame(Phylum = "Total", After_total)
  
  Before <- rbind(Before, Before_total)
  
  After <- rbind(After, After_total)
  
  
  Merged <- merge(Before, After, by = "Phylum", all = TRUE, sort = FALSE)
  
  Merged <- dplyr::mutate(Merged, PC_ASV_rem = round(100*ASVs_aft/ASVs_bef, 1), PC_ab_rem = round(100*total_ab_aft/total_ab_bef, 1))
  
  Merged <- dplyr::select(Merged, 1, 20, 21, 2, 11, 3, 12, 4, 13, 5, 14, 6, 15, 7, 16, 8, 17, 9, 18, 10, 19)
  
  # in case some phyla have been kicked out: you need to reorder like it was "Before" to have total at the last place
  Merged[is.na(Merged)] <- 0 # Phyla that have been removed are NA
  colnames(Merged)[20:21] <- c("tot_ab_bef", "tot_ab_aft")
  Merged[, 20] <- round(Merged[, 20])
  Merged[, 21] <- round(Merged[, 21])
  
  #Merged$tot_ab_bef <- round(Before$total_ab_bef)
  #Merged$tot_ab_aft <- round(After$total_ab_aft)
  
  out <- list(Tr_prev_vs_log10_ab = Tr_prev_vs_log10_ab,
              Tr_prev_vs_log10_ab_col = Tr_prev_vs_log10_ab_col,
              phylum_df_summary = Merged, df_ab_prev_filt)
  
}




#######################################
#### plot_canonical_correspondance
#######################################

plot_canonical_correspondance <- function(ps_ccpna, physeq, the_formula, group_var = group_var, second_ccp_variable = second_ccp_variable) {
        
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
        ps_scores <- vegan::scores(ps_ccpna)
        sites <- data.frame(ps_scores$sites)
        sites$Sample <- rownames(sites)
        sampleDF <- sample_data(physeq)
        sampleDF$Sample <- rownames(sampleDF)
        sites <- dplyr::left_join(sites, sampleDF, by = "Sample")
        
        SVs <- data.frame(ps_scores$species)
        SVs$SV_id <- rownames(SVs)
        
        # tax <- tax_table(ps)@.Data %>%
        #         data.frame(stringsAsFactors = FALSE)
        tax <- as.data.frame(tax_table(physeq), stringsAsFactors = FALSE)
        tax$SV_id <- rownames(tax)
        main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                         "Coriobacteriales", "Verrucomicrobiales")
        tax$Order[!(tax$Order %in% main_orders)] <- "Other"
        tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
        tax$SV_id <- rownames(tax)
        # tax$SV_id <- seq_len(ncol(otu_table(ps)))
        
        SVs <- dplyr::left_join(SVs, tax, by = "SV_id")
        
        evals_prop <- 100 * ps_ccpna$CCA$eig[1:2] / sum(ps_ccpna$CA$eig)
        Tr1 <- ggplot() +
                geom_point(data = SVs, aes(x = CCA1, y = CCA2, col = Order), size = 1.5) +
                geom_point(data = sites, aes_string(x = "CCA1", y = "CCA2", shape = group_var), size = 2, alpha = 0.8) +
                # geom_text_repel(data = species %>% filter(CCA2 < -0.80106),
                #                 aes(x = CCA1, y = CCA2, label = otu_id),
                #                 size = 1.5, segment.size = 0.1) +
                # facet_grid(. ~ Group) +
                guides(col = guide_legend(override.aes = list(size = 3))) +
                labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
                     y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
                scale_color_brewer(palette = "Set2") +
                # coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.33) + # why actually the 0.33?
                theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0))) +
                theme_bw()
        
        if (!is.null(second_ccp_variable)){
                Tr2 <- ggplot() +
                        geom_point(data = SVs, aes(x = CCA1, y = CCA2, col = Order), size = 1.5) +
                        geom_point(data = sites, aes_string(x = "CCA1", y = "CCA2", shape = second_ccp_variable), size = 2, alpha = 0.8) +
                        # geom_text_repel(data = species %>% filter(CCA2 < -0.80106),
                        #                 aes(x = CCA1, y = CCA2, label = otu_id),
                        #                 size = 1.5, segment.size = 0.1) +
                        # facet_grid(. ~ Group) +
                        guides(col = guide_legend(override.aes = list(size = 3))) +
                        labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
                             y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
                        scale_color_brewer(palette = "Set2") +
                        # coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.33) + # why actually the 0.33?
                        theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0))) +
                        theme_bw()
        } else {
                Tr2 <- NULL
        }
        
        list(Tr1 = Tr1, Tr2 = Tr2)

        
}


#######################################
#### plot_bar_own
#######################################
# based on plot_bar from phyloseq, difference, orders and colors the Samples based on group_var, and orders the fill based on abundance
# I guess inputs can be guessed on

plot_bar_own <- function(physeq, x = "Sample", y = "Abundance", group_var, fill = NULL,
                         color_sample_names = TRUE, facet_grid = NULL){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        
        if (!is.factor(sample_data(physeq)[[group_var]])) {sample_data(physeq)[[group_var]] <- as.factor(sample_data(physeq)[[group_var]])}
        
        if (is.null(fill)) { fill = "Phylum"}
        
        mdf <- psmelt(physeq)
        
        # order samples accoridng to levels
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        mdf$Sample <- factor(mdf$Sample, levels = LookUpDF$Sample, ordered = TRUE)
        # mdf$Sample <- factor(mdf$Sample, levels = c("A-15A", "A-5A", "A-2A", "A-1A", "B-15A", "B-5A", "B-2A", "B-1A"), ordered = TRUE)
        
        # order fill levels according to abundance over all samples
        mdf[, fill] <- as.character(mdf[, fill])
        mdf[is.na(mdf[, fill]), fill] <- "NA"
        sums <- group_by_(mdf, fill) %>% summarise(sum_abundance = sum(Abundance)) %>% arrange(sum_abundance)
        mdf[, fill] <- factor(mdf[, fill], levels = as.character(sums[[1]]), ordered = TRUE)
        
        if (length(levels(LookUpDF$Group)) <= 7 && color_sample_names){
                color_lookup <- data.frame(level = levels(LookUpDF$Group), color = cbPalette[2:(length(levels(LookUpDF$Group)) + 1)])
                #LookUpDF$Group[match(levels(mdf$Sample), LookUpDF$Sample)]
                colxaxis <- as.character(color_lookup$color[match(LookUpDF$Group[match(levels(mdf$Sample), LookUpDF$Sample)], color_lookup$level)])
        } else {
                colxaxis <- rep("black", nrow(LookUpDF))
        }
        
        
        if (length(levels(mdf[, fill])) <= 8) {
                fill_colors <- cbPalette[1:length(levels(mdf[, fill]))]
                names(fill_colors) <- rev(levels(mdf[, fill]))
        } else {
                fill_colors <- rev(viridis(length(levels(mdf[, fill]))))
                names(fill_colors) <- rev(levels(mdf[, fill]))
        }
        
        
        Tr <- ggplot(mdf, aes_string(x = x, y = y, fill = fill))
        Tr <- Tr + 
                geom_bar(stat = "identity", position = "stack") +
                theme_bw() +
                scale_fill_manual(values = fill_colors) +
                xlab("") +
                theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, colour = colxaxis))
        
        if (!is.null(facet_grid)) {
                Tr <- Tr + facet_wrap(~Sample.type)
        }
        
        
        return(Tr)
}

########
plot_bar_own1 <- function(physeq, x = "Sample", y = "Abundance", group_var, fill = NULL,
                         color_sample_names = TRUE, facet_grid = NULL){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#C110CE","#83BA03", "#02BAB3" )
  
  if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
  
  if (!is.factor(sample_data(physeq)[[group_var]])) {sample_data(physeq)[[group_var]] <- as.factor(sample_data(physeq)[[group_var]])}
  
  if (is.null(fill)) { fill = "Phylum"}
  
  mdf <- psmelt(physeq)
  
  # order samples accoridng to levels
  LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
  LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
  mdf$Sample <- factor(mdf$Sample, levels = LookUpDF$Sample, ordered = TRUE)
  # mdf$Sample <- factor(mdf$Sample, levels = c("A-15A", "A-5A", "A-2A", "A-1A", "B-15A", "B-5A", "B-2A", "B-1A"), ordered = TRUE)
  
  # order fill levels according to abundance over all samples
  mdf[, fill] <- as.character(mdf[, fill])
  mdf[is.na(mdf[, fill]), fill] <- "NA"
  sums <- group_by_(mdf, fill) %>% summarise(sum_abundance = sum(Abundance)) %>% arrange(sum_abundance)
  mdf[, fill] <- factor(mdf[, fill], levels = as.character(sums[[1]]), ordered = TRUE)
  
  if (length(levels(LookUpDF$Group)) <= 7 && color_sample_names){
    color_lookup <- data.frame(level = levels(LookUpDF$Group), color = cbPalette[2:(length(levels(LookUpDF$Group)) + 1)])
    #LookUpDF$Group[match(levels(mdf$Sample), LookUpDF$Sample)]
    colxaxis <- as.character(color_lookup$color[match(LookUpDF$Group[match(levels(mdf$Sample), LookUpDF$Sample)], color_lookup$level)])
  } else {
    colxaxis <- rep("black", nrow(LookUpDF))
  }
  
  
  if (length(levels(mdf[, fill])) <= 8) {
    fill_colors <- cbPalette[1:length(levels(mdf[, fill]))]
    names(fill_colors) <- rev(levels(mdf[, fill]))
  } else {
    fill_colors <- rev(viridis(length(levels(mdf[, fill]))))
    names(fill_colors) <- rev(levels(mdf[, fill]))
  }
  
  
  Tr <- ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  Tr <- Tr + 
    geom_bar(stat = "identity", position = "stack") +
    theme_bw() +
    scale_fill_manual(values = fill_colors) +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, colour = colxaxis))
  
  if (!is.null(facet_grid)) {
    Tr <- Tr + facet_wrap(~Sample.type)
  }
  
  
  return(Tr)
}
#######
##############
###Deseq_plot
############################################################################################
#logChangePlot Camila

# Deseq_plot=function(sigtab, alpha, label){
#   
#   myColors <-c("#F0E442", "#E69F00", "#56B4E9", "#009E73")
#   
#   phylum.levels<-c("Proteobacteria","Actinobacteria","Bacteroidetes","Firmicutes")
#   phylum.levels<-factor(phylum.levels)
#   names(myColors) <- levels(phylum.levels)
#   colScale <- scale_colour_manual(name = "Phylum",values = myColors)
# sigtab<-sigtab[which(sigtab$p_val_adj <= alpha), ]
# x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
# x = sort(x, TRUE)
# sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# # Genus order
# x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
# x = sort(x, TRUE)
# sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
# sigtab$name<-apply(sigtab, 1, get_pretty_taxon_name)
# plot_log<-ggplot(sigtab, aes(x=log2FoldChange, y=name, color=Phylum)) + geom_point(size=4) +geom_yline(yintercept=0, color="red")+colScale+ylab(NULL)
# if(label==TRUE){plot_log<- plot_log+geom_text(aes(label=Taxa),position = position_dodge(0.9),vjust = -1.5, size = 4)+ theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5))+theme_bw()}  else(plot_log<-plot_log+ theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust=0.5))+theme_bw() ) 
# return(plot_log)
# }


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the pweeks
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#prevanlence filter
prevalence_function=function(ps, percent){
  prevdf = apply(X = otu_table(ps), MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts to this data.frame
  prevdf = data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(ps), tax_table(ps))
  #Prevalence Filtering
  prevalenceThreshold = percent * nsamples(ps)
  prevalenceThreshold
  # Execute prevalence filter, using `prune_taxa()` function
  keepTaxa= rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
  ps_1 = prune_taxa(keepTaxa, ps)
  ps_1
  return(ps_1)
}

#Prevanlece plot
prevalencePlot=function(ps, percent){
  return_list<-list()
  prevdf = apply(X = otu_table(ps), MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts to this data.frame
  prevdf = data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(ps), tax_table(ps))
  #Prevalence Filtering
  prevalenceThreshold = percent * nsamples(ps)
  countThreshold=0.1*(sum(prevdf$TotalAbundance))
  # prevalenceThreshold
  
  plot1<-ggplot(prevdf, aes(x=TotalAbundance, y=Prevalence))+geom_point()+theme_bw() + geom_vline(xintercept = countThreshold)+ geom_hline(yintercept = prevalenceThreshold)+scale_x_log10()
  plot2<- ggplot(prevdf, aes(x=TotalAbundance, y=Prevalence))+geom_point()+theme_bw() +facet_grid(.~Phylum)+ geom_vline(xintercept = countThreshold)+ geom_hline(yintercept = prevalenceThreshold)+  scale_x_log10()
  return_list[[1]]<-plot1 
  return_list[[2]]<-plot2 
  names(return_list)<-c("all", "Phylum")
  return(return_list)
}

#Abundance filter
Abundance_filter = function(ps, percent){
  ######Abundance
  ps2<-transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
  Matrix<-as.data.frame(otu_table(ps2))
  Keep <- (colSums(Matrix)*100/(sum(colSums(Matrix)))) >= percent 
  Matrix_Filter <- Matrix[,Keep]
  keepTaxa = colnames(Matrix_Filter)
  ps3 = prune_taxa(keepTaxa, ps)
  return(ps3)
}

#Metadata and abundance table
info=function(ps){
  table_info<-as.data.frame(sample_data(ps))
  abundance_table<-as.data.frame(otu_table(ps))
  community_table_complete<-merge(table_info, abundance_table, by="row.names")
  row.names(community_table_complete)<-community_table_complete$Row.names
  community_table_complete<-community_table_complete[-1]
  return(community_table_complete)
}

#Phylogenetic tree
phylogenetic_tree=function(ps){
  
  seqs <- colnames(otu_table(ps))
  names(seqs) <- seqs # This propagates to the tip labels of the tree
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
  phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
  dm <- dist.ml(phang.align)
  treeNJ <- NJ(dm) 
  fit = pml(treeNJ, data=phang.align)
  
  fitGTR <- update(fit, k=4, inv=0.2)
  fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
  detach("package:phangorn", unload=TRUE)
  return(fitGTR)}

#PCoA
ordination1=function(ps,distance, colorby,title, shape, label){
  ord<- ordinate(ps, method = "PCoA", distance = distance)
  plot <- plot_ordination(ps, ord, color = colorby, shape=shape, label=label)
  plot <- plot + theme_bw() + ggtitle(title)  + geom_point(size = 2.5)  + ggtitle(title)
  return(plot)
}

#DeseQ2
deseqfunction<-function(ps_deseq, type, formula, contrast, group1,group2){
  Groupds = phyloseq_to_deseq2(ps_deseq, ~ Disease.Progression+Group_name)
  Groupds <- estimateSizeFactors(Groupds, type= type)
  Groupds = DESeq(Groupds, test="Wald")
  res = results(Groupds, cooksCutoff = FALSE, contrast = c(contrast, group1, group2))
  res_df<-as.data.frame(res)
  sigtab= res[which(res$padj <= 1), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_deseq)[rownames(sigtab), ], "matrix"))
  return(sigtab)
}

#Deseq_plot
Deseq_plot=function(sigtab, alpha, label){
  myColors <-c("#F0E442", "#E69F00", "#56B4E9", "#009E73")
  phylum.levels<-c("Proteobacteria","Actinobacteria","Bacteroidetes","Firmicutes")
  phylum.levels<-factor(phylum.levels)
  names(myColors) <- levels(phylum.levels)
  colScale <- scale_colour_manual(name = "Phylum",values = myColors)

  sigtab<-sigtab[which(sigtab$p_val_adj <= alpha), ]
  sigtab_2<-cbind(row.names(sigtab), sigtab)
  colnames(sigtab_2)[1]<-"Taxa"
  colnames(sigtab_2)[6]<-"sig"
  x = tapply(sigtab_2$log2FoldChange, sigtab_2$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtab_2$Phylum = factor(as.character(sigtab_2$Phylum), levels=names(x))
  # Genus order
  
  sigtab_2$name<- apply(sigtab_2, 1, get_pretty_taxon_name)
  x = tapply(sigtab_2$log2FoldChange, sigtab_2$name, function(x) max(x))
 # x = tapply(sigtab_2$log2FoldChange, sigtab_2$Genus, function(x) max(x))
  x = sort(x, TRUE)
  #sigtab_2$Genus = factor(as.character(sigtab_2$Genus), levels=names(x))
  sigtab_2$name = factor(as.character(sigtab_2$name), levels=names(x))
  plot_log<-ggplot(sigtab_2, aes(x=log2FoldChange, y=name, color=Phylum)) + geom_point(size=4)+colScale +geom_vline(xintercept=0, color="red")+ylab(NULL)
  if(label==TRUE){plot_log<- plot_log+geom_text(aes(label=Taxa),position = position_dodge(0.9),vjust = -1.5, size = 4)+theme_bw()}  else(plot_log<-plot_log+theme_bw())#+ theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5) ) 
  return(plot_log)
}

biplot_pcoa <- function (physeq, color = "Group", shape = NULL, axis1 = 1, axis2 = 2, 
                         show.taxa = TRUE, label = "Genus", repel = FALSE, show.legend = TRUE, 
                         custom_palette = NULL, dist_method=dist_method) 
{
  dist <- phyloseq::distance(physeq = physeq, method = "bray")
  ordination <- phyloseq::ordinate(physeq, method = "PCoA", 
                                   distance = dist)
  axes = c(axis1, axis2)
  DF <- phyloseq::plot_ordination(physeq, ordination, axes = axes, 
                                  type = "biplot", justDF = TRUE)
  x = colnames(DF)[1]
  y = colnames(DF)[2]
  ord_map = aes_string(x = x, y = y, color = color, shape = shape, 
                       na.rm = TRUE)
  label_map <- aes_string(x = x, y = y, label = label, na.rm = TRUE)
  if (length(extract_eigenvalue(ordination)[axes]) > 0) {
    eigvec = extract_eigenvalue(ordination)
    fracvar = eigvec[axes]/sum(eigvec)
    percvar = round(100 * fracvar, 1)
  }
  DF_Taxa <- DF[DF$id.type == "Taxa", ]
  DF_Samples <- DF[DF$id.type == "Samples", ]
  score <- fracvar[1] * (DF_Taxa[x]^2) + fracvar[2] * (DF_Taxa[y]^2)
  colnames(score) <- NULL
  DF_Taxa <- cbind(DF_Taxa, Score = score)
  rm(score)
  DF_Taxa <- dplyr::arrange(DF_Taxa, desc(Score))
  DF_Taxa <- head(DF_Taxa, 50)
  Tr <- ggplot(DF_Samples, ord_map) + geom_point(na.rm = TRUE, 
                                                 size = 2, show.legend = show.legend) + theme(aspect.ratio = 1)
  color_values <- get_my_palette_colors(custom_palette, color, 
                                        levels(DF_Samples[[color]]), offset = 0)
  Tr <- Tr + scale_color_manual(values = color_values)
  if (show.taxa) {
    if (repel) {
      Tr <- Tr + geom_text_repel(label_map, data = rm.na.phyloseq(DF_Taxa,
                                                                  label), size = 3, color = "black", vjust = "center",
                                 hjust = "middle", show.legend = FALSE, na.rm = TRUE)

    }
    else {
      Tr <- Tr + geom_text(label_map, data = rm.na.phyloseq(DF_Taxa,
                                                            label), check_overlap = TRUE, size = 3, color = "black",
                           vjust = "center", hjust = "middle", show.legend = FALSE,
                           na.rm = TRUE)

    }
  }
  strivar = paste0("PC", axes, "  [", percvar, "%]")
  Tr = Tr + xlab(strivar[1]) + ylab(strivar[2])
  return(Tr)
}

get_my_palette_colors <- function (pal, name = NULL, levels = NULL, offset = 0) 
{
  l_palette <- get_my_palette_name(pal, name)
  if (is.null(l_palette)) {
    warning(paste("Cannot find palette for", name))
    return(NULL)
  }
  if (!is.null(names(l_palette))) {
    if (identical(sort(names(l_palette)), sort(levels))) {
      l_colors = l_palette
    }
  }
  else {
    l_n = length(levels)
    l_colors = brewer.pal(offset + max(l_n, 3), l_palette)
    l_colors = l_colors[(offset + 1):(offset + l_n)]
    names(l_colors) = levels
  }
  return(l_colors)
}

extract_eigenvalue<-function (ordination) {ordination$values$Relative_eig}
 


get_my_palette_name <- function (pal, name) 
{
  if (is.null(name)) {
    stop("get_my_palette_name() requires name of factor")
  }
  l_palette = NULL
  if (!is.null(pal)) {
    l_palette = pal[[name]]
  }
  else {
    l_palette = "Set2"
    warning(paste0("Cannot find palette for ", name, "; using ", 
                   l_palette))
  }
  return(l_palette)
}

rm.na.phyloseq <- function (DF, key.var) 
{
  DF <- subset(DF, !is.na(eval(parse(text = key.var))))
  if (class(DF[, key.var]) == "factor") {
    DF[, key.var] <- factor(as(DF[, key.var], "character"))
  }
  return(DF)
}







#Deseq_plot_fun
Deseq_plot_fun=function(sigtab, alpha, label){
  
  sigtab<-sigtab[which(sigtab$p_val_adj <= alpha), ]
  sigtab_2<-cbind(row.names(sigtab), sigtab)
  colnames(sigtab_2)[1]<-"Taxa"
  colnames(sigtab_2)[6]<-"sig"
  x = tapply(sigtab_2$log2FoldChange, sigtab_2$Level_2, function(x) max(x))
  x = sort(x, TRUE)
  sigtab_2$Phylum = factor(as.character(sigtab_2$Level_2), levels=names(x))
  # Genus order
  
  sigtab_2$name<- apply(sigtab_2, 1, get_pretty_taxon_name_fun)
  x = tapply(sigtab_2$log2FoldChange, sigtab_2$name, function(x) max(x))
  # x = tapply(sigtab_2$log2FoldChange, sigtab_2$Genus, function(x) max(x))
  x = sort(x, TRUE)
  #sigtab_2$Genus = factor(as.character(sigtab_2$Genus), levels=names(x))
  sigtab_2$name = factor(as.character(sigtab_2$name), levels=names(x))
  plot_log<-ggplot(sigtab_2, aes(x=log2FoldChange, y=name, color=Level_2)) + geom_point(size=4)+geom_vline(xintercept=0, color="red")+ylab(NULL)
  if(label==TRUE){plot_log<- plot_log+geom_text(aes(label=Taxa),position = position_dodge(0.9),vjust = -1.5, size = 4)+theme_bw() }  else(plot_log<-plot_log+theme_bw())
  return(plot_log)
}

Deseq_plot_fun_path=function(sigtab, alpha, label){
  
  sigtab<-sigtab[which(sigtab$p_val_adj <= alpha), ]
  sigtab_2<-cbind(row.names(sigtab), sigtab)
  colnames(sigtab_2)[1]<-"Taxa"
  colnames(sigtab_2)[6]<-"sig"
  x = tapply(sigtab_2$log2FoldChange, sigtab_2$Pathway, function(x) max(x))
  x = sort(x, TRUE)
  sigtab_2$Label = factor(as.character(sigtab_2$Pathway), levels=names(x))
  # Genus order
  
  sigtab_2$name<- apply(sigtab_2, 1, get_pretty_taxon_name_path)
  x = tapply(sigtab_2$log2FoldChange, sigtab_2$name, function(x) max(x))
  # x = tapply(sigtab_2$log2FoldChange, sigtab_2$Genus, function(x) max(x))
  x = sort(x, TRUE)
  #sigtab_2$Genus = factor(as.character(sigtab_2$Genus), levels=names(x))
  sigtab_2$name = factor(as.character(sigtab_2$name), levels=names(x))
  plot_log<-ggplot(sigtab_2, aes(x=log2FoldChange, y=name, color= Pathway)) + geom_point(size=4)+geom_vline(xintercept=0, color="red")+ylab(NULL)
  if(label==TRUE){plot_log<- plot_log+geom_text(aes(label=Taxa),position = position_dodge(0.9),vjust = -1.5, size = 4)+theme_bw() }  else(plot_log<-plot_log+theme_bw())
  return(plot_log)
}




#volcano_plot
volcano_plot=function(sigtab){
  p<- ggplot(data = sigtab, aes(x = log2FoldChange, y = -log10(padj),  colour=Phylum)) +
    geom_point()+ theme_bw() +
    geom_vline(xintercept = c(-1,1),colour = "grey")+
    geom_hline(yintercept = -log10(0.05),colour = "grey")
  p<-p+geom_text(aes(x = log2FoldChange, y = -log10(padj),label=Genus),position = position_dodge(0.9),vjust = -1.5, size = 4, data=sigtab[(abs(sigtab$padj))<=0.1,] )
  return(p)
}

#Odds ratio
OR_Function=function(ps_odds, Healthy, Disease, alpha, strata){
  
  info_general<-info(ps_odds)
  info_matched_Healthy<-subset(info_general, info_general$Method == Healthy)
  info_matched_Disease<-subset(info_general, info_general$Method == Disease)
  
  #An odds ratio (OR) is a measure of association between an exposure and an outcome.
  #The OR represents the odds that an outcome will occur given a particular exposure, 
  #compared to the odds of the outcome occurring in the absence of that exposure. 
  
  #OR=1 Exposure does not affect odds of outcome
  #OR>1 Exposure associated with higher odds of outcome
  #OR<1 Exposure associated with lower odds of outcome
  
  # a	The number of individuals who both suffer from exposure and disease.
  # b	The number of individuals who suffer from exposure but are healthy.
  # c	The number of individuals who suffer from disesase but not exposed.
  # d	The number of individuals who neither suffered from exposure nor disease.
  
  #A point estimate of the study result represented by a black box. 
  #This black box also gives a representation of the size of the study. 
  #The bigger the box, the more participants in the study
  
  #The range of values within which you can be 95% certain the true value lies.
  
  ########################
  #       Disease  Healthy
  #          +     -
  # taxa  +  a |   b
  #       -------------
  #       - c  |   d
  ########################
  
  matrix_list<-list()
  metadata_odds<-data.frame()
  tabletext<-data.frame()
  return_list<-list()
  f=0
  l=0
  for(j in taxa_names(ps_odds))
  {
    f=1+f
    a<-apply(info_matched_Disease[j], 2, function(c)sum(c!=0)) 
    c<-apply(info_matched_Disease[j], 2, function(c)sum(c==0)) 
    
    b<-apply(info_matched_Healthy[j], 2, function(c)sum(c!=0))
    d<-apply(info_matched_Healthy[j], 2, function(c)sum(c==0))
    
    oddsR<-oddsratio(a, b, c, d, conf.level=0.95,  p.calc.by.independence=TRUE) 
    
    metadata_odds<-rbind(metadata_odds,data.frame(mean  = oddsR$estimate,lower = oddsR$conf.int[1], upper = oddsR$conf.int[2], p.value=oddsR$p.value[1]) )
    
    if(!is.na(tax_table(ps_odds)[j,7])){tabletext<-rbind(tabletext, data.frame(RSV= paste(j, tax_table(ps_odds)[j,6], tax_table(ps_odds)[j,7], sep = "-" ) ,OR  = round(oddsR$estimate,2), p.value = round(oddsR$p.value,4) ))
    }else{tabletext<-rbind(tabletext, data.frame(RSV= paste(j, tax_table(ps_odds)[j,6], sep = "-" ) ,OR  = round(oddsR$estimate,2),p.value = round(oddsR$p.value,4) ))}
    matrix_list[[f]]<-oddsR
    
  }
  
  tabletext<-subset(tabletext, tabletext$OR != 0 & is.nan(tabletext$OR) == FALSE & is.infinite(tabletext$OR)==FALSE)
  
  metadata_odds<-subset(metadata_odds, row.names(metadata_odds)%in%row.names(tabletext))
  tabletext<-tabletext[order(tabletext$p.value),]
  metadata_odds<-metadata_odds[order(metadata_odds$p.value),]
  metadata_odds<-subset(metadata_odds[-4], row.names(metadata_odds)%in%row.names(tabletext))
  
  tabletext<-as.matrix(tabletext)
  
  plot<-forestplot(tabletext, 
                   hrzl_lines = gpar(col="#444444"),
                   metadata_odds,new_page = TRUE,
                   #clip=c(lower,upper), 
                   xlog=TRUE,
                   col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
                   vertices = TRUE)
  
  return_list<-list("forest.Plot"=plot, "tabletext"=tabletext, "metadata_odds"= metadata_odds, "matrix"=matrix_list)
  return(return_list)
}

#Mantel.Haenszel
Mantel.Haenszel= function(ps_odds, Healthy, Disease, strata){
  
  # a	The number of individuals who both suffer from exposure and disease.
  # b	The number of individuals who suffer from exposure but are healthy.
  # c	The number of individuals who suffer from disesase but not exposed.
  # d	The number of individuals who neither suffered from exposure nor disease.
  
  ####################
  #       Disease  Healthy
  #          +    -
  # taxa  +  a |  b
  #       ---------------
  #       - c  |  d
  
  
  confounder<-as.data.frame(unique(sample_data(ps_odds)$Date.of.birth))
  confounder[]<-lapply(confounder, as.character)
  
  return_list<-list()
  
  for(i in 1:dim(confounder)[1]){
    
    l=confounder[i,]
    ps_cohort= prune_samples(sample_data(ps_odds)$Date.of.birth==l,ps_odds)
    
    info_general<-info(ps_cohort)
    info_matched_Healthy<-subset(info_general, info_general$Method == Healthy)
    info_matched_Disease<-subset(info_general, info_general$Method == Disease)
    
    f=0
    matrix_list<-list()
    
    
    for(j in taxa_names(ps_cohort))
    {
      f=1+f
      a<-apply(info_matched_Disease[j], 2, function(c)sum(c!=0)) 
      c<-apply(info_matched_Disease[j], 2, function(c)sum(c==0)) 
      
      b<-apply(info_matched_Healthy[j], 2, function(c)sum(c!=0))
      d<-apply(info_matched_Healthy[j], 2, function(c)sum(c==0)) 
      
      matrix <- matrix(c(a,b,c,d),nrow=2,byrow=TRUE) 
      colnames(matrix) <- c(Disease,Healthy)
      rownames(matrix) <- c("Presence","Non-Presence")
      matrix_list[[f]]<-matrix
    }
    names(matrix_list)<-taxa_names(ps_cohort)
    return_list[[i]]<-matrix_list
  }
  
  cmh_list<-list()
  table<-data.frame()
  for(k in 1:length(taxa_names(ps_cohort))){
    
    matrix.array <- array(c(return_list[[1]][[k]],return_list[[2]][[k]], return_list[[3]][[k]],return_list[[4]][[k]],return_list[[5]][[k]],return_list[[6]][[k]],return_list[[7]][[k]], return_list[[8]][[k]],return_list[[9]][[k]], return_list[[10]][[k]]),dim=c(2,2,10))
    cmh<-cmh.test(matrix.array) 
    mh<-mantelhaen.test(matrix.array, conf.level = 0.95) 
    cmh_list[[k]]<-cmh
    table<-rbind(table, data.frame(RSV= paste(tax_table(ps_cohort)[k,6], sep = "-" ),p.value =round((cmh$parameter[[3]]), 4) , Pooled.Odd.Ratio=round((cmh$parameter[[5]]),2), MH.Estimate  =round(( cmh$parameter[[4]]),2), 
                                   Ratio.of.level.1=cmh$parameter[[6]], 
                                   Ratio.of.level.2=cmh$parameter[[7]],
                                   Ratio.of.level.3=cmh$parameter[[8]],
                                   Ratio.of.level.4=cmh$parameter[[9]], 
                                   Ratio.of.level.5=cmh$parameter[[10]],
                                   Ratio.of.level.6=cmh$parameter[[11]], 
                                   Ratio.of.level.7=cmh$parameter[[12]],
                                   Ratio.of.level.8=cmh$parameter[[13]],
                                   Ratio.of.level.9=cmh$parameter[[14]],
                                   Ratio.of.level.10=cmh$parameter[[15]],
                                   lower=mh$conf.int[1], 
                                   upper= mh$conf.int[2]) )
  }
  names(cmh_list)<-taxa_names(ps_cohort)
  row.names(table)<-taxa_names(ps_cohort)
  
  #table$p.adj<-round(p.adjust(table$p.value, method = "fdr"),4)
  
  table_filter<-subset(table, is.nan(table$MH.Estimate) == FALSE & is.infinite(table$MH.Estimate)==FALSE)
  metadata_filter<-data.frame(mean  = table_filter$MH.Estimate,lower = table_filter$lower, upper =  table_filter$upper, p.value=table_filter$p.value)
  
  table_filter<-table_filter[order(table_filter$p.value),]
  metadata_filter<-metadata_filter[order(metadata_filter$p.value),]
  
  
  tabletext_filter<-table_filter[1:4]
  tabletext_filter<-as.matrix(tabletext_filter)
  metadata_filter<-metadata_filter[-4]
  
  forestplot(tabletext_filter,
             hrzl_lines = gpar(col="#444444"),
             metadata_filter,new_page = TRUE,
             #clip=c(lower,upper),
             xlog=TRUE,
             col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
             vertices = TRUE)
  return(table_filter)
}
#Subset
subset_ps=function(ps_prune, ps.to.subset, prevalence.percent, abundance.percent){
  prevalence_ps= prevalence(ps_prune, prevalence.percent)
  abundance_ps=Abundance_filter(prevalence_ps, abundance.percent)
  count_prune=prune_samples(sample_names(abundance_ps),ps.to.subset)
  count_ps=prune_taxa(taxa_names(abundance_ps), count_prune)
  abun_glom<-tax_glom(abundance_ps, "Genus")
  count_glom<-tax_glom(count_ps, "Genus")
  list_ps=list(abundance_ps=abundance_ps,count_ps=count_ps, abun_glom=abun_glom, count_glom=count_glom)
  return(list_ps)
}

#Multiple_phyloseqObject
phyloseqObjects=function(ps,glom){
  abundance_ps<-transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
  count_ps=ps
  #abun_glom<-tax_glom(abundance_ps, "Genus")
  count_glom<-tax_glom(count_ps, glom)
  abun_glom<-transform_sample_counts(count_glom, function(OTU) OTU/sum(OTU))
  list_ps=list(abundance_ps=abundance_ps,count_ps=count_ps, abun_glom=abun_glom, count_glom=count_glom)
  return(list_ps)
}
#taxa_levels
taxa_levels=function(ps, taxoomic_rank){
  return_list<-list()
  ps<-tax_glom(ps, taxoomic_rank)
  taxa_table_glom<-as.data.frame(tax_table(ps))
  info_ps<-info(ps)
  dim(info_ps)
  names<-as.character(unlist(taxa_table_glom[taxoomic_rank]))
  colnames(info_ps)[8:dim(info_ps)[2]]<-names
  info_psm <- melt(info_ps,id.vars=colnames(info_ps)[1:7])
  
  plot1<-ggplot(info_psm, aes(x = Sample.ID, y = value))+ geom_area(aes(colour=variable, fill=variable,  group=variable),  position = 'stack')  + theme_bw()
  plot2<-plot1+ ylab("Abundance")+xlab(NULL)+scale_fill_brewer(palette="PiYG")+ theme(axis.text.x = element_blank(), legend.position="bottom") + guides(colour=guide_legend(title=taxoomic_rank), fill=guide_legend(title=taxoomic_rank))
  
  sample_data(ps)$Sample.ID <- factor(sample_data(ps)$Sample.ID,  levels = sample_data(ps)$Sample.ID[order(sample_data(ps)$Subject)])
  plot3<-plot_bar(ps, x="Sample.ID", fill=taxoomic_rank) + facet_wrap(~Method, scales="free_x")
  
  sample_data(ps)$Sample.ID <- factor(sample_data(ps)$Sample.ID,  levels = sample_data(ps)$Sample.ID[order(sample_data(ps)$Method)])
  plot4<-plot_bar(ps, x="Sample.ID", fill=taxoomic_rank) + facet_wrap(~Subject, scales="free_x")+theme(legend.position="bottom")
  
  return_list[[1]]<-plot2
  return_list[[2]]<-plot3
  return_list[[3]]<-plot4
  
  Wilcoxon_list<-list()
  Wilcoxon_df<-data.frame()
  
  for(i in names){
    wilcoxon_test<-pairwise.wilcox.test(as.numeric(unlist(info_ps[i])), info_ps$Method, p.adjust.method= "fdr")
    Wilcoxon_list[[i]]<-wilcoxon_test
    a<-melt(wilcoxon_test$p.value)
    b<-t(a[-c(1:2)])
    row.names(b)<-i
    Wilcoxon_df<-rbind(Wilcoxon_df,b)
  }
  colnames(Wilcoxon_df)<-c("MP.vs.MachereyNagel", "Qiagen.vs.MachereyNagel", "SOP.INRA.vs.MachereyNagel", "Zymogen.vs.MachereyNagel","SSI.vs.MachereyNagel", "NA", "Qiagen.vs.MP", "SOP.INRA.vs.MP", "Zymogen.vs.MP","SSI.vs.MP", "NA", "NA", "SOP.INRA.vs.Qiagen", "Zymogen.vs.Qiagen","SSI.vs.Qiagen", "NA", "NA", "NA", "Zymogen.vs.SOP.INRA","SSI.vs.SOP.INRA", "NA","NA","NA","NA", "SSI.vs.Zymogen")
  Wilcoxon_df<-Wilcoxon_df[-c(6,11,12,16,17,18,21,22,23,24)]
  
  names_selected<-Wilcoxon_df[which(Wilcoxon_df$MP.vs.MachereyNagel <= 0.05 |Wilcoxon_df$Qiagen.vs.MachereyNagel <= 0.05 | Wilcoxon_df$SOP.INRA.vs.MachereyNagel <= 0.05 | Wilcoxon_df$Zymogen.vs.MachereyNagel <= 0.05 | Wilcoxon_df$SSI.vs.MachereyNagel <= 0.05 | Wilcoxon_df$Qiagen.vs.MP <= 0.05 | Wilcoxon_df$SOP.INRA.vs.MP <= 0.05 | Wilcoxon_df$Zymogen.vs.MP <= 0.05 | Wilcoxon_df$SSI.vs.MP <= 0.05| Wilcoxon_df$SOP.INRA.vs.Qiagen <= 0.05 | Wilcoxon_df$Zymogen.vs.Qiagen <= 0.05 | Wilcoxon_df$SSI.vs.Qiagen <= 0.05 | Wilcoxon_df$Zymogen.vs.SOP.INRA <= 0.05 | Wilcoxon_df$SSI.vs.SOP.INRA<= 0.05 | Wilcoxon_df$SSI.vs.Zymogen <= 0.05), ]
  names_selected<-round(names_selected[], digits = 3)
  return_list[[4]]<-Wilcoxon_df
  return_list[[5]]<-names_selected
  
  info_psm_sign<-info_psm[which(info_psm$variable %in% row.names(names_selected)), ]
  
  
  plot<-ggplot(info_psm_sign, aes(x = variable, y = sqrt(value), color=Method)) + geom_boxplot()+ theme_bw() + ylab("Abundance")+xlab(NULL)+theme(legend.position="top")+guides(colour=guide_legend(title=element_blank()))
  
  
  return_list[[6]]<-plot
  
  names(return_list)<-c("area_plot", "method_plot","Subject_plot","Wilcoxon_test", "Wilcoxon_test_significant", "Wilcoxon_plot")
  
  return(return_list)
}
##Prevalences by different extraction methods
prevalence_methods <- function(ps, prevalence = 2) {
  return.list<-list()
  seqtab<-as.data.frame(otu_table(ps))
  samdf<-as.data.frame(sample_data(ps))
  seqs<-getSequences(as.matrix(seqtab))
  seq_relation<-data.frame("Sequence"=seqs, "ASVs_No"=(paste("ASV",c(1:length(seqs)), sep="_")))
  colnames(seqtab)<-(paste("ASV",c(1:length(seqs)), sep="_"))
  
  info_methods<-merge(samdf,seqtab, by=0)
  info_methods<-info_methods[-c(1:6,8)]
  info_methods_2<-data.frame()
  
  for(i in 1:dim(info_methods)[1]){
    sample<-info_methods[i,]
    for(j in 2:dim(sample)[2]){
      value<-as.numeric(sample[j])
      if(value!=0){sample[j]<-as.character(sample$Method)}
    } 
    info_methods_2<-rbind(info_methods_2, sample)
  }
  
  
  info_methods_2<-info_methods_2[-1]
  FinalNumbersMethods<-data.frame()
  for (k in 1:dim(info_methods_2)[2] ){
    df<-info_methods_2[k]
    ASV_info<-data.frame(
      ASV_ID=colnames(info_methods_2[k]),
      Macherey_Nagel=length(which(info_methods_2[k]=="Macherey Nagel")),
      MP=length(which(info_methods_2[k]=="MP")),
      Qiagen=length(which(info_methods_2[k]=="Qiagen")),
      SOP_INRA=length(which(info_methods_2[k]=="SOP INRA")),
      SSI=length(which(info_methods_2[k]=="SSI")),
      Zymogen=length(which(info_methods_2[k]=="Zymogen")))
    ASV_info$NoKits<-colSums(t(ASV_info[-1])!= 0)
    FinalNumbersMethods<-rbind(FinalNumbersMethods,ASV_info )
    
  }
  
  row.names(FinalNumbersMethods)<-FinalNumbersMethods$ASV_ID
  FinalNumbersMethods<- FinalNumbersMethods[-1]
  FinalNumbersSeq <- data.frame(Sequence = seqs, InNumberSamples = colSums(seqtab != 0), TotalCounts = colSums(seqtab))
  FinalNumbersMethodsSeq<-merge(FinalNumbersSeq, FinalNumbersMethods, by=0)
  row.names(FinalNumbersMethodsSeq)<-FinalNumbersMethodsSeq$Row.names
  FinalNumbersMethodsSeq<-FinalNumbersMethodsSeq[-1]
  
  FinalNumbersPerMethod<-data.frame(
    Macherey_Nagel= sum(FinalNumbersMethods$Macherey_Nagel!= 0),
    MP=sum(FinalNumbersMethods$MP!= 0),
    Qiagen=sum(FinalNumbersMethods$Qiagen != 0),
    SOP_INRA=sum(FinalNumbersMethods$SOP_INRA != 0),
    SSI=sum(FinalNumbersMethods$SSI!= 0),
    Zymogen=sum(FinalNumbersMethods$Zymogen!= 0))
  
  FinalNumbersMethodsSeq <- group_by(FinalNumbersMethodsSeq,  NoKits)
  CountDistribution <- dplyr::summarise(FinalNumbersMethodsSeq, No_ASVs = n(), TotalCounts = sum(TotalCounts))
  
  SummaryCountsMethod <-data.frame(
    Macherey_Nagel= sum((subset(info_methods, info_methods$Method=="Macherey Nagel"))[-1]),
    MP=sum((subset(info_methods, info_methods$Method=="MP"))[-1]),
    Qiagen=sum((subset(info_methods, info_methods$Method=="Qiagen"))[-1]),
    SOP_INRA=sum((subset(info_methods, info_methods$Method=="SOP INRA"))[-1]),
    SSI=sum((subset(info_methods, info_methods$Method=="SSI"))[-1]),
    Zymogen=sum((subset(info_methods, info_methods$Method=="Zymogen"))[-1]))
  
  MethodsSummary<-merge(t(FinalNumbersPerMethod), t(SummaryCountsMethod), by=0)
  colnames(MethodsSummary)<-c("Method","NoASV", "TotalCounts")
  MethodsSummary$PC.ASV<-((MethodsSummary$NoASV*100)/sum(MethodsSummary$NoASV))
  MethodsSummary$PC.Counts<-((MethodsSummary$TotalCounts*100)/sum(MethodsSummary$TotalCounts))
  
  keepTaxa=as.character(FinalNumbersMethodsSeq$Sequence[(FinalNumbersMethodsSeq$NoKits >= prevalence)])
  ps_prevalence = prune_taxa(keepTaxa, ps)
  
  
  
  MethodSummary$Method<-as.character(MethodSummary$Method)
  
  ASV<- ggplot(MethodSummary, aes(Method, PC.ASV, fill=Method))+ geom_bar(stat="identity")+theme_bw()+ ylab("%ASVs")+xlab(NULL)+scale_fill_brewer(palette="BrBG")+ theme(legend.position="none") 
  
  Counts<- ggplot(MethodSummary, aes(Method, PC.Counts, fill=Method))+ geom_bar(stat="identity")+theme_bw()+ ylab("%Counts")+xlab(NULL)+scale_fill_brewer(palette="BrBG")+ theme(legend.position="none") 
  
  return.list[[1]]<-FinalNumbersMethodsSeq
  return.list[[2]]<-MethodsSummary
  return.list[[3]]<-ASV
  return.list[[4]]<-Counts
  return.list[[5]]<-ps_prevalence
  return(return.list)
  
}

#Metadata_plots
metadataPlots=function(metadata){
  
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
#rarefaction
ggrare <- function(physeq, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  ## Args:
  ## - physeq: phyloseq class object, from which abundance data are extracted
  ## - step: Step size for sample size in rarefaction curves
  ## - label: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  
  ##          variable to map to colors in the plot. This can be a sample
  ##          variable (among the set returned by
  
  ##
  ##          Finally, The color scheme is chosen automatically by
  ## - color: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - plot:  Logical, should the graphic be plotted.
  ## - parallel: should rarefaction be parallelized (using parallel framework)
  ## - se:    Default TRUE. Logical. Should standard errors be computed. 
  ## require vegan
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    #cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  ## Get sample data 
  if (!is.null(sample_data(physeq, FALSE))) {
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  ## Add, any custom-supplied plot-mapped variables
  if( length(color) > 1 ){
    data$color <- color
    names(data)[names(data)=="color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if( length(label) > 1 ){
    labels$label <- label
    names(labels)[names(labels)=="label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot(data = data, aes_string(x = "Size", y = ".S", group = "Sample", color = color))
  p <- p + labs(x = "Sample Size", y = "Richness")
  if (!is.null(label)) {
    p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                       size = 4, hjust = 0)
  }
  p <- p + geom_line()
  if (se) { ## add standard error if available
    p <- p + geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}

#Abundance_plot
# Abundance_plot=function(ps, taxonomic_rank){
#   return_list<-list()
#   ps<-tax_glom(ps, taxonomic_rank)
#   taxa_table_glom<-as.data.frame(tax_table(ps))
#   info_ps<-info(ps)
#   dim(info_ps)
#   names<-as.character(unlist(taxa_table_glom[taxonomic_rank]))
#   colnames(info_ps)[11:dim(info_ps)[2]]<-names
#   info_psm <- melt(info_ps,id.vars=colnames(info_ps)[1:10])
#  
# 
#     plot<-ggplot(info_psm, aes(x = variable, y = sqrt(value), color=Group_name)) + geom_boxplot()+ theme_bw() + ylab("Abundance")+xlab(NULL)+theme(legend.position="top")+guides(colour=guide_legend(title=element_blank()))+scale_fill_viridis(discrete=TRUE)
#   
# return(plot)
# }



taxa.fungi=function (taxa){
  for(i in 1:nrow(taxa)){
    if(!is.na(taxa[i,1])){
      if(taxa[i,1]=="Eukaryota.Fungi.Fungi"){taxa[i,1]<-"Eukaryota.Fungi"}   
      
      if(taxa[i,1]=="Eukaryota.Fungi"){
        if(is.na(taxa[i,2])==TRUE){taxa[i,2]<-"Fungi.sp"}
        if(is.na(taxa[i,3])==TRUE){taxa[i,3]<-"Fungi.sp"}
        if(is.na(taxa[i,4])==TRUE){taxa[i,4]<-"Fungi.sp"}
        if(is.na(taxa[i,5])==TRUE){taxa[i,5]<-"Fungi.sp"}
        if(is.na(taxa[i,6])==TRUE){taxa[i,6]<-"Fungi.sp"}}
    }}
  return(taxa) 
}



Eukaryota_Filter=function(ps, F_QualityStats, R_QualityStats, ReadSummary)
{
  return_list=list()
  ps_phylum<- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
  Kingdom<-as.data.frame(table(tax_table(ps_phylum)[, "Kingdom"], exclude = NULL))
  colnames(Kingdom)<-c("Kingdom", "#ASV")
  rare_max_total <- quantile(sample_sums(ps), probs = .25)
  step_size <- 200
  rarefaction_curves <- suppressWarnings(rarefaction_curve_own_fast(physeq = ps, group_var = "Subject", max_total = rare_max_total, step_size = step_size))
  #psRare<-rarefy_even_depth(ps, sample.size = 7653)
  
  ps_Fungi<- subset_taxa(ps, Kingdom=="Eukaryota.Fungi")
  Phylum_fungi<-as.data.frame(table(tax_table(ps_Fungi)[, "Phylum"], exclude = NULL))
  colnames(Phylum_fungi)<-c("Phylum", "#ASV")
  assignment_distribution <- get_assignment_distribution(tax_table(ps_Fungi))
  #sum(taxa_sums(ps_Fungi)==0) 
  filterList_Fungi<- plot_abundance_prev_filter(physeq = ps_Fungi, prevalence = 10)
  rare_max_total <- quantile(sample_sums(ps_Fungi), probs = .25)
  rarefaction_curves_f <- suppressWarnings(rarefaction_curve_own_fast(physeq = ps_Fungi, group_var = "Subject", max_total = rare_max_total, step_size = step_size))
  psRare_f<-rarefy_even_depth(ps_Fungi)
  
  ps_Eukaryota<- subset_taxa(ps, Kingdom=="Eukaryota")
  Phylum<-as.data.frame(table(tax_table(ps_Eukaryota)[, "Phylum"], exclude = NULL))
  colnames(Phylum)<-c("Phylum", "#ASV")
  assignment_distribution <- get_assignment_distribution(tax_table(ps_Eukaryota))
  #sum(taxa_sums(ps_Eukaryota)==0) 
  filterList_Eukaryota <- plot_abundance_prev_filter(physeq = ps_Eukaryota, prevalence = 10)
  
  
  AlphaDivTotal <- as.data.frame(estimate_richness(ps_Fungi, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "Fisher")))
  plotRichness_shannonR<-plot_richness(ps_Fungi, x="Subject", measures= "Shannon", color="Subject") +  theme_bw()+xlab(NULL)+ theme(legend.position="none")  +geom_boxplot() +scale_color_viridis(discrete=TRUE)
  
  plotRichness_ObservedR<-plot_richness(ps_Fungi, x="Subject", measures= "Observed", color="Subject") +  theme_bw()+xlab(NULL)+ theme(legend.position="none")  +geom_boxplot()+scale_color_viridis(discrete=TRUE)
  
  ps_object_Eukaryota<-phyloseqObjects(ps_Eukaryota, glom = "Phylum")
  ps_object_Fungi<-phyloseqObjects(ps_Fungi, glom ="Genus")
  
  pcoa=ordination(ps_object_Fungi$abundance_ps, "bray", "Subject", "PCoA Bray-Curtis")
  
  abundance_eukaryota<-plot_bar(ps_object_Eukaryota$abun_glom, x="New.ID.s", fill="Phylum")+ facet_wrap(~Subject, scales="free_x",ncol = 2)+scale_fill_viridis(discrete=TRUE)+theme_bw()+ theme(legend.position="bottom")
  
  abundance_Fungi<-plot_bar(ps_object_Fungi$abun_glom, x="New.ID.s", fill="Genus")+ facet_wrap(~Subject, scales="free_x",ncol = 2)+scale_fill_viridis(discrete=TRUE)+theme_bw()+ theme(legend.position="bottom")
  sample.names<-sample_names(ps)
  
  Tr <- QS_Median_OverviewPlot(QStatsList = F_QualityStats, SampleNames = sample.names)
  Tr <- Tr + 
    geom_vline(xintercept = Input$trimLeft[1], color = 'darkred', linetype = "dashed", size = .5) +
    geom_vline(xintercept = Input$truncLen[1], color = 'darkred', linetype = "dashed", size = .5) +
    ggtitle(paste("Median quality scores: FW reads. No Samples: ", length(sample.names), "; trimLeft: ", Input$trimLeft[1],
                  "; truncLen: ", Input$truncLen[1], sep = ""))
  
  TrR <- QS_Median_OverviewPlot(QStatsList = R_QualityStats, SampleNames = sample.names)
  TrR <- TrR + 
    geom_vline(xintercept = Input$trimLeft[2], color = 'darkred', linetype = "dashed", size = .5) +
    geom_vline(xintercept = Input$truncLen[2], color = 'darkred', linetype = "dashed", size = .5) +
    ggtitle(paste("Median quality scores: RV reads. No Samples: ", length(sample.names), "; trimLeft: ", Input$trimLeft[1],
                  "; truncLen: ", Input$truncLen[2], sep = ""))
  
  
  TrTR <- NoReads_StepsSimple(ReadSummary = ReadSummary, SampleNames = sample.names, sort = TRUE)
  TrTR <- TrTR + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
                       text = element_text(size=14))
  
  
  return_list[[1]]<-ps
  return_list[[2]]<-rarefaction_curves$Tr_richness_group
  return_list[[3]]<-rarefaction_curves_f$Tr_richness_col
  return_list[[4]]<-ps_phylum
  return_list[[5]]<-Kingdom
  return_list[[6]]<-ps_Fungi
  return_list[[7]]<-Phylum_fungi
  return_list[[8]]<-filterList_Fungi[[1]]
  return_list[[9]]<-filterList_Fungi[[2]]
  return_list[[10]]<-pcoa
  return_list[[11]]<-AlphaDivTotal
  return_list[[12]]<-plotRichness_ObservedR
  return_list[[13]]<-plotRichness_shannonR
  return_list[[14]]<-abundance_eukaryota
  return_list[[15]]<-abundance_Fungi
  return_list[[16]]<-ps_object_Fungi
  return_list[[17]]<-Tr
  return_list[[18]]<-TrR
  return_list[[19]]<-TrTR
  #return_list[[16]]<-ps_object_Fungi
  #return_list[[17]]<-abundance_Fungi
  
  names(return_list)<-c("No.Filter","RarePlot","Rareplot_fungi", "ps.Phylum", "kingdom_tab","ps_Fungi", "Phylum_tab_fungi", "filter_prevalence_fungi", "fliter_prevalence_Phylum_fungi", "pcoa", "AlphaDivTotal_fungi","AlphaDivTObserved_fungi","AlphaDivTShannon_fungi", "abundance_eukaryota", "abundance_Fungi", "ps_object_Fungi", "Qp_F", "QP_r", "TR")
  return(return_list)
}




###################### OLD FUNCTIONS #############################

# #######################################
# #### NoReads_Steps
# #######################################
# 
# NoReads_Steps <- function(QStatsList = NULL, NoFilteredReads = NULL, mergers = NULL, mergers.nochim = NULL, SampleNames, sort = TRUE) {
#         
#         # Input:
#         # QStatsList: list of QStats data frames, such as FW_QualityStats
#         # SampleNames: List of Sample names that must be names of the QStatsList, NoFilteredReads, mergers, and mergers.nochim
#         # sort: if the samples should be sorted in the plot based on NoReads
#         
#         if (all(c(is.null(QStatsList), is.null(NoFilteredReads), is.null(mergers), is.null(mergers.nochim)))) {
#                 stop("QStatsList, derepFs, mergers, mergers.nochim can not all be NULL")
#         }
#         
#         if (!all(SampleNames %in% c(names(QStatsList), names(NoFilteredReads), names(mergers), names(mergers.nochim)))) {
#                 stop("not all SampleNames were found in the given data")
#         }
#         
#         cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#         
#         ## Construct the plotting data frame
#         
#         df.all <- NULL
#         
#         if (!is.null(QStatsList) & any(SampleNames %in% names(QStatsList))) {
#                 
#                 QStatsList <- QStatsList[SampleNames[which(SampleNames %in% names(QStatsList))]]
#                 
#                 for (i in seq_along(QStatsList)) {
#                         QStatsList[[i]]$Sample <- names(QStatsList[i])
#                 }
#                 
#                 df.all <- do.call("rbind",QStatsList)
#                 df.all <- df.all[!duplicated(df.all$Sample), c("Sample", "NoReads")]
#                 df.all$Type <- "all"
#                 
#         }
#         
#         df.filtered <- NULL
#         
#         if (!is.null(NoFilteredReads) & any(SampleNames %in% names(NoFilteredReads))) {
#                 
#                 NoFilteredReads <- NoFilteredReads[SampleNames[which(SampleNames %in% names(NoFilteredReads))]]
#                 
#                 df.filtered <- data.frame(Sample = names(NoFilteredReads), NoReads = NoFilteredReads, Type = "filtered")
#                 
#         }
#         
#         df.merged <- NULL
#         
#         if (!is.null(mergers) & any(SampleNames %in% names(mergers))) {
#                 
#                 mergers <- mergers[SampleNames[which(SampleNames %in% names(mergers))]]
#                 
#                 df.merged <- data.frame(Sample = names(mergers), NoReads = 0, Type = "merged")
#                 for(i in seq_along(mergers)) {
#                         df.merged$NoReads[i] <- sum(mergers[[i]]$abundance)
#                 }
#                 
#         } 
#         
#         df.nochim <- NULL
#         
#         if (!is.null(mergers.nochim) & any(SampleNames %in% names(mergers.nochim))) {
#                 
#                 mergers.nochim <- mergers.nochim[SampleNames[which(SampleNames %in% names(mergers.nochim))]]
#                 
#                 df.nochim <- data.frame(Sample = names(mergers.nochim), NoReads = 0, Type = "nochim")
#                 for(i in seq_along(mergers.nochim)) {
#                         df.nochim$NoReads[i] <- sum(mergers.nochim[[i]]$abundance)
#                 }
#                 
#         }
#         
#         df.plot <- rbind(df.all, df.filtered, df.merged, df.nochim)
#         
#         if (sort) {
#                 # sort always by highest level: all, filtered, merged, nochim
#                 if (!is.null(df.all)) {
#                         df.all <- dplyr::arrange(df.all, desc(NoReads))
#                         LevelsWant <- as.character(df.all$Sample)
#                 } else if (!is.null(df.filtered)) {
#                         df.filtered <- dplyr::arrange(df.filtered, desc(NoReads))
#                         LevelsWant <- as.character(df.filtered$Sample)
#                 } else if (!is.null(df.merged)) {
#                         df.merged <- dplyr::arrange(df.merged, desc(NoReads))
#                         LevelsWant <- as.character(df.merged$Sample)
#                 } else if (!is.null(df.nochim)) {
#                         df.nochim <- dplyr::arrange(df.nochim, desc(NoReads))
#                         LevelsWant <- as.character(df.nochim$Sample)
#                 }
#                 
#                 
#                 df.plot$Sample <- factor(df.plot$Sample)
#                 for (i in seq_along(LevelsWant)) {
#                         df.plot$Sample <- relevel(factor(df.plot$Sample), ref = LevelsWant[i])
#                 }
#                 
#         }
#         
#         
#         ggplot(data = df.plot, aes(x = Sample, y = NoReads, color = Type))  + 
#                 geom_hline(yintercept = 10000, color = 'darkred', linetype = "dashed", size = .35) +
#                 geom_point(size = 2) +
#                 scale_color_manual(values = c(cbPalette[2:4],cbPalette[7])) +
#                 ylab("No Reads") + 
#                 xlab("") +
#                 theme_bw() + 
#                 theme(panel.grid.minor = element_blank(),
#                       panel.grid.major.y = element_blank(),
#                       panel.grid.major.x = element_line(color = "#999999", size = .15),
#                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
#                       legend.title = element_blank())
#         
# }




