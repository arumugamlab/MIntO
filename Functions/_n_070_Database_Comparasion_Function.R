##Database comparasion function

Seq_name=function(tax_table){
  tax_table<-as.data.frame(tax_table)
  seqs<-row.names(tax_table)
  SeqDF<-data.frame("ASVs_No"=(paste("ASV",c(1:length(seqs)), sep="_")), "Sequence"=seqs)
  row.names(SeqDF)<-SeqDF$Sequence
  taxa_table<-merge(SeqDF[-2],tax_table, by= 0, all=TRUE )
  row.names(taxa_table)<-taxa_table$ASVs_No
  colnames(taxa_table)[1]<-"sequence"
  taxa_table[]<-lapply(taxa_table, as.character)
  taxa_table[is.na(taxa_table)] <- "NA"
  return(taxa_table)
}
return.list<-list()
database_comp=function(taxadb1, taxadb2, dbName1, dbName2, warning = FALSE){
 
  try(library(dplyr),install.packages("dplyr"))
  try(library(stringr),install.packages("stringr"))
  library(stringdist)
  return.list<-list()
  
  taxa.db1.names<-Seq_name(taxadb1)
  taxa.db2.names<-Seq_name(taxadb2)
  
  identical<-inner_join(taxa.db1.names, taxa.db2.names, by = c("sequence", "ASVs_No", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
  taxa.db1.names_diff<-anti_join(taxa.db1.names, taxa.db2.names, by = c("sequence", "ASVs_No", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
  taxa.db2.names_diff<-anti_join(taxa.db2.names, taxa.db1.names, by = c("sequence", "ASVs_No", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
  
  taxa_comp<-taxa.db1.names_diff
  for(i in 1:nrow(taxa.db1.names_diff)){
    
    taxa1<-taxa.db1.names_diff[i,]
    taxa2<-taxa.db2.names_diff[i,]
    dif<-(which(!(taxa1==taxa2)))
    
    for(n in dif){taxa_comp[i,n]<- paste(taxa1[n],taxa2[n], sep = " - ") }
    taxa_comp[i,10]<-(colnames(taxa1)[dif[1]])
    
    if(length(grep("NA", taxa_comp[i,dif[1]]))>0){
      if(taxa1[n]== "NA"){taxa_comp[i,11]<-paste("NA_", dbName1,sep="") }else{taxa_comp[i,11]<-paste("NA_", dbName2,sep="")}
    }else{taxa_comp[i,11]<-"Diff"}
    
    if(taxa_comp[i,11]== "Diff"){
      dist<-stringdist(as.character(taxa1[dif[1]]), as.character(taxa2[dif[1]]), weight=c(1,1,1,1))
      taxa_comp[i,12]<-dist}else{taxa_comp[i,12]<-"NA"}
  }
  
  colnames(taxa_comp)[c(10,11,12)]<-c("Different.from", "NA", "Distance")
  
  misspelled<-subset(taxa_comp, as.numeric(unlist(taxa_comp[12]))<=3)
  N.A<-subset(taxa_comp, taxa_comp[12]=="NA")
  different<-subset(taxa_comp[-11], as.numeric(unlist(taxa_comp[12]))>=3)
  
  summary<-data.frame(identical=dim(identical)[1],different= dim(different)[1],misspelled=dim(misspelled)[1], N.A= dim(N.A)[1])
  return.list[[1]]<-identical  #identical taxonomic assignation in both databases
  return.list[[2]]<-different[-dim(different)[2]]  #Different taxonomic assignation
  return.list[[3]]<-misspelled[-dim(misspelled)[2]] #Possible misspelled
  return.list[[4]]<-N.A[-dim(N.A)[2]]        #N.A in one of the database but assigned in the other one 
  return.list[[5]]<-summary
  names(return.list)<-c("identical", "different", "misspelled", "N.A", "summary")
  
  # message("Identical taxonomic assignation in both databases: ", dim(identical)[1], "\n",
  #       "Different taxonomic assignation: ", dim(different)[1],  "\n",
  #       "Possible misspelled: ", dim(misspelled)[1], "\n",
  #       "N.A in one of the database but assigned in the other one: ", dim(N.A)[1]
  #       )
  print(summary)
  return(return.list)
  
}




#####alignment score#################
alignment_Score=function(taxa, database){
taxa<-taxa
fastaFile<-database
taxa[which(taxa[,7]=="sp."), 7] <- NA

# - make a DF with Genus Species of fastaFile -
NameList <- strsplit(names(fastaFile), split = ";")
Class <- lapply(NameList, `[`, 3) %>% unlist()
Order<-  lapply(NameList, `[`, 4) %>% unlist()
Family<-  lapply(NameList, `[`, 5) %>% unlist()
Genus <- lapply(NameList, `[`, 6) %>% unlist()
Species <- lapply(NameList, `[`, 7) %>% unlist()

DF <- data.frame(Index = 1:length(fastaFile), Class=Class, Order=Order,Family=Family, Genus = Genus, Species = Species, Width = width(fastaFile), Sequence = fastaFile)
DF$Genus <- as.character(DF$Genus)
DF$Species <- as.character(DF$Species)
DF$Species[DF$Species == "sp."] <- NA
# --


Last_NotNA_Position <- apply(taxa, 1, function(x){length(which(!is.na(x)))})

taxa_included <- as.matrix(taxa[which(Last_NotNA_Position > 2),])
Last_NotNA_Position <- apply(taxa_included, 1, function(x){length(which(!is.na(x)))})

#taxa_included <- taxa_included[1:9,]

mat <- nucleotideSubstitutionMatrix(match = 5, mismatch = -4, baseOnly = FALSE) # based on ftp://ftp.ncbi.nih.gov/blast/matrices/NUC.4.4 via ebi local alignment tool

ASVSequences <- vector(mode = "character", length = nrow(taxa_included))
tax_levels <- vector(mode = "character", length = nrow(taxa_included))
assigned_taxonomies <- vector(mode = "character", length = nrow(taxa_included))
NoSequences <- vector(mode = "numeric", length = nrow(taxa_included))
MeanPID1 <- vector(mode = "numeric", length = nrow(taxa_included))
MaxPID1<- vector(mode = "numeric", length = nrow(taxa_included))
SDPID1 <- vector(mode = "numeric", length = nrow(taxa_included))
MeanPID3 <- vector(mode = "numeric", length = nrow(taxa_included))
MaxPID3<- vector(mode = "numeric", length = nrow(taxa_included))
SDPID3 <- vector(mode = "numeric", length = nrow(taxa_included))


for (i in 1:nrow(taxa_included)){
        
        ASVSequence <- DNAString(rownames(taxa_included)[i])
        
        current_taxonomy <- colnames(taxa_included)[Last_NotNA_Position[i]]
        assigned_taxonomy <- taxa_included[i, current_taxonomy]
        
        relevantDBSequences <- fastaFile[DF$Index[DF[[current_taxonomy]] == assigned_taxonomy & !is.na(DF[[current_taxonomy]])]]
        
        PID1s <- vector(mode = "numeric", length = length(relevantDBSequences))
        PID3s <- vector(mode = "numeric", length = length(relevantDBSequences))
        
        for (e in 1:length(relevantDBSequences)) {
                
                currentDBSequence <- relevantDBSequences[e]
                palign <- pairwiseAlignment(pattern = ASVSequence, subject = currentDBSequence, type = "local", gapOpening = 10, gapExtension = 0.5, substitutionMatrix = mat)
                
                PID1s[e] <- pid(palign, type = "PID1")
                PID3s[e] <- pid(palign, type = "PID3")
                
                
        }
        
        ASVSequences[i] <- as.character(ASVSequence)
        tax_levels[i] <- current_taxonomy
        assigned_taxonomies[i] <- assigned_taxonomy
        NoSequences[i] <- length(relevantDBSequences)
        MeanPID1[i] <- mean(PID1s)
        MaxPID1[i]<-max(PID1s)
        SDPID1[i] <- sd(PID1s)
        MeanPID3[i] <- mean(PID3s)
        MaxPID3[i]<-max(PID3s)
        SDPID3[i] <- sd(PID3s)
        
}
ResultDF <- data.frame(ASVSequence = ASVSequences, tax_level = tax_levels, assigned_taxonomy = assigned_taxonomies, NoSequences = NoSequences, MeanPID1 = MeanPID1, MaxPID1=MaxPID1,  SDPID1 = SDPID1, MeanPID3 = MeanPID3, MaxPID3=MaxPID3,
                       SDPID3 = SDPID3)
return(ResultDF)
}



DB_overview=function(fastaFile, UNITE=FALSE){
        
        
        fastaFile<-fastaFile
        NameList <- strsplit(names(fastaFile), split = ";")
        if(UNITE){
                Kingdom <- lapply(NameList, `[`, 1) %>% unlist()
                Kingdom <-sapply(strsplit(Kingdom, "__"), `[`, 2)
                Genus <- lapply(NameList, `[`, 6) %>% unlist()
                Genus <-sapply(strsplit(Genus, "__"), `[`, 2)
                Phylum <- lapply(NameList, `[`, 2) %>% unlist()
                Phylum <-sapply(strsplit(Phylum, "__"), `[`, 2)
                Species <- lapply(NameList, `[`, 7) %>% unlist()
                Species <-sapply(strsplit(Species, "__"), `[`, 2)
                Species<-sapply(strsplit(Species, "_"), `[`, 2)
                Species[which(Species=="sp")]<-NA

        }else{
        Kingdom <- lapply(NameList, `[`, 1) %>% unlist()
        Phylum <- lapply(NameList, `[`, 2) %>% unlist()
        Class<- lapply(NameList, `[`, 3) %>% unlist()
        Order<- lapply(NameList, `[`, 4) %>% unlist()
        Family<- lapply(NameList, `[`, 4) %>% unlist()
        Genus <- lapply(NameList, `[`, 6) %>% unlist()
        Species <- lapply(NameList, `[`, 7) %>% unlist()
        Species[which(Species=="sp")]<-NA}
     

        
        DF <- data.frame(Index = 1:length(fastaFile), Kingdom=Kingdom, Phylum= Phylum,Class=Class, Order=Order, Family=Family,Genus = Genus, Species = Species, Width = width(fastaFile), Sequence = fastaFile)
        
        
        Taxonomy<-data.frame(Kingdom=length(unique(DF$Kingdom)), Phylum= length(unique(DF$Phylum)), Genus=length(unique(DF$Genus)), Species=length(unique(DF$Species)))
        
        Tr <- ggplot(DF, aes(x =Width))
        length_Phylum_plot<-Tr + geom_histogram(aes(color = Phylum, fill = Phylum), bins = 100,  alpha =0.4) + theme_minimal() +theme(legend.position = "top") # + scale_fill_manual(values = c("#00AFBB", "#E7B800","#FC4E07")) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) 
        
        length_Kingdom_plot<-Tr + geom_histogram(aes(color = Kingdom, fill = Kingdom), bins = 100,  alpha =0.4) + theme_minimal() +theme(legend.position = "top")
        
 ######
        
        sequencesBefore <- nrow(DF)
        
        IndexesCompleteDuplicates <- which(duplicated(DF[,-1]))
        
        DF_unique <- DF[-IndexesCompleteDuplicates,]
        
        
        # - This was just for later changing the ambiguous species annotations -
        ambiguous_species <- DF_unique[which(duplicated(DF_unique$Sequence) | duplicated(DF_unique$Sequence, fromLast = TRUE)),]
        # --
        
        Indexes_ambiguous_species <- which(duplicated(DF_unique$Sequence))
        
        DF_unique <- DF_unique[-Indexes_ambiguous_species,]
        lens <- width((fastaFile))
        Summary_DB<-data.frame(nSequences=sequencesBefore, 
                               Duplicate= length(IndexesCompleteDuplicates), 
                               Ambiguous_species=length(Indexes_ambiguous_species), 
                               min_length=min(lens), 
                               max_length=max(lens), 
                               Mean_length=round(mean(lens, na.rm=TRUE),2),
                               Median_length=median(lens, na.rm=TRUE),
                               sd_length=round(sd(lens, na.rm=TRUE),3))
        # - remove duplicates (i.e. entries with identical sequences) -
        fastaFile <- unique(fastaFile)
        
        
        # - illustrate how many sequences and unique species are there for each genus -
        DF_uniqueSpecies <- DF_unique[!duplicated(DF_unique$Species), ]
        NoSpeciesBefore <- summarise(group_by(DF_uniqueSpecies, Genus), No_Species = n()) %>% complete(Genus, fill = list(No_Species = 0))
        NoSequences <- summarise(group_by(DF_unique, Genus), No_Sequences = n()) %>% complete(Genus, fill = list(No_Sequences = 0))
        
        NoSpeciesBefore <- merge(NoSequences, NoSpeciesBefore, by = "Genus") %>% arrange(Genus)
        NoSpeciesBeforeShow <- rbind(NoSpeciesBefore, data.frame("Genus" = "Total", "No_Sequences" = sum(NoSpeciesBefore$No_Sequences),
                                                                 "No_Species" = sum(NoSpeciesBefore$No_Species)))
        
        
        
        NoSequences_Species <- summarise(group_by(DF_unique, Species), No_Sequences = n()) %>% complete(Species, fill = list(No_Sequences = 0))
        NoSequences_Species<-NoSequences_Species[order(-NoSequences_Species$No_Sequences),]
        
        ressult.list<-list(DF=DF, 
                           length_Phylum_plot=length_Phylum_plot,
                           length_Kingdom_plot=length_Kingdom_plot,
                           ambiguous_species=ambiguous_species, 
                           Summary_DB=Summary_DB,
                           DF_uniqueSpecies=DF_uniqueSpecies,
                           NoSpeciesBeforeShow=NoSpeciesBeforeShow,
                           NoSequences_Species=NoSequences_Species,
                           fastaFileUnique=fastaFile, 
                           Taxonomy=Taxonomy)
        return(ressult.list)
        
}

##ASV_length
ASV_length_Plots=function(seqtab.nochim, taxa){
        seqtab.nochim<-seqtab.nochim
        taxa<-taxa
        length1<-as.data.frame(nchar(getSequences(seqtab.nochim)))
        colnames(length1)<-"length_ASV"
        row.names(length1)<-colnames(seqtab.nochim)
        
        length<-merge(length1, taxa, by="row.names")
        length[3:9] <- lapply(length[3:9], as.character)
        
        length2<-merge(length1, taxa, by="row.names")
        length2[3:9] <- lapply(length[3:9], as.character)
        
        length[is.na(length)] <-"NA"
        a <- ggplot(length, aes(x = length_ASV))
        length_kingom_plot<-a + geom_histogram(aes(color = Kingdom, fill = Kingdom), bins = 80,  alpha =0.4) + theme_minimal() +theme(legend.position = "top")  + scale_fill_manual(values = c("#00AFBB", "#E7B800","#FC4E07")) + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) 
        
        
        length_fungi<-subset(length2,length2$Kingdom =="Fungi")
        
        b <- ggplot(length_fungi, aes(x = length_ASV))
        length_fungi_plot<-b + geom_histogram(aes(color = Phylum, fill = Phylum), bins = 80,  alpha =0.4) + theme_minimal()+theme(legend.position = "top")  + scale_fill_manual(values = wes_palette(n = 5, name = "Cavalcanti1"))  + scale_color_manual(values = wes_palette(n = 5, name = "Cavalcanti1"))
        
       
       Kingdom_tab<- table(taxa[,"Kingdom"])
       Phylum_tab<-as.data.frame(table(length_fungi[,"Phylum"]))
        Summary_taxonomy<-data.frame(No.Kingdom=dim(length_fungi)[1],
                                      No.Phylum=length(which(!is.na(length_fungi$Phylum))), Phylum_NA=length(which(is.na(length_fungi$Phylum))),
                                      No.Genus=length(which(!is.na(length_fungi$Genus))), Genus_NA=length(which(is.na(length_fungi$Genus))),
                                      No.Species=length(which(!is.na(length_fungi$Species))), Species_NA=length(which(is.na(length_fungi$Species))) )        
        
        Summary_taxonomy<-as.data.frame(t(Summary_taxonomy))
        Summary_taxonomy$Taxonomy<-row.names(Summary_taxonomy)
        colnames(Summary_taxonomy)<-c("No.ASVs", "Taxonomy")
        Summary_taxonomy<-Summary_taxonomy %>%select(Taxonomy,No.ASVs )
        
        
        Summary_taxonomy2<-data.frame(Fungi=dim(length_fungi)[1],
                                     Phylum=length(which(!is.na(length_fungi$Phylum))),
                                     Genus=length(which(!is.na(length_fungi$Genus))),
                                     Species=length(which(!is.na(length_fungi$Species))))        
        
        Summary_taxonomy2<-as.data.frame(t(Summary_taxonomy2))
        Summary_taxonomy2$Taxonomy<-row.names(Summary_taxonomy2)
        colnames(Summary_taxonomy2)<-c("No.ASVs", "Taxonomy")
        Summary_taxonomy2<-Summary_taxonomy2 %>%select(Taxonomy,No.ASVs )
        
        output<-list(kingdom=length_kingom_plot, Phylum=length_fungi_plot,Kingdom_tab=Kingdom_tab, Phylum_tab=Phylum_tab, Summary_taxonomy=Summary_taxonomy, Summary_taxonomy2=Summary_taxonomy2 )
        return(output)
        
        }

Blast_taxonomy=function(Blast_tab=Blast_tab, 
                        Phylum_threshold=Phylum_threshold, 
                        Class_threshold=Class_threshold,
                        Order_threshold=Order_threshold,  
                        Family_threshold=Family_threshold,  
                        Genus_threshold=Genus_threshold,   
                        Species_threshold=Species_threshold){
        
        #- add PI---        
        
        BLAST$PI2<-round((BLAST$Pidentity*(BLAST$alignment.length)))/BLAST$length # hoping nominator is no_matches
        BLAST$PI2<-(BLAST$Pidentity*(BLAST$alignment.length-BLAST$mismatches))/BLAST$length      
        #---
        
        # - Split Taxonomy -
        
        taxaNames<-as.character(BLAST$subject)
        Taxonomy <- strsplit((taxaNames), split = ";")
        BLAST$Kingdom <- lapply(Taxonomy, `[`, 1) %>% unlist()
        BLAST$Phylum <- lapply(Taxonomy, `[`, 2) %>% unlist()
        BLAST$Class <- lapply(Taxonomy, `[`, 3) %>% unlist()
        BLAST$Order <- lapply(Taxonomy, `[`, 4) %>% unlist()
        BLAST$Family<-lapply(Taxonomy, `[`, 5) %>% unlist()
        BLAST$Genus <- lapply(Taxonomy, `[`, 6) %>% unlist()
        BLAST$Species <- lapply(Taxonomy, `[`, 7) %>% unlist()
        
        BLAST_tax<-BLAST %>% select(query, seq, Pidentity,alignment.length,length, PI2, Kingdom:Species)
        
        BLAST_tax$Phylum[which(BLAST_tax$Phylum=="NA")]<-NA
        BLAST_tax$Class[which(BLAST_tax$Class=="NA")]<-NA
        BLAST_tax$Order[which(BLAST_tax$Order=="NA")]<-NA
        BLAST_tax$Family[which(BLAST_tax$Family=="NA")]<-NA
        BLAST_tax$Genus[which(BLAST_tax$Genus=="NA")]<-NA
        BLAST_tax$Species[which(BLAST_tax$Species=="NA")]<-NA
        # --
        
        
        # - phylum level threshold we have to decide on mamasita
        BLAST_tax_filter = dplyr::filter(BLAST_tax, PI2 >= Phylum_threshold) 
        summary_tax<-summarise(group_by(BLAST_tax_filter, query), hits = n())
        # --
        
        tax_blast<-data.frame()
        tax_blast_Amb_genus<-data.frame()
        
        
        for(i in summary_tax$query ){
                
                ASV <-dplyr::filter(BLAST_tax_filter,query==i)
                
                # - Kingdom assignment -
                
                summary_ASV<-summarise(group_by(ASV, Kingdom), hits = n(), MaxPID=max(PI2), MedianPID = median(PI2)) %>% arrange(desc(MaxPID), hits)
                
                Kingdom_assignment = summary_ASV$Kingdom[1]
                
                # - phylum assignment -
                ASV_phylum <-dplyr::filter(ASV, Kingdom == Kingdom_assignment, PI2 >= Phylum_threshold)
                
                summary_ASV<-summarise(group_by(ASV_phylum, Phylum), hits = n(), MaxPID=max(PI2), MedianPID = median(PI2)) %>% arrange(desc(MaxPID), hits)
                
                phylum_assignment = summary_ASV$Phylum[1]
                # --
                
                # - class assignment -
                ASV_class <-dplyr::filter(ASV, Phylum == phylum_assignment, PI2 >= Class_threshold)
                summary_ASV<-summarise(group_by(ASV_class, Class), hits = n(), MaxPID=max(PI2), MedianPID = median(PI2)) %>% arrange(desc(MaxPID), hits)
                
                class_assignment = summary_ASV$Class[1]
                # --
                
                
                # - order assignment -
                ASV_order <-dplyr::filter(ASV_class, Class == class_assignment, PI2 >= Order_threshold)
                
                summary_ASV<-summarise(group_by(ASV_order, Order), hits = n(), MaxPID=max(PI2), MedianPID = median(PI2)) %>% arrange(desc(MaxPID), hits)
                
                order_assignment = summary_ASV$Order[1]
                # --
                
                # - family assignment -
                ASV_family <-dplyr::filter(ASV_order, Order == order_assignment, PI2 >= Family_threshold)
                
                summary_ASV<-summarise(group_by(ASV_family, Family), hits = n(), MaxPID=max(PI2), MedianPID = median(PI2)) %>% arrange(desc(MaxPID), hits)
                
                family_assignment = summary_ASV$Family[1]
                
                
                # --
                
                
                # - genus assignment -
                ASV_genus <-dplyr::filter(ASV_family, Family == family_assignment, PI2 >= Genus_threshold)
                
                summary_ASV<-summarise(group_by(ASV_genus, Genus), hits = n(), MaxPID=max(PI2), MedianPID = median(PI2)) %>% arrange(desc(MaxPID), hits)
                
                genus_assignment = summary_ASV$Genus[1]
                
                if(dim(summary_ASV)[1]>1){
                        genus_assignment_Amb<-paste(summary_ASV$Genus, collapse = '/')
                        Genus_PID_Amb <- paste(round(summary_ASV$MaxPID), collapse = '/')
                        
                        ASV_tax_amb<-data.frame(query=i, Sequence= ASV$seq[1], Genus_PID_Amb=Genus_PID_Amb, Kingdom= Kingdom_assignment, Phylum=phylum_assignment, Class=class_assignment, Order=order_assignment, Family=family_assignment, Genus=genus_assignment_Amb, Species=NA)
                        tax_blast_Amb_genus<- rbind(tax_blast_Amb_genus,ASV_tax_amb) 
                }
                # --
                
                
                # - species assignment -
                ASV_species <- dplyr::filter(ASV_genus, Genus == genus_assignment, PI2 >= Species_threshold)
                
                if(dim(ASV_species)[1]!=0)
                {summary_ASV<-summarise(group_by(ASV_species, Species), hits = n(), MaxPID=max(PI2), MedianPID = median(PI2)) %>% arrange(desc(MaxPID), hits)
                
                species_assignment <- paste(summary_ASV$Species, collapse = '/')
                species_PID <- paste(round(summary_ASV$MaxPID), collapse = '/')}else{species_assignment<-NA
                species_PID<-NA
                }
                
                # --
                
                # - Final Taxonomy -
                ASV_tax<-data.frame(query=i, Sequence= ASV$seq[1], species_PID=species_PID, Kingdom= Kingdom_assignment, Phylum=phylum_assignment, Class=class_assignment, Order=order_assignment, Family=family_assignment, Genus=genus_assignment, Species=species_assignment)
                
                tax_blast<-rbind(tax_blast,ASV_tax)
                
                
        }
        
        taxa_blast_final<-tax_blast
        row.names(taxa_blast_final)<-taxa_blast_final$Sequence
        
        taxa_blast_final<-taxa_blast_final[order(taxa_blast_final$query),]
        taxa_blast_final<-taxa_blast_final[-c(1:3)]
        
        
        
        
        
        result_list<- list(tax_blast_info=tax_blast, tax_blast_Amb_genus=tax_blast_Amb_genus, BLAST_Tax_complete= BLAST_tax, taxa_blast=taxa_blast_final)
        
}



database_comp_blast=function(taxadb1, taxadb2, dbName1, dbName2, warning = FALSE){
        
        try(library(dplyr),install.packages("dplyr"))
        try(library(stringr),install.packages("stringr"))
        library(stringdist)
        return.list<-list()
        
        taxa.db1.names<-taxadb1
        taxa.db2.names<-taxadb2
        
        identical<-inner_join(taxa.db1.names, taxa.db2.names, by = c("sequence", "ASVs_No", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
        taxa.db1.names_diff<-anti_join(taxa.db1.names, taxa.db2.names, by = c("sequence", "ASVs_No", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
        taxa.db2.names_diff<-anti_join(taxa.db2.names, taxa.db1.names, by = c("sequence", "ASVs_No", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
        
        taxa_comp<-taxa.db1.names_diff
        for(i in 1:nrow(taxa.db1.names_diff)){
                
                taxa1<-taxa.db1.names_diff[i,]
                taxa2<-taxa.db2.names_diff[i,]
                dif<-(which(!(taxa1==taxa2)))
                
                for(n in dif){taxa_comp[i,n]<- paste(taxa1[n],taxa2[n], sep = " - ") }
                taxa_comp[i,10]<-(colnames(taxa1)[dif[1]])
                
                if(length(grep("NA", taxa_comp[i,dif[1]]))>0){
                        if(taxa1[n]== "NA"){taxa_comp[i,11]<-paste("NA_", dbName1,sep="") }else{taxa_comp[i,11]<-paste("NA_", dbName2,sep="")}
                }else{taxa_comp[i,11]<-"Diff"}
                
                if(taxa_comp[i,11]== "Diff"){
                        dist<-stringdist(as.character(taxa1[dif[1]]), as.character(taxa2[dif[1]]), weight=c(1,1,1,1))
                        taxa_comp[i,12]<-dist}else{taxa_comp[i,12]<-"NA"}
        }
        
        colnames(taxa_comp)[c(10,11,12)]<-c("Different.from", "NA", "Distance")
        
        misspelled<-subset(taxa_comp, as.numeric(unlist(taxa_comp[12]))<=3)
        N.A<-subset(taxa_comp, taxa_comp[12]=="NA")
        different<-subset(taxa_comp[-11], as.numeric(unlist(taxa_comp[12]))>=3)
        
        summary<-data.frame(identical=dim(identical)[1],different= dim(different)[1],misspelled=dim(misspelled)[1], N.A= dim(N.A)[1])
        return.list[[1]]<-identical  #identical taxonomic assignation in both databases
        return.list[[2]]<-different[-dim(different)[2]]  #Different taxonomic assignation
        return.list[[3]]<-misspelled[-dim(misspelled)[2]] #Possible misspelled
        return.list[[4]]<-N.A[-dim(N.A)[2]]        #N.A in one of the database but assigned in the other one 
        return.list[[5]]<-summary
        names(return.list)<-c("identical", "different", "misspelled", "N.A", "summary")
        
        # message("Identical taxonomic assignation in both databases: ", dim(identical)[1], "\n",
        #       "Different taxonomic assignation: ", dim(different)[1],  "\n",
        #       "Possible misspelled: ", dim(misspelled)[1], "\n",
        #       "N.A in one of the database but assigned in the other one: ", dim(N.A)[1]
        #       )
        print(summary)
        return(return.list)
        
}






ASV_summary=function(ASV_ID, blast_tax_complete){
        
    i<-ASV_ID
    BLAST_tax_filter <- blast_tax_complete         
    ASV <-dplyr::filter(BLAST_tax_filter,query==i)
    Phylum_threshold=50
    Class_threshold=60
    Order_threshold=70
    Family_threshold=80
    Genus_threshold=90
    Species_threshold=97            
                # - Kingdom assignment -
                
                summary_ASV_kingdom<-summarise(group_by(ASV, Kingdom), hits = n(), MaxPID=max(PI2), MedianPID = median(PI2)) %>% arrange(desc(MaxPID), hits)
                
                Kingdom_assignment = summary_ASV_kingdom$Kingdom[1]
                
                # - phylum assignment -
                ASV_phylum <-dplyr::filter(ASV, Kingdom == Kingdom_assignment, PI2 >= Phylum_threshold)
                
                summary_ASV_phylum<-summarise(group_by(ASV_phylum, Phylum), hits = n(), MaxPID=max(PI2), MedianPID = median(PI2)) %>% arrange(desc(MaxPID), hits)
                
                phylum_assignment = summary_ASV_phylum$Phylum[1]
                # --
                
                # - class assignment -
                ASV_class <-dplyr::filter(ASV, Phylum == phylum_assignment, PI2 >= Class_threshold)
                summary_ASV_Class<-summarise(group_by(ASV_class, Class), hits = n(), MaxPID=max(PI2), MedianPID = median(PI2)) %>% arrange(desc(MaxPID), hits)
                
                class_assignment = summary_ASV_Class$Class[1]
                # --
                
                
                # - order assignment -
                ASV_order <-dplyr::filter(ASV_class, Class == class_assignment, PI2 >= Order_threshold)
                
                summary_ASV_order<-summarise(group_by(ASV_order, Order), hits = n(), MaxPID=max(PI2), MedianPID = median(PI2)) %>% arrange(desc(MaxPID), hits)
                
                order_assignment = summary_ASV_order$Order[1]
                # --
                
                # - family assignment -
                ASV_family <-dplyr::filter(ASV_order, Order == order_assignment, PI2 >= Family_threshold)
                
                summary_ASV_family<-summarise(group_by(ASV_family, Family), hits = n(), MaxPID=max(PI2), MedianPID = median(PI2)) %>% arrange(desc(MaxPID), hits)
                
                family_assignment = summary_ASV_family$Family[1]
                
                
                # --
                
                
                # - genus assignment -
                ASV_genus <-dplyr::filter(ASV_family, Family == family_assignment, PI2 >= Genus_threshold)
                
                summary_ASV_genus<-summarise(group_by(ASV_genus, Genus), hits = n(), MaxPID=max(PI2), MedianPID = median(PI2)) %>% arrange(desc(MaxPID), hits)
                
                genus_assignment = summary_ASV_genus$Genus[1]
                
                
                # - species assignment -
                ASV_species <- dplyr::filter(ASV_genus, Genus == genus_assignment, PI2 >= Species_threshold)
                
              
                summary_ASV_species<-summarise(group_by(ASV_species, Species), hits = n(), MaxPID=max(PI2), MedianPID = median(PI2)) %>% arrange(desc(MaxPID), hits)

        result_list<- list(ASV=ASV,
                           summary_ASV_kingdom=summary_ASV_kingdom, 
                           summary_ASV_phylum=summary_ASV_phylum,
                           summary_ASV_Class=summary_ASV_Class,
                           summary_ASV_order=summary_ASV_order,
                           summary_ASV_family=summary_ASV_family,
                           summary_ASV_genus=summary_ASV_genus, 
                           summary_ASV_species=summary_ASV_species)
        return(result_list)
        
}


## Primer_Check

PrimerCheck=function(FastaFile,FPrimer, RPrimer ){
        
        #Functions
        functionpath <-functionpath <- "~/Dada_Pipel-master/Functions"
        source(file.path(functionpath, "Functions_qPCRDesign.R"))    
        #--
        
        fastaFile<-FastaFile
        
        NameList <- strsplit(names(fastaFile), split = ";")
        Genus <- lapply(NameList, `[`, 6) %>% unlist()
        Species <- lapply(NameList, `[`, 7) %>% unlist()
        
        DF <- data.frame(Index = 1:length(fastaFile), Genus = Genus, Species = Species, Width = width(fastaFile), Sequence = fastaFile)
        
        #--
        
        ## Test your primers: To how many sequences is there a perfect match
        
        # PrimerF<- "GCATCGATGAAGAACGCAGC"
        # PrimerR <- "TCCTCCGCTTATTGATATGC"
        
        
        ForwardPrimer <- c(PrimerF = FPrimer)
        ReversePrimer <- c(PrimerR = RPrimer)
        
        PrimerToTest <- DNAString(ForwardPrimer)
        HitsForward <- vcountPattern(pattern = PrimerToTest, subject = fastaFile)
        
        PrimerToTest <- reverseComplement(DNAString(ReversePrimer)) # NB: need reverse complement here
        HitsReverse <- vcountPattern(pattern = PrimerToTest, subject = fastaFile)
        
        
        
        SummaryDF <- cbind(DF, data.frame(Forward = HitsForward, Reverse = HitsReverse))
        
        SummaryDF$Amplicon <- 0
        SummaryDF$Amplicon[SummaryDF$Forward > 0 & SummaryDF$Reverse > 0] <- 1
        
        SummaryDF_sum <- group_by(SummaryDF, Genus) %>% summarise(Total = n(), Forward = sum(Forward), Reverse = sum(Reverse), Amplicon = sum(Amplicon)) %>% arrange(desc(Amplicon))
        
        
        
        # print(xtable(SummaryDF_sum, align = "|c|c|c|c|c|", digits = 0), include.rownames = FALSE)
        
        
        
        # - This plot might actually be ok useful to get an overview of the data base -
        
        # SummaryDF_sum <- filter(SummaryDF_sum, Forward > 0 | Reverse > 0 | Amplicon > 0)
        
        SummaryDF_sum$Genus <- factor(SummaryDF_sum$Genus, levels = rev(SummaryDF_sum$Genus), ordered = T)
        
        SummaryDF_sum_PC <- mutate(SummaryDF_sum, Total = Total/sum(Total), Forward = Forward/sum(Forward), Reverse = Reverse/sum(Reverse), Amplicon = Amplicon/sum(Amplicon))
        
        SummaryDF_sum <- rbind(SummaryDF_sum, SummaryDF_sum_PC)
        SummaryDF_sum$Type <- "sequences"
        SummaryDF_sum$Type[(nrow(SummaryDF_sum_PC)+1):(2*nrow(SummaryDF_sum_PC))] <- "percentage"
        
        SummaryDFPlot <- gather(SummaryDF_sum, key = "Primer", value = "Value", Total:Amplicon)
        
        SummaryDFPlot$Primer <- factor(SummaryDFPlot$Primer, levels = c("Total", "Forward", "Reverse", "Amplicon"), ordered = T)
        SummaryDFPlot$Type <- factor(SummaryDFPlot$Type, levels = c("sequences", "percentage"), ordered = T)
        
        Genus <- levels(SummaryDF_sum$Genus)
        if (length(Genus) < 16) {
                GenusColors <- rev(QuantColors15[1:length(Genus)])       
        } else if (length(Genus) < 19) {
                GenusColors <- tol18rainbow[1:length(Genus)]
        } else {
                GenusColors <- viridis(length(Genus))
        }
        
        names(GenusColors) <- Genus 
        
        
        Tr <- ggplot(SummaryDFPlot, aes(x = Primer, y = Value, fill = Genus))
        Tr <- Tr +
                geom_bar(stat = "identity", position = "stack") +
                facet_wrap(~ Type, scales = "free_y") +
                theme_bw() +
                scale_fill_manual("", values = GenusColors) +
                xlab("") +
                ylab("Number sequences colored by Genus") +
                theme(legend.position = "none")
        
        output<-list(plot_Primer=Tr, SummaryDFPlot=SummaryDFPlot, SummaryDF_sum=SummaryDF_sum)
        
      return(output)  
        
}
