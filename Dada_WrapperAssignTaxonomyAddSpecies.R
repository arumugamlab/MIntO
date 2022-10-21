
working_dir <- "data/continuous_culture/"

#install.packages('phangorn')
load(paste0(working_dir,"Dada_Data/DenoisedData.RData"))
functpath = "Functions/"
source(file.path(functpath, "Dada_TaxonomyFunctions.R"))

# Silva 
assignTaxonomyaddSpecies(seqtab = seqtab.nochim, 
                         minBoot = 80, #50,
                         allowMultiple = TRUE,
                         PathToRefs = "Silva_DB/",
                         RefDataBase = "silva_nr99_v138.1_train_set.fa.gz", #"silva_nr99_v138.1_wSpecies_train_set.fa.gz", # using silva training set with or without species names
                         SpeciesDB = "silva_species_assignment_v138.1.fa.gz",
                         PathToSave =paste0(working_dir,"Dada_Silva_80/"), #Dada_Silva_sp_80/ # using silva training set with or without species names
                         tryRC = FALSE)

# RDP 
assignTaxonomyaddSpecies(seqtab = seqtab.nochim, 
                         minBoot = 80, #50,
                         allowMultiple = TRUE,
                         PathToRefs = "RDP_DB/",
                         RefDataBase = "rdp_train_set_18.fa.gz",
                         SpeciesDB = "rdp_species_assignment_18.fa.gz",
                         PathToSave =paste0(working_dir,"Dada_RDP_80/"), 
                         tryRC = FALSE)


# Silva
load(paste0(working_dir,"Dada_Silva_80/Taxonomy.RData"))
taxa.species_silva_df <- as.data.frame(taxa.species)
taxa.species_silva_difficile <- rownames(taxa.species_silva_df[which(taxa.species_silva_df$Species == 'difficile'),])

# RDP
load(paste0(working_dir,"Dada_RDP_80/Taxonomy.RData"))
taxa.species_rdp_df <- as.data.frame(taxa.species)
taxa.species_rdp_scindens <- rownames(taxa.species_rdp_df[which(taxa.species_rdp_df$Species == 'scindens'),])

# View 16S sequence matching diffcile or scindens annotations
View(taxa.species_rdp_df[which(rownames(taxa.species_rdp_df) == taxa.species_silva_difficile),])
View(taxa.species_silva_df[which(rownames(taxa.species_silva_df) == taxa.species_rdp_scindens),])
