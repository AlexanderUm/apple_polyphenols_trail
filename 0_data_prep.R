#-------------------------------------------------------------------------------
# Load objects 
#-------------------------------------------------------------------------------
load("prm.Rdata")

set.seed(prm.ls$general$seed)

PrevalenceTables <- list()


#-------------------------------------------------------------------------------
# Variables
#-------------------------------------------------------------------------------
# Path
SeqPath <- prm.ls[["data_prep"]][["seq_path"]]

MetaPath <- prm.ls[["data_prep"]][["meta_path"]]

DirOut <- prm.ls[["general"]][["DirOut"]]

# Columns 
ColSeqID <- prm.ls[["general"]][["SeqID_col"]]

ColSampleID <- prm.ls[["general"]][["SampleID_col"]]

ColGroup <- prm.ls[["general"]][["Group_col"]]

ColsToFactors <- prm.ls[["data_prep"]][["Meta_Cols_to_factors"]]


# Filtering
ReadsMin <- prm.ls[["data_prep"]][["min_reads_tax"]]

RareDepth <- prm.ls[["data_prep"]][["rare_depth"]]

TaxaKingdomKeep <- prm.ls[["data_prep"]][["tax_king_keep"]]

TaxaToRemove <- prm.ls[["data_prep"]][["tax_kick"]]

TaxaLevelsGlom <- prm.ls[["data_prep"]][["tax_levels"]]

SamplesToKick <- prm.ls$data_prep$samples_kick


#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------
source("R/phy_shorten_tax_names.R")
source("R/phy_norm.R")
source("R/phy_prevalence_calc.R")

#-------------------------------------------------------------------------------
# Data import
#-------------------------------------------------------------------------------
# Import data into phyloseq 
Ps0 <- qza_to_phyloseq(features = paste0(SeqPath, "/asv_table.qza"), 
                        tree = paste0(SeqPath, "/tree/rooted-tree.qza"), 
                        taxonomy = paste0(SeqPath, "/taxonomy_07.qza"))

# Read in metadata and select samples with highest reads count 
MetaData <- read.csv(MetaPath) 
              
PsMeta <- MetaData %>% 
                mutate(across(all_of(ColsToFactors), 
                              function(x){factor(x, levels = unique(x))}), 
                       ColNames = .data[[ColSeqID]], 
                       SampleID = paste0("S", SampleID)) %>% 
                mutate(TimeNumeric = case_match(TimePoint, 1~0, 2~6, 3~12)) %>% 
                mutate(TimeFactor = paste0("Week ", TimeNumeric)) %>% 
                mutate(TimeFactor = factor(TimeFactor, 
                                           levels = unique(TimeFactor))) %>% 
                column_to_rownames(var = "ColNames") %>% 
                arrange(match(.data[[ColSeqID]], sample_names(Ps0))) %>% 
                mutate(nReads = sample_sums(Ps0)) %>% 
                slice(which.max(nReads), .by = all_of(ColSampleID))

if(!is.null(SamplesToKick)) {
  
  for (i in names(SamplesToKick)) {
    
    PsMeta <- PsMeta %>% 
                filter(!.data[[i]] %in% SamplesToKick[[i]])
    
  }
  
}


# Add metadata to phyloseq object
Ps1 <- prune_samples(rownames(PsMeta), Ps0)

PsMeta <- PsMeta[sample_names(Ps1), ]

sample_data(Ps1) <- PsMeta

# Write out metadata 
write.csv(PsMeta, 
          file = paste0(DirOut, "/samples_data.csv"))

#-------------------------------------------------------------------------------
# Prune taxa 
#-------------------------------------------------------------------------------
# Taxa with less than 10 reads in total   
Ps1 <- prune_taxa(taxa_sums(Ps1) > ReadsMin, Ps1)

# Keep ASVs
Ps1 <- prune_taxa(tax_table(Ps1)[, "Kingdom"] %in% TaxaKingdomKeep, Ps1)

# Kick ASVs
for(i in 1:length(TaxaToRemove)) {
  
  InstTaxaToRemove <- TaxaToRemove[i]
  
  Ps1 <- prune_taxa(!tax_table(Ps1)[, names(InstTaxaToRemove)] %in% InstTaxaToRemove, Ps1)
  
}


#-------------------------------------------------------------------------------
# Tax glom and count transformation
#-------------------------------------------------------------------------------
PsLs <- list()

for(i in TaxaLevelsGlom)  {
  
  if(i == "ASV") {InstPs <- Ps1
  
  # Adjust taxa names
  taxa_names(InstPs) <- phy_shorten_tax_names(InstPs) %>% 
                              as.data.frame() %>% 
                              setNames("feature") %>% 
                              mutate(feature2 = if(n() > 1) {
                                 paste0(feature, "__asv", row_number())} else {
                                  paste0(feature, "__asv")}, .by = "feature") %>% 
                              pull(feature2)
  
  } else { 
    
    InstPs <- tax_glom(Ps1, i, NArm = FALSE) 
  
    # Adjust taxa names
    taxa_names(InstPs) <- phy_shorten_tax_names(InstPs) %>% 
                                as.data.frame() %>% 
                                setNames("feature") %>% 
                                mutate(feature2 = if(n() > 1) {
                                  paste0(feature, "__", 
                                         tolower(str_sub(i, 1, 1)), 
                                         row_number())} else {
                                  feature}, .by = "feature") %>% 
                                pull(feature2)
  }
  
 
  
  # Calculate prevalence per taxa 
  PrevalenceTables[[i]] <- phy_prevalence_calc(InstPs, 
                                    per_group_cols = ColGroup)
  
  # Transform count 
  InstPsRare <- rarefy_even_depth(InstPs, 
                               sample.size = RareDepth,
                               rngseed = prm.ls[["general"]][["seed"]])
  
  InstPsRareLog <- transform_sample_counts(InstPsRare, function(x) {log2(x + 1)})
  
  PsPrecent <- transform_sample_counts(InstPs, function(x) {x/sum(x)*100})
  
  PsLs[[i]] <- list("Raw" = InstPs, 
                    "Rare" = InstPsRare, 
                    "Rare_log2" = InstPsRareLog,
                    "CSS" = phy_css_norm(InstPs), 
                    "CLR" = phy_clr_norm(InstPs), 
                    "TSS" = phy_tss_norm(InstPs, log2_transform = FALSE), 
                    "Rel" = PsPrecent)
}


#-------------------------------------------------------------------------------
# Color schema 
#-------------------------------------------------------------------------------
# Groups color
ColorSchLs <- list()
  
ColorSchLs[[ColGroup]] <- brewer.pal(7, "Dark2")[1:length(levels(PsMeta[[ColGroup]]))] %>% 
                                    setNames(levels(PsMeta[[ColGroup]]))
  
ColorSchLs[["TimeFactor"]] <- brewer.pal(7, "Set1")[1:length(levels(PsMeta[["TimeFactor"]]))] %>% 
                                   setNames(levels(PsMeta[["TimeFactor"]]))


#-------------------------------------------------------------------------------
# Write data
#-------------------------------------------------------------------------------
dir.create(paste0(DirOut, "/supp"), recursive = TRUE, showWarnings = FALSE)

save(list = c("PsLs", "PsMeta", "MetaData", 
              "ColorSchLs", "PrevalenceTables"),
     file = paste0(DirOut, "/supp/0_data.RData"))

rm(list = ls())
