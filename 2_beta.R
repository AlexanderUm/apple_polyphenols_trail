################################################################################
# Beta diversity
################################################################################
#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
load("prm.Rdata")

set.seed(prm.ls$general$seed)

DirOut <- prm.ls[["general"]][["DirOut"]]


# Metadata columns 
ColGroup <- prm.ls[["general"]][["Group_col"]]

ColTimeNumeric <- prm.ls[["general"]][["TimeNumeric_col"]]

ColParticipantID <- prm.ls[["general"]][["PraticipantID_col"]]


# Phyloseq
TaxaLvl <- prm.ls[["beta"]][["Tax_lvl"]]

Norm <- prm.ls[["beta"]][["Norm"]]


# Beta
BetaDistances <- prm.ls[["beta"]][["distances"]]

# Adonis
Nperm <- prm.ls[["beta"]][["n_perm"]]

AdonisFormula <- prm.ls[["beta"]][["AdonisFormula"]]


# Directories
dir.create(paste0(DirOut, "/beta/plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(DirOut, "/beta/tabs"), recursive = TRUE, showWarnings = FALSE)


# Custom functions
source("R/phy_dists_ls.R")
source("R/plot_PCoAs.R")
source("R/rda_extract_for_plot.R")
source("R/plot_extracted_rda_data.R")


# Data 
load(paste0(DirOut, "/supp/0_data.Rdata"))

# Objects 
ResBetaLs <- list()


#-------------------------------------------------------------------------------
# Calculate distances 
#-------------------------------------------------------------------------------
Par <- expand.grid("TaxLvl" = TaxaLvl, 
                   "Norm" = Norm, 
                   stringsAsFactors = FALSE)

DistsLs <- list()

for(i in 1:nrow(Par)) {
  
  PsInst <- PsLs[[Par$TaxLvl[i]]][[Par$Norm[i]]]
  
  InstDist <- phy_dist_ls(PsInst, dists = BetaDistances) %>% 
                setNames(names(BetaDistances))
  
  DistsLs[[Par$TaxLvl[i]]][[Par$Norm[i]]] <- InstDist
  
}


#-------------------------------------------------------------------------------
# Test significance with Adonis2
#-------------------------------------------------------------------------------
Par <- expand.grid("TaxLvl" = TaxaLvl, 
                   "Norm" = Norm, 
                   "Form" = AdonisFormula,
                   "Dist" = names(BetaDistances),
                    stringsAsFactors = FALSE)

AdonisRes <- NULL

for(i in 1:nrow(Par)) {
  
  InstDist <- DistsLs[[Par$TaxLvl[i]]][[Par$Norm[i]]][[Par$Dist[i]]]
  
  InstForm <- paste0("InstDist ~ ", Par$Form[i])
  
  AdonisRes <-  adonis2(formula = as.formula(InstForm), 
                        data = PsMeta, 
                        by = "terms", 
                        permutations = Nperm, 
                        parallel = 4) %>% 
                  tidy() %>% 
                  mutate(TaxaLvl = Par$TaxLvl[i], 
                         Normalization = Par$Norm[i], 
                         Distance = Par$Dist[i], 
                         Formula = InstForm, 
                         Permuations = Nperm) %>% 
                  add_case() %>% 
                  bind_rows(AdonisRes, .) %>% 
                  suppressWarnings()
  
}

write.csv(AdonisRes, 
          file = paste0(DirOut, "/beta/tabs/AdonisResults.csv"), 
          row.names = FALSE, na = "")


#-------------------------------------------------------------------------------
# Plot results 
#-------------------------------------------------------------------------------
Par <- expand.grid("TaxLvl" = TaxaLvl, 
                   "Norm" = Norm, 
                   stringsAsFactors = FALSE)

OrdPlots <- list()

for(i in 1:nrow(Par)) {

  #-------------------------------------------------------------------------------
  # PCoA
  #-------------------------------------------------------------------------------
  PCoA_plot <- plot_PCoAs(list_of_dists = DistsLs[[Par$TaxLvl[i]]][[Par$Norm[i]]], 
                           metadata = PsMeta, 
                           group_col = "Treatment", 
                           connect_obervations_col = "ParticipantID", 
                           order_connect_col = "TimeFactor",
                           custom_colors = ColorSchLs$Treatment, 
                           shape_col = "TimeFactor") 
  
  #-------------------------------------------------------------------------------
  # dbRDA
  #-------------------------------------------------------------------------------
  # Custom function: Make RDA data  
  InstRdaData <- rda_extract_for_plot(dists_ls = DistsLs[[Par$TaxLvl[i]]][[Par$Norm[i]]], 
                                      metadata = PsMeta, 
                                      form = AdonisFormula)
  
  # Plot 
  dbRDA_plot <- plot_extracted_rda_data(extracted_data = InstRdaData, 
                                         connect_observations = TRUE, 
                                         add_elepses = TRUE, 
                                         group_col = "Treatment", 
                                         observ_time_group_col = "ParticipantID", 
                                         time_col = "TimeFactor", 
                                         color_vec = ColorSchLs$Treatment)
  
  for(j in c("PCoA_plot", "dbRDA_plot")) {
    
    p <- get(j)
    
    ggsave(filename = paste0(DirOut, "/beta/plots/", 
                             j, "-",
                             Par$TaxLvl[i], "_", 
                             Par$Norm[i], ".svg"), 
           plot = p$Comb, 
           width = 8, height = 5.5)
    
    OrdPlots[[Par$TaxLvl[i]]][[Par$Norm[i]]][[j]] <- p
    
  }

}

