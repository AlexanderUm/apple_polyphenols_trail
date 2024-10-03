################################################################################
# DA
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
TaxaLvl <- prm.ls[["DA"]][["Tax_lvl"]]

Norm <- prm.ls[["DA"]][["Norm"]]

PlotNorm <- prm.ls[["DA"]][["PlotNorm"]]

# DA 

# Directories
dir.create(paste0(DirOut, "/DA/plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(DirOut, "/DA/tabs"), recursive = TRUE, showWarnings = FALSE)

# Data 
load(paste0(DirOut, "/supp/0_data.Rdata"))

# Objects 
ResDaLs <- list()


#-------------------------------------------------------------------------------
# DA with LinDA
#-------------------------------------------------------------------------------
Par <- expand.grid("TaxLvl" = TaxaLvl, 
                   "Norm" = Norm, 
                   stringsAsFactors = FALSE)

ResLinDf <- NULL

for(i in 1:nrow(Par)) {
  
  PsInst <- PsLs[[Par$TaxLvl[i]]][[Par$Norm[i]]]  
  
  # Convert to mstat object
  mstat <- mStat_convert_phyloseq_to_data_obj(PsInst)
  
  # Change to class 
  mstat$meta.dat <- PsMeta
  
  LinRes <- generate_taxa_trend_test_long(
                        data.obj = mstat,
                        subject.var = ColParticipantID,
                        time.var = ColTimeNumeric,
                        group.var = ColGroup,
                        prev.filter = prm.ls$DA$Min_prev,
                        feature.level = "original",
                        feature.dat.type = "count")
  
  ResLinDf <- LinRes$original[[1]] %>% 
                      bind_cols(Par[i, ]) %>% 
                      mutate(Method ="LD_tr", 
                             Prevalence_CutOff = prm.ls$DA$Min_prev) %>% 
                      bind_rows(ResLinDf, .)
                    
}

ResLinDf <- ResLinDf %>% 
                arrange(P.Value)

write.csv(ResLinDf, 
          file = paste0(DirOut, "/DA/tabs/LinDAComb.csv"), 
          row.names = FALSE, na = "")


#-------------------------------------------------------------------------------
# Spaghetti plots
#-------------------------------------------------------------------------------
Par <- expand.grid("TaxLvl" = TaxaLvl, 
                   "Norm" = PlotNorm, 
                   stringsAsFactors = FALSE)

for(i in 1:nrow(Par)) {
  
  PathIst <- paste0(DirOut, "/DA/plots/", Par$TaxLvl[i], "/", Par$Norm[i])
  
  dir.create(PathIst, showWarnings = FALSE, recursive = TRUE)
  
  TaxToPlot <- ResLinDf %>% 
                  filter(TaxLvl == Par$TaxLvl[i]) %>% 
                  filter(Adjusted.P.Value <= prm.ls$DA$Qval_cutoff)
  
  if(nrow(TaxToPlot) < 1) {
    
    TaxToPlot <- ResLinDf %>% 
                  filter(TaxLvl == Par$TaxLvl[i]) %>% 
                  head(n=20)
    
  }
  
  TaxToPlotVec <- TaxToPlot$Variable
  
  
  OtuInst <- PsLs[[Par$TaxLvl[i]]][[Par$Norm[i]]] %>% 
              otu_table() %>% 
              t() %>% 
              as.data.frame() %>% 
              bind_cols(PsMeta)
  
  for(j in TaxToPlotVec) {
    
    SpagPlot <- ggplot(OtuInst, aes(y = .data[[j]], 
                        x = .data[[ColTimeNumeric]], 
                        color = .data[[ColGroup]])) + 
                  geom_point(size = 2) +
                  geom_line(aes(group = .data[[ColParticipantID]])) +
                  scale_x_continuous(breaks = sort(unique(OtuInst[[ColTimeNumeric]]))) + 
                  theme_bw() + 
                  scale_color_manual(values = ColorSchLs$Treatment) + 
                  ggtitle(gsub("__", " ", j)) + 
                  xlab("Time (Weeks)") + 
                  ylab(paste0("Abundance (", Par$Norm[i], ")"))
    
    ggsave(filename = paste0(PathIst, "/", j, ".svg"), 
           plot = SpagPlot, width = 4.5, height = 3)
    
  }
  
  #-----------------------------------------------------------------------------
  # Write out OTU and Taxa tables
  write.csv(OtuInst, 
            file = paste0(DirOut, "/DA/tabs/", 
                          Par$TaxLvl[i], 
                          "_abundance_", Par$Norm[i], ".csv"))
  
  PsLs[[Par$TaxLvl[i]]][[Par$Norm[i]]] %>% 
                tax_table() %>% 
                as.matrix() %>% 
                as.data.frame() %>% 
    write.csv(., file = paste0(DirOut, "/DA/tabs/", 
                                Par$TaxLvl[i], 
                                "_taxonomy.csv"))
}
