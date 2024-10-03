################################################################################
# List of parameter that will be used 
################################################################################
prm.ls <- list()

dir.create("out/supp", recursive = TRUE, showWarnings = FALSE)

#-------------------------------------------------------------------------------
# General 
#-------------------------------------------------------------------------------
prm.ls[["general"]] <- list("seed" = 39653, 
                            "Group_col" = "Treatment", 
                            "SeqID_col" = "SeqID",
                            "SampleID_col" = "SampleID", 
                            "TimeNumeric_col" = "TimeNumeric", 
                            "PraticipantID_col" = "ParticipantID",
                            "n_core" = 4, 
                            "DirOut" = "out") 


#-------------------------------------------------------------------------------
# Data preparation (0_data_prep.R)
#-------------------------------------------------------------------------------
prm.ls[["data_prep"]] <- list("seq_path" = "data/qiime/",
                              "meta_path" = "data/meta/meta_data.csv", 
                              "Meta_Cols_to_factors" = c("Treatment", 
                                                         "ParticipantID"),
                              "samples_kick" = list("SeqID" = c("102S40")), # NULL if keep all
                              "min_reads_tax" = 10, 
                              "rare_depth" = 10000, 
                              "tax_king_keep" = c("d__Bacteria", "d__Archaea"), 
                              "tax_kick" = c("Genus" = "Mitochondria", 
                                             "Genus" = "Chloroplast"), 
                              "tax_levels" = c("ASV", "Genus", "Family", "Phylum"))


#-------------------------------------------------------------------------------
# Composition 
#-------------------------------------------------------------------------------
prm.ls[["composition"]] <- list("Tax_lvl" = c("ASV", "Genus", "Phylum"),
                                "Norm" = "Rel",
                                "ColGroup" = "Treatment",
                                "out_dir_path" = "out/alpha", 
                                "MinPrev" = 0.2, 
                                "PlotTopN" = c("ASV" = 25, 
                                               "Genus" = 25, 
                                               "Phylum" = 25), 
                                "ShowSampNames" = FALSE,
                                "DirOut" = "out/compostion")


#-------------------------------------------------------------------------------
# Alpha 
#-------------------------------------------------------------------------------
prm.ls[["alpha"]] <- list("Tax_lvl" = "ASV",
                          "Norm" = "Rare",
                          "measures" = c("Observed", "Shannon",
                                         "InvSimpson", "PhyloDiverity"), 
                          "NcolPlot" = 1, 
                          "out_dir_path" = "out/alpha")


#-------------------------------------------------------------------------------
# Beta
#-------------------------------------------------------------------------------
prm.ls[["beta"]] <- list("Tax_lvl" = c("ASV", "Genus"),
                         "Norm" = "CSS",
                         "distances" = c("Unweighted UniFrac" = "unifrac", 
                                         "Weighted UniFrac" = "wunifrac", 
                                         "Jaccard" = "jaccard", 
                                         "Bray-Curtis" = "bray"),
                         "n_perm" = 999, 
                         "AdonisFormula" = "TimeNumeric*Treatment")


#-------------------------------------------------------------------------------
# DA
#-------------------------------------------------------------------------------
prm.ls[["DA"]] <- list("Tax_lvl" = c("ASV", "Genus"),
                         "Norm" = c("Raw"),
                         "PlotNorm" = c("TSS"),
                         "Pval_cutoff" = 0, # 0 means will not be used for filtering
                         "Qval_cutoff" = 0.05, 
                         "Min_prev" = 0.25, 
                         "Main_cols" = "Diet",
                         "Maas_fix_effect" = c("Batch", "Diet"), 
                         "Maas_method" = "LM", 
                         "Maas_pcorr_method" = "BH", 
                         "Maas_Reff" = "Diet,glucose", 
                         "Maas_norm" = "NONE",
                         "Maas_tansform" = "NONE",
                         "heat_width" = 12,
                         "out_dir_path" = "out/DA")


#-------------------------------------------------------------------------------
lib.to.load <- c("phyloseq", "tidyverse", "qiime2R", 
                 "metagenomeSeq", "RColorBrewer", "broom", 
                 "ggsignif", "vegan", "usedist", "cowplot", "ape", 
                 "Maaslin2", "ComplexHeatmap", "circlize", "MicrobiomeStat")

for (i in lib.to.load) {library(i, character.only = TRUE)}


#-------------------------------------------------------------------------------
# Write out parameters file 
#-------------------------------------------------------------------------------
save(list = c("prm.ls"), file = "prm.Rdata")

rm(list = ls())
