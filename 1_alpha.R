################################################################################
# Alpha diversity
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
TaxaLvl <- prm.ls[["alpha"]][["Tax_lvl"]]

Norm <- prm.ls[["alpha"]][["Norm"]]


# Alpha 
AlphaIdexes <- prm.ls[["alpha"]][["measures"]]

# Plot 
NcolPlot <- prm.ls[["alpha"]][["NcolPlot"]]

# Directories
dir.create(paste0(DirOut, "/alpha/plots"), recursive = TRUE, showWarnings = FALSE)

dir.create(paste0(DirOut, "/alpha/tabs"), recursive = TRUE, showWarnings = FALSE)

# Custom functions
source("R/phy_alpha.R")

# Data 
load(paste0(DirOut, "/supp/0_data.Rdata"))

# Objects 
ResAlphaLs <- list()

PsInst <- PsLs[[TaxaLvl]][[Norm]]


################################################################################
# Calculate alpha diversity
#-------------------------------------------------------------------------------
AlphaDf <- phy_alpha(PsInst, 
                      measures = AlphaIdexes) %>% 
              bind_cols(., PsMeta) 

AlphaDfLong <- AlphaDf %>% 
                  pivot_longer(cols = all_of(AlphaIdexes))

write.csv(AlphaDf, 
          file = paste0(DirOut, "/alpha/tabs/alpha_diversity_per_sample.csv"), 
          row.names = FALSE)


#-------------------------------------------------------------------------------
# Calculate LMM using LinDA
#-------------------------------------------------------------------------------
# Format alpha diversity data
MstatAlpha <- list()

for(i in 1:length(AlphaIdexes)) {
  
  MstatAlpha[[i]] <- AlphaDf %>% 
                        select(AlphaIdexes[i])
  
}

# Mstat object
MstatObj <- mStat_convert_phyloseq_to_data_obj(PsInst)

# Add original metadata 
MstatObj$meta.dat <- PsMeta


AlphaMstat <- generate_alpha_trend_test_long(
                        data.obj = MstatObj, 
                        alpha.obj = MstatAlpha,
                        alpha.name = AlphaIdexes,
                        time.var = ColTimeNumeric, 
                        subject.var = ColParticipantID,
                        group.var = ColGroup)

# Write out results 
for(i in names(AlphaMstat)) {
  
  write.csv(AlphaMstat[[i]], 
            file = paste0(DirOut, 
                          "/alpha/tabs/LinDa_trend--", 
                          i, ".csv"), row.names = FALSE)
  
}


#-------------------------------------------------------------------------------
# Plot data as a spaghetti plot 
#-------------------------------------------------------------------------------
# Shorten digits function 
short_digt <- function(number, Ndigits = 3) {
  
  DigOut <- ifelse(round(number, digits = Ndigits) == 0, 
                     paste0("<0.", paste(rep("0", Ndigits-1), collapse = ""), "1"), 
                     paste0("=", sprintf(paste0("%.", Ndigits, "f"), 
                                         round(number, digits = Ndigits))))
  
  return(DigOut)
  
}

# Extract statistical information 
StatTextDf <- NULL

for(i in names(AlphaMstat)) {
  
    StatTextDf <- AlphaMstat[[i]] %>% 
                  mutate(Text = paste0("P", short_digt(P.Value, 3), 
                                       " [Est", short_digt(Estimate, 2), "]")) %>% 
                  .[nrow(AlphaMstat[[i]]), "Text"] %>% 
                  mutate(name = i, 
                         y = max(AlphaDfLong[AlphaDfLong$name == i, "value"]) + 
                           (max(AlphaDfLong[AlphaDfLong$name == i, "value"]) - 
                              min(AlphaDfLong[AlphaDfLong$name == i, "value"]))*0.2, 
                         x = min(AlphaDfLong[AlphaDfLong$name == i, ColTimeNumeric])) %>% 
                  bind_rows(StatTextDf, .)
  
}

AlphaPlot <- ggplot(AlphaDfLong, 
                  aes(x = .data[[ColTimeNumeric]], 
                      group = .data[[ColParticipantID]], 
                      color = .data[[ColGroup]])) +
              geom_line(aes(y = value), 
                        alpha = 0.75, linewidth = 1) + 
              geom_point(aes(y = value, 
                             shape = .data[[ColGroup]]), 
                         size = 2) +
              geom_label(data = StatTextDf, 
                         aes(label = Text, y = y, x = x), 
                         inherit.aes = FALSE, hjust = 0) +
              facet_grid(c("name"), scales = "free") +
              theme_bw() + 
              scale_color_manual(values = ColorSchLs[[ColGroup]]) + 
              theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.title.y = element_blank()) + 
              xlab("Week") + 
              scale_x_continuous(breaks = sort(unique(AlphaDfLong[[ColTimeNumeric]])))

AlphaPlotLs <- list("plot" = AlphaPlot, 
                    "w" = NcolPlot*3.5 + 1.5, 
                    "h" = ceiling(length(AlphaIdexes)/NcolPlot)*2.5)
  
ggsave(filename = paste0(DirOut, "/alpha/plots/alpha_spag.png"), 
       plot = AlphaPlotLs$plot, 
       width = AlphaPlotLs$w, 
       height = AlphaPlotLs$h, 
       dpi = 600)
  

ResAlphaLs[["Tables"]][["LinDA"]] <- AlphaMstat

ResAlphaLs[["Plots"]][["Spagg"]] <- AlphaPlotLs


#-------------------------------------------------------------------------------
save(list = c("ResAlphaLs"), 
     file = paste0(DirOut, "/supp/res_alpha.Rdata"))

# Clean environment 
rm(list = ls())
gc()
