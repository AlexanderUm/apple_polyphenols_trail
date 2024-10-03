
alpha_kw_or_wil <- function(alpha_df, 
                            meta_data, 
                            col_group, 
                            plot_ncol, 
                            sig_text_size = 3) {
  
  require("tidyverse")
  require("FSA")
  
  # Checks  
  if(!identical(rownames(alpha_df), rownames(meta_data))) {
    
    stop("Row names in alpha diversity table and metadata are not identical.")
    
  }
  
  if(length(levels(meta_data[[col_group]])) < 2) { 
    
    stop("Group column of hte meta data has less then two level. No statistical test can be perforemed.")
    
    }

  
  # Combine data
  DataComb <- bind_cols(alpha_df, meta_data)
  
  DataCombLong <-  DataComb %>% 
                      pivot_longer(cols = all_of(colnames(alpha_df)), 
                                   names_to = "Index")
  
  #-----------------------------------------------------------------------------
  # Statistical testing 
  #-----------------------------------------------------------------------------
  StatResComb <- NULL
  
  for(iInd in colnames(alpha_df)) { 
    
    InstForm <- paste0(iInd, "~", col_group)
    
    # Wilcox test
    if(length(levels(DataComb[[col_group]])) == 2) {
      
      InstRes <- wilcox.test(as.formula(InstForm), DataComb) %>% 
                        tidy() 
      
    }
    
    if(length(levels(DataComb[[col_group]])) > 2) {
      
      KwRes <- kruskal.test(as.formula(InstForm), DataComb) 
      
      InstRes <- dunnTest(as.formula(InstForm), DataComb)$res  %>% 
                        mutate(Test = "Dunn", 
                               p.value = P.adj) %>% 
                        add_row(Comparison = "Overall", 
                                Test = "Kruskal-Wallis", 
                                p.value = KwRes$p.value,
                                .before = 1)
    }
    
    StatResComb <- InstRes %>% 
                    mutate(Index = iInd) %>% 
                    rbind(StatResComb, .)
  }
  
  #-----------------------------------------------------------------------------
  # Plot results
  #-----------------------------------------------------------------------------
  # Base plot
  BasePlot <- ggplot(DataCombLong, 
                   aes(y = value, x = .data[[col_group]])) + 
                    geom_jitter(aes(colour = .data[[col_group]]), 
                                height = 0, 
                                width = 0.15, 
                                alpha = 0.5, 
                                na.rm = TRUE) +
                    geom_violin(fill = NA) + 
                    facet_grid(c("Index"), scales = "free") + 
                    theme_bw() + 
                    theme(axis.line = element_line(color='black'),
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), 
                          axis.title.y = element_blank())
  
  # Plot annotation dataframe
  MinMaxDf <- DataCombLong %>% 
                  reframe(y_min = min(value), 
                          y_max = max(value), 
                          .by = Index)
  
  StatResCombF <- StatResComb %>% 
                      mutate(p_short = round(p.value, 3)) %>% 
                      mutate(p_text = ifelse(p_short == 0, 
                                             "P<0.001", 
                                             paste0("P=", sprintf("%.3f", 
                                                                  p_short)))) %>% 
                      left_join(., MinMaxDf, by = "Index")
  
  # Plot in correspondence with number of groups 
  # If 2 groups add wilcoxon results 
  if(length(levels(DataComb[[col_group]])) == 2) {
    
    SigDf <- StatResCombF %>% 
                mutate(Start = levels(DataComb[[col_group]])[1], 
                       End = levels(DataComb[[col_group]])[2], 
                       y_sig = y_max + (y_max - y_min) * 0.2, 
                       y_inv_point = y_max + (y_max - y_min) * 0.4)
    
    FinalPlot <- BasePlot +  
                geom_signif(data = SigDf,
                            aes(xmin = Start,
                                xmax = End,
                                annotations = p_text,
                                y_position = y_sig),
                            textsize = sig_text_size, 
                            vjust = -0.2,
                            manual = TRUE, 
                            margin_top = 1) + 
                geom_point(data = SigDf,
                           aes(x = End, 
                               y = y_inv_point), 
                           x=NA, na.rm = TRUE) 
  }
  
  # More than two levels (Kruskal-Wallis)
  if(length(levels(DataComb[[col_group]])) > 2) {
    
    KwAnnot <- StatResCombF %>% 
                filter(Test == "Kruskal-Wallis") %>% 
                mutate(Text = paste0(Test, ": ", p_text), 
                       x_text = levels(DataComb[[col_group]])[1], 
                       y_text = y_max + (y_max - y_min) * 0.25, 
                       y_inv_point = y_max + (y_max - y_min) * 0.4)
    
    FinalPlot <-  BasePlot + 
                    geom_text(data = KwAnnot, 
                              aes(x = x_text, 
                                  y = y_text, 
                                  label = Text), 
                              hjust = 0.1, 
                              size = sig_text_size) + 
                    geom_point(data = KwAnnot,
                               aes(x = x_text, 
                                   y = y_inv_point), 
                               x=NA) 
                  
    if(any(KwAnnot[["p.value"]] <= 0.05)) {
      
      # Data frame for significance levels 
      SigDunnDf <- StatResCombF %>% 
                      filter(Test == "Dunn", 
                             p.value <= 0.05) %>% 
                      mutate(Start = str_split(.$Comparison, 
                                               " - ", simplify = TRUE)[, 1], 
                             End = str_split(.$Comparison, 
                                             " - ", simplify = TRUE)[, 2]) %>% 
                      group_by(Index) %>% 
                      mutate(y_sig = y_max + ((y_max - y_min)*(1:n()*0.35))) %>% 
                      mutate(y_inv_point = y_max + ((y_max - y_min)*((n()+2)*0.35)), 
                             y_text = y_max + ((y_max - y_min)*((n()+1)*0.45))) %>% 
                      ungroup()
      
      # Adjust y text location taking into account significance bars 
      KwAnnot <- left_join(KwAnnot, SigDunnDf[c("Index", "y_text")], 
                           by = "Index", 
                           multiple = "first") %>% 
                  mutate(y_text = ifelse(is.na(y_text.y), 
                                         y_text.x, 
                                         y_text.y))
      
      # Final plot 
      FinalPlot <- BasePlot +  
                    geom_text(data = KwAnnot, 
                              aes(x = x_text, 
                                  y = y_text, 
                                  label = Text), 
                              hjust = 0.1, 
                              size = sig_text_size) +
                    geom_signif(data = SigDunnDf,
                                aes(xmin = Start,
                                    xmax = End,
                                    annotations = p_text,
                                    y_position = y_sig),
                                textsize = sig_text_size, 
                                vjust = -0.2,
                                manual = TRUE, 
                                margin_top = 1) + 
                    geom_point(data = SigDunnDf,
                               aes(x = End, 
                                   y = y_inv_point), 
                               x=NA) 
    } 
  }
  
  ResOut <- list("Stat_tab" = StatResComb, 
                 "Plot" = FinalPlot)
  
  return(ResOut)
  
}