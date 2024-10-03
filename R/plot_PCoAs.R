#-------------------------------------------------------------------------------
# Plot PCoAs plots 
#-------------------------------------------------------------------------------
plot_PCoAs <- function(list_of_dists, 
                       metadata, 
                       group_col, 
                       connect_obervations_col = NULL,
                       order_connect_col = NULL,
                       shape_col = NULL,
                       add_ellipsis = TRUE, 
                       custom_colors = NULL, 
                       ncols = 2) {
  
  require(ape)
  require(tidyverse)
  require(cowplot)
  
  pcoa.plots.ls <- NULL
  
  for (i.dist in names(list_of_dists)) { 
    
    dist.inst <- list_of_dists[[i.dist]]
    
    # Make PCoA object
    a.pcoa <- ape::pcoa(dist.inst)
    
    axis.pcoa <- a.pcoa$vectors[, 1:2] %>%
                    as.data.frame() %>%
                    bind_cols(., metadata)
    
    if(!is.null(custom_colors)) { 
      
      axis.pcoa <- axis.pcoa %>% 
                      arrange(across(all_of(c(connect_obervations_col, 
                                            order_connect_col))))
      
      }
    
    # Extract percentage
    var.c <- round((a.pcoa$values$Relative_eig[1:2]*100), 1)
    
    # Base plot 
    pcoa.p <- ggplot(axis.pcoa,
                     aes(x = Axis.1, 
                         y = Axis.2)) +
                    ggtitle(paste0("Distance: ", i.dist)) +
                    coord_fixed() +
                    theme_bw() +
                    xlab(paste0(colnames(axis.pcoa)[1], " [", var.c[1], "%]")) +
                    ylab(paste0(colnames(axis.pcoa)[2], " [", var.c[2], "%]")) +
                    theme(plot.title = element_text(size=11, face="italic"), 
                          panel.grid = element_blank())
      
    if(add_ellipsis) {
      
      pcoa.p <- pcoa.p +
                  stat_ellipse(geom = "polygon", 
                               aes(fill = .data[[group_col]]), 
                               alpha = 0.1, level = 0.90, size =2)
    }
    
    if(!is.null(connect_obervations_col)) { 
      
      pcoa.p <- pcoa.p +
        geom_path(aes(group = .data[[connect_obervations_col]], 
                      color = .data[[group_col]]), 
                  alpha = 0.75)
      
    }
  
    # Add point layer 
    if(!is.null(custom_colors)) {
      
        pcoa.p <- pcoa.p +
                    geom_point(aes(color = .data[[group_col]], 
                                   shape = .data[[shape_col]]),
                               alpha = 0.75, size = 2)
    } else {
      
      pcoa.p <- pcoa.p +
                  geom_point(aes(color = .data[[group_col]]),
                             alpha = 0.75, size = 2)
    }

    
    if(!is.null(custom_colors)) {
      
      pcoa.p <- pcoa.p +
                  scale_color_manual(values = custom_colors) +
                  scale_fill_manual(values = custom_colors)
    }
    
    pcoa.plots.ls[[i.dist]] <- pcoa.p
    
  }
  
  
  #-------------------------------------------------------------------------------
  # Combine plots 
  #-------------------------------------------------------------------------------
  p.ls.inst.f <- lapply(pcoa.plots.ls, 
                        function(x){x+theme(legend.position="none")})
  
  legd.inst <- get_legend(pcoa.plots.ls[[1]])
  
  p.grid.1 <- plot_grid(plotlist = p.ls.inst.f, 
                        ncol = 2)
  
  p.grid.2 <- plot_grid(p.grid.1, 
                        legd.inst, 
                        ncol = ncols, 
                        rel_widths = c(0.8, .2))
  
  out.ls <- list("Comb" = p.grid.2, 
                 "ind" = pcoa.plots.ls)
  
  return(out.ls)
  
}
