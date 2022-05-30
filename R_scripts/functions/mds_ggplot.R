mds_ggplot = function(DGEL_input,
                      group_labels,
                      col_palette,
                      top_genes=500) {
  mds = plotMDS(DGEL_input, top=top_genes,plot=FALSE)
  dim1_label = paste0(mds$axislabel," 1 (",(mds$var.explained[1] * 100 ) %>% round(1)," %)")
  dim2_label = paste0(mds$axislabel," 2 (",(mds$var.explained[2] * 100 ) %>% round(1)," %)")
  
  mds_data = data.frame(dim1 = mds$x, dim2 = mds$y, group = DGEL_input$samples$group) %>% 
    mutate(treatment = recode(group, !!!group_labels))

  mds_plot = ggplot(mds_data, aes(x=dim1, y=dim2, color = treatment)) + 
    geom_point() +
    scale_color_manual(values = col_palette) +
    labs(x = dim1_label,
         y = dim2_label) +
    geom_dl(aes(label = treatment),
            method = "smart.grid") +
    coord_cartesian(clip = "off") +
    theme(aspect.ratio = 1,
          legend.position = "left",
          legend.text.align = 0,
          legend.spacing.x = unit(0.5, "mm"),
          legend.title = element_blank(),
          legend.key = element_blank())
  return(mds_plot)
  }