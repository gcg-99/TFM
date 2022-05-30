count_cpm = function(data_path = "./bash_output/ChIPseq_mapping/",
                     treatment = ".sorted.bam",
                     control = ".sorted.bam",
                     experiment_names = c(str_extract(treatment, "ChIP\\d+|HB\\d+"), str_extract(control, "Input\\d+|control")),
                     tilling_window = tilling_window,
                     fit_peaks = ChIP_enrich_results$fit_peaks
                     
                     )  { 
  
  
  # Count the reads and normalize them
  chip_counts = summarizeOverlaps(features = tilling_window, 
    reads = c(file.path(data_path, treatment), file.path(data_path, control))
  )
  chip_counts = assays(chip_counts)[[1]]
  chip_cpm = t(t(chip_counts)*(1000000/colSums(chip_counts))) %>% as.data.frame()
  
  
  
  # Plot ChIP vs Input signal
  
  # Find enriched tilling windows
  chip_enriched_regions = countOverlaps(tilling_window, fit_peaks) > 0
  
  # ChIP2 vs Input2 signal plot
  chip_cpm$enriched_regions = as.numeric(chip_enriched_regions) %>%
    str_replace(., '1', "Enriched region") %>% str_replace(., '0', "Non-enriched region")
  colnames(chip_cpm) = c("ChIP", "Input", "enriched_regions")
  
  chip_enrichment_plot = ggplot(data = chip_cpm,
                                aes(x = Input, y = ChIP, color = enriched_regions), subset = .(label == "Enriched region")) +
    geom_point(alpha=0.5) + geom_abline(slope = 1) + xlab(paste(experiment_names[2], "CPM")) + ylab(paste(experiment_names[1], "CPM")) + 
    ggtitle(paste(experiment_names[1], "vs", experiment_names[2])) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title=element_blank(),
          legend.key=element_rect(colour = NA, fill = NA)) + scale_color_manual(values=c('#4472C4','gray'))
  
  return(list("ChIP_CPM_counts" = chip_cpm, "CPM_plot" = chip_enrichment_plot))
  
  }

