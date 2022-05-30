enrichR_call = function(data_path = "./",
                       treatment = ".sorted.bam",
                       control = ".sorted.bam",
                       experiment_names = c(str_extract(treatment, "ChIP\\d+|HB\\d+"), str_extract(control, "Input\\d+|control")),
                       genome = tilling_window,
                       filter_qval = 0.05,
                       countConfiguration = countConfigSingleEnd()) { 
  
  
  # Enrichment calling (https://rdrr.io/bioc/normr/man/normR-enrichR.html)
  
  chip_fit = enrichR(
    
    # ChIP file
    treatment   = file.path(data_path, treatment),
    
    # control file
    control     = file.path(data_path, control),
    
    # genome version
    genome      = genome,
    verbose     = FALSE,
    
    # window size for counting
    countConfig = countConfiguration)
  
  summary(chip_fit)
  
  
  
  # Assign statistical significance and calculate enrichment to ChIP enriched regions
  
  # extract peaks location
  chip_fit_peaks = getRanges(chip_fit)
  
  # annotate peaks with the q-values and enrichment
  chip_fit_peaks$qvalue = getQvalues(chip_fit)
  chip_fit_peaks$enrichment = getEnrichment(chip_fit)
  
  # select enriched peaks
  chip_fit_peaks = subset(chip_fit_peaks, !is.na(component))
  chip_fit_peaks = subset(chip_fit_peaks, qvalue < filter_qval)
  chip_fit_peaks = chip_fit_peaks[order(chip_fit_peaks$qvalue)] # order the peaks based on the q-val
  
  # collapse nearby enriched regions
  chip_fit_peaks = reduce(chip_fit_peaks)
  
  
  # Percentage of reads in ChIP peaks
  
  # extract counts from the fit object per window
  chip_fit_counts = data.frame(getCounts(chip_fit)) 
  
  # change the column names of the data.frame
  if (colnames(chip_fit_counts) == c("treatment", "control")) {
    colnames(chip_fit_counts) = experiment_names
  }
  else {
    colnames(chip_fit_counts) = rev(experiment_names)
  }
  
  
  # extract the q-value corresponding to each bin
  chip_fit_counts$qvalue = getQvalues(chip_fit)
  
  # Define which regions are peaks using a q value cutoff
  chip_fit_counts$enriched[is.na(chip_fit_counts$qvalue)]  = 'Not Peak'
  chip_fit_counts$enriched[chip_fit_counts$qvalue > 0.05]  = 'Not Peak'
  chip_fit_counts$enriched[chip_fit_counts$qvalue <= 0.05] = 'Peak'
  
  # remove the q value column
  chip_fit_counts$qvalue = NULL 
  
  # reshape the data.frame into a long format
  chip_fit_counts_df = tidyr::pivot_longer(
    data = chip_fit_counts, 
    cols = -enriched,
    names_to = 'experiment',
    values_to = 'counts'
  )
  
  # sum the number of reads in the Peak and Not Peak regions
  chip_fit_counts_df = group_by(.data = chip_fit_counts_df, experiment, enriched)
  chip_fit_counts_df = summarize(.data = chip_fit_counts_df, num_of_reads = sum(counts))
  
  # calculate the percentage of reads
  chip_fit_counts_df = group_by(.data = chip_fit_counts_df, experiment)
  chip_fit_counts_df = mutate(.data = chip_fit_counts_df, total=sum(num_of_reads))
  chip_fit_counts_df$percentage = with(chip_fit_counts_df, round(num_of_reads/total,2))
  
  chip_fit_counts_df
  
  
  
  # Plot percentage of reads in called peaks
  
  # Order bars of DAP fit counts so they're in the same order than ChIP's (for better plot comparisons)
  if (chip_fit_counts_df$experiment %>% grepl("control", .) %>% sum()>0) {
    chip_per_reads_plot = ggplot(chip_fit_counts_df, aes(x = fct_rev(fct_reorder(experiment, percentage)),  y = percentage,  fill = enriched)) +
      geom_bar(stat='identity', position='dodge') + 
      geom_text(label=paste0(round(chip_fit_counts_df$percentage*100, 2), "%"),
                position = position_dodge(width = 1),
                vjust = -0.5, size = 3.5) + 
      theme(axis.text = element_text(size=10, face='bold'),
            axis.title = element_text(size=12,face="bold"),
            plot.title = element_text(hjust = 0.5), 
            legend.title=element_blank(), 
            legend.key=element_rect(colour = NA, fill = NA)) + 
      xlab('Experiment') + ylab('% of reads in region') +
      ggtitle(paste("Percentage of", experiment_names, "reads in called peaks")) +
      scale_fill_manual(values=c('gray','#4472C4'))
    
  }
  else {
    chip_per_reads_plot = ggplot(chip_fit_counts_df, aes(x = experiment,  y = percentage,  fill = enriched)) +
      geom_bar(stat='identity', position='dodge') + 
      geom_text(label=paste0(round(chip_fit_counts_df$percentage*100, 2), "%"),
                position = position_dodge(width = 1),
                vjust = -0.5, size = 3.5) + 
      theme(axis.text = element_text(size=10, face='bold'),
            axis.title = element_text(size=12,face="bold"),
            plot.title = element_text(hjust = 0.5), 
            legend.title=element_blank(), 
            legend.key=element_rect(colour = NA, fill = NA)) + 
      xlab('Experiment') + ylab('% of reads in region') +
      ggtitle(paste("Percentage of", experiment_names, "reads in called peaks")) +
      scale_fill_manual(values=c('gray','#4472C4'))
      }
  
  chip_per_reads_plot
  
  return(list("NormRFit_object" = chip_fit, 
                "fit_peaks" = chip_fit_peaks, 
                "fit_counts" = chip_fit_counts, 
                "fit_counts_df" = chip_fit_counts_df, 
                "reads_plot" = chip_per_reads_plot))
  
  }
