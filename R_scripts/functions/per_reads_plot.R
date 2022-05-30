per_reads_plot = function(filtered_peak_features,
                          all_features,
                          reads = c("peak file", "input file"),
                          mode = "IntersectionStrict",
                          experiment_name = "",
                          ggplot_title = paste0("Percentage of ", experiment_name, " reads in called peaks")) { 
  
  # Count reads of TF peaks
  peaks_counts = summarizeOverlaps(
    features = filtered_peak_features, 
    reads = reads,
    mode = mode)
  peaks_counts = assays(peaks_counts)[[1]]
  
  TF_peak_counts = peaks_counts[,1] %>% sum()
  Input_peak_counts = peaks_counts[,2] %>% sum()
  
  
  
  # Count all reads of ChIP-seq
  total_counts = summarizeOverlaps(
    features = all_features, 
    reads = reads,
    mode = mode)
  total_counts = assays(total_counts)[[1]]
  
  TF_total_peak_counts = total_counts[,1] %>% sum()
  Input_total_peak_counts = total_counts[,2] %>% sum()
  
  
  
  # Count matrix of ChIP3 peaks
  counts_df = as.data.frame(cbind(rbind(TF_peak_counts, TF_total_peak_counts-TF_peak_counts, 
                                        Input_peak_counts, Input_total_peak_counts-Input_peak_counts),
                                  rbind(TF_total_peak_counts, TF_total_peak_counts,
                                        Input_total_peak_counts, Input_total_peak_counts)))
  colnames(counts_df) = c("Reads", "Total")
  counts_df = counts_df %>% as.tibble() %>% mutate(counts_df, experiment = c("ChiP", "ChiP", "Input", "Input"),
                                                   enriched = c("Peak", "Not Peak", "Peak", "Not Peak"),
                                                   percentage = Reads/Total)
  
  
  # Plot the percentage of reads in peaks
  per_reads_plot = ggplot(counts_df, aes(x = experiment,  y = percentage,  fill = enriched)) +
    geom_bar(stat='identity', position='dodge') + 
    geom_text(label=paste0(round(counts_df$percentage*100, 2), "%"),
              position = position_dodge(width = 1),
              vjust = -0.5, size = 3.5) + 
    theme(axis.text = element_text(size=10, face='bold'),
          axis.title = element_text(size=12,face="bold"),
          plot.title = element_text(hjust = 0.5), 
          legend.title=element_blank(), 
          legend.key=element_rect(colour = NA, fill = NA)) + 
    xlab('Experiment') + ylab('% of reads in region') +
    ggtitle(ggplot_title) +
    scale_fill_manual(values=c('gray','red')) + theme_update(axis.text.x = element_text(color = "black"), 
                                                             axis.text.y = element_text(color = "black"),
                                                             axis.ticks = element_line(color = "black"),
                                                             panel.border = element_rect(color="black", fill = NA),
                                                             panel.grid.major = element_blank(), 
                                                             panel.grid.minor = element_blank(),
                                                             panel.background = element_blank(),
                                                             legend.background = element_blank(),
                                                             plot.background = element_rect(fill = "transparent", color = NA))
  
  
  
  return(list(counts_df, per_reads_plot))
  
}