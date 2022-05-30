# find_tests function -----------------------------------------------------
# wrapper function to retrieve tests in data easily
find_tests = function(edger_merged_results) {
  edger_merged_results %>% colnames %>% str_extract("^logFC_.*") %>% discard(is.na) %>% gsub("logFC_","",.)
}


# volcano_plot function ---------------------------------------------------
# edger_merged_results is a data.frame from edgeR, with contrast names appended to logFC and FDR colnames (e.g., logFC_hb21i_vs_gus)
# comparison is the contrast used to calculate fold change (use find_tests to find possible comparisons) (e.g., "hb21i_vs_gus")

volcano_plot = function(edger_merged_results,
                        comparison = NULL,
                        gene_id_column = "Geneid",
                        logfdr_threshold = -log10(0.05),
                        logfc_threshold = log2(2),
                        top_genes = 5,
                        col_up = "red",
                        col_dn = "blue",
                        col_ns = "gray") {
  if(is.null(comparison)){stop("Please select a comparison")}

  # obtain column name of logFC and FDR data depending on comparison name
  logfc_column = paste0("logFC_",comparison)
  fdr_column = paste0("FDR_",comparison)
  
  # manipulate edger object to specify which points are colored
  vp_input = edger_merged_results %>% as.data.frame %>% 
    dplyr::select(all_of(c(gene_id_column,logfc_column,fdr_column))) %>% 
    mutate(log10_fdr = -log10(get(fdr_column))) %>% 
    mutate(vp_color = case_when(
      (log10_fdr > logfdr_threshold) & (get(logfc_column) > logfc_threshold) ~ "sig_up",
      (log10_fdr > logfdr_threshold) & (get(logfc_column) < -logfc_threshold) ~ "sig_dn",
      TRUE ~ "non_s"))
  
  # generate labels for "top_genes" up and down regulated
  top_up = vp_input %>% filter(vp_color == "sig_up") %>% top_n(top_genes,get(logfc_column)) %>% pull(Geneid)
  top_dn = vp_input %>% filter(vp_color == "sig_dn") %>% top_n(top_genes,-get(logfc_column)) %>% pull(Geneid)
  gene_symbol = suppressMessages(mapIds(org.At.tair.db,
                                        column = "SYMBOL",
                                        keys = pull(vp_input,get(gene_id_column)),
                                        keytype = "TAIR",
                                        multiVals = "first"))
  
  vp_input = vp_input %>% 
    mutate(gene_symbol = coalesce(gene_symbol,get(gene_id_column))) %>% 
    mutate(gene_label = case_when(
      (get(gene_id_column) %in% c(top_up,top_dn)) ~ gene_symbol,
      TRUE ~ ""))
  
  # label for count of up and down genes
  label_up = vp_input %>% filter(vp_color == "sig_up") %>% dplyr::count() %>% as.numeric() %>% paste("Up:", .)
  label_dn = vp_input %>% filter(vp_color == "sig_dn") %>% dplyr::count() %>% as.numeric() %>% paste("Down:", .)
  
  # generate the plot
  vplot = ggplot(vp_input, aes(x=get(logfc_column),y=log10_fdr,color=vp_color))+
    geom_point(shape=16, alpha = 0.5) +
    scale_x_continuous(limits = symmetric_limits) +
    scale_color_manual(values = c("sig_dn" = col_dn,
                                  "sig_up" = col_up,
                                  "non_s" = col_ns)) +
    geom_text_repel(aes(label=gene_label),
                    max.overlaps = Inf,
                    min.segment.length = 0,
                    #box.padding = 5,
                    #point.padding = 5,
                    show.legend = FALSE) +
    annotate(geom="text",
             label = label_up, 
             color = col_up,
             x = Inf, y = -Inf,
             hjust = 1.1, vjust = -1) +
    annotate(geom="text",
             label = label_dn, 
             color = col_dn,
             x = -Inf, y = -Inf,
             hjust = -0.1, vjust = -1) +
    labs(title = comparison,
         x = "log2(FC)",
         y = "-log10(adj. p-val)") +
    coord_cartesian(clip = "off") +
    theme(legend.position="none",
          plot.title = element_text(hjust = 0.5))
  
  return(vplot)
}