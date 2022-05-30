batch_enrichedheatmaps = function(input_list,
                                  ext5 = 3000,
                                  ext3 = 1000,
                                  tss = NULL,
                                  gtf = "./input/Arabidopsis_thaliana.TAIR10.52.gtf.gz",
                                  value_column = NULL,
                                  mean_mode = "w0",
                                  w = 25,
                                  background = 0,
                                  smooth = TRUE,
                                  anno_color = "black",
                                  axis_param = list(side = "left"),
                                  heatmap_col = c("white", anno_color),
                                  heatmap_col_breaks = c(0, 1),
                                  use_raster = TRUE, 
                                  raster_resize_mat = TRUE, 
                                  raster_quality = 1,
                                  show_heatmap_legend = FALSE, 
                                  top_annotation = NULL) {
  
  if (!is(input_list, "CompressedGRangesList")) {
    stop("Input is not in CompressedGRangesList format")
  }
  
  if(is.null(tss)) {
    print(paste("TSS coordinates not provided, generating TSS coordinates from:", gtf))
    tss = plyranges::read_gff(gtf) %>% filter(type == "gene") %>% 
      GenomicRanges::resize(width=1, fix='start') %>% 
      plyranges::select(gene_id) %>% `names<-`(.$gene_id) %>% sort
  }
  
  peaks_in_input_list = input_list %>% sapply(length) %>% unname
  print(paste(c("Generating normalized matrix for", length(input_list), "GRanges with", peaks_in_input_list, "peaks"), collapse =" "))
  norm_matrix_list = lapply(input_list, EnrichedHeatmap::normalizeToMatrix,
                            target = tss,
                            extend = c(ext5,ext3),
                            verbose = FALSE,
                            value_column = value_column,
                            mean_mode = mean_mode, 
                            w = w, 
                            background = background, 
                            smooth = smooth)
  
  print("Generating EnrichedHeatmap plots")
  if (is.null(top_annotation)) {
    top_anno = ComplexHeatmap::HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = anno_color),
                                                                       axis_param = axis_param))
  } else {
    top_anno = top_annotation
  }
  
  # Create color function list only if value_column is provided, otherwise use 0 to 1 scale (pressence/absence) for every plot
  if (!is.null(value_column)) {
    quantile_list = lapply(norm_matrix_list, quantile, probs = heatmap_col_breaks)
    col_fun_list = mapply(
      norm_matrix_list,
      FUN = circlize::colorRamp2,
      breaks = quantile_list,
      MoreArgs = list(colors = heatmap_col)
    )
  } else {
    col_fun_list = rep(list(heatmap_col), length(norm_matrix_list))
  }
  order_list = lapply(norm_matrix_list, function(x) order(enriched_score(x), decreasing = TRUE))
  
  enrich_heatmap_list = mapply(
    EnrichedHeatmap::EnrichedHeatmap,
    norm_matrix_list,
    column_title = names(norm_matrix_list),
    name = names(norm_matrix_list),
    col = col_fun_list,
    row_order = order_list,
    MoreArgs = list(
      use_raster = use_raster,
      raster_resize_mat = raster_resize_mat,
      raster_quality = raster_quality,
      show_heatmap_legend = show_heatmap_legend,
      top_annotation = top_anno
    )
  )
  
  print("Generating HeatmapList")
  ht_list = NULL
  for (i in seq_along(enrich_heatmap_list)) {
    ht_list = ht_list + enrich_heatmap_list[[i]]
  }
  
  return(ht_list)
  
}