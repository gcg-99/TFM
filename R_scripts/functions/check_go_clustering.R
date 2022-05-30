check_go_clustering = function(go_id) {
  # if(!exists("athGO")){
  #   athGO = godata("org.At.tair.db", ont="BP")
  # }
  # GO_anc = as.list(GOBPANCESTOR)
  # GO_def = suppressMessages(AnnotationDbi::select(GO.db,
  #                                                 keys=(keys(GO.db)),
  #                                                 columns=c("GOID","TERM","DEFINITION"),
  #                                                 keytype="GOID"))
  
  ss_rel = mgoSim(go_id, go_id, semData=athGO, measure="Rel", combine=NULL)
  GO_hclust = (1 - ss_rel) %>% as.dist %>%
    hclust(method = "ward.D2")
  
  GO_hclust_n = length(GO_hclust$order)
  clustrange = as.clustrange(GO_hclust, diss = (1 - ss_rel), ncluster = GO_hclust_n)
  
  crange_pbc = clustrange$stats$PBC
  max_crange_pbc = find_local_max(crange_pbc)
  
  heatmap_list = NULL
  i = 1
  
  for(k in max_crange_pbc){
    test = cutree(GO_hclust, k=k) %>% as.factor
    
    df = data.frame("id" = rownames(ss_rel),
                    "cluster" = test)
    
    repr_GO = GO_anc[go_id] %>% stack() %>% `colnames<-`(c("ancestral", "go_id")) %>% 
      left_join(dplyr::select(df, c("cluster", "id")), by = setNames("id", "go_id")) %>% 
      group_by(cluster, ancestral) %>% add_count() %>% ungroup %>% 
      mutate(IC = GO_semsim@IC[ancestral]) %>% 
      mutate(importance = (n*IC^3)) %>% 
      group_by(cluster) %>% 
      dplyr::slice(which.max(importance)) %>% 
      left_join(dplyr::select(GO_def, c("GOID","TERM")), by = setNames("GOID", "ancestral")) %>% 
      dplyr::select(c("ancestral","cluster","importance","TERM"))
    
    df2 = df %>% tibble %>% left_join(dplyr::select(repr_GO, c("cluster", "TERM")), by ="cluster") %>% 
      rename("anc_term" = "TERM") %>% mutate(cluster = factor(cluster, levels = levels(test))) %>% 
      mutate(anc_term = factor(anc_term, levels = unique(.$anc_term[order(.$cluster)]) ))
    
    anno_df = df2[,c("id","anc_term")] %>% column_to_rownames("id")
    
    set.seed(1)
    col_vector = distinctColorPalette(levels(anno_df$anc_term) %>% length)
    custom_colors = setNames(col_vector, levels(anno_df$anc_term))
    
    # Heatmap parameters (annotation, row labels, color function for matrix)
    ha = HeatmapAnnotation(df = anno_df,
                           which = "row",
                           annotation_label = paste0("GO Clusters (k = ",k,")"),
                           show_annotation_name = FALSE,
                           col = list(anc_term = custom_colors),
                           annotation_legend_param = list(labels_gp = gpar(fontsize = 6),
                                                          grid_height = unit(2, "mm"),
                                                          grid_width = unit(2, "mm")))
    
    row_labs = structure(paste(go_id, GO_data$Description), names = paste0(go_id))
    col_fun = colorRamp2(c(0,1), c("white", "red"))
    
    term_order = unlist(lapply(levels(test), function(le) {
      l = test == le
      if(sum(l) <= 1) {
        return(which(l))
      } else {
        mm = ss_rel[l, l, drop = FALSE]
        which(l)[hclust(stats::dist(mm))$order]
      }
    }))
    
    # Generate Heatmap
    GO_heatmap = Heatmap(ss_rel,
                         name = paste0("k_",k),
                         row_order = term_order, column_order = term_order,
                         col = col_fun,
                         row_labels = row_labs[rownames(ss_rel)],
                         row_names_max_width = max_text_width(row_labs),
                         row_names_gp = gpar(fontsize = 4),
                         right_annotation = ha,
                         show_column_names = FALSE,
                         show_heatmap_legend = FALSE)

    heatmap_list[[i]] = GO_heatmap
    i = i+1
    }
  
  # ht_list = NULL
  # for (i in seq_along(heatmap_list)) {
  #   ht_list = ht_list + heatmap_list[[i]]
  # }
  # 
  return(heatmap_list)
}