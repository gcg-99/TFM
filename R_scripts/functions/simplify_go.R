simplify_go = function(go_input,
                       GO_sem_sim = NULL,
                       go_data = "org.At.tair.db",
                       ont = "BP",
                       GO_def = NULL,
                       clust_method = "ward.D2")
  {
  if (is(go_input, "enrichResult")) {
    #print("Input is clusterProfiler::enrichGO enrichResult format")
    go_id = as.data.frame(go_input) %>% dplyr::pull(unique(ID))
  } else {
      go_id = go_input
    }
  
  # Create databases for GOs if not given -----------------------------------
  if(is.null(GO_sem_sim)) {
    print(paste("GO semantic similarity database not provided, generating it from:", go_data, "with", ont, "data"))
    GO_sem_sim = GOSemSim::godata("org.At.tair.db", ont="BP")
  }
  
  GO_anc = as.list(GOBPANCESTOR)
  
  if(is.null(GO_def)) {
  GO_def = suppressMessages(AnnotationDbi::select(GO.db,
                                                  keys=(keys(GO.db)),
                                                  columns=c("GOID","TERM","DEFINITION"),
                                                  keytype="GOID"))
  }

  # Extract GO ids and cluster based on Semantic Similarity using binary_cut
  #go_id = as.data.frame(go_input) %>% dplyr::pull(unique(ID))
  ss_rel = GOSemSim::mgoSim(go_id, go_id, semData = GO_sem_sim, measure = "Rel", combine = NULL)
  set.seed(888) # for clustering reproducibility
  bcut_ss_rel = simplifyEnrichment::binary_cut(ss_rel, 
                                               try_all_partition_fun = TRUE
                                               ) %>% 
    factor(., levels = c(setdiff(names(sort(table(cl), decreasing = TRUE)), 0), 0))
  
  # Order terms according to binary_cut cluster partition (to give to complexHeatmap instead of dendrogram)
  term_order = unlist(lapply(levels(bcut_ss_rel), function(le) {
    l = bcut_ss_rel == le
    if(sum(l) <= 1) {
      return(which(l))
    } else {
      mm = ss_rel[l, l, drop = FALSE]
      which(l)[hclust(stats::dist(mm))$order]
    }
  }))
  
  # Make dataframe with GO and cluster number
  df = data.frame("id" = rownames(ss_rel),
                  "cluster" = bcut_ss_rel)
  
  # Obtain representative GO term for each cluster (to do: improve "importance" definition)
  repr_GO = GO_anc[go_id] %>% stack() %>% `colnames<-`(c("ancestral", "go_id")) %>% 
    left_join(dplyr::select(df, c("cluster", "id")), by = setNames("id", "go_id")) %>% 
    group_by(cluster, ancestral) %>% add_count() %>% ungroup %>% 
    mutate(IC = GO_semsim@IC[ancestral]) %>% 
    mutate(importance = (n*IC^2)) %>% 
    group_by(cluster) %>% 
    dplyr::slice(which.max(importance)) %>% 
    left_join(dplyr::select(GO_def, c("GOID","TERM")), by = setNames("GOID", "ancestral")) %>% 
    dplyr::select(c("ancestral","cluster","importance","TERM"))
  
  # Generate annotation with representative GO term
  df2 = df %>% tibble %>% left_join(dplyr::select(repr_GO, c("cluster", "TERM")), by ="cluster") %>% 
    rename("anc_term" = "TERM") %>% mutate(cluster = factor(cluster, levels = levels(bcut_ss_rel))) %>% 
    mutate(anc_term = factor(anc_term, levels = unique(.$anc_term[order(.$cluster)]) ))
  
  anno_df = df2[,c("id","anc_term")] %>% column_to_rownames("id")
  
  # Create distinct color palette for annotation
  set.seed(1)
  col_vector = distinctColorPalette(levels(anno_df$anc_term) %>% length)
  custom_colors = setNames(col_vector, levels(anno_df$anc_term))
  
  # Heatmap parameters (annotation, row labels, color function for matrix)
  ha = HeatmapAnnotation(df = anno_df,
                         which = "row",
                         annotation_label = "GO Cluster",
                         show_annotation_name = FALSE,
                         col = list(anc_term = custom_colors))
  
  row_labs = structure(paste(go_id, GO_data$Description), names = paste0(go_id))
  col_fun = colorRamp2(c(0,1), c("white", "red"))
  
  # Generate Heatmap
  GO_heatmap = Heatmap(ss_rel,
                       name = "Similarity",
                       row_order = term_order, column_order = term_order,
                       col = col_fun,
                       row_labels = row_labs[rownames(ss_rel)],
                       row_names_max_width = max_text_width(row_labs),
                       right_annotation = ha,
                       show_column_names = FALSE)
  
  return(GO_heatmap)
}