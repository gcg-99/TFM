get_direct_targets = function(peak_file,
                              de_file,
                              target_region = NULL,
                              gtf = "./input/Arabidopsis_thaliana.TAIR10.52.gtf.gz",
                              up_region = 3000,
                              down_region = 1000,
                              gene_id_column_in_de_file = "Geneid",
                              gene_id_column_in_target_region = "gene_id",
                              peak_name_column = "peak_name",
                              output = c("vector", "filtered"),
                              exp_name = "")
{
  output = match.arg(output)
  
  if (is.null(target_region)) {
    print(paste0("target_region coordinates not provided, generating them from: ", gtf, 
                 " (", up_region," bp upstream and ", down_region, " bp downstream of tss"))
    
    target_region = plyranges::read_gff(gtf) %>% plyranges::filter(type == "gene") %>%
      plyranges::anchor_5p() %>% plyranges::mutate(width = 1) %>% plyranges::stretch((2*up_region)-1) %>%
      plyranges::anchor_5p() %>% plyranges::stretch(down_region-up_region) %>% 
      `mcols<-`(value = list("gene_id" = .$gene_id))
  }
  
  peaks_in_genes = plyranges::find_overlaps(peak_file, target_region) %>% 
    plyranges::mutate(combined_name = paste0(get(gene_id_column_in_target_region),"_",get(peak_name_column))) %>% 
    `names<-`(.$combined_name)
  
  # To do: add distance to corresponding tss
  # split peaks_in_genes, obtain corresponding tss to the gene, calculate distance, and merge into granges
  
  gene_ids_of_bound_genes = peaks_in_genes %>% as.data.frame() %>% dplyr::pull(get(gene_id_column_in_target_region))
  
  direct_targets = de_file %>% dplyr::filter(get(gene_id_column_in_de_file) %in% gene_ids_of_bound_genes)
  
  n_direct_targets = length(row.names(direct_targets))
  n_de = length(row.names(de_file))
  
  print(paste0(n_direct_targets, " out of ", n_de, 
               " differentially expressed genes are potential direct targets of ", exp_name, " (" ,
               round(100*n_direct_targets/n_de, 2),"%)"))
  
  direct_target_ids = direct_targets %>% pull(get(gene_id_column_in_de_file))
  
  if(output == "vector") {
    return(direct_target_ids)
  } else {
    return(list("direct_targets" = direct_targets, "peaks_in_genes" = peaks_in_genes))
  }
}