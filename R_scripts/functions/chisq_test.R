chisq_test = function(de_file,
                      peak_file,
                      logfc_threshold = 1,
                      target_region = NULL,
                      total_genes = NULL,
                      only_nuclear = TRUE,
                      gtf = "./input/Arabidopsis_thaliana.TAIR10.52.gtf.gz",
                      up_region = 3000,
                      down_region = 1000,
                      verbose = FALSE,
                      gene_id_column_in_de_file = "Geneid",
                      gene_id_column_in_target_region = "gene_id",
                      peak_name_column = "peak_name") {
  
  # create "target_region" to determine if a gene is bound if not provided
  if (is.null(target_region)) {
    print(paste0("target_region coordinates not provided, generating them from: ", gtf, 
                 " (", up_region," bp upstream and ", down_region, " bp downstream of tss"))
    
    target_region = plyranges::read_gff(gtf) %>% plyranges::filter(type == "gene") %>%
      plyranges::anchor_5p() %>% plyranges::mutate(width = 1) %>% plyranges::stretch((2*up_region)-1) %>%
      plyranges::anchor_5p() %>% plyranges::stretch(down_region-up_region) %>% 
      `mcols<-`(value = list("gene_id" = .$gene_id))
  }

  # obtain number of total genes if not provided
  if (is.null(total_genes)) {
    total_genes = plyranges::read_gff(gtf) %>% filter(type == "gene") %>% 
      as.data.frame() %>% pull(gene_id)
  }
  
  # take nuclear genes located in a chromosome
  if (only_nuclear == TRUE) {
    total_genes = total_genes %>% str_subset("AT\\d")
  }

  # calculate which genes are bound (have a peak in their target_region)
  bound_genes = plyranges::find_overlaps(peak_file, target_region) %>% as.data.frame() %>% 
    pull(get(gene_id_column_in_target_region))
  
  # retrieve which genes are differentially expressed according to logfc_threshold and the direction of expression
  logfc_column = de_file %>% colnames() %>% str_subset("^logFC")
  
  diff_exp_genes = de_file %>% as.data.frame() %>% pull(get(gene_id_column_in_de_file))
  up_genes = de_file %>% filter(get(logfc_column) > logfc_threshold) %>% pull(Geneid)
  dn_genes = de_file %>% filter(get(logfc_column) < -logfc_threshold) %>% pull(Geneid)
  
  # for debugging
  if(verbose == TRUE) {
    print(paste(length(total_genes),"total genes"))
    print(paste(length(diff_exp_genes),"differentially expressed genes"))
    print(paste(length(up_genes),"up-regulated genes"))
    print(paste(length(dn_genes),"down-regulated genes"))
    }
  
  # generate table of total genes with factors for each condition (e.g., differentially expressed, etc.)
  chi_input = data.frame("gene_id" = total_genes) %>% 
    mutate(bound = factor(ifelse(gene_id %in% bound_genes, "bound", "not_bound"),
                          levels = c("bound", "not_bound")),
           diff_exp = factor(ifelse(gene_id %in% diff_exp_genes, "differentially_expressed", "not_differentially_expressed"),
                             levels = c("differentially_expressed", "not_differentially_expressed")),
           up = factor(ifelse(gene_id %in% up_genes, "up", "not"),
                       levels = c("up", "not")),
           dn = factor(ifelse(gene_id %in% dn_genes, "dn", "not"),
                       levels = c("dn", "not")))
  
  # generate tables for chi-sq test (see https://statsandr.com/blog/chi-square-test-of-independence-in-r/)
  table_list = list(table_de = table(chi_input$bound, chi_input$diff_exp),
                    table_up = table(chi_input$bound, chi_input$up),
                    table_dn = table(chi_input$bound, chi_input$dn))
  
  # apply chisq.test function to each table and extract useful information (p-value, expected and observed values)
  test_list = lapply(table_list, chisq.test)
  
  expected = sapply(test_list, function(x) x$expected[1,1])
  observed = sapply(test_list, function(x) x$observed[1,1])
  p_value = sapply(test_list, function(x) x$p.value)
  significance = symnum(p_value, 
                        corr = FALSE, 
                        na = FALSE, 
                        cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                        symbols = c("***", "**", "*", ".", " ")) %>% as.vector()
  
  # generate output data.frame with useful information
  test_result = data.frame(gene_expression = c("Differentially expressed",
                                               "Up-regulated",
                                               "Down-regulated"),
                           expected,
                           observed,
                           p_value, 
                           significance,
                           row.names = NULL) %>% 
    mutate(direction = case_when(observed > expected ~ "enriched",
                                 observed < expected ~ "depleted",
                                 observed == expected ~ "same"),
           .after = observed) %>% 
    mutate(fold_over_expected = round(observed/expected, 2),
           .after = observed)
  
  return(c(table_list, test_result))
}