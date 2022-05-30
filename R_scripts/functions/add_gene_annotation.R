add_gene_annotation = function(edgeR_output,
                               gene_id_column = "Geneid",
                               tidy_output = TRUE) {
  # Load prerequisite libraries and generate ensemble plants biomaRt object 
  if(!exists("ensembl")){
    ensembl = biomaRt::useMart(biomart="plants_mart",
                               dataset="athaliana_eg_gene",
                               host="https://plants.ensembl.org")
    }
  edgeR_output = as.data.frame(edgeR_output)
  gene_list = edgeR_output %>% dplyr::pull(gene_id_column)
  print(paste0(length(gene_list)," genes detected"))
  
  # Get gene symbol from org.At.tair.db (separated by comma if there are multiple symbols)
  print("Retrieving gene symbol")
  gene_symbol = suppressMessages(AnnotationDbi::mapIds(org.At.tair.db,
                                                       column = "SYMBOL",
                                                       keys = gene_list,
                                                       keytype = "TAIR",
                                                       multiVals = "list")) %>%
    map_chr(., ~str_c(.x, collapse=", "))
  
  # Get long gene description from org.At.tair.db (only one description is retrieved)
  print("Retrieving long gene description")
  long_description = suppressMessages(AnnotationDbi::mapIds(org.At.tair.db,
                                                            column = "GENENAME",
                                                            keys = gene_list,
                                                            keytype = "TAIR",
                                                            multiVals = "first"))
  
  # Get short gene description from ensemble plants and remove [Source:IDXXX]
  print("Retrieving short gene description")
  short_description = biomaRt::getBM(attributes = c("tair_locus","description"),
                                     filters = "tair_locus",
                                     values = gene_list,
                                     mart = ensembl) %>%
    dplyr::rename(short_description = description) %>%
    dplyr::mutate(short_description = gsub(" \\[Source.*","",short_description))
  
  # Write output table, if no gene symbol, add gene id
  print("Generating output table")
  annotated_edgeR_output = edgeR_output %>%
    mutate(gene_symbol = dplyr::coalesce(gene_symbol,get(gene_id_column)), .after=all_of(gene_id_column),
           long_description = long_description) %>%
    dplyr::left_join(short_description, by = setNames("tair_locus",all_of(gene_id_column))) %>%
    dplyr::relocate(short_description, .after = gene_symbol) %>%
    { if(tidy_output == TRUE) dplyr::select(., -c(Chr,Start,End,Strand,Length)) else . }
  # Return output
  return(annotated_edgeR_output)
}
