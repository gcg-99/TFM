get_consensus_peaks_GRangesList = function(granges_list,
                                           peak_names,
                                           rep_number = 2) {
  
  peak_coverage = GenomicRanges::GRangesList(granges_list) %>% GenomicRanges::coverage()
  
  covered_ranges = IRanges::slice(peak_coverage, lower=rep_number, rangesOnly=T) %>% 
    GenomicRanges::GRanges() %>% 
    plyranges::mutate(peak_name = paste0("peak_",rep(1:length(.)))) #%>% 
  #plyranges::mutate(name = gene_id) %>% 
  #plyranges::mutate(type = rep("gene",length(.)))
  
  return(covered_ranges)
}