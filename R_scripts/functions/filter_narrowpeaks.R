filter_narrowpeaks = function(file_to_filter,
                              qvalue_threshold = 0.05,
                              peak_min_width = 10,
                              peak_max_width = Inf,
                              filtered_file) {
  
  read_narrowpeaks(file_to_filter) %>% filter((qValue > -log10(qvalue_threshold) & 
                                                  width > peak_min_width &
                                                  width < peak_max_width)) %>%
    write_narrowpeaks(., file = filtered_file)
}
