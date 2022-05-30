get_goids = function(GOenrich) {
  GOenrich %>% as.data.frame() %>% pull(ID)
}