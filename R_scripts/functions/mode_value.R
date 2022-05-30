mode_value = function(x) {
  ux = unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}