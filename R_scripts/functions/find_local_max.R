find_local_max <- function(x,
                           partial = TRUE,
                           decreasing = FALSE) {
  if (decreasing) {
    if (partial) {
      which(diff(c(FALSE, diff(x) > 0, TRUE)) > 0)
    } else {
      which(diff(diff(x) > 0) > 0) + 1
    }
  } else {
    if (partial) {
      which(diff(c(TRUE, diff(x) >= 0, FALSE)) < 0)
    } else {
      which(diff(diff(x) >= 0) < 0) + 1
    }
  }
}
