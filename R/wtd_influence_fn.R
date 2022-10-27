wtd_influence_fn <- function(a, y, e_hat, mu_hat, w_hat, truncation_pt = NULL) {
  if (is.null(truncation_pt)) {
    n <- length(y)
    truncation_pt <- 5/sqrt(n)/log(n)
  }
  e_hat <- pmin(e_hat, 1 - truncation_pt)
  e_hat <- pmax(e_hat, truncation_pt)
  w_hat*(a*y/e_hat + (1 - a/e_hat)*mu_hat)
}
