eqwt_influence_fn <- function(a, y, e_hat, w_hat = 1, mu_hat, truncation_pt = mean(a)/20, adaptive = FALSE) {
  e_hat <- pmin(e_hat, 1 - truncation_pt)
  e_hat <- pmax(e_hat, truncation_pt)
  if (adaptive) {
    k <- (1-e_hat)/e_hat*a/e_hat
    pi_hat <- mean(k)
    n_hat <- sum(a/e_hat)
    nn <- length(y)
    adj_term <- 1/pi_hat*mean((y - mu_hat)*k)*(1 - n_hat/nn)
  } else adj_term <- 0
  w_hat*(a*y/e_hat + (1 - a/e_hat)*mu_hat + adj_term)
}

