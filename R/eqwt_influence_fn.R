eff_influence_fn <- function(a, y, e_hat, mu_hat, truncation_pt = mean(a)/20, adaptive = FALSE) {
  e_hat <- pmin(e_hat, 1 - truncation_pt)
  e_hat <- pmax(e_hat, truncation_pt)
  if (adaptive) {
    k <- (1-e_hat)/e_hat*a/e_hat
    pi_hat <- mean(k)
    n_hat <- sum(a/e_hat)
    nn <- length(y)
    adj_term <- 1/pi_hat*mean((y - mu_hat)*k)*(1 - n_hat/nn)
  } else adj_term <- 0
  a*y/e_hat + (1 - a/e_hat)*mu_hat + adj_term
}

influence_fn_given_w <- function(a, y, w, w_star, z, e_fn, mu_fn, truncation_pt = mean(a)/20, adaptive = FALSE) {
  eff_influence_fn(a = a*(w == w_star),
                   y = y,
                   e_hat = e_fn(w_star, z),
                   mu_hat = mu_fn(w_star, z),
                   truncation_pt = truncation_pt,
                   adaptive = adaptive)
}

eqwt_influence_fn <- function(a, y, w, support_W = unique(w), z, e_fn, mu_fn, truncation_pt = mean(a)/20, adaptive = FALSE) {
  dplyr::bind_cols(purrr::map(support_W, ~influence_fn_given_w(a, y, w, ., z, e_fn, mu_fn, truncation_pt, adaptive)))
}
