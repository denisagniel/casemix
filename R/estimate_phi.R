estimate_phi <- function(data, a, y, truncation_pt = 1e-7, adaptive = FALSE) {
  avals <- unique(pull(data, !!a))
  mu_nms <- paste0('mu_', avals)
  e_nms <- paste0('e_', avals)
  if (!all(mu_nms %in% colnames(data))) stop("mean functions have not been estimated yet. please see `estimate_mu` function.")
  if (!all(e_nms %in% colnames(data))) stop("propensity scores have not been estimated yet. please see `estimate_e` function.")
  updated_data <- data
  for (aa in avals) {
    this_mu <- paste0('mu_', aa)
    this_e <- paste0('e_', aa)
    this_phi <- paste0('phi_', aa)
    this_a <- paste0('a_', aa)
    updated_data <- mutate(updated_data,
                           !!this_a := 1*(!!sym(a) == aa),
                           !!this_phi := wtd_influence_fn(a = !!sym(this_a),
                                                          y = !!sym(y),
                                                          e_hat = !!sym(this_e),
                                                          mu_hat = !!sym(this_mu),
                                                          w_hat = 1,
                                                          truncation_pt = truncation_pt,
                                                          adaptive = adaptive))
  }
  updated_data
}
