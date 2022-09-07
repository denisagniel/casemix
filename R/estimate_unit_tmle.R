#' TMLE estimator for unit quality.
#'
#' @param data data.frame containing the information to analyze
#' @param a string identifying the column in `data` that denotes the units of interest
#' @param aval string identifying the particular unit being estimated
#' @param y string identifying the column in `data` that denotes the outcome of interest
#' @param w string identifying the column in `data` that denotes the equity weight
#' @param e string identifying the column in `data` that denotes the propensity score
#' @param mu string identifying the column in `data` that denotes the mean function
#'
#' @export

estimate_unit_tmle <- function(data, a, aval, y, w = NULL, e, mu, truncation_pt = 1e-12) {
  if (is.null(w)) {
    data <- dplyr::mutate(data, w = 1)
    w <- 'w'
  }
  if (is.null(truncation_pt)) {
    n <- nrow(data)
    truncation_pt <- 5/sqrt(n)/log(n)
  }
  if_data <- dplyr::mutate(data,
                           inf_fn = wtd_influence_fn(a = 1*(!!rlang::sym(a) == aval),
                                                     y = !!rlang::sym(y),
                                                     e_hat = !!rlang::sym(e),
                                                     w_hat = !!rlang::sym(w),
                                                     mu_hat = !!rlang::sym(mu),
                                                     truncation_pt = truncation_pt,
                                                     adaptive = FALSE
                           ))
  tmle_ds <- dplyr::mutate(if_data,
                           h = 1*(!!rlang::sym(a) == aval)/pmax(!!rlang::sym(e), truncation_pt)*!!rlang::sym(w),
                           mu_h = pmin(pmax(truncation_pt, !!rlang::sym(mu)), 1 - truncation_pt))
  tmle_fm <- as.formula(glue::glue('{y} ~ h + offset(qlogis(mu_h))'))
  if (any(round(pull(tmle_ds, sym(w)), 9)) < 0) browser()
  tmle_fit <- glm(tmle_fm, data = tmle_ds, family = binomial)
  eps <- coef(tmle_fit)[2]
  tmle_ds <- dplyr::mutate(tmle_ds,
                           muhat_star = plogis(qlogis(mu_h) + eps*h))

  out <- dplyr::summarise(tmle_ds,
                   !!a := aval,
                   plugin_est = mean(inf_fn),
                   tmle_est = mean(muhat_star*!!rlang::sym(w)),
                   est_se = sqrt(var(inf_fn)/dplyr::n()))
  # if (out$est_se > 0.1) browser()
  out
}
