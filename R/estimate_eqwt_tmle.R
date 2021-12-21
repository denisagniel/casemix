#' TMLE estimator for unit quality.
#'
#' @param data data.frame containing the information to analyze
#' @param a string identifying the column in `data` that denotes the binary variable indicating the unit of interest
#' @param y string identifying the column in `data` that denotes the outcome of interest
#' @param w string identifying the column in `data` that denotes the equity weight
#' @param e string identifying the column in `data` that denotes the propensity score
#' @param mu string identifying the column in `data` that denotes the mean function
#' @param truncation_pt a number in (0,1) to which the propensity score is truncated
#' @param adaptive logical flag of whether to use the adaptive normalization of Khan and Ugander (2021)
#'
#' @export

estimate_eqwt_tmle <- function(data, a, y, w = NULL, e, mu, truncation_pt, adaptive = FALSE) {
  if (is.null(w)) {
    dplyr::mutate(data, w = 1)
    w <- 'w'
  }
  if_data <- dplyr::mutate(data,
                    inf_fn = eqwt_influence_fn(a = !!rlang::sym(a),
                                               y = !!rlang::sym(y),
                                               e_hat = !!rlang::sym(e),
                                               w_hat = !!rlang::sym(w),
                                               mu_hat = !!rlang::sym(mu),
                                               truncation_pt = truncation_pt,
                                               adaptive = adaptive
                    ))
  tmle_ds <- dplyr::mutate(if_data,
                    !!e := pmin(pmax(e, truncation_pt), 1-truncation_pt),
                    h = !!rlang::sym(a)/!!rlang::sym(e))
  tmle_fm <- as.formula(glue::glue('{y} ~ -1 + offset(qlogis({mu})) + h'))
  tmle_fit <- glm(tmle_fm, data = tmle_ds, family = binomial)
  eps <- coef(tmle_fit)
  tmle_ds <- dplyr::mutate(tmle_ds,
                           muhat_star = plogis(qlogis(!!rlang::sym(mu)) + eps*h))

  dplyr::summarise(tmle_ds,
            plugin_est = mean(inf_fn),
            tmle_est = mean(muhat_star),
            est_var = var(inf_fn)/dplyr::n())
}
