#' Plug-in estimator for unit quality.
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
estimate_eqwt_plugin <- function(data, a, y, w = NULL, e, mu, truncation_pt, adaptive = FALSE) {
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
  dplyr::summarise(if_data,
            plugin_est = mean(inf_fn),
            est_var = var(inf_fn)/dplyr::n())
}
