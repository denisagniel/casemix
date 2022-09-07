#' TMLE for unit quality.
#'
#' @param data data.frame containing the information to analyze
#' @param folds optional string identifying the column in `data` that denotes the folds for cross-fitting
#' @param id string identifying the column in `data` that denotes the ID variable
#' @param a string identifying the column in `data` that denotes the units of interest
#' @param y string identifying the column in `data` that denotes the outcome of interest
#' @param covars string vector identifying the columns in `data` that denote the control characteristics
#' @param truncation_pt a number in (0,1) to which the propensity score is truncated (from above and below); default is 5/sqrt(n)/log(n)
#' @param adaptive logical flag of whether to use the adaptive normalization of Khan and Ugander (2021)
#' @param K optional integer identifying the number of folds for cross-fitting; required if folds = NULL
#' @param lrnr mlr3 learner object which will be used to estimate the mean function and propensity score
#' @param separate_e logical flag for whether propensity scores for each unit should be estimated separately or in a big multinomial model
#' @param separate_mu logical flag for whether mean functions for each unit should be estimated separately or in a big joint model
#' @param calibrate_e logical flag for whether propensity scores should be calibrated after fitting
#' @param calibrate_mu logical flag for whether mean functions should be calibrated after fitting
#'
#' @export
popwt_tmle <- function(data,
                       folds = NULL,
                       id,
                       a, y,
                       covars,
                       truncation_pt = 1e-7,
                       adaptive = FALSE,
                       K = 2,
                       lrnr = lrn('classif.ranger'),
                       separate_e = FALSE,
                       separate_mu = FALSE,
                       calibrate_e = TRUE,
                       calibrate_mu = TRUE) {
  if (is.null(folds)) {
    ds <- make_folds(data, a, K)
    folds <- 'fold'
  } else ds <- data


  ds <- ds %>%
    mutate(wt = rep(1, nrow(ds)))
  out <- estimate_wtd_tmle(ds, folds = folds,
                             id = id,
                             wt = 'wt',
                           wvars = covars,
                           zvars = NULL,
                             a = a,
                             y = y,
                             lrnr = lrnr,
                             separate_e = separate_e,
                             separate_mu = separate_mu,
                             calibrate_e = calibrate_e,
                             calibrate_mu = calibrate_mu,
                             truncation_pt = truncation_pt)
  shrink_estimates(out, 'tmle_est', 'est_se')
}
