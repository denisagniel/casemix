#' Plug-in estimator for equity-weighted unit quality.
#'
#' @param data data.frame containing the information to analyze
#' @param folds optional string identifying the column in `data` that denotes the folds for cross-fitting
#' @param id string identifying the column in `data` that denotes the ID variable
#' @param a string identifying the column in `data` that denotes the units of interest
#' @param y string identifying the column in `data` that denotes the outcome of interest
#' @param wvars string vector identifying the columns in `data` that denote the protected characteristics
#' @param zvars string vector identifying the columns in `data` that denote the control (non-protected) characteristics
#' @param truncation_pt a number in (0,1) to which the propensity score is truncated (from above and below); default is 5/sqrt(n)/log(n)
#' @param adaptive logical flag of whether to use the adaptive normalization of Khan and Ugander (2021)
#' @param K optional integer identifying the number of folds for cross-fitting; required if folds = NULL
#' @param lrnr mlr3 learner object which will be used to estimate the mean function and propensity score
#' @param lrnr_e mlr3 learner object which will be used to estimate the propensity score (ignored if `lrnr` is specified)
#' @param lrnr_mu mlr3 learner object which will be used to estimate the mean function (ignored if `lrnr` is specified)
#' @param separate_e logical flag for whether propensity scores for each unit should be estimated separately or in a big multinomial model
#' @param separate_mu logical flag for whether mean functions for each unit should be estimated separately or in a big joint model
#' @param epsilon positive scalar that indicates the amount that the optimization of equity balance constraints is allowed to deviate from the required constraints
#' @param calibrate_e logical flag for whether propensity scores should be calibrated after fitting
#' @param calibrate_mu logical flag for whether mean functions should be calibrated after fitting
#'
#' @export
#'
#' @import mlr3 mlr3learners mlr3extralearners dplyr purrr rlang
#' @importFrom tibble tibble
#' @importFrom glue glue
#' @importFrom progressr progressor
#' @importFrom stringr str_remove
#' @importFrom ebnm ebnm_normal
#'
eqwt_plugin <- function(data,
                        folds = NULL,
                        id,
                        a, y, wvars, zvars,
                        truncation_pt = NULL,
                        adaptive = FALSE,
                        K = 2,
                        lrnr = lrn('classif.ranger'),
                        lrnr_e = NULL,
                        lrnr_mu = NULL,
                        separate_e = FALSE,
                        separate_mu = FALSE,
                        epsilon = 1e-12,
                        calibrate_e = TRUE,
                        calibrate_mu = TRUE) {
  if (is.null(folds)) {
    ds <- make_folds(data, a, K)
    folds <- 'fold'
  } else ds <- data


  ds <- estimate_equity_wts(ds, folds = folds,
                            id = id,
                            wvars = wvars,
                            zvars = zvars,
                            a = a,
                            y = y,
                            lrnr = lrnr,
                            lrnr_e = lrnr_e,
                            lrnr_mu = lrnr_mu,
                            separate_e = separate_e,
                            separate_mu = separate_mu,
                            epsilon = epsilon,
                            calibrate_e = calibrate_e,
                            calibrate_mu = calibrate_mu,
                            sub_k = K)

  eqwt_ds <- unique(select(ds, all_of(wvars), equity_wt))
  out <- estimate_wtd_plugin(ds, folds = folds,
                             id = id,
                             wt = 'equity_wt',
                             wvars = wvars,
                             zvars = zvars,
                             a = a,
                             y = y,
                             lrnr = lrnr,
                             lrnr_e = lrnr_e,
                             lrnr_mu = lrnr_mu,
                             separate_e = separate_e,
                             separate_mu = separate_mu,
                             calibrate_e = calibrate_e,
                             calibrate_mu = calibrate_mu,
                             truncation_pt = truncation_pt)
  mutate(out,
         wt_ds = list(eqwt_ds))
}
