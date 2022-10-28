#' Targeted minimum loss-based estimation for unit quality using sample splitting, standardized to a population satsifying certain equity properties.
#'
#' This function is used for estimating many causal effects, all standardized to a common population distribution that satisfies certain equity properties, with a primary application to quality measurement, where many units (e.g., doctors, hospitals, schools) are measured repeatedly (e.g., on different patients), and we want to know what the average outcome would be for each unit, if they were measured on a popuation that had equal representation for observations with certain characteristics.
#'
#' @param data data.frame containing the information to analyze
#' @param folds optional string identifying the column in `data` that denotes the folds for cross-fitting
#' @param id string identifying the column in `data` that denotes the ID variable
#' @param a string identifying the column in `data` that denotes the units of interest
#' @param y string identifying the column in `data` that denotes the outcome of interest
#' @param wvars string vector identifying the columns in `data` that denote the protected characteristics
#' @param zvars string vector identifying the columns in `data` that denote the control (non-protected) characteristics
#' @param truncation_pt a number in (0,1) to which the propensity score is truncated (from above and below); default is 5/sqrt(n)/log(n)
#' @param K optional integer identifying the number of folds for cross-fitting; required if folds = NULL
#' @param lrnr mlr3 learner object which will be used to estimate the mean function and propensity score
#' @param lrnr_e mlr3 learner object which will be used to estimate the propensity score (ignored if `lrnr` is specified)
#' @param lrnr_mu mlr3 learner object which will be used to estimate the mean function (ignored if `lrnr` is specified)
#' @param separate_e logical flag for whether propensity scores for each unit should be estimated separately or in a big multinomial model
#' @param separate_mu logical flag for whether mean functions for each unit should be estimated separately or in a big joint model
#' @param epsilon positive scalar that indicates the amount that the optimization of equity balance constraints is allowed to deviate from the required constraints
#' @param calibrate_e logical flag for whether propensity scores should be calibrated after fitting
#' @param calibrate_mu logical flag for whether mean functions should be calibrated after fitting
#' @param calibration_grps tuning parameter for calibration, lower numbers induce more smoothness, with the default being 500
#' @param condition_on a string indicating a variable within which to estimate conditional quality estimates.
#'
#' @param verbose logical flag for whether to print informative message about progress.
#'
#' @return A `tibble` with the following columns:\itemize{
#'   \item a column with the same name as the argument \code{a} which indicates the unit.
#'   \item \code{tmle_est}: the TMLE estimate.
#'   \item \code{tmle_se}: the estimated standard error for the TMLE.
#'   \item \code{shrinkage_est}: an estimate that is shrunk using empirical bayes via the \code{ebnm} package.
#'   \item \code{shrinkage_est_se}: the estimated standard error for the shrinkage estimate.
#'   \item \code{reliability}: the reliability of the shrinkage estimate.
#' }
#' @export
#'
#' @import mlr3 mlr3learners dplyr rlang
#' @importFrom tibble tibble
#' @importFrom glue glue
#' @importFrom progressr progressor
#' @importFrom stringr str_remove
#' @importFrom ebnm ebnm_normal
#' @importFrom limSolve ldp
#' @importFrom quadprog solve.QP
#' @importFrom purrr pmap_df map_df map reduce
#'
eqwt_tmle <- function(data,
                      folds = NULL,
                      id,
                      a, y, wvars, zvars,
                      truncation_pt = NULL,
                      K = 2,
                      lrnr = lrn('classif.ranger'),
                      lrnr_e = NULL,
                      lrnr_mu = NULL,
                      separate_e = FALSE,
                      separate_mu = FALSE,
                      epsilon = 1e-12,
                      calibrate_e = TRUE,
                      calibrate_mu = TRUE,
                      calibration_grps = 500,
                      condition_on = NULL,
                      verbose = TRUE) {
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
                            calibration_grps = calibration_grps,
                            condition_on = condition_on,
                            sub_k = K,
                            verbose = verbose)

  eqwt_ds <- unique(select(ds, all_of(wvars), equity_wt))

  out <- estimate_wtd_tmle(ds, folds = folds,
                             id = id,
                             wvars = wvars,
                             zvars = zvars,
                           wt = 'equity_wt',
                             a = a,
                             y = y,
                           lrnr = lrnr,
                           lrnr_e = lrnr_e,
                           lrnr_mu = lrnr_mu,
                             separate_e = separate_e,
                             separate_mu = separate_mu,
                             calibrate_e = calibrate_e,
                             calibrate_mu = calibrate_mu,
                           calibration_grps = calibration_grps,
                           condition_on = condition_on,
                             truncation_pt = truncation_pt,
                           verbose = verbose)



  ##########################
  ## calculate reliability and shrunken estimators
  mutate(out,
         wt_ds = list(eqwt_ds))
}
