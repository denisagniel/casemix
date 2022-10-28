#' Doubly-robust estimator for unit quality, using plug-in influence-function-based estimation and sample splitting, as in double machine learning.
#'
#' This function is used for estimating many causal effects, all standardized to a common population distribution, with a primary application to quality measurement, where many units (e.g., doctors, hospitals, schools) are measured repeatedly (e.g., on different patients), and we want to know what the average outcome would be for each unit, if they were measured on the whole popuation (e.g., if each doctor saw all patients).
#'
#' @param data data.frame containing the information to analyze
#' @param folds optional string identifying the column in `data` that denotes the folds for cross-fitting
#' @param id string identifying the column in `data` that denotes the ID variable
#' @param a string identifying the column in `data` that denotes the units of interest
#' @param y string identifying the column in `data` that denotes the outcome of interest
#' @param covars string vector identifying the columns in `data` that denote the characteristics to control for confounding
#' @param truncation_pt a number in (0,1) to which the propensity score is truncated (from above and below); default is 5/sqrt(n)/log(n)
#' @param K optional integer identifying the number of folds for cross-fitting; required if folds = NULL
#' @param lrnr mlr3 learner object which will be used to estimate the mean function and propensity score
#' @param lrnr_e mlr3 learner object which will be used to estimate the propensity score (ignored if `lrnr` is specified)
#' @param lrnr_mu mlr3 learner object which will be used to estimate the mean function (ignored if `lrnr` is specified)
#' @param separate_e logical flag for whether propensity scores for each unit should be estimated separately or in a big multinomial model
#' @param separate_mu logical flag for whether mean functions for each unit should be estimated separately or in a big joint model
#' @param calibrate_e logical flag for whether propensity scores should be calibrated after fitting
#' @param calibrate_mu logical flag for whether mean functions should be calibrated after fitting
#' @param calibration_grps tuning parameter for calibration, lower numbers induce more smoothness, with the default being 500
#' @param infl_fn_only logical flag for whether to return a tibble with the influence function for each observation and each unit (which can be used to calculate contrasts between units).
#' @param condition_on a string indicating a variable within which to estimate conditional quality estimates.
#' @param contrasts tibble or data.frame encoding pre-specified contrasts between units to return
#' @param verbose logical flag for whether to print informative messages about progress.
#' @return If \code{infl_fn_only = FALSE}, a `tibble` with the following columns:\itemize{
#'   \item a column with the same name as the argument \code{a} which indicates the unit.
#'   \item \code{plugin_est}: the plug-in estimate.
#'   \item \code{plugin_se}: the estimated standard error for the plug-in.
#'   \item \code{shrinkage_est}: an estimate that is shrunk using empirical bayes via the \code{ebnm} package.
#'   \item \code{shrinkage_est_se}: the estimated standard error for the shrinkage estimate.
#'   \item \code{reliability}: the reliability of the shrinkage estimate.
#' }. Else if \code{infl_fn_only = TRUE}, an augmented version of \code{data}, with influence functions for each unit denoted by columns named \code{phi_} followed by the unit name.
#' @export
#'
#' @import mlr3 mlr3learners dplyr rlang
#' @importFrom tibble tibble
#' @importFrom glue glue
#' @importFrom progressr progressor
#' @importFrom stringr str_remove
#' @importFrom ebnm ebnm_normal
#' @importFrom purrr pmap_df map_df map reduce
#'
popwt_plugin <- function(data,
                        folds = NULL,
                        id,
                        a, y, covars,
                        truncation_pt = NULL,
                        K = 2,
                        lrnr = NULL,
                        lrnr_e = NULL,
                        lrnr_mu = NULL,
                        separate_e = FALSE,
                        separate_mu = FALSE,
                        calibrate_e = TRUE,
                        calibrate_mu = TRUE,
                        calibration_grps = 500,
                        infl_fn_only = FALSE,
                        condition_on = NULL,
                        contrasts = NULL,
                        verbose = TRUE) {
  if (is.null(folds)) {
    ds <- make_folds(data, a, K)
    folds <- 'fold'
  } else ds <- data

  ds <- ds %>%
    mutate(wt = rep(1, nrow(ds)),
           !!sym(a) := as.character(!!sym(a)))
  out <- estimate_wtd_plugin(ds, folds = folds,
                             id = id,
                             wt = 'wt',
                             wvars = covars,
                             zvars = NULL,
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
                             truncation_pt = truncation_pt,
                             infl_fn_only = infl_fn_only,
                             condition_on = condition_on,
                             contrasts = contrasts,
                             verbose = verbose)
  out
}
