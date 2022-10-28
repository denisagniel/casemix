#' Workhorse function for TMLE. Intended to be wrapped by `popwt_tmle` and `eqwt_tmle`.
#'
#' @param data data.frame containing the information to analyze
#' @param folds optional string identifying the column in `data` that denotes the folds for cross-fitting
#' @param id string identifying the column in `data` that denotes the ID variable
#' @param a string identifying the column in `data` that denotes the units of interest
#' @param y string identifying the column in `data` that denotes the outcome of interest
#' @param wvars string vector identifying the columns in `data` that denote the protected characteristics
#' @param zvars string vector identifying the columns in `data` that denote the control (non-protected) characteristics
#' @param wt string vector identifying the column in `data` that denotes the observation weights
#' @param truncation_pt a number in (0,1) to which the propensity score is truncated (from above and below)
#' @param lrnr mlr3 learner object which will be used to estimate the mean function and propensity score
#' @param lrnr_e mlr3 learner object which will be used to estimate the propensity score (ignored if `lrnr` is specified)
#' @param lrnr_mu mlr3 learner object which will be used to estimate the mean function (ignored if `lrnr` is specified)
#' @param separate_e logical flag for whether propensity scores for each unit should be estimated separately or in a big multinomial model
#' @param separate_mu logical flag for whether mean functions for each unit should be estimated separately or in a big joint model
#' @param calibrate_e logical flag for whether propensity scores should be calibrated after fitting
#' @param calibrate_mu logical flag for whether mean functions should be calibrated after fitting
#' @param calibration_grps tuning parameter for calibration, lower numbers induce more smoothness, with the default being 500
#' @param condition_on a string indicating a variable within which to estimate conditional quality estimates.
#' @param verbose logical flag for whether to print informative messages about progress.
#'
estimate_wtd_tmle <- function(ds,
                                folds,
                                id,
                                a, y, wvars, zvars,
                                wt,
                                truncation_pt = 1e-7,
                              lrnr = lrn('classif.ranger'),
                              lrnr_e = NULL,
                              lrnr_mu = NULL,
                                separate_e = TRUE,
                                separate_mu = TRUE,
                                calibrate_e = FALSE,
                                calibrate_mu = FALSE,
                              calibration_grps = 500,
                              condition_on = NULL,
                              verbose = FALSE) {
  ####################
  ## estimate nuisance functions
  if (!is.null(lrnr)) {
    if (!is.null(lrnr_mu) | !is.null(lrnr_e)) stop('If `lrnr` is non-null, `lrnr_mu` and `lrnr_e` should be NULL.')
    ds <- left_join(ds,
                    estimate_e(ds, folds, id, c(wvars, zvars), a, lrnr, separate = separate_e, calibrate = calibrate_e, calibration_grps = calibration_grps, verbose = verbose), by = id)
    ds <- left_join(ds,
                    estimate_mu(ds, folds, id, c(wvars, zvars), y, a, lrnr, separate = separate_mu, calibrate = calibrate_mu, calibration_grps = calibration_grps, verbose = verbose), by = id)
  } else {
    ds <- left_join(ds,
                    estimate_e(ds, folds, id, c(wvars, zvars), a, lrnr_e, separate = separate_e, calibrate = calibrate_e, calibration_grps = calibration_grps, verbose = verbose), by = id)
    ds <- left_join(ds,
                    estimate_mu(ds, folds, id, c(wvars, zvars), y, a, lrnr_mu, separate = separate_mu, calibrate = calibrate_mu, calibration_grps = calibration_grps, verbose = verbose), by = id)
  }

  avals <- unique(pull(ds, all_of(a)))
  out <- map_df(avals, ~estimate_unit_tmle(data = ds,
                                           a = a,
                                           aval = .,
                                           y = y,
                                           w = wt,
                                           e = glue('e_{.}'),
                                           mu = glue('mu_{.}'),
                                           truncation_pt = truncation_pt,
                                           condition_on = condition_on)
  )


  ##########################
  ## calculate reliability and shrunken estimators
  shrink_estimates(out, 'tmle_est', 'tmle_se')
}
