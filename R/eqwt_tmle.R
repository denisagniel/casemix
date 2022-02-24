#' TMLE estimator for unit quality.
#'
#' @param data data.frame containing the information to analyze
#' @param folds optional string identifying the column in `data` that denotes the folds for cross-fitting
#' @param id string identifying the column in `data` that denotes the ID variable
#' @param a string identifying the column in `data` that denotes the units of interest
#' @param y string identifying the column in `data` that denotes the outcome of interest
#' @param wvars string vector identifying the columns in `data` that denote the protected characteristics
#' @param zvars string vector identifying the columns in `data` that denote the control (non-protected) characteristics
#' @param K optional integer identifying the number of folds for cross-fitting; required if folds = NULL
#' @param lrnr mlr3 learner object which will be used to estimate the mean function and propensity
#' @param epsilon positive scalar that indicates the amount that the optimization is allowed to deviate from the required constraints
#'
#' @export
#'
eqwt_tmle <- function(data,
                                 folds = NULL,
                                 id,
                                 a, y, wvars, zvars,
                                 K = 2,
                                 lrnr = lrn('classif.ranger'),
                                 epsilon = 1e-12) {
  if (is.null(folds)) {
    ds <- make_folds(data, a, K)
    folds <- 'fold'
  } else ds <- data
  ds <- left_join(ds,
                  estimate_e(ds, folds, id, c(wvars, zvars), a, lrnr))
  ds <- left_join(ds,
                  estimate_mu(ds, folds, id, c(wvars, zvars), y, a, lrnr))
  ds <- estimate_phi(ds, a, y)
  avals <- ds %>% pull(!!sym(a)) %>% unique
  ds <- estimate_equity_wts(ds, avals = avals, wvars = wvars, epsilon = epsilon)
  eqwt_ds <- ds %>% select(all_of(wvars), equity_wt, alt_wt1, alt_wt2) %>% unique

  out <- map_df(avals, ~estimate_eqwt_tmle(data = ds,
                                    a = a,
                                    aval = .,
                                    y = y,
                                    w = 'equity_wt',
                                    e = glue('e_{.}'),
                                    mu = glue('mu_{.}'))
  )


  ##########################
  ## calculate reliability and shrunken estimators
  out <- shrink_estimates(out, 'tmle_est', 'est_se')
  mutate(out,
         wt_ds = list(eqwt_ds))
}
