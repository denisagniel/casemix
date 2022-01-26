#' Plug-in estimator for unit quality.
#'
#' @param data data.frame containing the information to analyze
#' @param folds optional string identifying the column in `data` that denotes the folds for cross-fitting
#' @param id string identifying the column in `data` that denotes the ID variable
#' @param a string identifying the column in `data` that denotes the binary variable indicating the unit of interest
#' @param y string identifying the column in `data` that denotes the outcome of interest
#' @param wvars string vector identifying the columns in `data` that denote the protected characteristics
#' @param zvars string vector identifying the columns in `data` that denote the control (non-protected) characteristics
#' @param truncation_pt a number in (0,1) to which the propensity score is truncated
#' @param adaptive logical flag of whether to use the adaptive normalization of Khan and Ugander (2021)
#' @param K optional integer identifying the number of folds for cross-fitting; required if folds = NULL
#' @param lrnr mlr3 learner object which will be used to estimate the mean function and propensity
#' @param epsilon positive scalar that indicates the amount that the optimization is allowed to deviate from the required constraints
#'
#' @export
estimate_eqwt_plugin <- function(data,
                                 folds = NULL,
                                 id,
                                 a, y, wvars, zvars,
                                 truncation_pt = 1e-7,
                                 adaptive = FALSE,
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
  ds <- estimate_equity_wts(ds, avals, wvars, epsilon)
  phis <- paste0('phi_', avals)
  ests <- dplyr::summarise_at(ds,
                      all_of(phis),
            list(pluginest = ~mean(.*equity_wt),
            se = ~sqrt(var(.*equity_wt)/dplyr::n())))
  pivot_longer(ests, cols = everything()) %>%
    separate(name, into = c('phi', a, 'type'), sep = '_') %>%
    select(-phi) %>%
    pivot_wider(id_cols = a, names_from = type, values_from = value)

}
