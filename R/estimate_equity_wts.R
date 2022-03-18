#' Estimate weights that obey marginal balance constraints on protected characteristics.
#'
#' These weights ensure that each characteristic represented by a variable in `wvars` is balanced in the weighted population.
#'
#' @param data data.frame containing the information to analyze
#' @param a optional string identifying the column in `data` that denotes the units of interest; if not provided, must provide `avals`
#' @param avals optional string identifying unique values of `a` found in the data; if not provided, must provide `a`
#' @param wvars string vector identifying the columns in `data` that denote the protected characteristics
#' @param epsilon positive scalar that indicates the amount that the optimization is allowed to deviate from the required constraints
#'
#'
#' @export
estimate_equity_wts <- function(data, folds, id, wvars, zvars, a, y, lrnr, separate_e = TRUE, separate_mu = TRUE, epsilon = 1e-12, tune = FALSE, evals = 20, callibrate_e = TRUE, callibrate_mu = TRUE, sub_k = 2, truncation_pt = 1e-7) {
  all_folds <- pull(data, folds)
  unique_folds <- unique(all_folds)

  map_df(unique_folds,
                        ~fold_wt(data,
                                 folds,
                                 out_fold = .,
                                 id,
                                 wvars,
                                 zvars,
                                 a,
                                 y,
                                 lrnr = lrnr,
                                 separate_e,
                                 separate_mu,
                                 epsilon,
                                 tune,
                                 evals,
                                 callibrate_e,
                                 callibrate_mu,
                                 K = sub_k))
}

fold_wt <- function(data,
                    folds,
                    out_fold,
                    id,
                    wvars,
                    zvars,
                    a,
                    y,
                    lrnr,
                    separate_e,
                    separate_mu,
                    epsilon,
                    tune,
                    evals,
                    callibrate_e,
                    callibrate_mu,
                    K = 2) {
  train_data <- data %>%
    filter(!!sym(folds) != out_fold)
  out_data <- data %>%
    filter(!!sym(folds) == out_fold)
  #######################
  ## create sub-folds on training data
  train_data <- make_folds(train_data, unit_id = a, K = K, fold_name = 'subfold')


  #######################
  ## estimate mu and e
  train_data <- left_join(train_data,
                  estimate_e(train_data, 'subfold', id, c(wvars, zvars), a, lrnr, separate = separate_e, tune = tune, callibrate = callibrate_e), by = id)
  train_data <- left_join(train_data,
                  estimate_mu(train_data, 'subfold', id, c(wvars, zvars), y, a, lrnr, separate = separate_mu, tune = tune, callibrate = callibrate_mu), by = id)
  train_data <- estimate_phi(train_data, a, y)
  avals <- unique(dplyr::pull(train_data, !!sym(a)))
  train_data <- calculate_wts(train_data, avals = avals, wvars = wvars, epsilon = epsilon)
  eqwt_ds <- train_data %>% select(all_of(wvars), equity_wt) %>% unique
  out_data %>%
    inner_join(eqwt_ds, by = wvars)
}
