#' Estimate weights that obey marginal balance constraints on protected characteristics.
#'
#' These weights ensure that each characteristic represented by a variable in `wvars` is balanced in the weighted population.
#'
#' @param data data.frame containing the information to analyze
#' @param folds optional string identifying the column in `data` that denotes the folds for cross-fitting
#' @param id string identifying the column in `data` that denotes the ID variable
#' @param wvars string vector identifying the columns in `data` that denote the protected characteristics
#' @param zvars string vector identifying the columns in `data` that denote the control (non-protected) characteristics
#' @param a optional string identifying the column in `data` that denotes the units of interest (i.e., the thing that is being ranked in terms of quality, like health care providers or schools)
#' @param y string identifying the column in `data` that denotes the outcome of interest
#' @param lrnr mlr3 learner object which will be used to estimate the mean function and propensity score
#' @param lrnr_e mlr3 learner object which will be used to estimate the propensity score (ignored if `lrnr` is specified)
#' @param lrnr_mu mlr3 learner object which will be used to estimate the mean function (ignored if `lrnr` is specified)
#' @param separate_e logical flag for whether propensity scores for each unit should be estimated separately or in a big multinomial model
#' @param separate_mu logical flag for whether mean functions for each unit should be estimated separately or in a big joint model
#' @param epsilon positive scalar that indicates the amount that the optimization of equity balance constraints is allowed to deviate from the required constraints
#' @param calibrate_e logical flag for whether propensity scores should be calibrated after fitting
#' @param calibrate_mu logical flag for whether mean functions should be calibrated after fitting
#' @param calibration_grps tuning parameter for calibration, lower numbers induce more smoothness, with the default being 500
#' @param sub_k number of sub-folds to make for doing local cross-fitting for weight estimation only
#' @param verbose logical flag for whether to print informative message about progress.
#'
#' @export
estimate_equity_wts <- function(data, folds, id, wvars, zvars, a, y, lrnr, lrnr_e = NULL, lrnr_mu = NULL, separate_e = TRUE, separate_mu = TRUE, epsilon = 1e-12, calibrate_e = TRUE, calibrate_mu = TRUE, calibration_grps = 500, sub_k = 2, verbose = TRUE) {
  ##############################
  ## find the folds in the data so that weight estimation can be done on separate folds
  #################################
  all_folds <- pull(data, folds)
  unique_folds <- unique(all_folds)

  #######################
  ## run the workhorse function on each fold
  #####################
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
                                 lrnr_e = lrnr_e,
                                 lrnr_mu = lrnr_mu,
                                 separate_e,
                                 separate_mu,
                                 epsilon,
                                 calibrate_e,
                                 calibrate_mu,
                                 calibration_grps = calibration_grps,
                                 K = sub_k,
                                 verbose = verbose))
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
                    lrnr_e,
                    lrnr_mu,
                    separate_e,
                    separate_mu,
                    epsilon,
                    calibrate_e,
                    calibrate_mu,
                    calibration_grps,
                    K = 2,
                    verbose = TRUE) {

  train_data <- data %>%
    filter(!!sym(folds) != out_fold)
  out_data <- data %>%
    filter(!!sym(folds) == out_fold)

  #######################
  ## create sub-folds on training data
  ##   these subfolds are only used for fitting preliminary propensity scores and mean functions for weight estimation only
  ##   in some sense, these don't need to be as *correct*
  train_data <- make_folds(train_data, unit_id = a, K = K, fold_name = 'subfold')


  #######################
  ## estimate mu (mean functions) and e (propensity scores) on sub-folds
  if (!is.null(lrnr)) {
    if (!is.null(lrnr_mu) | !is.null(lrnr_e)) stop('If `lrnr` is non-null, `lrnr_mu` and `lrnr_e` should be NULL.')
    train_data <- left_join(train_data,
                            estimate_e(train_data, 'subfold', id, c(wvars, zvars), a, lrnr, separate = separate_e, calibrate = calibrate_e, calibration_grps = calibration_grps, verbose = verbose), by = id)
    train_data <- left_join(train_data,
                            estimate_mu(train_data, 'subfold', id, c(wvars, zvars), y, a, lrnr, separate = separate_mu, calibrate = calibrate_mu, calibration_grps = calibration_grps, verbose = verbose), by = id)
  } else {
    train_data <- left_join(train_data,
                            estimate_e(train_data, 'subfold', id, c(wvars, zvars), a, lrnr_e, separate = separate_e, calibrate = calibrate_e, calibration_grps = calibration_grps, verbose = verbose), by = id)
    train_data <- left_join(train_data,
                            estimate_mu(train_data, 'subfold', id, c(wvars, zvars), y, a, lrnr_mu, separate = separate_mu, calibrate = calibrate_mu, calibration_grps = calibration_grps, verbose = verbose), by = id)
  }


  #######################
  ## add influence function to training data
  train_data <- estimate_phi(train_data, a, y)

  ######################
  ## calculate optimal weights based on influence function
  avals <- unique(dplyr::pull(train_data, !!sym(a)))
  train_data <- calculate_wts(train_data, avals = avals, wvars = wvars, epsilon = epsilon)

  ######################
  ## subset data down and output weights on the output (non-training) data
  eqwt_ds <- train_data %>% select(all_of(wvars), equity_wt) %>% unique
  out_data %>%
    inner_join(eqwt_ds, by = wvars)
}
