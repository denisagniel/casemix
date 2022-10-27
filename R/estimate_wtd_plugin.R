
#' Workhorse function for plug-in estimator. Intended to be wrapped by `popwt_plugin` and `eqwt_plugin`.
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
#' @param infl_fn_only logical flag for whether to return a tibble with the influence function for each observation and each unit (which can be used to calculate contrasts between units).
#' @param contrasts tibble or data.frame encoding pre-specified contrasts between units to return
#' @param condition_on a string indicating a variable within which to estimate conditional quality estimates.
#' @param verbose logical flag for whether to print informative messages about progress.
#'
estimate_wtd_plugin <- function(ds,
                        folds,
                        id,
                        a, y, wvars, zvars,
                        wt,
                        truncation_pt = 1e-7,
                        adaptive = FALSE,
                        lrnr = lrn('classif.ranger'),
                        lrnr_e = NULL,
                        lrnr_mu = NULL,
                        separate_e = TRUE,
                        separate_mu = TRUE,
                        calibrate_e = FALSE,
                        calibrate_mu = FALSE,
                        calibration_grps = 500,
                        infl_fn_only = FALSE,
                        contrasts = NULL,
                        condition_on = NULL,
                        verbose = TRUE) {

  ####################
  ## estimate nuisance functions
  if (!is.null(lrnr)) {
    if (!is.null(lrnr_mu) | !is.null(lrnr_e)) stop('If `lrnr` is non-null, `lrnr_mu` and `lrnr_e` should be NULL.')
    ds <- left_join(ds,
                    estimate_e(ds, folds, id, c(wvars, zvars), a, lrnr, separate = separate_e, calibrate = calibrate_e, calibration_grps = calibration_grps, verbose = verbose), by = id)
    ds <- left_join(ds,
                    estimate_mu(ds, folds, id, c(wvars, zvars), y, a, lrnr, separate = separate_mu, calibrate = calibrate_mu, calibration_grps = calibration_grps,verbose = verbose), by = id)
  } else {
    ds <- left_join(ds,
                    estimate_e(ds, folds, id, c(wvars, zvars), a, lrnr_e, separate = separate_e, calibrate = calibrate_e, calibration_grps = calibration_grps, verbose = verbose), by = id)
    ds <- left_join(ds,
                    estimate_mu(ds, folds, id, c(wvars, zvars), y, a, lrnr_mu, separate = separate_mu, calibrate = calibrate_mu, calibration_grps = calibration_grps, verbose = verbose), by = id)
  }

  ds <- estimate_phi(ds, a, y, truncation_pt = truncation_pt)
  # browser()
  avals <- unique(dplyr::pull(ds, !!sym(a)))
  phis <- paste0('phi_', avals)
  if (infl_fn_only) {
    return(ds)
  } else {
    if (is.null(condition_on)) {
      ests <- dplyr::summarise_at(ds,
                                  all_of(phis),
                                  list(pluginest = ~mean(.*!!sym(wt)),
                                       pluginse = ~sqrt(var(.*!!sym(wt))/dplyr::n())))

      out <- tidyr::pivot_longer(ests, cols = everything()) %>%
        tidyr::separate(name, into = c('phi', a, 'type'), sep = '_') %>%
        dplyr::select(-phi) %>%
        tidyr::pivot_wider(id_cols = a, names_from = type, values_from = value) %>%
        dplyr::rename(plugin_est = pluginest,
                      plugin_se = pluginse)

      if (!is.null(contrasts)) {
        browser()
          # wide_phi <- dplyr::select(ds, all_of(c(id, phis)))
          # long_phi <- tidyr::pivot_longer(wide_phi, names_to = 'phi_nm', values_to = 'inf_fn', cols = phis)
          # long_phi <- dplyr::mutate(long_phi, !!sym(a) := stringr::str_remove(phi_nm, 'phi_'))
          # contrast_ds <- dplyr::inner_join(long_phi, contrasts) %>%
          #   group_by(phi_nm) %>%
          #   mutate(phi_var = var(inf_fn)) %>%
          #   ungroup %>%
          #   mutate(phi_wt = phi_var/sum(phi_var))
          #
          # contrast_res <- contrast_ds %>%
          #   select(-id, -phi_nm, -unit) %>%
          #   summarise_at(vars(-phi_sd, -inf_fn), ~
      }
    } else {
      ds <- group_by(ds, !!sym(condition_on))
      ests <- dplyr::summarise_at(ds,
                                  all_of(phis),
                                  list(pluginest = ~mean(.*!!sym(wt)),
                                       pluginse = ~sqrt(var(.*!!sym(wt))/dplyr::n())))
      out <- tidyr::pivot_longer(ests, cols = colnames(ests)[!colnames(ests) == condition_on]) %>%
        tidyr::separate(name, into = c('phi', a, 'type'), sep = '_') %>%
        dplyr::select(-phi) %>%
        tidyr::pivot_wider(id_cols = c(a, condition_on), names_from = type, values_from = value) %>%
        dplyr::rename(plugin_est = pluginest,
                      plugin_se = pluginse)
    }





    ##########################
    ## calculate reliability and shrunken estimators
    shrink_estimates(out, 'plugin_est', 'plugin_se')
  }
}
