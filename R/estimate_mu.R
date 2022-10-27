estimate_mu <- function(data, folds, id, x, y, a, lrnr, task_name = 'mu', separate = TRUE, calibrate = TRUE, calibration_grps = 500, verbose = FALSE) {
  if (lrnr$predict_type != 'prob') lrnr$predict_type <- 'prob'
  data <- dplyr::mutate(data, row_id = 1:nrow(data)) ## create an internal id with known characteristics
  data <- dplyr::mutate_if(data, is.character, as.factor) ## character vectors aren't typically allowed

  if (!separate) {
    if (verbose) print('Fitting single mean function model for all units. Progress bar advances as prediction occurs for each unit on a given fold - one progress bar per fold.')
    ############################
    ## identify the task for mlr3
    xy_dat <- select(data, any_of(c(a, x, y)))
    if (!inherits(dplyr::pull(xy_dat, !!y), 'factor')) {
      xy_dat <- dplyr::mutate_at(xy_dat, vars(y), as.factor)
    }
    this_task <- mlr3::as_task_classif(xy_dat, target = y, id = task_name, twoclass = TRUE)

    ####################
    ## estimate on K-1 folds and predict on the other fold
    all_folds <- dplyr::pull(data, folds)
    unique_folds <- unique(all_folds)
    predictions <- purrr::map_df(unique_folds,
                          ~learn_fold_mu(task = this_task,
                                           train_ids = which(all_folds != .),
                                           test_ids = which(all_folds == .),
                                           lrnr = lrnr,
                                           a = a,
                                         calibrate = calibrate,
                                         calibration_grps = calibration_grps,
                                         verbose = verbose))
  } else {
    avals <- unique(dplyr::pull(data, !!rlang::sym(a)))
    if (verbose) {
      print('Fitting separate mean function model for each unit. Progress bar advances as prediction occurs for each unit on a given fold - one progress bar per fold.')
      progressr::with_progress({
      p <- progressr::progressor(steps = length(avals))
      pred_js <- purrr::map(avals, function(aa) {
        p(message = glue::glue('Estimating outcome function for unit {aa}'))

        ############################
        ## identify the task for mlr3
      xy_dat <- select(data, any_of(c(y, x)))
      if (!inherits(dplyr::pull(xy_dat, !!y), 'factor')) {
        xy_dat <- dplyr::mutate_at(xy_dat, vars(y), as.factor)
      }
      this_task <- mlr3::as_task_classif(xy_dat, target = y, id = paste0(task_name, '_', aa))

      ####################
      ## estimate on K-1 folds and predict on the other fold
      all_folds <- dplyr::pull(data, folds)
      unique_folds <- unique(all_folds)
      all_as <- dplyr::pull(data, !!rlang::sym(a))
      purrr::map_df(unique_folds,
             ~learn_fold_mu_separate(task = this_task,
                                    train_ids = which(all_folds != . & all_as == aa),
                                    test_ids = which(all_folds == .),
                                    lrnr = lrnr,
                                    a = aa,
                                    calibrate = calibrate,
                                    calibration_grps = calibration_grps))
    })
    })
    } else {
      pred_js <- purrr::map(avals, function(aa) {

        ############################
        ## identify the task for mlr3
        xy_dat <- select(data, any_of(c(y, x)))
        if (!inherits(dplyr::pull(xy_dat, !!y), 'factor')) {
          xy_dat <- dplyr::mutate_at(xy_dat, vars(y), as.factor)
        }
        this_task <- mlr3::as_task_classif(xy_dat, target = y, id = paste0(task_name, '_', aa))

        ####################
        ## estimate on K-1 folds and predict on the other fold
        all_folds <- dplyr::pull(data, folds)
        unique_folds <- unique(all_folds)
        all_as <- dplyr::pull(data, !!rlang::sym(a))
        purrr::map_df(unique_folds,
                      ~learn_fold_mu_separate(task = this_task,
                                              train_ids = which(all_folds != . & all_as == aa),
                                              test_ids = which(all_folds == .),
                                              lrnr = lrnr,
                                              a = aa,
                                              calibrate = calibrate,
                                              calibration_grps = calibration_grps))
      })
    }
    ############################
    ## put everything into a single dataset
    predictions <- purrr::reduce(pred_js, dplyr::inner_join, by = 'row_id')
  }

  ####################
  ## merge back with original ids
  ids <- select(data, !!id, row_id)
  out_ds <- inner_join(ids, predictions, by = 'row_id') %>%
    select(-row_id)

  out_ds
}


learn_fold_mu <- function(task, train_ids, test_ids, lrnr, a, avals = NULL, calibrate, calibration_grps = 500, verbose = FALSE) {
  ### train the model on training data and obtain predicted values on test data
  lrnr$train(task, row_ids = train_ids)
  if (is.null(avals)) avals <- unique(dplyr::pull(task$data(), !!rlang::sym(a)))
  if (verbose) {
    progressr::with_progress({
      p <- progressr::progressor(steps = length(avals))
      mus <- purrr::map(avals, function(aa) {
        p(message = glue::glue('Estimating outcome function for unit {aa}'))
        predict_unit_mu(lrnr, task, task$data(), test_ids, train_ids, a, aa, calibrate = calibrate, calibration_grps = calibration_grps)
      })
    })
  } else {
      mus <- purrr::map(avals, function(aa) {
        predict_unit_mu(lrnr, task, task$data(), test_ids, train_ids, a, aa, calibrate = calibrate, calibration_grps = calibration_grps)
      })
  }

  purrr::reduce(mus, dplyr::inner_join, by = 'row_id')
}

predict_unit_mu <- function(lrnr, task, data, test_ids, train_ids, a, aval, calibrate, calibration_grps) {
  ## collect ids
  id_data <- data %>%
    dplyr::mutate(row_id = 1:nrow(data))
  ##############################
  ## set the unit identifier to `aval` because we want to obtain the predicted value for each observation in the world where they belong to unit `aval`
  a_ds <- dplyr::mutate(data, !!a := aval)
  ## collect predicted values for all observations and put in a dataset
  predicted_vals <- lrnr$predict_newdata(task, newdata = a_ds)
  pred_ds <- tibble::tibble(row_id = 1:nrow(data), !!glue::glue('mu_{aval}') := predicted_vals$prob[,2])

  ######################
  ## get predicted vals for test/held out data
  test_pred_ds <- pred_ds %>%
    dplyr::filter(row_id %in% test_ids)
  if (calibrate) {
    #####################
    ## if calibrating, get predicted values from training data
    train_pred_ds <- pred_ds %>%
      dplyr::filter(row_id %in% train_ids)
    #####################
    ## only calibrate on data that truly belongs to unit `aval`
    a_ds <- id_data %>%
      dplyr::filter(!!sym(a) == aval)
    ## calibrate:
    cal_pred <- calibrate_prediction(train_pred_ds, a_ds, test_pred_ds, pred_nm = glue('mu_{aval}'), task_y = task$target_names, target_val = task$class_names[2], grps = calibration_grps)
    return(cal_pred)
  } else return(test_pred_ds) ## if no calibration, just pass back predicted values from the test data
}

learn_fold_mu_separate <- function(task, train_ids, test_ids, lrnr, aa, calibrate, calibration_grps) {
  #####################
  #### train the model on training data and obtain predicted values on test data
  lrnr$train(task, row_ids = train_ids)
  # browser()
  predicted_vals <- lrnr$predict(task, row_ids = test_ids)
  pred_ds <- tibble::tibble(row_id = test_ids,
         !!glue::glue('mu_{aa}') := predicted_vals$prob[,2])
  if (calibrate) {
    # browser()
    train_predicted_vals <- lrnr$predict(task, row_ids = train_ids)
    train_pred_ds <- tibble::tibble(row_id = train_ids,
                              !!glue::glue('mu_{aa}') := train_predicted_vals$prob[,2])
    join_ds <- task$data() %>% transmute(row_id = 1:nrow(task$data()), !!rlang::sym(task$target_names))
    cal_pred <- calibrate_prediction(train_pred_ds, join_ds, pred_ds, pred_nm = glue('mu_{aa}'), task_y = task$target_names, target_val = task$class_names[2], grps = calibration_grps)
    return(cal_pred)
  } else return(pred_ds)
}
