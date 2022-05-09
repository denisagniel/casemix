estimate_mu <- function(data, folds, id, x, y, a, lrnr, task_name = 'mu', separate = TRUE, tune = TRUE, evals = 20, calibrate = TRUE, verbose = TRUE) {
  if (lrnr$predict_type != 'prob') lrnr$predict_type <- 'prob'
  data <- dplyr::mutate(data, row_id = 1:nrow(data)) ## create an internal id with known characteristics
  data <- dplyr::mutate_if(data, is.character, as.factor) ## character vectors aren't typically allowed

  if (!separate) {
    if (verbose) print('Fitting single mean function model for all units.')
    ############################
    ## identify the task for mlr3
    xy_dat <- select(data, any_of(c(a, x, y)))
    if (!inherits(dplyr::pull(xy_dat, !!y), 'factor')) {
      xy_dat <- dplyr::mutate_at(xy_dat, vars(y), as.factor)
    }
    this_task <- mlr3::as_task_classif(xy_dat, target = y, id = task_name, twoclass = TRUE)

    ##########################
    ## this is a first pass at this - we may want to be more intentional about how the tuning is done
    if (tune) {
      if (lrnr$id == 'classif.kknn') {
        ts <- lts('classif.kknn.rbv2')
        tlrn <- ts$get_learner()
      } else if (lrnr$id == 'classif.ranger') {

        tlrn <- lrn('classif.ranger',
                    num.trees = to_tune(lower = 1, upper = 2000),
                    replace = to_tune(),
                    sample.fraction = to_tune(lower = 0.1, upper = 1),
                    min.node.size = to_tune(lower = 1, upper = 100),
                    splitrule = to_tune(),
                    num.random.splits = to_tune(lower = 1, upper = 100),
                    mtry = to_tune(lower = 1, upper = length(x)),
                    max.depth = to_tune(lower = 1, upper = 10))
      }
      else {
        ts <- lts(lrnr)
        tlrn <- lrnr
      }
      if (tlrn$predict_type != 'prob') tlrn$predict_type <- 'prob'
      instance = tune(
        method = "random_search",
        task = this_task,
        learner = tlrn,
        resampling = rsmp("holdout"),
        measure = msr("classif.mbrier"),
        term_evals = evals
      )
      lrnr$param_set$values = instance$result_learner_param_vals
    }

    ####################
    ## estimate on K-1 folds and predict on the other fold
    all_folds <- dplyr::pull(data, folds)
    unique_folds <- unique(all_folds)
    # predictions <- furrr::future_map_dfr(unique_folds,
    predictions <- purrr::map_df(unique_folds,
                          ~learn_fold_mu(task = this_task,
                                           train_ids = which(all_folds != .),
                                           test_ids = which(all_folds == .),
                                           lrnr = lrnr,
                                           a = a,
                                         calibrate = calibrate))
    # if (callibrate) {
    #   predictions <- callibrate_mu(predictions, data, a, y, folds)
    # }
  } else {
    if (verbose) print('Fitting separate mean function model for each unit.')
    avals <- unique(dplyr::pull(data, !!rlang::sym(a)))
    #
    progressr::with_progress({
      p <- progressr::progressor(steps = length(avals))
      pred_js <- furrr::future_map(avals, function(aa) {
        p()
      # ds_j <- filter(data, !!rlang::sym(a) == aa)
      xy_dat <- select(data, any_of(c(y, x)))
      if (!inherits(dplyr::pull(xy_dat, !!y), 'factor')) {
        xy_dat <- dplyr::mutate_at(xy_dat, vars(y), as.factor)
      }
      this_task <- mlr3::as_task_classif(xy_dat, target = y, id = paste0(task_name, '_', aa))

      all_folds <- dplyr::pull(data, folds)
      unique_folds <- unique(all_folds)
      all_as <- dplyr::pull(data, !!rlang::sym(a))
      purrr::map_df(unique_folds,
             ~learn_fold_mu_separate(task = this_task,
                                    train_ids = which(all_folds != . & all_as == aa),
                                    test_ids = which(all_folds == .),
                                    lrnr = lrnr,
                                    a = aa,
                                    calibrate = calibrate))
    })
    })
    predictions <- purrr::reduce(pred_js, dplyr::inner_join, by = 'row_id')
  }
  ids <- select(data, !!id, row_id)
  out_ds <- inner_join(ids, predictions, by = 'row_id') %>%
    select(-row_id)

  out_ds
  # select(predictions, !!id, starts_with('mu'))
}


learn_fold_mu <- function(task, train_ids, test_ids, lrnr, a, avals = NULL, calibrate) {
  lrnr$train(task, row_ids = train_ids)
  # browser()
  if (is.null(avals)) avals <- unique(dplyr::pull(task$data(), !!rlang::sym(a)))
  # mus <- furrr::future_map(avals, ~predict_unit_mu(lrnr, task, task$data(test_ids), test_ids, a, .))
  progressr::with_progress({

    p <- progressr::progressor(steps = length(avals))
    es <- furrr::future_map(avals, function(aa) {
      p()
      predict_unit_mu(lrnr, task, task$data(), test_ids, train_ids, a, aa, calibrate = calibrate)
    })
  })
  purrr::reduce(mus, dplyr::inner_join, by = 'row_id')
}

predict_unit_mu <- function(lrnr, task, data, test_ids, train_ids, a, aval, calibrate) {
  # browser()
  id_data <- data %>%
    dplyr::mutate(row_id = 1:nrow(data))
  a_ds <- dplyr::mutate(data, !!a := aval)
  predicted_vals <- lrnr$predict_newdata(task, newdata = a_ds)
  pred_ds <- tibble::tibble(row_id = 1:nrow(data), !!glue::glue('mu_{aval}') := predicted_vals$prob[,2])
  train_pred_ds <- pred_ds %>%
    dplyr::filter(row_id %in% train_ids)
  test_pred_ds <- pred_ds %>%
    dplyr::filter(row_id %in% test_ids)
  if (calibrate) {
    a_ds <- id_data %>%
      dplyr::filter(!!sym(a) == aval)
    cal_pred <- calibrate_prediction(train_pred_ds, a_ds, test_pred_ds, pred_nm = glue('mu_{aval}'), task_y = task$target_names, target_val = task$class_names[2])
    return(cal_pred)
  } else return(pred_ds)
}

learn_fold_mu_separate <- function(task, train_ids, test_ids, lrnr, aa, calibrate) {
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
    join_ds <- task$data() %>% transmute(row_id = 1:nrow(task$data()), y)
    cal_pred <- calibrate_prediction(train_pred_ds, join_ds, pred_ds, pred_nm = glue('mu_{aa}'), task_y = task$target_names, target_val = task$class_names[2])
    return(cal_pred)
  } else return(pred_ds)
}
