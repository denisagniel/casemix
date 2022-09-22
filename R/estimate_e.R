estimate_e <- function(data, folds, id, x, a, lrnr, task_name = 'e', separate = TRUE, evals = 20, calibrate = TRUE, verbose = FALSE) {
  if (lrnr$predict_type != 'prob') lrnr$predict_type <- 'prob'
  data <- dplyr::mutate(data, row_id = 1:nrow(data))
  data <- dplyr::mutate_if(data, is.character, as.factor)

  if (!separate) {
    if (verbose) print('Fitting single propensity score model for all units.')
    xy_dat <- select(data, any_of(c(a, x)))
    if (!inherits(dplyr::pull(xy_dat, !!a), 'factor')) {
      xy_dat <- dplyr::mutate_at(xy_dat, vars(a), as.factor)
    }
    this_task <- mlr3::as_task_classif(xy_dat, target = a, id = task_name)

    # if (tune) {
    #   if (lrnr$id == 'classif.kknn') {
    #     ts <- lts('classif.kknn.rbv2')
    #     tlrn <- ts$get_learner()
    #   } else if (lrnr$id == 'classif.ranger') {
    #     tlrn <- lrn('classif.ranger',
    #                 num.trees = to_tune(lower = 1, upper = 2000),
    #                 replace = to_tune(),
    #                 sample.fraction = to_tune(lower = 0.1, upper = 1),
    #                 min.node.size = to_tune(lower = 1, upper = 100),
    #                 splitrule = to_tune(),
    #                 num.random.splits = to_tune(lower = 1, upper = 100),
    #                 mtry = to_tune(lower = 1, upper = length(x)),
    #                 max.depth = to_tune(lower = 1, upper = 10))
    #   }
    #   else {
    #     ts <- lts(lrnr)
    #     tlrn <- lrnr
    #   }
    #   if (tlrn$predict_type != 'prob') tlrn$predict_type <- 'prob'
    #   instance = tune(
    #     method = "random_search",
    #     task = this_task,
    #     learner = tlrn,
    #     resampling = rsmp("holdout"),
    #     measure = msr("classif.mbrier"),
    #     term_evals = evals
    #   )
    #   lrnr$param_set$values = instance$result_learner_param_vals
    # }

    all_folds <- dplyr::pull(data, folds)
    unique_folds <- unique(all_folds)
    # predictions <- furrr::future_map_dfr(unique_folds,
    predictions <- purrr::map_df(unique_folds,
                          ~learn_fold_e(task = this_task,
                                        train_ids = which(all_folds != .),
                                        test_ids = which(all_folds == .),
                                        lrnr = lrnr,
                                        a = a,
                                        calibrate = calibrate))
  } else {
    avals <- unique(dplyr::pull(data, !!rlang::sym(a)))
    if (verbose) print('Fitting separate propensity score model for each unit.')
    # progressr::with_progress({
      p <- progressr::progressor(steps = length(avals))
      # pred_js <- furrr::future_map(avals, function(aa) {

        pred_js <- purrr::map(avals, function(aa) {
          p()
        ds_j <- dplyr::mutate(data, a_j = as.factor(ifelse(!!rlang::sym(a) == aa, TRUE, FALSE)))
        xy_dat <- select(ds_j, any_of(c('a_j', x)))
        this_task <- mlr3::as_task_classif(xy_dat, target = 'a_j', id = paste0(task_name, '_', aa))

        all_folds <- dplyr::pull(ds_j, folds)
        unique_folds <- unique(all_folds)
        purrr::map_df(unique_folds,
                      ~learn_fold_e_separate(task = this_task,
                                             train_ids = which(all_folds != .),
                                             test_ids = which(all_folds == .),
                                             lrnr = lrnr,
                                             a = aa,
                                             calibrate = calibrate))
      })
    # })

    predictions <- purrr::reduce(pred_js, inner_join, by = 'row_id')
  }

  # if (calibrate) {
  #   predictions <- callibrate_e(predictions, data, a, folds)
  # }
  out_ds <- inner_join(data, predictions, by = 'row_id')

  select(out_ds, !!id, starts_with('e_'))
}


learn_fold_e <- function(task, train_ids, test_ids, lrnr, a, calibrate, verbose = FALSE) {
  lrnr$train(task, row_ids = train_ids)
  # browser()
  predicted_vals <- lrnr$predict(task, row_ids = test_ids)
  avals <- unique(dplyr::pull(task$data(), !!a))
  progressr::with_progress({

    p <- progressr::progressor(steps = length(avals))
  es <- purrr::map(avals, function(aa) {
    p()
    tibble::tibble(row_id = test_ids,
                           !!glue::glue('e_{aa}') := predicted_vals$prob[,aa])
  })
  })
  train_predicted_vals <- lrnr$predict(task, row_ids = train_ids)
  train_es <- purrr::map(avals, ~tibble::tibble(row_id = train_ids,
                                         !!glue::glue('e_{.}') := train_predicted_vals$prob[,.]))
  if (calibrate) {
    if (verbose) print('Calibrating propensity scores on one of the folds.')
    join_ds <- task$data() %>% transmute(row_id = 1:nrow(task$data()), !!sym(a))
    cal_es <- map(1:length(es), function(i) {
      e_ds <- es[[i]]
      train_e_ds <- train_es[[i]]
      aa <- stringr::str_remove(colnames(e_ds)[2], 'e_')

      join_a <- join_ds %>%
        dplyr::mutate(a_j = 1*(!!sym(a) == aa)) %>%
        select(row_id, a_j)
      cal_pred <- calibrate_prediction(train_e_ds, join_a, e_ds, pred_nm = glue('e_{aa}'), task_y = 'a_j', target_val = 1)
    })
  }
  purrr::reduce(es, dplyr::inner_join, by = 'row_id')
  # reduce(es, inner_join)
}

learn_fold_e_separate <- function(task, train_ids, test_ids, lrnr, aa, calibrate = TRUE, verbose = FALSE) {
  # browser()
  lrnr$train(task, row_ids = train_ids)
  test_predicted_vals <- lrnr$predict(task, row_ids = test_ids)
  test_pred_ds <- tibble::tibble(row_id = test_ids,
                           !!glue::glue('e_{aa}') := test_predicted_vals$prob[,'TRUE'])
  if (calibrate) {
    if (verbose) print('Calibrating propensity scores.')
    train_preds <- lrnr$predict(task, row_ids = test_ids)
    train_pred_ds <- tibble::tibble(row_id = test_ids,
                                    !!glue::glue('e_{aa}') := train_preds$prob[,'TRUE'])
    join_ds <- task$data() %>% transmute(row_id = 1:nrow(task$data()), a_j)

    cal_pred <- calibrate_prediction(train_pred_ds, join_ds, test_pred_ds, pred_nm = glue('e_{aa}'), task_y = task$target_names, target_val = task$class_names[2])
    return(cal_pred)
  } else return(test_pred_ds)


}
