estimate_e <- function(data, folds, id, x, a, lrnr, task_name = 'e', separate = TRUE, tune = TRUE, evals = 20, callibrate = TRUE) {
  if (lrnr$predict_type != 'prob') lrnr$predict_type <- 'prob'
  data <- mutate(data, row_id = 1:nrow(data))
  data <- mutate_if(data, is.character, as.factor)

  if (!separate) {
    xy_dat <- select(data, any_of(c(a, x)))
    if (!inherits(pull(xy_dat, !!a), 'factor')) {
      xy_dat <- mutate_at(xy_dat, vars(a), as.factor)
    }
    this_task <- mlr3::as_task_classif(xy_dat, target = a, id = task_name)

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

    all_folds <- pull(data, folds)
    unique_folds <- unique(all_folds)
    predictions <- map_df(unique_folds,
                          ~learn_fold_e(task = this_task,
                                        train_ids = which(all_folds != .),
                                        test_ids = which(all_folds == .),
                                        lrnr = lrnr,
                                        a = a))
  } else {
    avals <- unique(pull(data, !!sym(a)))
    pred_js <- map(avals, function(aa) {
      ds_j <- mutate(data, a_j = as.factor(ifelse(!!sym(a) == aa, TRUE, FALSE)))
      xy_dat <- select(ds_j, any_of(c('a_j', x)))
      this_task <- mlr3::as_task_classif(xy_dat, target = 'a_j', id = paste0(task_name, '_', aa))

      all_folds <- pull(ds_j, folds)
      unique_folds <- unique(all_folds)
      map_df(unique_folds,
                       ~learn_fold_e_separate(task = this_task,
                                              train_ids = which(all_folds != .),
                                              test_ids = which(all_folds == .),
                                              lrnr = lrnr,
                                              a = aa))
    })
    predictions <- purrr::reduce(pred_js, inner_join)
  }

  if (callibrate) {
    predictions <- callibrate_e(predictions, data, a, folds)
  }
  out_ds <- inner_join(data, predictions, by = 'row_id')

  select(out_ds, !!id, starts_with('e'))
}


learn_fold_e <- function(task, train_ids, test_ids, lrnr, a) {
  lrnr$train(task, row_ids = train_ids)
  predicted_vals <- lrnr$predict(task, row_ids = test_ids)
  avals <- unique(pull(task$data(), !!a))
  es <- map(avals, ~tibble(row_id = test_ids,
                           !!glue('e_{.}') := predicted_vals$prob[,.]))
  reduce(es, inner_join)
}

learn_fold_e_separate <- function(task, train_ids, test_ids, lrnr, aa) {
  lrnr$train(task, row_ids = train_ids)
  predicted_vals <- lrnr$predict(task, row_ids = test_ids)
  tibble(row_id = test_ids,
                           !!glue('e_{aa}') := predicted_vals$prob[,'TRUE'])
}
