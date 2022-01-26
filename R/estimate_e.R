estimate_e <- function(data, folds, id, x, a, lrnr, task_name = 'e') {
  if (lrnr$predict_type != 'prob') lrnr$predict_type <- 'prob'
  data <- mutate(data, row_id = 1:nrow(data))

  xy_dat <- select(data, any_of(c(a, x)))
  if (!inherits(pull(xy_dat, !!a), 'factor')) {
    xy_dat <- mutate_at(xy_dat, vars(a), as.factor)
  }
  this_task <- mlr3::as_task_classif(xy_dat, target = a, id = task_name)

  all_folds <- pull(data, folds)
  unique_folds <- unique(all_folds)
  predictions <- map_df(unique_folds,
                        ~learn_fold_e(task = this_task,
                                       train_ids = which(all_folds != .),
                                       test_ids = which(all_folds == .),
                                       lrnr = lrnr,
                                      a = a))
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

