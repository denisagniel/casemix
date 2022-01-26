estimate_mu <- function(data, folds, id, x, y, lrnr) {
  if (lrnr$predict_type != 'prob') lrnr$predict_type <- 'prob'
  data <- mutate(data, row_id = 1:nrow(data))

  xy_dat <- select(data, any_of(c(x, y)))
  if (!inherits(pull(xy_dat, !!y), 'factor')) {
    xy_dat <- mutate_at(xy_dat, vars(y), as.factor)
  }
  this_task <- mlr3::as_task_classif(xy_dat, target = y, id = task_name)

  all_folds <- pull(data, folds)
  unique_folds <- unique(all_folds)
  incl <- pull(data, include_train)
  predictions <- map_df(unique_folds,
                        ~learn_fold_prob(task = this_task,
                                         train_ids = which(all_folds != . & incl),
                                         test_ids = which(all_folds == .),
                                         lrnr = lrnr))
  # out_ds <- inner_join(data, predictions, by = 'row_id')

  select(predictions, !!id, starts_with('mu'))
}


learn_fold_prob <- function(task, train_ids, test_ids, lrnr, a, avals = NULL) {
  lrnr$train(task, row_ids = train_ids)
  if (is.null(avals)) avals <- task$data %>% pull(!!sym(a)) %>% unique
  mus <- map(avals, ~predict_unit(lrnr, task, task$data(test_ids), a, .))
  browser()
  reduce(mus, full_join)
}

predict_unit <- function(lrnr, task, data, a, aval) {
  new_ds <- data %>% mutate(!!a := aval)
  predicted_vals <- lrnr$predict_newdata(task, newdata = new_ds)
  mutate(data, !!glue('mu_{aval}') := predicted_vals[,2])
}
