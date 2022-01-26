estimate_mu <- function(data, folds, id, x, y, a, lrnr, task_name = 'mu') {
  if (lrnr$predict_type != 'prob') lrnr$predict_type <- 'prob'
  data <- mutate(data, row_id = 1:nrow(data))

  xy_dat <- select(data, any_of(c(a, x, y)))
  if (!inherits(pull(xy_dat, !!y), 'factor')) {
    xy_dat <- mutate_at(xy_dat, vars(y), as.factor)
  }
  this_task <- mlr3::as_task_classif(xy_dat, target = y, id = task_name)

  all_folds <- pull(data, folds)
  unique_folds <- unique(all_folds)
  predictions <- map_df(unique_folds,
                        ~learn_fold_mu(task = this_task,
                                         train_ids = which(all_folds != .),
                                         test_ids = which(all_folds == .),
                                         lrnr = lrnr,
                                         a = a))
  # out_ds <- inner_join(data, predictions, by = 'row_id')

  select(predictions, !!id, starts_with('mu'))
}


learn_fold_mu <- function(task, train_ids, test_ids, lrnr, a, avals = NULL) {
  lrnr$train(task, row_ids = train_ids)
  if (is.null(avals)) avals <- task$data() %>% pull(!!sym(a)) %>% unique
  mus <- map(avals, ~predict_unit_mu(lrnr, task, task$data(test_ids), test_ids, a, .))
  reduce(mus, inner_join)
}

predict_unit_mu <- function(lrnr, task, data, ids, a, aval) {
  new_ds <- data %>% mutate(!!a := aval)
  predicted_vals <- lrnr$predict_newdata(task, newdata = new_ds)
  tibble(id = ids, !!glue('mu_{aval}') := predicted_vals$prob[,2])
}
