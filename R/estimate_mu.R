estimate_mu <- function(data, folds, id, x, y, a, lrnr, task_name = 'mu', separate = TRUE, tune = TRUE, evals = 20, callibrate = TRUE) {

  if (lrnr$predict_type != 'prob') lrnr$predict_type <- 'prob'
  data <- mutate(data, row_id = 1:nrow(data))
  data <- mutate_if(data, is.character, as.factor)

  xy_dat <- select(data, any_of(c(a, x, y)))
  if (!inherits(pull(xy_dat, !!y), 'factor')) {
    xy_dat <- mutate_at(xy_dat, vars(y), as.factor)
  }
  this_task <- mlr3::as_task_classif(xy_dat, target = y, id = task_name, twoclass = TRUE)

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
                        ~learn_fold_mu(task = this_task,
                                         train_ids = which(all_folds != .),
                                         test_ids = which(all_folds == .),
                                         lrnr = lrnr,
                                         a = a))
  if (callibrate) {
    predictions <- callibrate_mu(predictions, data, a, y, folds)
  }
  ids <- select(data, !!id, row_id)
  out_ds <- inner_join(ids, predictions, by = 'row_id') %>%
    select(-row_id)

  out_ds
  # select(predictions, !!id, starts_with('mu'))
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
  tibble(row_id = ids, !!glue('mu_{aval}') := predicted_vals$prob[,2])
}
