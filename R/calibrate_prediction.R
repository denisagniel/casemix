calibrate_prediction <- function(pred_ds, join_ds, test_ds, pred_nm, task_y, target_val, grps = 500) {
  pred_ds <- pred_ds %>%
    mutate(pred_grp = Hmisc::cut2(!!sym(pred_nm), g = grps),
           predbar = !!sym(pred_nm))
  test_ds <- test_ds %>%
    mutate(predbar = !!sym(pred_nm))
  cal_ds <- pred_ds %>%
    inner_join(join_ds, by = 'row_id')
  grp_sum <- cal_ds %>%
    group_by(pred_grp) %>%
    summarise(ybar = mean(!!sym(task_y) == target_val),
              predbar = mean(!!sym(pred_nm)))
  if (nrow(grp_sum) <= 10) {
    test_ds <- test_ds %>%
      inner_join(grp_sum) %>%
      rename(yhat = ybar)
  } else {
    mod_fit <- try(mgcv::gam(ybar ~ s(predbar, k = 10), data = grp_sum, family = binomial), silent = TRUE)
    if (class(mod_fit) == 'try-error') return(test_ds %>% transmute(row_id, !!pred_nm := !!sym(pred_nm)))

    test_ds <- test_ds %>%
      mutate(yhat = predict(mod_fit, newdata = test_ds, type = 'response'))
  }


  test_ds %>%
    transmute(row_id, !!pred_nm := yhat)
}
