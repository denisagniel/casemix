calibrate_prediction <- function(pred_ds, join_ds, test_ds, pred_nm, task_y, target_val, grps = 500) {
  ### break the predictions up into groups of similar predicted values in both the test and training data
  pred_ds <- pred_ds %>%
    mutate(pred_grp = cut(!!sym(pred_nm), breaks = grps))
  if (any(is.na(pred_ds$pred_grp)) | any(pred_ds$pred_grp == Inf) | any(pred_ds$pred_grp == -Inf)) browser()
  test_ds <- test_ds

  ### join up the dataset with ids
  cal_ds <- pred_ds %>%
    inner_join(join_ds, by = 'row_id')

  ## calculate the average predicted value and the average observed value within the groups
  grp_sum <- cal_ds %>%
    group_by(pred_grp) %>%
    summarise(ybar = mean(!!sym(task_y) == target_val),
              predbar = mean(!!sym(pred_nm)))
  if (nrow(grp_sum) <= 20) {
    warning("Too few predicted values to calibrate. Using uncalibrated estimates.")
    ## if there aren't more than 20 predicted probabilities, don't do anything, return the original dataset
    return(test_ds %>%
             transmute(row_id, !!pred_nm := !!sym(pred_nm)))
  } else {
    ## otherwise, fit a gam model regressing the observed outcomes on the predicted values
    mod_fit <- try(mgcv::gam(ybar ~ s(predbar, k = 20), data = grp_sum, family = quasibinomial, control = gam.control(maxit = 1000)), silent = TRUE)
    if ('try-error' %in% class(mod_fit)) {
      warning("Calibration model failed. Using uncalibrated estimates.")
      return(test_ds %>% transmute(row_id, !!pred_nm := !!sym(pred_nm)))
    }
    test_ds <- test_ds %>%
      mutate(yhat = predict(mod_fit, newdata = test_ds, type = 'response'))
  }


  test_ds %>%
    transmute(row_id, !!pred_nm := yhat)
}
