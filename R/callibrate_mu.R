callibrate_mu <- function(predictions, data, a, y, folds) {
  a_ds <- predictions %>%
    inner_join(data %>% select(row_id, !!folds, !!a, !!y)) %>%
    mutate(!!folds := as.factor(!!sym(folds)))
  avals <- a_ds %>% pull(!!a) %>% unique
  new_predix <- map(avals, function(aa) {
    fit <- glm(as.formula(glue('{y} ~ mu_{aa}*{folds}')), data = a_ds %>% filter(!!sym(a) == aa), family = binomial)
    transmute(predictions, row_id = row_id,
              !!glue('mu_{aa}') := predict(fit, newdata = a_ds, type = 'response'))
  })
  out <- purrr::reduce(new_predix, inner_join)
  out
}

