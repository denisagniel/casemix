callibrate_e <- function(predictions, data, a, folds) {
  a_ds <- predictions %>%
    inner_join(data %>% select(row_id, !!folds, !!a)) %>%
    mutate(!!folds := as.factor(!!sym(folds)))
  avals <- a_ds %>% pull(!!a) %>% unique
  new_predix <- map(avals, function(aa) {
    fit <- glm(as.formula(glue('1*({a} == {aa}) ~ e_{aa}*{folds}')), data = a_ds, family = binomial)
    transmute(predictions, row_id = row_id,
              !!glue('e_{aa}') := predict(fit, newdata = a_ds, type = 'response'))
  })
  purrr::reduce(new_predix, inner_join)
}
