#' Shrink and obtain reliabilities for a given set of estimates of unit quality
#'
#' @param data data.frame containing the information to analyze
#' @param ests string identifying the column in `data` that denotes the estimates
#' @param ses string identifying the column in `data` that denotes the standard errors of the estimates
#' @param est_nm name for shrunken estimates in the returned data.frame
#' @param rel_nm name for reliability in the returned data.frame
#'
#'
#' @export

shrink_estimates <- function(data, ests, ses, est_nm = 'shrinkage_est', rel_nm = 'reliability', ...) {
  these_ests <- data %>% pull(!!sym(ests))
  these_ses <- data %>% pull(!!sym(ses))
  # shrink_obj <- FEShR::fe_shrink(as.vector(these_ests), M = map(these_ses^2, as.matrix), centering = 'gen', ...)
  res <- ebnm::ebnm_normal(these_ests, s = these_ses, mode = 'estimate')
  data %>%
    mutate(!!est_nm := res$posterior$mean,
           !!glue('{est_nm}_se') := res$posterior$sd,
           !!rel_nm := res$fitted_g$sd^2/(res$fitted_g$sd^2 + (!!sym(glue('{est_nm}_se')))^2))
}
