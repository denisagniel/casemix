#' Shrink and obtain reliabilities for a given set of estimates of unit quality
#'
#' @param data data.frame containing the information to analyze
#' @param ests string identifying the column in `data` that denotes the estimates
#' @param ses string identifying the column in `data` that denotes the standard errors of the estimates
#' @param est_nm name for shrunken estimates in the returned data.frame
#' @param lam_nm name for lambda parameter in the returned data.frame
#' @param rel_nm name for reliability in the returned data.frame
#'
#'
#' @export

shrink_estimates <- function(data, ests, ses, est_nm = 'shrinkage_est', lam_nm = 'lambda', rel_nm = 'reliability', ...) {
  these_ests <- data %>% pull(!!sym(ests))
  these_ses <- data %>% pull(!!sym(ses))
  shrink_obj <- FEShR::fe_shrink(as.vector(these_ests), M = map(these_ses^2, as.matrix), centering = 'gen', ...)
  data %>%
    mutate(!!est_nm := as.vector(shrink_obj$thetahat),
           !!lam_nm := as.vector(shrink_obj$Lambda_opt),
           !!rel_nm := !!sym(lam_nm)/(!!sym(lam_nm) + (!!sym(ses))^2))
}
