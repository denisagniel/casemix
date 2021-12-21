make_folds <- function(data, unit_id, K) {
  data <- dplyr::group_by(data, !!rlang::sym(unit_id))
  data <- dplyr::mutate(data, fold = sample(1:K, size = dplyr::n(), replace = TRUE))
  data <- dplyr::ungroup(data)
  data
}
