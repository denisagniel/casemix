#' Split dataset into folds.
#'
#' @param data data.frame containing the information to analyze
#' @param unit_id string identifying the identifier for units in `data`
#' @param K number of folds
#'
#'
#' @export
make_folds <- function(data, unit_id, K, fold_name = 'fold') {
  data <- dplyr::group_by(data, !!rlang::sym(unit_id))
  data <- dplyr::mutate(data, !!fold_name := sample(1:K, size = dplyr::n(), replace = TRUE))
  data <- dplyr::ungroup(data)
  data
}
