#' Split dataset into folds.
#'
#' @param data data.frame containing the information to analyze
#' @param unit_id string identifying the identifier for units in `data`
#' @param K number of folds
#'
#'
#' @export
make_folds <- function(data, unit_id, K, fold_name = 'fold') {
  if (fold_name %in% colnames(data)) {
    warning(glue::glue('The column {fold_name} was already present in the dataset and is being overwritten.'))
  }
  dplyr::mutate(data, !!fold_name := sample(1:K, size = dplyr::n(), replace = TRUE))
}
