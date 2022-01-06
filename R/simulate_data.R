#' Simulate quality data
#'
#' @param p number of units (e.g., providers)
#' @param n vector of length `p` specifying number of per-unit measurements (e.g., number of patients per provider)
#' @param pi_x vector of length `p` specifying the unit-specific probability of binary covariate `X1`
#' @param mu_x vector of length `p` specifying the unit-specific mean of a normal (sd = 1) covariate `X2`
#' @param qq vector of length `p` specifying the true unit-specific quality
#'
#' @export

simulate_data <- function(p, n, pi_x, mu_x, qq){

  # specify unit-level parameters
  unit_ds <- dplyr::tibble(
    id = 1:p,
    n,
    pi_x,
    mu_x,
    qq
  )

  # compute patient-level outcomes for each unit
  sim_inds <- function(n, pi_x, mu_x, id, qq) {
    tibble(id = id,
           qq = qq,
           x1 = rbinom(n, prob = pi_x, size = 1),
           x2 = rnorm(n, mean = mu_x),
           eta = plogis(x1 - x2 - x1 * x2 + qq + 2*pnorm(qq)*x2),
           y = rbinom(n, prob = eta, size = 1))
  }

  ind_ds <- purrr::pmap_df(unit_ds, sim_inds)

  # compute true theta's
  theta_0 <- map_df(1:p, function(i) {
    this_q <- ind_ds %>%
      filter(id == i) %>%
      pull(qq) %>%
      unique
    ind_ds %>%
      summarise(id = i,
                theta_0 = mean(plogis(x1 - x2 - x1*x2 + this_q + 2*pnorm(this_q)*x2)))
  })

  vals <- list(data = ind_ds, theta = theta_0)
  return(vals)
}
