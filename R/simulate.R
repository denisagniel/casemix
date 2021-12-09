#' Simulate individual level data
#' @export
sim_inds <- function(n, pi_x, mu_x, id, qq) {
  tibble(id = id,
         qq = qq,
         x1 = rbinom(n, prob = pi_x, size = 1),
         x2 = rnorm(n, mean = mu_x),
         eta = plogis(x1 - x2 - x1*x2 + qq + 2*pnorm(qq)*x2),
         y = rbinom(n, prob = eta, size = 1))
}
