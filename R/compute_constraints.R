compute_constraints <- function(data, wvars, epsilon = 0) {
  #####################
  ## figure out how many unique values there are for each protected characteristic
  k_j <- summarise_at(data, vars(all_of(wvars)), ~length(unique(.)))

  #####################
  ## calculate the prob. mass function
  grp_ds <- group_by_at(data, wvars)
  count_ds <- summarise(grp_ds, n = n())
  count_ds <- ungroup(count_ds)
  f_ds <- mutate(count_ds, f = n/sum(n))

  #############################
  ## these constraints enforce marginal balance
  wxfW <- map(wvars, function(w) {
    ####################
    ## find the distinct values of the protected characteristic
    distinct_w <- distinct(f_ds, !!sym(w))
    these_vals <- pull(distinct_w, !!sym(w))

    ####################
    ## give them names
    these_nms <- glue::glue('{w}_{these_vals}')

    #####################
    ## put the binary variables indicating the unique values of the protected characteristic into the pmf dataset
    for (i in 1:nrow(distinct_w)) {
      distinct_w <- mutate(distinct_w, !!these_nms[i] := 1*(!!sym(w) == these_vals[i]))
    }
    wfW <- left_join(f_ds, distinct_w)

    ########################
    ## multiply the binary indicators by the pmf
    mutate_at(wfW, these_nms, ~.*f)
  })
  wxfW_ds <- reduce(wxfW, full_join)
  wxfW_ds <- select(wxfW_ds, all_of(wvars), f, starts_with(wvars))
  ##########################
  ## create a matrix that identifies w_j * f_W(w)
  half_Amat <- select(wxfW_ds, starts_with(wvars))
  half_Amat <- select(half_Amat, -any_of(wvars))

  pp <- nrow(half_Amat)
  ##########################
  ## set up the constraint matrix
  Amat <- as.matrix(cbind(
    wxfW_ds$f, ### this column ensures that the reweighted density remains a density and sums to approx. 1 (\geq 1 - epsilon)
    -wxfW_ds$f, ### this column ensures that the reweighted density remains a density and sums to approx. 1 (< 1 + epsilon)
    half_Amat, ### this column ensures that the \sum \xi_w w_j f_W(w) \geq k_j - epsilon
    -half_Amat, ### this column ensures that the -\sum \xi_w w_j f_W(w) \geq -k_j - epsilon or alternatively that
                ###   \sum \xi_w w_j f_W(w) < k_j + epsilon
    diag(pp) ### these columns ensure that all weights are \geq 0
    ))
  ks <- unlist(map(k_j, function(k) rep(1/k, k)))
  bvec <- c(1 - epsilon, -1 - epsilon, ks - epsilon, -ks - epsilon, rep(0, pp))


  return(list(Amat = Amat, bvec = bvec, pmf_ds = wxfW_ds))
}
