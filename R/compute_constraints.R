compute_constraints <- function(data, wvars, epsilon = 0) {
  k_j <- summarise_at(data, vars(all_of(wvars)), ~length(unique(.)))
  grp_ds <- group_by_at(ds, wvars)
  count_ds <- summarise(grp_ds, n = n())
  count_ds <- ungroup(count_ds)
  f_ds <- mutate(count_ds, f = n/sum(n))

  wxfW <- map(wvars, function(w) {
    distinct_w <- distinct(f_ds, !!sym(w))
    these_vals <- pull(distinct_w, !!sym(w))
    these_nms <- glue::glue('{w}_{these_vals}')
    for (i in 1:nrow(distinct_w)) {
      distinct_w <- mutate(distinct_w, !!these_nms[i] := 1*(!!sym(w) == these_vals[i]))
    }
    wfW <- left_join(f_ds, distinct_w)
    mutate_at(wfW, these_nms, ~.*f)
  })
  wxfW_ds <- reduce(wxfW, full_join)
  wxfW_ds <- select(wxfW_ds, all_of(wvars), f, starts_with(wvars))
  half_Amat <- select(wxfW_ds, starts_with(wvars))
  half_Amat <- select(half_Amat, -any_of(wvars))
  pp <- nrow(half_Amat)
  Amat <- as.matrix(cbind(half_Amat, -half_Amat, rep(1, pp), -rep(1, pp), diag(pp)))
  ks <- unlist(map(k_j, function(k) rep(1/k, k)))
  bvec <- c(ks - epsilon, -ks - epsilon, -1, 0, rep(0, pp))
  return(list(Amat = Amat, bvec = bvec))
}
