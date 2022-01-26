estimate_equity_wts <- function(data, avals, wvars, epsilon = 1e-6) {
  constraints_list <- compute_constraints(data, wvars, epsilon)
  phis <- paste0('phi_', avals)
  phi_reformat <- reformat_phis(data, phis, wvars, constraints_list$pmf_ds)
  qp_sln <- quadprog::solve.QP(Dmat = diag(phi_reformat$phi_w),
                               Amat = constraints_list$Amat,
                               bvec = constraints_list$bvec,
                               dvec = rep(0, nrow(phi_reformat)))
  xi_ds <- select(phi_reformat, all_of(wvars))
  xi_ds <- mutate(xi_ds, equity_wt = qp_sln$solution)
  left_join(data, xi_ds)
}

reformat_phis <- function(data, phis, wvars, pmf_ds, phi_fn = sum) {
  grp_ds <- group_by_at(data, all_of(wvars))
  summarise_phi_ds <- summarise_at(grp_ds, all_of(phis), ~sum(.^2))
  rowwise_phi_ds <- rowwise(ungroup(summarise_phi_ds))
  out_ds <- mutate(rowwise_phi_ds, phi_w = phi_fn(c_across(phis)))

  ##############################
  ## make sure the structure of this matrix matches up with the structure used to construct the constraints
  out_wvars <- select(out_ds, all_of(wvars))
  pmf_ds_wvars <- select(pmf_ds, all_of(wvars))
  if (!all(out_wvars == pmf_ds_wvars)) {
    # browser()
    stop('the construction of the phis and the constraints do not match.')
  }
  ungroup(out_ds)
}
