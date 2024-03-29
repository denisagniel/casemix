% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_equity_wts.R
\name{estimate_equity_wts}
\alias{estimate_equity_wts}
\title{Estimate weights that obey marginal balance constraints on protected characteristics.}
\usage{
estimate_equity_wts(
  data,
  folds,
  id,
  wvars,
  zvars,
  a,
  y,
  lrnr,
  lrnr_e = NULL,
  lrnr_mu = NULL,
  separate_e = TRUE,
  separate_mu = TRUE,
  epsilon = 1e-12,
  calibrate_e = TRUE,
  calibrate_mu = TRUE,
  calibration_grps = 500,
  sub_k = 2,
  verbose = TRUE
)
}
\arguments{
\item{data}{data.frame containing the information to analyze}

\item{folds}{optional string identifying the column in `data` that denotes the folds for cross-fitting}

\item{id}{string identifying the column in `data` that denotes the ID variable}

\item{wvars}{string vector identifying the columns in `data` that denote the protected characteristics}

\item{zvars}{string vector identifying the columns in `data` that denote the control (non-protected) characteristics}

\item{a}{optional string identifying the column in `data` that denotes the units of interest (i.e., the thing that is being ranked in terms of quality, like health care providers or schools)}

\item{y}{string identifying the column in `data` that denotes the outcome of interest}

\item{lrnr}{mlr3 learner object which will be used to estimate the mean function and propensity score}

\item{lrnr_e}{mlr3 learner object which will be used to estimate the propensity score (ignored if `lrnr` is specified)}

\item{lrnr_mu}{mlr3 learner object which will be used to estimate the mean function (ignored if `lrnr` is specified)}

\item{separate_e}{logical flag for whether propensity scores for each unit should be estimated separately or in a big multinomial model}

\item{separate_mu}{logical flag for whether mean functions for each unit should be estimated separately or in a big joint model}

\item{epsilon}{positive scalar that indicates the amount that the optimization of equity balance constraints is allowed to deviate from the required constraints}

\item{calibrate_e}{logical flag for whether propensity scores should be calibrated after fitting}

\item{calibrate_mu}{logical flag for whether mean functions should be calibrated after fitting}

\item{calibration_grps}{tuning parameter for calibration, lower numbers induce more smoothness, with the default being 500}

\item{sub_k}{number of sub-folds to make for doing local cross-fitting for weight estimation only}

\item{verbose}{logical flag for whether to print informative message about progress.}
}
\description{
These weights ensure that each characteristic represented by a variable in `wvars` is balanced in the weighted population.
}
