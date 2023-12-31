% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Simulation Functions.R
\name{accelerated.em}
\alias{accelerated.em}
\title{Extrapolating Conditional Expectations to Accelerate EM Procedures}
\usage{
accelerated.em(
  data,
  par,
  e.step,
  m.step,
  mix.step,
  loglik,
  iters,
  alpha.bound = 100,
  backtrack = 0.5,
  max.backtrack = 0,
  pre.em = 10,
  skip.accel = 1,
  check.convergence = TRUE,
  convergence.eps = 1e-08,
  convergence.iters = 100
)
}
\arguments{
\item{data}{The data used to fit the model, usually but not always a data
matrix or data frame. This object is passed directly to the e.step, m.step,
and loglik functions and is not used by the acceleration procedure itself.}

\item{par}{List of initial parameters and any other parameters or values that
may be needed. This object is passed directly to the e.step, m.step,
and loglik functions and is not used by the acceleration procedure itself.
Specifically, the conditional expectations should be save into par
throughout all steps.}

\item{e.step}{A function that performs the E-step of an EM procedure. Takes
two arguments: par and data. Will be called as part of the accelerated EM
procedure.}

\item{m.step}{A function that performs the (C)M-step of an EM procedure. Also
takes two arguments: par and data.}

\item{mix.step}{A function that takes two pars, one before and one after, and
extrapolates the conditional expectations contained into a new par.}

\item{loglik}{A function that calculates the log-likelihood given two
arguments par and data. This should return a single numeric value, and should
throw an error if parameters are invalid.}

\item{iters}{The maximum number of iterations to attempt of the procedure.}

\item{alpha.bound}{The upper bound for the extrapolation factor. This can be
used to mitigate overly aggressive extrapolations that result in invalid
conditional expectations.}

\item{backtrack}{The backtracking factor between 0 and 1 for when an
extrapolated step fails to produce valid parameters, as determined by the
loglik function throwing an error.}

\item{max.backtrack}{The maximum number of times to attempt backtracking
before falling back to performing a regular EM update step.}

\item{pre.em}{The number of iterations of regular EM to run to warm up
potentially very bad initial parameter values.}

\item{skip.accel}{The number of iterations between which to attempt a
conditionally extrapolated accelerated EM update. If this value is equal to
one, then every iteration will attempt an accelerated EM.}

\item{check.convergence}{Whether or not to check convergence of the
log-likelihood.}

\item{convergence.eps}{Threshold for consecutive log-likelihood values.}

\item{convergence.iters}{Number of successive log-likelihood differences
below `convergence.eps` required to declare convergence.}
}
\description{
Performs extrapolation on the conditional expectations to accelerate the
given EM procedure.
}
\details{
The defaults are generally what is used in the manuscript.
}
