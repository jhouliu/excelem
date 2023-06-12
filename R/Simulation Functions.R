# library(Rfast);library(ggplot2);library(dplyr);library(SQUAREM)

# Helper Functions ----

#' Multivariate Normal Density function
#'
#' @description
#' Calculates the (log) multivariate normal density function for a matrix of
#' observations. Based on the Rfast version except it does a Cholesky
#' decomposition on the covariance sigma, and re-uses it to compute the
#' log-determinant instead of starting from scratch. Will throw an error if the
#' Cholesky decomposition fails.
dmvnorm = function (x, mu, sigma, logged = FALSE) {
  chol             = cholesky(sigma)
  if (any(is.nan(diag(chol, names = FALSE)))) stop("Bad Sigma")
  quat             = -0.5 * mahala(x, mu, chol, ischol = TRUE)
  pow              = length(mu)/2
  logcon           = pow * 1.83787706640935 + sum(log(diag(chol, names = FALSE)))
  den              = quat - logcon
  if (!logged) den = exp(den)
  den
}

aitken.acceleration = function(x) {(x[3] - x[2]) / (x[2] - x[1])}
aitken.factor = function(x) (x[2] - x[1]) / (2 * x[2] - x[1] - x[3])
tr = function(X) if (length(X)==1) X else sum(diag(X, names = FALSE))
logdet = function(X) determinant.matrix(X)$modulus

#' Wrapper for SQUAREM package function
#'
#' @param data As in accelerated.em().
#' @param par As in accelerated.em().
#' @param e.step As in accelerated.em().
#' @param m.step As in accelerated.em().
#' @param loglik As in accelerated.em().
#' @param pre.em As in accelerated.em().
#' @param wrapper A function that converts the parameters, often of a list form,
#' to the vectorized representation expected by SQUAREM. This function should
#' try to preserve alternate metadata by assigning the original par to an
#' attribute of the parameter vector.
#' @param unwrapper The inverse of the wrapper function.
#' @param ... Further arguments to be passed to SQUAREM control.
#'
#' @details
#' The wrapper and unwrapper functions permit re-use of the e.step, m.step, and
#' loglik functions used in the extrapolated conditional expectation
#' acceleration procedure. If dedicated e.step, m.step, and loglik functions
#' conforming to the forms expected by SQUAREM are available, set the wrapper
#' and unwrapper functions to be the identity.
#'
squarem.wrapper = function(data, par, e.step, m.step, loglik, wrapper, unwrapper, pre.em = 100, ...) {
  info = list()
  par$memory.env = environment()
  fit = squarem(
    wrapper(par),
    data,
    fixptfn = function(par, data) {wrapper(m.step(e.step(unwrapper(par), data), data))},
    objfn = function(par, data) {
      par = unwrapper(par)
      ll = loglik(par, data)
      assign("info", list(get("info", par$memory.env), c(ll, unname(proc.time()["elapsed"]))), par$memory.env)
      return(ll)
    },
    control = list(minimize = FALSE, ...)
  )
  info = matrix(unlist(info), ncol = 2, byrow = TRUE)
  return(list(
    par = unwrapper(fit$par),
    loglik.trace = info[,1],
    alpha.trace = rep(0, nrow(info)),
    timing = c(0, diff(info[,2])),
    update = rep("squarem", nrow(info)),
    converged = fit$convergence
  ))
}

# Main Framework ----

#' Extrapolating Conditional Expectations to Accelerate EM Procedures
#'
#' @description
#' Performs extrapolation on the conditional expectations to accelerate the
#' given EM procedure.
#'
#' @param data The data used to fit the model, usually but not always a data
#' matrix or data frame. This object is passed directly to the e.step, m.step,
#' and loglik functions and is not used by the acceleration procedure itself.
#' @param par List of initial parameters and any other parameters or values that
#' may be needed. This object is passed directly to the e.step, m.step,
#' and loglik functions and is not used by the acceleration procedure itself.
#' Specifically, the conditional expectations should be save into par
#' throughout all steps.
#' @param e.step A function that performs the E-step of an EM procedure. Takes
#' two arguments: par and data. Will be called as part of the accelerated EM
#' procedure.
#' @param m.step A function that performs the (C)M-step of an EM procedure. Also
#' takes two arguments: par and data.
#' @param mix.step A function that takes two pars, one before and one after, and
#' extrapolates the conditional expectations contained into a new par.
#' @param loglik A function that calculates the log-likelihood given two
#' arguments par and data. This should return a single numeric value, and should
#' throw an error if parameters are invalid.
#' @param iters The maximum number of iterations to attempt of the procedure.
#' @param alpha.bound The upper bound for the extrapolation factor. This can be
#' used to mitigate overly aggressive extrapolations that result in invalid
#' conditional expectations.
#' @param backtrack The backtracking factor between 0 and 1 for when an
#' extrapolated step fails to produce valid parameters, as determined by the
#' loglik function throwing an error.
#' @param max.backtrack The maximum number of times to attempt backtracking
#' before falling back to performing a regular EM update step.
#' @param pre.em The number of iterations of regular EM to run to warm up
#' potentially very bad initial parameter values.
#' @param skip.accel The number of iterations between which to attempt a
#' conditionally extrapolated accelerated EM update. If this value is equal to
#' one, then every iteration will attempt an accelerated EM.
#' @param check.convergence Whether or not to check convergence of the
#' log-likelihood.
#' @param convergence.eps Threshold for consecutive log-likelihood values.
#' @param convergence.iters Number of successive log-likelihood differences
#' below `convergence.eps` required to declare convergence.
#'
#' @details
#' The defaults are generally what is used in the manuscript.
#'
accelerated.em = function(
    data, par, e.step, m.step, mix.step, loglik, iters,
    alpha.bound = 100, backtrack = 0.5, max.backtrack = 0,
    pre.em = 10, skip.accel = 1,
    check.convergence = TRUE, convergence.eps = 1e-8, convergence.iters = 100
) {
  # Ensure configuration parameters are valid
  stopifnot(iters >= 0)
  stopifnot(backtrack >= 0)
  stopifnot(backtrack < 1)
  stopifnot(alpha.bound >= 0)
  stopifnot(pre.em >= 3)
  stopifnot(max.backtrack >= 0)

  # Internal helper
  backtrack.scale = c(backtrack^(0:max.backtrack), 0)
  step = function(alpha, min.loglik) {
    for (a in alpha * backtrack.scale) {
      ll = tryCatch({
        new.par = m.step(mix.step(par, old.par, a), data)
        loglik(new.par, data)
      }, error = function(e) NaN)
      if (isTRUE(is.finite(ll) && ll > min.loglik)) break
    }
    return(list(par = new.par, ll = ll, alpha = a))
  }

  timing = alpha.trace = loglik.trace = numeric(iters) * NA
  update = rep(NA_character_, iters)
  converged = FALSE

  # Main loop
  last.aitken = 1
  start.time = proc.time()["elapsed"]
  for (iter in 1:iters) {
    if (iter >= pre.em && iter %% skip.accel == 0) {
      # Exceeds pre.em and every skip.accel iter, do aitken acceleration
      aitken.value = min(aitken.factor(loglik.trace[iter - 3:1]), alpha.bound)
      if (is.finite(aitken.value) && aitken.value > 0) { # Valid alpha value
        stepped = step(aitken.value, loglik.trace[iter - 1])
        update[iter] = if (stepped$alpha > 0) "aitken" else "aitken_fail"
      } else { # Invalid alpha value, use last alpha
        stepped = step(last.aitken, loglik.trace[iter - 1])
        update[iter] = if (stepped$alpha > 0) "old" else "old_fail"
      }
      new.par = stepped$par
      ll = stepped$ll
      alpha = stepped$alpha
    } else {
      # During other iterations, only do em
      new.par = m.step(par, data)
      ll = loglik(new.par, data)
      alpha = 0
      update[iter] = if (iter <= pre.em) "pre-em" else "em"
    }

    old.par = par
    par = e.step(new.par, data)
    loglik.trace[iter] = ll

    if (iter > 1 && !isTRUE(loglik.trace[iter] > loglik.trace[iter - 1])) {
      # Update went poorly, undo update and do regular EM
      update[iter] = "very_bad"
      loglik.trace[iter] = loglik(par <- e.step(m.step(old.par, data), data), data)
    }
    timing[iter] = unname(proc.time()["elapsed"])
    alpha.trace[iter] = alpha
    loglik.conv.check = tail(loglik.trace[1:iter], convergence.iters + 1)
    if (check.convergence && (iter > convergence.iters + 1) && all(diff(loglik.conv.check) <= convergence.eps) &&
      isTRUE(all(loglik.conv.check == loglik.conv.check[1]) || cor(head(loglik.conv.check, -1), tail(loglik.conv.check, -1)) <= 0)) {
      loglik.trace = as.numeric(na.omit(loglik.trace))
      alpha.trace = as.numeric(na.omit(alpha.trace))
      timing = as.numeric(na.omit(timing))
      update = as.character(na.omit(update))
      converged = TRUE
      break
    }
  }

  return(list(
    par = par,
    loglik.trace = loglik.trace,
    alpha.trace = alpha.trace,
    timing = diff(c(start.time, timing)),
    update = update,
    converged = converged
  ))
}

# Factor Analyzer ----

loglik.FA = function(par, data) {
  Sigma = par$BBt + diag(par$D, par$p)
  chol = cholesky(Sigma)
  if (any(is.nan(diag(chol, names = FALSE)))) return(NaN)
  invSigma = chol2inv(chol)
  return(-sum(log(diag(chol))) - 0.5 * sum(par$Cww/par$n * invSigma))
}

e.step.FA = function(par, data) {
  invDB = par$B / par$D
  BtinvDB = crossprod(par$B, invDB)
  par$gamma = invDB - (invDB %*% spdinv(diag(1, par$q, names = FALSE) + BtinvDB) %*% BtinvDB)
  par$Omega = diag(1, par$q, names = FALSE) - crossprod(par$gamma, par$B)

  par$E.Cwa = par$Cww %*% par$gamma
  par$E.Caa = crossprod(par$gamma, par$E.Cwa) + par$n * par$Omega
  return(par)
}

m.step.FA = function(par, data) {
  chol.E.Caa = cholesky(par$E.Caa)
  if (any(is.nan(chol.E.Caa))) return(par)
  par$B = par$E.Cwa %*% chol2inv(chol.E.Caa)
  par$BBt = tcrossprod(par$B)
  par$D = (diag(par$Cww, names = FALSE) - mahala(par$E.Cwa, rep(0, par$q), chol.E.Caa, ischol = TRUE)) / par$n
  # if (any(par$D <= 0)) warning("Heywood case")
  return(par)
}

cm.step.FA = function(par, data) {
  chol.E.Caa = cholesky(par$E.Caa)
  if (any(is.nan(diag(chol.E.Caa, names = FALSE)))) return(par)
  par$B = par$E.Cwa %*% chol2inv(chol.E.Caa)
  par$BBt = tcrossprod(par$B)

  obj = function(D) {
    obj.par = par
    obj.par$D = exp(D)
    return(-loglik.FA(obj.par, data))
  }
  par$D = exp(nlm(obj, log(par$D))$estimate)
  return(par)
}

mix.step.FA = function(par, old.par, alpha) {
  new.par = par
  new.par$E.Cwa = (alpha + 1) * par$E.Cwa - alpha * old.par$E.Cwa
  new.par$E.Caa = (alpha + 1) * par$E.Caa - alpha * old.par$E.Caa
  return(new.par)
}

init.par.FA = function(data, q, method = "PCA") {
  p = ncol(data)
  n = nrow(data)
  par = list()
  par$p = p
  par$q = q
  par$n = n

  if (method == "FA") {
    FA = factanal(data, q)
    par$B = FA$loadings
    par$D = FA$uniquenesses
  } else if (method == "PCA") {
    PCA = prcomp(data)
    par$B = PCA$rotation[, 1:q]
    par$D = PCA$sdev^2
  } else {
    par$B = matrix(0, p, q)
    par$D = rep(1, p)
  }

  par$BBt = tcrossprod(par$B)
  par$mu = colMeans(data)
  centred.w = t(t(data) - par$mu)
  Cww = crossprod(centred.w)
  par$Cww = Cww
  return(e.step.FA(par, data))
}

gen.data.FA = function(n, p, q, sd = 1, seed = sample.int(2147483647L, 1)) {
  set.seed(seed)
  latent.factors = matrix(runif(n * q, -1, 1), n, q)
  loading.matrix = mvtnorm::rmvnorm(p, rep(0, q), diag(1, q))
  sigma = rep(sd, length.out = p)
  noise = mvtnorm::rmvnorm(n, rep(0, p), diag(sigma))
  observed.data = tcrossprod(latent.factors, loading.matrix) + noise
  observed.data = eachrow(observed.data, colMeans(observed.data), "-")
  return(observed.data)
}

# GMM ----

loglik.GMM = function(par, data) {
  if (any(par$pi <= 0)) return(NaN)
  tryCatch({
    logwts = vapply(1:par$G, function(g) {
      dmvnorm(data, par$mu[[g]], par$sigma[[g]], log = TRUE) + log(par$pi[g])
    }, numeric(par$n))
    return(sum(log(rowsums(exp(logwts)))))
  }, error = function(e) NaN)
}

e.step.GMM = function(par, data) {
  par$logwts = vapply(1:par$G, function(g) {
    dmvnorm(data, par$mu[[g]], par$sigma[[g]], log = TRUE) + log(par$pi[[g]])
  }, numeric(par$n))
  wts = exp(par$logwts - rowMaxs(par$logwts, TRUE)); wts = wts / rowsums(wts)
  par$zhat = wts
  return(par)
}

m.step.GMM = function(par, data) {
  z.total = colsums(par$zhat)
  if (any(par$z.total <= 0)) return(par)
  par$pi = z.total / par$n
  w = t(par$zhat) / z.total
  par$mu = asplit(w %*% data, 1)
  par$sigma = lapply(1:par$G, function(g) {
    centred.data = eachrow(data, par$mu[[g]], "-")
    crossprod(centred.data * w[g,], centred.data)
  })
  return(par)
}

mix.step.GMM = function(par, old.par, alpha) {
  new.par = par
  new.par$zhat = (alpha + 1) * par$zhat - alpha * old.par$zhat
  return(new.par)
}

init.par.GMM = function(data, G) {
  km = kmeans(data, G)
  par = list()
  par$n = nrow(data)
  par$d = ncol(data)
  par$G = G
  par$zhat = matrix(0, par$n, par$G)
  par$zhat[cbind(1:par$n, km$cluster)] = 1
  par = m.step.GMM(par, data)
  return(par)
}

gen.data.GMM = function(n, d, sigma = 1, seed = sample.int(2147483647L, 1)) {
  set.seed(seed)
  nn = as.numeric(rmultinom(1, n, rep(1, 2^d)))
  mus = expand.grid(rep(list(c(-1,1)), d))
  data = do.call(rbind, lapply(1:nrow(mus), function(g) {
    mvtnorm::rmvnorm(nn[g], as.numeric(mus[g,]), diag(sigma, d))
  }))
  return(data)
}

# Variance Components 3rd re-write (assume R = identity) ----

loglik.VC3 = function(par, data) {
  ll = vapply(1:par$groups, function(j) {
    inv.sigmaj = diag(1/par$sigmasq, par$nj[j]) - par$shortcut1[[j]] / par$sigmasq
    sum(dmvnorm.inv.sigma(data$y[[j]], par$Xbeta[[j]], inv.sigmaj, logged = TRUE))
  }, 1)
  return(sum(ll))
}

e.step.VC3 = function(par, data) {
  invD = solve(par$D)
  par$E = lapply(1:par$groups, function(j) {
    E.b = solve(par$ZtZ[[j]] + par$sigmasq * invD,
                crossprod(data$Z[[j]], data$y[[j]] - par$Xbeta[[j]]))
    V.b = solve(par$ZtZ[[j]] / par$sigmasq + invD)
    E.bb = V.b + tcrossprod(E.b)
    return(list(E.b = E.b, E.bb = E.bb))
  })
  return(par)
}

m.step.VC3 = function(par, data) {
  new.beta = as.numeric(par$invsumXtX %*% Reduce(`+`, lapply(1:par$groups, function(j) {
    crossprod(data$X[[j]], data$y[[j]] - data$Z[[j]] %*% par$E[[j]]$E.b)
  })))

  new.D = Reduce(`+`, lapply(par$E, `[[`, "E.bb")) / par$groups

  new.sigmasq = (1/par$n) * sum(vapply(1:par$groups, function(j) {
    ej = data$y[[j]] - par$Xbeta[[j]] - data$Z[[j]] %*% par$E[[j]]$E.b
    return(tr(par$sigmasq * par$shortcut1[[j]]) + sum(ej^2))
  }, 1))

  par$beta = new.beta
  par$Xbeta = lapply(data$X, function(X) as.numeric(X %*% new.beta))
  par$D = new.D
  par$invD = spdinv(new.D)
  par$sigmasq = new.sigmasq
  par$shortcut1 = lapply(1:par$groups, function(j) {
    data$Z[[j]] %*% tcrossprod(spdinv(par$sigmasq * par$invD + par$ZtZ[[j]]), data$Z[[j]])
  })

  return(par)
}

cm1.step.VC3 = function(par, data) {
  new.D = Reduce(`+`, lapply(par$E, `[[`, "E.bb")) / par$groups
  new.sigmasq = (1/par$n) * sum(vapply(1:par$groups, function(j) {
    ej = data$y[[j]] - par$Xbeta[[j]] - data$Z[[j]] %*% par$E[[j]]$E.b
    return(tr(par$sigmasq * par$shortcut1[[j]]) + sum(ej^2))
  }, 1))
  par$D = new.D
  par$invD = spdinv(new.D)
  par$sigmasq = new.sigmasq

  sigmasqinvD = new.sigmasq * par$invD
  intermediate = lapply(par$ZtZ, function(ZtZ) spdinv(sigmasqinvD + ZtZ))
  new.beta = as.numeric(solve(
    par$sumXtX - Reduce(`+`, lapply(1:par$groups, function(j) par$XtZ[[j]] %*% tcrossprod(intermediate[[j]], par$XtZ[[j]]))),
    par$sumXty - Reduce(`+`, lapply(1:par$groups, function(j) par$XtZ[[j]] %*% intermediate[[j]] %*% par$Zty[[j]]))
  ))
  par$beta = new.beta
  par$Xbeta = lapply(data$X, function(X) as.numeric(X %*% new.beta))
  par$shortcut1 = lapply(1:par$groups, function(j) data$Z[[j]] %*% tcrossprod(spdinv(sigmasqinvD + par$ZtZ[[j]]), data$Z[[j]]))
  return(par)
}

cm2.step.VC3 = function(par, data) {
  new.D = Reduce(`+`, lapply(par$E, `[[`, "E.bb")) / par$groups
  par$D = new.D
  par$invD = spdinv(new.D)
  sigmasqinvD = par$sigmasq * par$invD

  intermediate = lapply(par$ZtZ, function(ZtZ) spdinv(sigmasqinvD + ZtZ))
  new.beta = as.numeric(solve(
    par$sumXtX - Reduce(`+`, lapply(1:par$groups, function(j) par$XtZ[[j]] %*% tcrossprod(intermediate[[j]], par$XtZ[[j]]))),
    par$sumXty - Reduce(`+`, lapply(1:par$groups, function(j) par$XtZ[[j]] %*% intermediate[[j]] %*% par$Zty[[j]]))
  ))
  par$beta = new.beta
  new.beta = par$beta
  par$Xbeta = lapply(data$X, function(X) as.numeric(X %*% new.beta))

  obj = function(sigmasq) {
    par$sigmasq = exp(sigmasq)
    sigmasqinvD = par$sigmasq * par$invD
    par$shortcut1 = lapply(1:par$groups, function(j) data$Z[[j]] %*% tcrossprod(spdinv(sigmasqinvD + par$ZtZ[[j]]), data$Z[[j]]))
    return(-loglik.VC3(par, data))
  }
  par$sigmasq = exp(nlm(obj, log(par$sigmasq))$estimate)
  sigmasqinvD = par$sigmasq * par$invD
  par$shortcut1 = lapply(1:par$groups, function(j) data$Z[[j]] %*% tcrossprod(spdinv(sigmasqinvD + par$ZtZ[[j]]), data$Z[[j]]))
  return(par)
}

mix.step.VC3 = function(par, old.par, alpha) {
  new.par = par
  new.par$E = lapply(1:par$groups, function(j) {
    E.b = (alpha + 1) * par$E[[j]]$E.b - alpha * old.par$E[[j]]$E.b
    E.bb = (alpha + 1) * par$E[[j]]$E.bb - alpha * old.par$E[[j]]$E.bb
    return(list(E.b = E.b, E.bb = E.bb))
  })
  return(new.par)
}

init.par.VC3 = function(data) {
  par = list()
  # constants
  par$groups = length(data$y)
  par$nj = sapply(data$y, length)
  par$n = sum(par$nj)
  par$p = ncol(data$X[[1]])
  par$q = ncol(data$Z[[1]])
  # precompute constant quantities
  par$XtX = lapply(data$X, crossprod)
  par$ZtZ = lapply(data$Z, crossprod)
  par$XtZ = lapply(1:par$groups, function(j) crossprod(data$X[[j]], data$Z[[j]]))
  par$Zty = lapply(1:par$groups, function(j) crossprod(data$Z[[j]], data$y[[j]]))
  par$Xty = lapply(1:par$groups, function(j) crossprod(data$X[[j]], data$y[[j]]))
  par$sumXtX = Reduce(`+`, par$XtX)
  par$sumXty = Reduce(`+`, par$Xty)
  par$invsumXtX = spdinv(par$sumXtX)
  # parameters
  par$beta = rep(0, ncol(data$X[[1]]))
  par$sigmasq = 1
  par$D = diag(1, ncol(data$Z[[1]]))
  par$invD = solve(par$D)
  par$Xbeta = lapply(data$X, function(X) as.numeric(X %*% par$beta))
  par$shortcut1 = lapply(1:par$groups, function(j) {
    data$Z[[j]] %*% spdinv(par$sigmasq * par$invD + par$ZtZ[[j]]) %*% t(data$Z[[j]])
  })
  par = e.step.VC3(par, data)
  return(par)
}

gen.data.VC3 = function(groups, nj, p, q, noise.sigmasq = 1, seed = sample.int(2147483647L, 1)) {
  set.seed(seed)
  data.X = lapply(1:groups, function(j) matrix(runif(nj[j] * p, -1, 1), nj[j], p))
  data.Z = lapply(1:groups, function(j) matrix(runif(nj[j] * q, -1, 1), nj[j], q))
  data.beta = rnorm(p)
  data.D = rWishart(1, q + 1, diag(q))[,,1]
  data.b = lapply(1:groups, function(j) mvtnorm::rmvnorm(1, rep(0, q), data.D)[1,])
  data.epsilon = lapply(1:groups, function(j) rnorm(nj[j], 0, sqrt(noise.sigmasq)))
  data.y = lapply(1:groups, function(j) as.numeric(data.X[[j]] %*% data.beta + data.Z[[j]] %*% data.b[[j]] + data.epsilon[[j]]))
  list(y = data.y, X = data.X, Z = data.Z)
}

# SQUAREM adapter ----

## Factor Analysis ----

fixpt.par.wrap.FA = function(par) {
  theta = c(
    par$B,
    par$D
  )
  attr(theta, "stuff") = par
  return(theta)
}

fixpt.par.unwrap.FA = function(par) {
  new.par = attr(par, "stuff")
  new.par$B = matrix(par[1:(new.par$p * new.par$q)], new.par$p, new.par$q)
  new.par$D = par[new.par$p * new.par$q + 1:new.par$p]
  return(new.par)
}

fixpt.par.wrap.2.FA = function(par) {
  theta = c(
    par$B,
    log(par$D)
  )
  attr(theta, "stuff") = par
  return(theta)
}

fixpt.par.unwrap.2.FA = function(par) {
  new.par = attr(par, "stuff")
  new.par$B = matrix(par[1:(new.par$p * new.par$q)], new.par$p, new.par$q)
  new.par$D = exp(par[new.par$p * new.par$q + 1:new.par$p])
  return(new.par)
}

## GMM ----

fixpt.par.wrap.GMM = function(par) {
  theta = c(
    par$pi,
    unlist(par$mu),
    unlist(lapply(par$sigma, function(x) x[upper.tri(x, diag = TRUE)]))
  )
  attr(theta, "stuff") = par
  return(theta)
}

fixpt.par.unwrap.GMM = function(par) {
  new.par = attr(par, "stuff")
  new.par$pi = par[1:new.par$G]
  new.par$mu = asplit(matrix(par[new.par$G + 1:(new.par$d * new.par$G)], nrow = new.par$d, ncol = new.par$G), 2)
  new.par$sigma = lapply(1:new.par$G, function(i) {
    S = matrix(0, new.par$d, new.par$d)
    S[upper.tri(S, diag = TRUE)] = par[new.par$G + new.par$d * new.par$G + (i - 1) * new.par$d * (new.par$d + 1) / 2 + 1:(new.par$d * (new.par$d + 1) / 2)]
    S[lower.tri(S, diag = FALSE)] = t(S)[lower.tri(S, diag = FALSE)]
    return(S)
  })
  return(new.par)
}

fixpt.par.wrap.2.GMM = function(par) {
  theta = c(
    log(par$pi[-1]) - log(par$pi[1]),
    unlist(par$mu),
    unlist(lapply(par$sigma, function(x) {
      x = cholesky(x)
      x[upper.tri(x, diag = TRUE)]
    }))
  )
  attr(theta, "stuff") = par
  return(theta)
}

fixpt.par.unwrap.2.GMM = function(par) {
  new.par = attr(par, "stuff")
  par = c(0, as.numeric(par))
  new.par$pi = exp(par[1:new.par$G])
  new.par$pi = new.par$pi / sum(new.par$pi)
  new.par$mu = asplit(matrix(par[new.par$G + 1:(new.par$d * new.par$G)], nrow = new.par$d, ncol = new.par$G), 2)
  new.par$sigma = lapply(1:new.par$G, function(i) {
    S = matrix(0, new.par$d, new.par$d)
    S[upper.tri(S, diag = TRUE)] = par[new.par$G + new.par$d * new.par$G + (i - 1) * new.par$d * (new.par$d + 1) / 2 + 1:(new.par$d * (new.par$d + 1) / 2)]
    S = crossprod(S)
    return(S)
  })
  return(new.par)
}

## Variance Components ----

fixpt.par.wrap.VC3 = function(par) {
  theta = c(
    par$beta,
    par$sigmasq,
    par$D[upper.tri(par$D, diag = TRUE)]
  )
  attr(theta, "stuff") = par
  return(theta)
}

fixpt.par.unwrap.VC3 = function(par) {
  new.par = attr(par, "stuff")
  new.par$beta = par[1:new.par$p]
  new.par$sigmasq = par[new.par$p + 1]
  new.par$D[upper.tri(new.par$D, diag = TRUE)] = par[-(1:(new.par$p +  1))]
  new.par$D[lower.tri(new.par$D, diag = FALSE)] = t(new.par$D)[lower.tri(new.par$D, diag = FALSE)]
  return(new.par)
}

fixpt.par.wrap.2.VC3 = function(par) {
  chol = cholesky(par$D)
  theta = c(
    par$beta,
    log(par$sigmasq),
    chol[upper.tri(chol, diag = TRUE)]
  )
  attr(theta, "stuff") = par
  return(theta)
}

fixpt.par.unwrap.2.VC3 = function(par) {
  new.par = attr(par, "stuff")
  new.par$beta = par[1:new.par$p]
  new.par$sigmasq = exp(par[new.par$p + 1])
  chol = new.par$D * 0
  chol[upper.tri(chol, diag = TRUE)] = par[-(1:(new.par$p +  1))]
  new.par$D = crossprod(chol)
  return(new.par)
}
