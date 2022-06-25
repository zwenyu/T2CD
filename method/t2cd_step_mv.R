require(arfima)
require(fracdiff)
require(splines)
require(ggplot2)
require(mvtnorm)
source('./method/bsWLS.R')
source('./helperfunction/helperfn.R')

# helper functions from univariate method
get.m <- function(x.2, dfrac, wt = rep(1, length(x.2))) {
  wt <- c(wt)
  diff.wt <- diffseries_keepmean(wt, dfrac)
  diff.wx <- diffseries_keepmean(wt*x.2, dfrac)
  sum(wt*diff.wx*diff.wt)/sum(wt*diff.wt^2)
}

get.sd <- function(x.2, dfrac, mean, wt = rep(1, length(x.2))) {
  diff_p = t(diffseries_keepmean(t(wt*(x.2 - mean)), dfrac))
  sqrt(sum(wt*diff_p^2)/sum(wt))
}

loglik_res_step = function(par, x.2, wt = rep(1, length(x.2))){
  
  dfrac = par[1]
  
  if (length(par) == 1) {
    mean <- get.m(x.2 = x.2, dfrac = dfrac, wt = wt)
    sd <- get.sd(x.2 = x.2, dfrac = dfrac, mean = mean, wt = wt)
  } else if (length(par) == 2) {
    mean <- par[2]
    sd <- get.sd(x.2 = x.2, dfrac = dfrac, mean = mean, wt = wt)
  } else if (length(par) == 3) {
    mean <- par[2]
    sd <- par[3]
  }
  
  logL = loglik_norm(x.2 = x.2, dfrac = dfrac, mean = mean, sd = sd, wt = wt)
  
  return(logL)
}

loglik_t_res_step = function(par, x.2, wt = rep(1, length(x.2))) {
  dfrac = par[1]
  mean = par[2]
  sd = par[3]
  df = par[4]
  logL = loglik_t(x.2 = x.2, dfrac = dfrac, mean = mean, sd = sd, df = df, wt = wt)
  
  return(logL)
}

loglik_t = function(x.2, dfrac, mean, sd, df, wt = rep(1, length(x.2))) {
  diff_p = c(diffseries_keepmean(wt*(x.2-mean), dfrac))
  sum(wt*ldt_mc(x = diff_p, mean = 0, sd = sd, df = df))
}

loglik_norm = function(x.2, dfrac, mean, sd, wt = rep(1, length(x.2))) {
  diff_p = c(diffseries_keepmean(wt*(x.2-mean), dfrac))
  sum(wt*ldnorm_mc(x = diff_p, mean = 0, sd = sd))
}

ldt_mc <- function(x, mean, sd, df) {
  # scale <- sd*sqrt((df - 2)/df)
  # dt((x - mean)/scale, df = df, log = TRUE) - log(scale)
  -log(pi)/2 - log(df - 2)/2 + log((gamma((df + 1)/2)/gamma(df/2))) +
    -(df + 1)*log(1 + (x - mean)^2/(sd^2*(df - 2)))/2 - log(sd^2)/2
}

ldnorm_mc <- function(x, mean, sd) {
  -log(2*pi*sd^2)/2 - (x - mean)^2/(2*sd^2)
}

fit_res_step = function(t.1, x.1, x.2, deg, seqby, resd.seqby, idx, N,
                        t_resd = FALSE, start.2 = NULL, prof = TRUE,
                        ll.1 = ll.1, cp.only = FALSE){
  
  x <- c(x.1, x.2)
  ll.1 <- c(ll.1, rep(0, length(x) - length(ll.1)))
  
  if (!cp.only) {
    if (!t_resd | (t_resd & is.null(start.2))) {
      
      if (is.null(start.2)) {
        start.norm <- 0
        if (!prof) {
          start.norm <- c(start.norm, mean(x.2), sd(x.2))
        }
      } else {
        start.norm <- start.2$d
        if (!prof) {
          start.norm <- c(start.norm, start.2$m, start.2$sd)
        }
      }
      
      optim.2 = optim(par = start.norm,
                      fn = loglik_res, method = "L-BFGS-B",
                      x = x, sigmoid = FALSE, ll.1 = ll.1,
                      t_resd = FALSE, idx = idx,
                      control = list("fnscale" = -1))
      fit_d = optim.2$par[1]
      if (length(optim.2$par) == 1) {
        fit_m <- get.m(x.2 = x.2, dfrac = fit_d)
        fit_sd <- get.sd(x.2 = x.2, dfrac = fit_d, mean = fit_m)
      } else if (length(optim.2$par) == 2) {
        fit_m <- optim.2$par[2]
        fit_sd <- get.sd(x.2 = x.2, dfrac = fit_d, mean = fit_m)
      } else {
        fit_m <- optim.2$par[2]
        fit_sd <- abs(optim.2$par[3])
      }
      ll = optim.2$value
      fit_df <- NA
    }
    
    if (t_resd) {
      
      if (is.null(start.2)) {
        start.t <- c(fit_d, fit_m, fit_sd, 3)
      } else {
        start.t <- c(start.2$d, start.2$m, start.2$sd, start.2$df)
      }
      
      opt <- optim(start.t, loglik_res, x = x,
                   lower = c(-Inf, -Inf, -Inf, 2 + 10^(-14)), upper = c(Inf, Inf, Inf, 300),
                   sigmoid = FALSE, ll.1 = ll.1,
                   t_resd = TRUE, idx = idx,
                   control = list("fnscale" = -1), method = "L-BFGS-B")
      ll = opt$value
      fit_d <- opt$par[1]
      fit_m <- opt$par[2]
      fit_sd <- abs(opt$par[3])
      fit_df <- opt$par[4]
    }
  } else {
    start.normt <- start.2$d
    if (!prof | t_resd) {
      start.normt <- c(start.normt, start.2$m, start.2$sd)
    }
    if (t_resd) {
      start.normt <- c(start.normt, start.2$df)
    }
    ll <- loglik_res(start.normt, x = x,
                     sigmoid = FALSE, ll.1 = ll.1,
                     t_resd = t_resd, idx = idx)
    fit_d <- start.2$d
    fit_m <- start.2$m
    fit_sd <- start.2$sd
    fit_df <- start.2$df
  }
  return(list(fit_M=ll,fit_d=fit_d,fit_m=fit_m, fit_sd=fit_sd, fit_df = fit_df))
}

t2cd_step = function(dat, t.max = 72, tau.range = c(10, 50), deg = 3,
                     seqby = 1, resd.seqby = 5,
                     use_scale = TRUE, t_resd = FALSE, prof = TRUE){
  # dat: input time series
  # t.max: cutoff for time considered
  # tau.range candidate change point range
  # deg: degree for B-spline
  # segby, resd.seqby: interval between knots
  # use_scale: if true, scale time series
  res1 = search_dtau_step(dat, t.max, tau.range, deg,
                          seqby = seqby, resd.seqby = resd.seqby,
                          use_scale = use_scale, t_resd = t_resd, prof = prof)
  return(res1)
}

search_dtau_step = function(dat, t.max = 72, tau.range = c(10, 50), deg = 3,
                            seqby = 1, resd.seqby = 5,
                            use_scale = TRUE, t_resd = FALSE, prof = TRUE,
                            cp.only = FALSE,
                            fix.d = NULL, fix.m = NULL, fix.sd = NULL,
                            fix.df = NA) {
  # select data below t.max
  if (is.na(t.max)){
    t.max = max(dat$tim, na.rm = T)
  }
  
  tim.ind = !is.na(dat$tim) & dat$tim <= t.max
  t.maxidx = which(tim.ind == T)
  
  res = dat$res[t.maxidx]
  if (use_scale){
    res_mean = scale(res[t.maxidx], center = F) # scaling
  }else{
    res_mean = matrix(res[t.maxidx], ncol = 1)
  }
  tim = dat$tim[t.maxidx]
  N = length(res_mean)
  
  # initialize result vectors
  tau.idx = which(tim >= tau.range[1] & tim <= tau.range[2])
  
  # iterate through each tau, return log-likelihood
  foreach.tau <- vector("list", length = length(tau.idx))
  for (j in 1:length(tau.idx)) {
    # cat("j=", j, "\n")
    idx = tau.idx[j]
    
    # optimize for polynomial component
    x.1 = res_mean[1:idx]
    t.1 = tim[1:idx]
    n.1 = length(x.1)
    
    fit1 = refitWLS(t.1, x.1, deg = deg, seqby = seqby, resd.seqby = resd.seqby)
    resd1.1 = x.1 - fit1$fit.vals
    var.resd1.1 = fit1$var.resd
    ll.1 = dnorm(resd1.1, log = TRUE, sd = sqrt(var.resd1.1))
    
    # optimize for ARFIMA
    x.2 = res_mean[(idx+1):N]
    
    if (!cp.only) {
      if (j == 1) {
        start.2 <- NULL
      } else {
        start.2 <- list("d" = foreach.tau[[j-1]]$d,
                        "m" = foreach.tau[[j-1]]$m,
                        "sd" = foreach.tau[[j-1]]$sd,
                        "df" = foreach.tau[[j-1]]$df)
      }
    } else {
      start.2 <- list("d" = fix.d, "m" = fix.m,
                      "sd" = fix.sd, "df" = fix.df)
    }
    
    fit_res = tryCatch(fit_res_step(t.1 = t.1,  x.1 = x.1, x.2 = x.2, deg = deg,
                                    seqby = seqby, resd.seqby = resd.seqby, idx = idx,
                                    N = N, t_resd = t_resd,
                                    start.2 = start.2, prof = prof, ll.1 = ll.1,
                                    cp.only = cp.only),
                       error = function(e) {return(NA)})
    if ((!t_resd & any(is.na(fit_res[1:4]))) | (t_resd & any(is.na(fit_res)))) {
      foreach.tau[[j]] <- list("M" = -Inf,
                               "d" = NA, "m" = NA, "df" = NA,
                               "sd" = NA,
                               "fit.vals" = rep(NA, length(fit1$fit.vals)),
                               "var.resd" = rep(NA, length(fit1$var.resd)))
    }else{
      foreach.tau[[j]] <- list("M" = fit_res$fit_M,
                               "d" = fit_res$fit_d,
                               "m" = fit_res$fit_m,
                               "sd" = fit_res$fit_sd,
                               "df" = fit_res$fit_df,
                               "fit.vals" = fit1$fit.vals,
                               "var.resd" = fit1$var.resd)
    }
    
  }
  
  M <- unlist(lapply(foreach.tau, function(x) {x$M}))
  d <- unlist(lapply(foreach.tau, function(x) {x$d}))
  m <- unlist(lapply(foreach.tau, function(x) {x$m}))
  sd <- unlist(lapply(foreach.tau, function(x) {x$sd}))
  df <- unlist(lapply(foreach.tau, function(x) {x$df}))
  
  # tau and d at maximum log-likelihood
  M_df = data.frame(tau = tim[tau.idx], M = M, d = d, m = m, sd = sd, df = df)
  max.idx = which.max(M)
  max.tau = tim[tau.idx[max.idx]]
  max.d = d[max.idx]
  max.m = m[max.idx]
  max.sd = sd[max.idx]
  max.df = df[max.idx]
  fit.vals <- foreach.tau[[max.idx]]$fit.vals
  var.resd <- foreach.tau[[max.idx]]$var.resd
  
  retlist <- list(M_df = M_df, res = res, tim = tim, tau.idx = tau.idx,
                  tau = max.tau, d = max.d, m = max.m, sd = max.sd, df = max.df,
                  idx = tau.idx[max.idx], logL = M[max.idx],
                  fit.vals = fit.vals, var.resd = var.resd,
                  ll.1 = ll.1)
  
  return(retlist)
}

# helper functions for multivariate loglikelihood
loglik_res_step_mv = function(param, x, idx, ll.1,
                              t_resd = FALSE){
  dfrac = param[1]
  if (t_resd) {
    p <- nrow(x)
    mean <- param[1 + 1:p]
    sd <- param[1 + p + 1:p]
    df <- param[1 + 2*p + 1:p]
  }
  
  logL <- 0
  for (k in 1:nrow(x)) {
    if (!t_resd) {
      par <- dfrac
    } else {
      par <- c(dfrac, mean[k], sd[k], df[k])
    }
    N.k <- max(which(!is.na(ll.1[k, ])))
    ll <- loglik_res(par = par, x = x[k, 1:N.k],
                     sigmoid = FALSE, ll.1 = ll.1[k, 1:N.k],
                     t_resd = t_resd, idx = idx[k])
    
    logL <- logL + ll
  }
  
  return(logL)
}

# multivariate implementation for T2CD-step
# grid search for d and tau
t2cd_step_mv = function(dat, t.max = 72, tau.range = c(10, 50), deg = 3,
                        seqby = 1, resd.seqby = 5,
                        use_scale = TRUE,
                        maxiter=10, tol=1e-6, t_resd = FALSE){
  # dat: input time series
  # t.max: cutoff for time considered
  # tau.range candidate change point range
  # deg: degree for B-spline
  # seqby, resd.seqby: interval between knots
  # use_scale: if true, scale time series
  # maxiter: maximum number of iterations solving for tau and FI parameters
  # tol: stops iterating if difference between the most recent objective values is below tol
  
  if (!is.matrix(dat$res)) {
    dat$res <- matrix(dat$res, nrow = 1, ncol = length(dat$res))
  }
  if (!is.matrix(dat$tim)) {
    dat$tim <- matrix(dat$tim, nrow = 1, ncol = length(dat$tim))
  }
  
  # select data below t.max
  if (is.na(t.max)){
    t.max = min(apply(dat$tim, 1, max, na.rm = T))
  }
  
  max.tim.ind =
    max(apply(dat$tim, 1, function(x) {max(which(!is.na(x) & x <= t.max))}))
  
  res = matrix(dat$res[, 1:max.tim.ind], nrow = nrow(dat$res))
  tim = matrix(dat$tim[, 1:max.tim.ind], nrow = nrow(dat$tim))
  
  res[!is.na(tim) & tim > t.max] <- NA
  tim[!is.na(tim) & tim > t.max] <- NA
  if (use_scale){
    res_mean = t(scale(t(res), center = F)) # scaling
  }else{
    if (!is.matrix(res)) {
      res_mean = matrix(res, ncol = 1)
    } else {
      res_mean = res
    }
  }
  N = ncol(res_mean)
  p = nrow(res_mean)
  
  # alternate between optimizing for tau and FI parameters
  iter_k = 1
  iter_flag = TRUE
  neglogL_prev = Inf
  init.d = init.m = init.sd = init.df =
    fix.d = fix.m = fix.sd = fix.df =
    init.tau = init.idx = tau = idx = rep(NA, p)
  hist.d = hist.neglogL = c()
  
  while (iter_flag){
    # retain all data in regime 2
    reg2 = matrix(NA, nrow = p, ncol = N)
    fit.vals <- var.resd <- ll.1 <- matrix(NA, nrow = p, ncol = N)
    for (k in 1:p){
      cat("k=", k, "\n")
      if (iter_k==1){
        res_k = t2cd_step(list(res=res_mean[k,], tim=tim[k,]),
                          t.max, tau.range, deg,
                          seqby = seqby, resd.seqby = resd.seqby,
                          use_scale = FALSE, t_resd = t_resd)
        init.d[k] = res_k$d
        init.m[k] = res_k$m
        init.sd[k] = res_k$sd
        init.df[k] = res_k$df
        init.tau[k] = res_k$tau
        init.idx[k] = res_k$idx
        ll.1[k, 1:length(na.omit(res_mean[k,]))] <-
          c(res_k$ll.1, rep(0, length(na.omit(res_mean[k,])) - length(res_k$ll.1)))
        fit.vals[k, 1:length(res_k$fit.vals)] <- res_k$fit.vals
        var.resd[k, 1:length(res_k$var.resd)] <- res_k$var.resd
      }else{
        res_k = search_dtau_step(list(res=res_mean[k,], tim=tim[k,]), t.max,
                                 tau.range, deg,
                                 seqby = seqby, resd.seqby = resd.seqby,
                                 use_scale = FALSE,
                                 fix.d = d_current, fix.m = m_current[k],
                                 fix.sd = sd_current[k], fix.df = df_current[k],
                                 t_resd = t_resd, cp.only = TRUE)
        fix.d[k] = res_k$d
        fix.m[k] = res_k$m
        fix.sd[k] = res_k$sd
        fix.df[k] = res_k$df
        tau[k] = res_k$tau
        idx[k] = res_k$idx
        ll.1[k, 1:length(na.omit(res_mean[k,]))] <-
          c(res_k$ll.1, rep(0, length(na.omit(res_mean[k,])) - length(res_k$ll.1)))
        fit.vals[k, 1:length(res_k$fit.vals)] <- res_k$fit.vals
        var.resd[k, 1:length(res_k$var.resd)] <- res_k$var.resd
      }
      reg2[k, (res_k$idx+1):N] = res_mean[k, (res_k$idx+1):N]
    }
    
    if (iter_k==1){
      d_current = mean(init.d)
      m_current = init.m
      sd_current = init.sd
      df_current = init.df
      tau = init.tau
      idx = init.idx
    }
    
    if (p > 1) {
      x = reg2
      # optimizing over d, m (fixed taus, I think)
      if (!t_resd) {
        optim_params = optim(par = d_current,
                             fn = loglik_res_step_mv, method = "Brent", lower = -10,
                             upper = 10, x = x, idx = idx, ll.1 = ll.1,
                             t_resd = t_resd, control = list("fnscale" = -1))
        d_current = optim_params$par
        df_current <- sd_current <- m_current <- rep(NA, nrow(res))
        for (l in 1:length(m_current)) {
          m_current[l] <-  get.m(x.2 = na.omit(x[l, (idx[l] + 1):N]), dfrac = d_current)
          sd_current[l] <-  get.sd(x.2 = na.omit(x[l, (idx[l] + 1):N]), dfrac = d_current,
                                   mean = m_current[l])
        }
        neglogL_current = -optim_params$value
      } else {
        optim_params = optim(par = c(d_current, m_current, sd_current, df_current),
                             fn = loglik_res_step_mv, method = "L-BFGS-B",
                             lower = c(-Inf, rep(c(-Inf, -Inf, 2 + 10^(-14)), each = length(m_current))),
                             upper = c(Inf, rep(c(Inf, Inf, 300), each = length(m_current))),
                             x = x, idx = idx, ll.1 = ll.1,
                             t_resd = t_resd, control = list("fnscale" = -1))
        d_current = optim_params$par[1]
        m_current = optim_params$par[1 + 1:length(m_current)]
        sd_current = optim_params$par[(1 + length(m_current)) + 1:length(m_current)]
        df_current = optim_params$par[(1 + 2*length(m_current)) + 1:length(m_current)]
        neglogL_current = -optim_params$value
      }
      
      # update iter_flag
      if (iter_k==maxiter | abs(neglogL_prev-neglogL_current)<tol){
        iter_flag = FALSE
      }else{
        iter_k = iter_k + 1
      }
      # update
      hist.d = c(hist.d, d_current)
      hist.neglogL = c(hist.neglogL, neglogL_current)
      neglogL_prev = neglogL_current
    } else {
      iter_flag = FALSE
    }
    
  }
  
  return(list(res = res, res_mean = res_mean, tim = tim,
              tau = tau, idx = idx, d = d_current,
              m = m_current, sd = sd_current, df = df_current,
              univ_tau = init.tau, univ_idx = init.idx, univ_d = init.d,
              hist.d = hist.d, hist.neglogL = hist.neglogL,
              iter_k = iter_k,
              fit.vals = fit.vals, var.resd = var.resd))
}

# plot sequences and fitted lines
plot_t2cd_step_mv = function(results, tau.range = c(10, 50),
                             use_scale = TRUE, return_plot = FALSE) {
  res = results$res
  res_mean = results$res_mean
  tim = results$tim
  tau.idx = results$tau.idx
  N = ncol(res)
  
  # select optimal parameters
  opt_d = results$d
  opt_tau = results$tau
  opt_idx = results$idx
  m <- results$m
  
  fit.vals <- var.resd1 <- var.resd1.1 <- mu <- fit1 <-
    matrix(nrow = nrow(res_mean), ncol = ncol(res_mean))
  
  for (k in 1:nrow(res)) {
    ### fitted values for first regime
    fit1[k, 1:opt_idx[k]] = results$fit.vals[k, 1:opt_idx[k]]
    var.resd1.1[k, 1:opt_idx[k]] = results$var.resd[k, 1:opt_idx[k]]
    N.k <- length(na.omit(res_mean[k, ]))
    x.2 = res_mean[k, (opt_idx[k]+1):N.k]
    mu[k, 1:N.k] = c(fit1[k, 1:opt_idx[k]],
                     trend_fi(s = x.2, eff.d = opt_d,  mu = m[k]))
    
    # Back on original scale
    if (use_scale){
      fit.vals[k, ] = mu[k, ]*(attributes(res_mean)$'scaled:scale'[k])
      var.resd1[k, 1:opt_idx[k]] =
        var.resd1.1[k, 1:opt_idx[k]]*(attributes(res_mean)$'scaled:scale'[k])^2
    }else{
      fit.vals[k, ] = mu[k, ]
      var.resd1[k, 1:opt_idx[k]] = var.resd1.1[k, 1:opt_idx[k]]
    }
  }
  
  # plotting
  if (return_plot){
    plot(tim[1, ], res[1, ],
         ylim = c(min(c(res[1, ], fit.vals[1, ])),
                  max(c(res[1, ], fit.vals[1, ]))), type = 'l',
         main = paste('Values fitted with d: ', round(opt_d,3), ' tau: ', round(opt_tau,3)),
         xlab = 'Time (hour)', ylab = 'Resistance (ohm)')
    if (is.na(opt_tau)){
      lines(tim[1, ], fit.vals[1, ], col = "blue", lwd = 1)
    }else{
      opt_tau.idx = which(tim[1, ] == opt_tau)
      lines(tim[1, 1:opt_idx], fit.vals[1, 1:opt_idx], col = "blue", lwd = 1)
      lines(tim[1, (opt_idx+1):N], fit.vals[1, (opt_idx+1):N], col = "green", lwd = 1)
      abline(v = opt_tau, lty = 2, col = "red")
    }
    abline(v = tau.range, lty = 1, col = "red")
  }
  
  return(list(fit.vals = fit.vals,
              opt_idx = opt_idx, N = N,
              scaling = attributes(res_mean)$'scaled:scale',
              var.resd1 = var.resd1))
}

# parametric bootstrap using outputs from t2cd_step and plot.t2cd_step
bootstrap_sample_step_mv = function(results, plot_results, seed = 0){
  
  set.seed(seed)
  res = results$res
  tim = results$tim
  N = ncol(res)
  opt_idx = results$idx
  
  samp <- matrix(nrow = nrow(res), ncol = ncol(res))
  
  for (k in 1:nrow(res)) {
    # regime 1
    fit.vals1 = plot_results$fit.vals[k, 1:opt_idx[k]]
    var.resd1 = plot_results$var.resd1[k, 1:opt_idx[k]]
    noise1 = rnorm(opt_idx[k], 0, sqrt(var.resd1))
    
    # regime 2
    opt_d = results$d
    m = results$m[k]
    sd = results$sd[k]
    df <- results$df[k]
    sim = sim_fi(N-opt_idx[k], opt_d, mu = m, sig = sd, df = df)
    seq_fi = sim$s
    
    samp[k, ] = c(fit.vals1 + noise1, seq_fi*plot_results$scaling[k])
  }
  
  return(list(res=samp, tim=tim))
}