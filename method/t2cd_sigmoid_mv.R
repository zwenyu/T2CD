require(arfima)
require(fracdiff)
require(splines)
require(ggplot2)
require(mvtnorm)
source('./method/bsWLS.R')
source('./helperfunction/helperfn.R')

# helper functions from univariate method
loglik_res_sigmoid_alpha = function(param,
                                    tim_cp = NULL, tau.idx = NULL,
                                    x, ll.1,
                                    pen = FALSE, C = NULL, t_resd = FALSE,
                                    rest) {
  loglik_res(param = c(param, rest), tim_cp = tim_cp, tau.idx = tau.idx,
             x = x, ll.1 = ll.1, pen = pen, C = C, t_resd = t_resd,
             sigmoid = TRUE)
}

loglik_res_sigmoid_rest = function(param,
                                   tim_cp = NULL, tau.idx = NULL,
                                   x, ll.1,
                                   pen = FALSE, C = NULL, t_resd = FALSE,
                                   alpha) {
  loglik_res(param = c(alpha, param), tim_cp = tim_cp, tau.idx = tau.idx,
             x = x, ll.1 = ll.1, pen = pen, C = C, t_resd = t_resd,
             sigmoid = TRUE)
}

get.wt <- function(alpha, tim_cp, tau.idx, N) {
  alpha0 <- alpha[1]
  alpha1 <- alpha[2]
  wt_cp = sigmoid(alpha0+alpha1*tim_cp) # 0 to 1
  wt = c(rep(wt_cp[1], tau.idx[1]-1),
         wt_cp,
         rep(ifelse(wt_cp[length(wt_cp)] != 0, wt_cp[length(wt_cp)], 1),
             N-tau.idx[length(tau.idx)]))
  return(wt)
}

loglik_res = function(param, sigmoid = TRUE,
                      tim_cp = NULL, tau.idx = NULL,
                      x, ll.1,
                      pen = FALSE, C = NULL, t_resd = FALSE,
                      idx = NULL) {
  
  if (sigmoid) {
    alpha0 = param[1]
    alpha1 = param[2]
    
    # weights
    wt <- get.wt(alpha = c(alpha0, alpha1), tim_cp = tim_cp, tau.idx = tau.idx,
                 N = length(x))
    k <- 2
  } else {
    wt <- c(rep(0, idx), rep(1, (length(x) - idx)))
    k <- 0
  }
  dfrac = param[k + 1]
  
  if (t_resd) {
    df = param[length(param)]
  }
  
  if (((length(param) - k) == 1 & !t_resd)) {
    m <- get.m(x.2 = x, dfrac = dfrac, wt = wt)
    sd <- get.sd(x.2 = x, dfrac = dfrac, mean = m, wt = wt)
  } else if ((length(param) - k) == 2 & !t_resd) {
    m <- param[k + 2]
    sd <- get.sd(x.2 = x, dfrac = dfrac, mean = m, wt = wt)
  } else if (((length(param) - k) == 3) & !t_resd | t_resd) {
    m <- param[k + 2]
    sd <- param[k + 3]
  }
  
  logL = sum((1-wt)*ll.1) +
    ifelse(!t_resd,
           loglik_res_step(par = c(dfrac, m, sd), x.2 = x, wt = wt),
           loglik_t_res_step(par = c(dfrac, m, sd, df), x.2 = x, wt = wt)) +
    ifelse(pen, C*sum(wt[length(wt)] - wt[1]), 0)
  
  return(logL)
}

t2cd_sigmoid = function(dat, t.max = 72, tau.range = c(10, 50),
                        init.tau = c(15, 30, 45), deg = 3, C = 1000,
                        seqby = 1, resd.seqby = 5, use_scale = TRUE,
                        t_resd = FALSE, prof = TRUE, several.init = TRUE){
  # dat: input time series
  # t.max: cutoff for time considered
  # tau.range candidate change point range
  # init.tau: candidate taus to initialize learning
  # deg: degree for B-spline
  # C: regularization coefficient
  # segby, resd.seqby: interval between knots
  # use_scale: if true, scale time series
  res1 = search_dtau_sigmoid(dat, t.max, tau.range, init.tau, deg, C,
                             seqby = seqby, resd.seqby = resd.seqby, use_scale = use_scale,
                             t_resd = t_resd, prof = prof, several.init = several.init)
  # increase C if change point not in candidate range
  multiplier = 2
  while (is.na(res1$tau)){
    C_new = multiplier*C
    res1 = search_dtau_sigmoid(dat, t.max, tau.range, init.tau, deg, C_new,
                               seqby = seqby, resd.seqby = resd.seqby, use_scale = use_scale,
                               t_resd = t_resd, prof = prof, several.init = several.init)
    multiplier = multiplier + 1
  }
  return(res1)
}

sigmoid = function(a){return(1/(1+exp(-a)))}
softmax = function(a,b){
  return(exp(a)/(exp(a) + exp(b)))
}

search_dtau_sigmoid = function(dat, t.max = 72, tau.range = c(10, 50),
                               init.tau = c(15, 30, 45), deg = 3, C = 1000,
                               seqby = 1, resd.seqby = 5, use_scale = T,
                               t_resd = FALSE, prof = TRUE, several.init = T){
  # select data below t.max
  if (is.na(t.max)){
    t.max = min(apply(dat$tim, 1, max, na.rm = T))
  }
  
  tim.ind = !is.na(dat$tim) & dat$tim <= t.max
  t.maxidx = which(apply(tim.ind, 2, all) == T)
  
  res = matrix(dat$res[, t.maxidx], nrow = nrow(dat$res))
  tim = matrix(dat$tim[, t.maxidx], nrow = nrow(dat$tim))
  if (use_scale){
    res_mean = t(scale(t(res), center = F)) # scaling
  }else{
    res_mean = res
  }
  N = ncol(res_mean)
  
  # points in change range
  tau.idx = which(apply(tim >= tau.range[1] & tim <= tau.range[2], 2, all))
  tim_cp = matrix(tim[,tau.idx], nrow = nrow(tim))
  
  # optimize for spline component
  fit1 = refitWLS(tim, res_mean, deg = deg,
                  seqby = seqby, resd.seqby = resd.seqby)
  resd1 = res_mean - fit1$fit.vals
  var.resd1 = fit1$var.resd
  ll.1.i = dnorm(resd1, log = TRUE, sd = sqrt(var.resd1))
  
  x = res_mean
  lastN = N
  
  opt_logL = -Inf
  if (several.init) {
    for (tau_i in init.tau){
      tau_i.idx = which.min(apply(tim <= tau_i, 2, all) == T)
      start <- c(-tau_i, 1, 0)
      lower <- rep(-Inf, length(start))
      upper <- rep(Inf, length(start))
      if (!prof & !t_resd) {
        start <- c(start, mean(x[tau_i.idx:lastN]), sd(x[tau_i.idx:lastN]))
        lower <- c(lower, rep(-Inf, 2))
        upper <- c(upper, rep(Inf, 2))
      }
      if (t_resd) {
        start <- c(start, mean(x[tau_i.idx:lastN]), sd(x[tau_i.idx:lastN]), 3)
        lower <- c(lower, -Inf, -Inf, 2 + 10^(-14))
        upper <- c(upper, Inf, Inf, 300)
      }
      optim_i = optim(par = start,
                      fn = loglik_res, method = "L-BFGS-B",
                      tim_cp = tim_cp, tau.idx = tau.idx, x = x,
                      ll.1 = ll.1.i, pen = TRUE, C = C,
                      control = list("fnscale" = -1), t_resd = t_resd,
                      lower = lower, upper = upper, sigmoid = TRUE)
      logL_i = loglik_res(optim_i$par,
                          tim_cp = tim_cp, tau.idx = tau.idx,
                          x = x, ll.1 = ll.1.i, t_resd = t_resd,
                          sigmoid = TRUE)
      if (logL_i > opt_logL){
        opt_logL = logL_i
        optim_params = optim_i
        fit.vals <- fit1$fit.vals
        var.resd <- var.resd1
        ll.1 = ll.1.i
      }
    }
    opt_param = optim_params$par
  } else {
    
    tau_i = 50
    tau_i.idx = which.min(apply(tim <= tau_i, 2, all) == T)
    
    start <- c(0)
    lower <- rep(-Inf, length(start))
    upper <- rep(Inf, length(start))
    if (!prof & !t_resd) {
      start <- c(start, mean(x[tau_i.idx:lastN]), sd(x[tau_i.idx:lastN]))
      lower <- c(lower, rep(-Inf, 2))
      upper <- c(upper, rep(Inf, 2))
    }
    if (t_resd) {
      start <- c(start, mean(x[tau_i.idx:lastN]), sd(x[tau_i.idx:lastN]), 3)
      lower <- c(lower, -Inf, -Inf, 2 + 10^(-14))
      upper <- c(upper, Inf, Inf, 300)
    }
    
    initial <- optim(par = start,
                     fn = loglik_res, method = "L-BFGS-B",
                     x = x, idx = tau_i.idx,
                     ll.1 = ll.1,
                     control = list("fnscale" = -1), t_resd = t_resd,
                     lower = lower, upper = upper, sigmoid = FALSE)
    
    
    start <- c(-tau_i, 1)
    lower.alpha <- c(-Inf, -Inf)
    upper.alpha <- c(Inf, Inf)
    
    optim_alpha = optim(par = start,
                        fn = loglik_res_sigmoid_alpha, method = "L-BFGS-B",
                        tim_cp = tim_cp, tau.idx = tau.idx, x = x,
                        ll.1 = ll.1, pen = TRUE, C = C,
                        control = list("fnscale" = -1), t_resd = t_resd,
                        lower = lower.alpha, upper = upper.alpha, rest = initial$par)
    
    optim_rest = optim(par = initial$par,
                       fn = loglik_res_sigmoid_rest, method = "L-BFGS-B",
                       tim_cp = tim_cp, tau.idx = tau.idx, x = x,
                       ll.1 = ll.1, pen = TRUE, C = C,
                       control = list("fnscale" = -1), t_resd = t_resd,
                       lower = lower, upper = upper, alpha = optim_alpha$par)
    
    logL_i = loglik_res(c(optim_alpha$par, optim_rest$par),
                        tim_cp = tim_cp, tau.idx = tau.idx,
                        x = x, ll.1 = ll.1, t_resd = t_resd,
                        sigmoid = TRUE)
    while (abs(logL_i - opt_logL) > 10^(-3)) {
      
      opt_logL <- logL_i
      
      optim_alpha = optim(par = optim_alpha$par,
                          fn = loglik_res_sigmoid_alpha, method = "L-BFGS-B",
                          tim_cp = tim_cp, tau.idx = tau.idx, x = x,
                          ll.1 = ll.1, pen = TRUE, C = C,
                          control = list("fnscale" = -1), t_resd = t_resd,
                          lower = lower.alpha, upper = upper.alpha, rest = optim_rest$par)
      
      optim_rest = optim(par = optim_rest$par,
                         fn = loglik_res_sigmoid_rest, method = "L-BFGS-B",
                         tim_cp = tim_cp, tau.idx = tau.idx, x = x,
                         ll.1 = ll.1, pen = TRUE, C = C,
                         control = list("fnscale" = -1), t_resd = t_resd,
                         lower = lower, upper = upper, alpha = optim_alpha$par)
      
      logL_i = loglik_res(c(optim_alpha$par, optim_rest$par),
                          tim_cp = tim_cp, tau.idx = tau.idx,
                          x = x, ll.1 = ll.1, t_resd = t_resd,
                          sigmoid = TRUE)
      
      
      
    }
    opt_param = c(optim_alpha$par, optim_rest$par)
  }
  
  alpha0 = opt_param[1]
  alpha1 = opt_param[2]
  opt_d = opt_param[3]
  
  
  # weights
  wt <- get.wt(alpha = c(alpha0, alpha1), tim_cp = tim_cp, tau.idx = tau.idx,
               N = length(x))
  
  opt_tau.idx = which(wt>=0.5, arr.ind = TRUE)[1]-1
  opt_tau = tim[opt_tau.idx[1]]
  
  # range of change location
  opt_taurange1.idx = which(wt>=0.1, arr.ind = TRUE)[1]-1
  opt_taurange2.idx = which(wt>=0.9, arr.ind = TRUE)[1]-1
  opt_taurange1 = tim[opt_taurange1.idx[1]]
  opt_taurange2 = tim[opt_taurange2.idx[1]]
  opt_taurange1[is.na(opt_taurange1)] = tau.range[1]
  opt_taurange2[is.na(opt_taurange2)] = tau.range[2]
  
  if (t_resd) {
    opt_m <- opt_param[4]
    opt_sd <- abs(opt_param[5])
    opt_df <- opt_param[6]
  } else {
    if (prof) {
      opt_m <- get.m(x.2 = x, dfrac = opt_d, wt = wt)
      opt_sd <- get.sd(x.2 = x, dfrac = opt_d, mean = opt_m, wt = wt)
    } else if (!prof) {
      opt_m <- opt_param[4]
      opt_sd <- abs(opt_param[5])
    }
    opt_df = NA
  }
  
  return(list(res = res, tim = tim, tau.idx = tau.idx, m = opt_m,
              sd = opt_sd, df = opt_df,
              d = opt_d, tau = opt_tau, idx = opt_tau.idx,
              tau.range1 = opt_taurange1, tau.range2 = opt_taurange2,
              param = opt_param, logL = opt_logL,
              fit.vals = fit.vals, var.resd = var.resd,
              ll.1 = ll.1))
}

# helper functions for multivariate loglikelihood
loglik_res_sigmoid_mv = function(param, tim_cp, tau.idx, p, x,
                                 ll.1.mat, alpha0, alpha1,
                                 C = NULL, pen = FALSE, t_resd = FALSE){
  dfrac = param[1]
  if (t_resd) {
    m <- param[1 + 1:p]
    sd <- param[1 + p + 1:p]
    df <- param[1 + 2*p + 1:p]
  }
  
  logL = 0
  for (k in 1:p) {
    parev <- c(alpha0[k], alpha1[k], dfrac)
    if (t_resd) {
      parev <- c(parev, m[k], sd[k], df[k])
    }
    logL = logL +
      loglik_res(param = parev, sigmoid = TRUE,
                 tim_cp = na.omit(tim_cp[k, ]), tau.idx = na.omit(tau.idx[k, ]),
                 x = na.omit(x[k, ]), ll.1 = na.omit(ll.1.mat[k, ]),
                 pen = pen, t_resd = t_resd, C = C)
  }
  
  return(logL)
}

# multivariate implementation for T2CD-sigmoid
# optimize the likelihood for d and tau
# option to initialize tau at multiple indices
# loglikelihood, penalty to enforce tau within tau.range
t2cd_sigmoid_mv = function(dat, t.max = 72, tau.range = c(10, 50),
                           init.tau = c(15, 30, 45), deg = 3, C = 1000,
                           seqby = 1, resd.seqby = 5, use_scale = TRUE, t_resd = FALSE){
  # dat: input time series
  # t.max: cutoff for time considered
  # tau.range candidate change point range
  # init.tau: candidate taus to initialize learning
  # deg: degree for B-spline
  # C: regularization coefficient
  # seqby, resd.seqby: interval between knots
  # use_scale: if true, scale time series
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
    res_mean = res
  }
  N = ncol(res)
  p = nrow(res)
  x = res_mean
  
  # points in change range
  tim_cp = tim
  tim_cp[!is.na(tim) & (tim < tau.range[1] | tim > tau.range[2])] <- NA
  
  min.tau.idx <- apply(tim_cp, 1, function(x) {min(which(!is.na(x)))})
  max.tau.idx <- apply(tim_cp, 1, function(x) {max(which(!is.na(x)))})
  
  tau.idx <- matrix(NA, nrow = nrow(tim), ncol = max(max.tau.idx - min.tau.idx) + 1)
  for (k in 1:p) {
    tau.idx[k, 1:(max.tau.idx[k] - min.tau.idx[k] + 1)] <- min.tau.idx[k]:max.tau.idx[k]
  }
  
  tim_cp = tim_cp[, !apply(tim_cp, 2, function(x) {sum(is.na(x))}) == nrow(tim_cp), drop = FALSE]
  
  # determine range of d to search and find initializing parameters
  if (!t_resd) {
    init.param = matrix(NA, p, 3)
  } else {
    init.param = matrix(NA, p, 6)
  }
  init.d = rep(NA, p)
  fit.vals <- var.resd <- ll.1.mat <- matrix(NA, ncol = N, nrow = p)
  for (k in 1:p){
    res_k = t2cd_sigmoid(list(res=matrix(x[k,], 1), tim=matrix(tim[k,], 1)),
                         t.max, tau.range, init.tau, deg, C,
                         seqby = seqby, resd.seqby = resd.seqby,
                         use_scale = FALSE, t_resd = t_resd)
    init.param[k,] = res_k$par
    init.d[k] = res_k$d
    ll.1.mat[k, 1:length(res_k$ll.1)] <- res_k$ll.1
    fit.vals[k, 1:length(res_k$fit.vals)] <- res_k$fit.vals
    var.resd[k, 1:length(res_k$var.resd)] <- res_k$var.resd
  }
  
  alpha0 = init.param[,1]
  alpha1 = init.param[,2]
  univ_d = init.param[,3]
  opt_d = mean(univ_d)
  
  if (p > 1) {
    start <- c(mean(init.d))
    lower <- -Inf
    upper <- Inf
    if (t_resd) {
      start <- c(start, init.param[, 4], init.param[, 5], init.param[, 6])
      lower <- c(lower, rep(-Inf, 2*p), rep(2 + 10^(-14), p))
      upper <- c(upper, rep(Inf, 2*p), rep(300, p))
    }
    
    optim_params = optim(par = start,
                         fn = loglik_res_sigmoid_mv, method = "L-BFGS-B",
                         lower = lower, upper = upper,
                         tim_cp = tim_cp, tau.idx = tau.idx, p = p, x = x,
                         ll.1.mat = ll.1.mat, C = C,
                         alpha0 = alpha0, alpha1 = alpha1, pen = TRUE,
                         control = list("fnscale" = -1), t_resd = t_resd)
    opt_param = c(alpha0, alpha1, optim_params$par)
    opt_logL = loglik_res_sigmoid_mv(optim_params$par,
                                     tim_cp = tim_cp, tau.idx = tau.idx,
                                     p = p, x = x,
                                     ll.1.mat = ll.1.mat, alpha0 = alpha0,
                                     alpha1 = alpha1, t_resd = t_resd)
    opt_d = optim_params$par[1]
  } else {
    opt_param = res_k$param
    opt_logL = res_k$logL
  }
  # weights
  wt <- matrix(NA, nrow = p, ncol = N)
  for (k in 1:p) {
    wt[k, 1:length(na.omit(x[k, ]))] <-
      get.wt(alpha = c(alpha0[k], alpha1[k]),
             tim_cp = na.omit(tim_cp[k, ]),
             tau.idx = na.omit(tau.idx[k, ]),
             N = length(na.omit(x[k, ])))
  }
  
  df <- sd <- m <- rep(NA, nrow(x))
  
  for (k in 1:length(m)) {
    if (!t_resd) {
      m[k] <- get.m(x.2 = na.omit(x[k, ]), dfrac = opt_d, wt = na.omit(wt[k, ]))
      sd[k] <- get.sd(x.2 = na.omit(x[k, ]), dfrac = opt_d, mean = m[k],
                      wt = na.omit(wt[k, ]))
    } else {
      m[k] <- opt_param[2*p + 1 + k]
      sd[k] <- abs(opt_param[2*p + 1 + p + k])
      df[k] <- opt_param[2*p + 1 + 2*p + k]
    }
  }
  
  
  opt_tau.idx = apply(wt, 1, function(x){return(which(x>=0.5, arr.ind = TRUE)[1]-1)})
  opt_tau = c()
  for (k in 1:p){
    opt_tau = c(opt_tau, tim[k,opt_tau.idx[k]])
  }
  
  return(list(res = res, res_mean = res_mean, tim = tim,
              tim_cp = tim_cp,
              tau.idx = tau.idx,
              idx = opt_tau.idx,
              d = opt_d, univ_d = univ_d, tau = opt_tau, param = opt_param,
              logL = opt_logL,
              m = m, sd = sd, df = df,
              fit.vals = fit.vals, var.resd = var.resd))
}

# plot sequences and fitted lines
plot_t2cd_sigmoid_mv = function(results, tau.range = c(10, 50),
                                return_plot = TRUE, use_scale = TRUE){
  res = results$res
  tim = results$tim
  tim_cp = results$tim_cp
  tau.idx = results$tau.idx
  res_mean = results$res_mean # scaling
  N = ncol(res_mean)
  p = nrow(res_mean)
  
  fit1 <- results$fit.vals
  var.resd1 <- results$var.resd
  
  # select optimal parameters
  opt_d = results$d
  opt_tau = results$tau
  opt_param = results$param
  
  # fitted values
  alpha0 = opt_param[1:p]
  alpha1 = opt_param[p + 1:p]
  m = results$m
  
  # weights
  wt <- matrix(nrow = p, ncol = N)
  for (k in 1:p) {
    N.k <- length(na.omit(res[k, ]))
    wt[k, 1:N.k] <-
      get.wt(alpha = c(alpha0[k], alpha1[k]),
             tim_cp = na.omit(tim_cp[k, ]),
             tau.idx = na.omit(tau.idx[k, ]),
             N = N.k)
  }
  
  # update variables if using original or first difference
  fit.vals <- mu.2 <- diff_p <- matrix(nrow = nrow(res), ncol = ncol(res))
  for (k in 1:p) {
    N.k <- length(na.omit(res[k, ]))
    x.2 = res_mean[k, ]
    d = opt_d
    diff_p[k, 1:N.k] = t(diffseries_keepmean(t(na.omit(wt[k, ]*(x.2-m[k]))), d))
    
    mu.2[k, ] = wt[k, ]*(x.2-m[k]) - diff_p[k, ]
    fit.vals[k, ] = (mu.2[k, ] + wt[k, ]*m[k]) +(1-wt[k, ])*fit1[k, ]
    if (use_scale) {
      fit.vals[k, ] = fit.vals[k, ]*(attributes(res_mean)$'scaled:scale'[k])
      var.resd1[k, ] = var.resd1[k, ]*attributes(res_mean)$'scaled:scale'[k]^2
    }
    
  }
  
  # plotting
  if (return_plot){
    plot(tim[1, ], res[1, ],
         ylim = c(min(c(res[1, ], fit.vals[1, ])), max(c(res[1, ], fit.vals[1, ]))), type = 'l',
         main = paste('Values fitted with d: ', round(opt_d,3), ' tau: ', round(mean(opt_tau),3)),
         xlab = 'Time (hour)', ylab = 'Resistance (ohm)')
    
    if (is.na(opt_tau[1])){
      lines(tim[1, ], fit.vals[1, ], col = "blue", lwd = 1)
    }else{
      opt_tau.idx = which(tim == opt_tau[1])
      lines(tim[1, 1:opt_tau.idx],
            fit.vals[1, 1:opt_tau.idx], col = "blue", lwd = 1)
      lines(tim[1, (opt_tau.idx):ncol(fit.vals)],
            fit.vals[1, (opt_tau.idx):ncol(fit.vals)], col = "green", lwd = 1)
      abline(v = opt_tau[1], lty = 2, col = "red")
    }
    abline(v = tau.range, lty = 1, col = "red")
  }
  
  return(list(fit.vals = fit.vals,
              var.resd1 = var.resd1,
              wt = wt, scaling = attributes(res_mean)$'scaled:scale'))
}

# parametric bootstrap using outputs from t2cd_step and plot.t2cd_step
bootstrap_sample_sigmoid_mv = function(results, plot_results, seed = 0){
  
  set.seed(seed)
  res = results$res
  tim = results$tim
  N <- ncol(res)
  p <- nrow(res)
  
  opt_d = results$d
  
  samp <- matrix(nrow = p, ncol = N)
  
  for (k in 1:p) {
    # regime 1
    fit.vals1 = plot_results$fit.vals[k, ]
    var.resd1 = plot_results$var.resd1[k, ]
    noise1 = rnorm(N, 0, sqrt(var.resd1))
    
    # regime 2
    
    m = results$m[k]
    sd = results$sd[k]
    df = results$df[k]
    wt = plot_results$wt[k, ]
    sim = sim_fi(N, opt_d, mu = m, sig = sd, df = df)
    seq_fi = sim$s
    
    samp[k, ] = c((1-wt)*(fit.vals1 + noise1) + wt*(seq_fi)*plot_results$scaling[k])
  }
  
  return(list(res=samp, tim=tim))
}