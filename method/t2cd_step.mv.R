require(arfima)
require(fracdiff)
require(splines)
require(ggplot2)
require(mvtnorm)
source('./method/bsWLS.R')
source('./method/t2cd_step.R')
source('./helperfunction/logL_fima.R')

# multivariate implementation for T2CD-step
# grid search for d and tau
t2cd_step.mv = function(dat, t.max = 72, tau.range = c(10, 50), deg = 3,
                        seqby = 1, resd.seqby = 5, 
                        use_arf = TRUE, use_scale = TRUE,
                        maxiter=10, tol=1e-6){
  # dat: input time series
  # t.max: cutoff for time considered
  # tau.range candidate change point range
  # deg: degree for B-spline  
  # seqby, resd.seqby: interval between knots  
  # use_arf: if true, use arfima estimates from arfima package
  # use_scale: if true, scale time series
  # maxiter: maximum number of iterations solving for tau and FI parameters
  # tol: stops iterating if difference between the most recent objective values is below tol
  
  # select data below t.max
  if (is.na(t.max)){
    t.max = min(apply(dat$tim, 1, max, na.rm = T))
  }
  
  tim.ind = !is.na(dat$tim) & dat$tim <= t.max
  t.maxidx = which(apply(tim.ind, 2, all) == T)
  
  res = matrix(dat$res[, t.maxidx], nrow = nrow(dat$res))
  tim = matrix(dat$tim[, t.maxidx], nrow = nrow(dat$tim))
  N = ncol(res)
  p = nrow(res)
  
  # points in change range
  tau.idx = which(apply(tim >= tau.range[1] & tim <= tau.range[2], 2, all))
  
  # alternate between optimizing for tau and FI parameters
  iter_k = 1
  iter_flag = TRUE
  neglogL_prev = Inf
  init.d = init.m = fix.d = fix.m = init.tau = init.idx = tau = idx = rep(NA, p)
  hist.d = hist.neglogL = c()
  
  while (iter_flag){
    # retain all data in regime 2
    reg2 = matrix(NA, nrow = p, ncol = N)
    for (k in 1:p){
      if (iter_k==1){
        res_k = t2cd_step(list(res=res[k,], tim=tim[k,]), t.max, tau.range, deg,
                          seqby = seqby, resd.seqby = resd.seqby,
                          use_arf = use_arf, use_scale = use_scale) 
        init.d[k] = res_k$d
        init.m[k] = res_k$m
        init.tau[k] = res_k$tau
        init.idx[k] = res_k$idx        
      }else{
        res_k = search_tau_step(list(res=res[k,], tim=tim[k,]), t.max, tau.range, deg,
                          seqby = seqby, resd.seqby = resd.seqby,
                          use_arf = use_arf, use_scale = use_scale,
                          fix.d = d_current, fix.m = m_current[k])         
        fix.d[k] = res_k$d
        fix.m[k] = res_k$m
        tau[k] = res_k$tau
        idx[k] = res_k$idx        
      }
      reg2[k, (res_k$idx+1):N] = res[k, (res_k$idx+1):N]
    }
    
    # preprocessing
    dflag = 'original'
    if (use_scale){
      reg2 = t(scale(t(reg2), center = F))
    }
    x.2 = reg2
    if (iter_k==1){
      init.d.2 = init.d  
      init.m.2 = init.m
    }else{
      init.d.2 = fix.d
      init.m.2 = fix.m
    }
    
    # helper functions for loglikelihood
    negloglik = function(param){
      m = param[1:p]
      dfrac = param[p+1]
      
      n.2 = rowSums(!is.na(x.2))
      diff_p = t(diffseries_keepmean(t(x.2-m), dfrac))
      
      diff_p[diff_p==0] = NA
      neglogL = sum(n.2*log(rowMeans(diff_p^2, na.rm = T)))
      
      return(neglogL)
    }
    
    # optimizing
    optim_params = optim(par = c(init.m.2, mean(init.d.2)),
                         fn = negloglik, method = "BFGS")
    d_current = optim_params$par[p + 1]
    m_current = optim_params$par[1:p]
    neglogL_current = negloglik(optim_params$par)
    
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
    
  }
  
  return(list(res = res, tim = tim, tau.idx = tau.idx,
              tau = tau, idx = idx, d = d_current, 
              univ_tau = init.tau, univ_idx = init.idx, univ_d = init.d, 
              hist.d = hist.d, hist.neglogL = hist.neglogL,
              iter_k = iter_k, dflag = dflag))
}


# grid search for tau given FI parameters
# return the likelihood for each combination
search_tau_step = function(dat, t.max = 72, tau.range = c(10, 50), deg = 3,
                           dflag = 'original', seqby = 1, resd.seqby = 5,
                           use_arf = FALSE, use_scale = TRUE,
                           fix.d = NULL, fix.m = NULL){
  stopifnot(use_arf==FALSE)
  
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
  M = c()
  d = c()
  m = c()
  
  # iterate through each tau, return log-likelihood  
  for (j in 1:length(tau.idx)){
    tau_j = tau.idx[j]
    
    # optimize for polynomial component
    x.1 = res_mean[1:tau_j]
    t.1 = tim[1:tau_j]
    n.1 = length(x.1)
    
    fit1 = refitWLS(t.1, x.1, deg = deg, seqby = seqby, resd.seqby = resd.seqby)
    resd1.1 = x.1 - fit1$fit.vals
    var.resd1.1 = fit1$var.resd
    ll.1 = sum(dnorm(resd1.1, log = TRUE, sd = sqrt(var.resd1.1)))
    
    # optimize for ARFIMA
    if (dflag == 'original'){
      x.2 = res_mean[(tau_j+1):N]
    }else{
      x.2 = diff(res_mean[(tau_j):N], 1)
    }
    # helper functions for loglikelihood
    negloglik = function(param){
      m = param[1]
      dfrac = param[2]
      if (dfrac<=-0.5 | dfrac>=1.5){
        return(1e+10)
      }
      
      diff_p = c(diffseries_keepmean(matrix(x.2-m, ncol = 1), dfrac))
      
      neglogL = sum(diff_p^2)
      
      return(neglogL)
    }
    loglik = function(param){
      m = param[1]
      dfrac = param[2]
      
      diff_p = c(diffseries_keepmean(matrix(x.2-m, ncol = 1), dfrac))
      
      logL = sum(dnorm(diff_p, log = TRUE, sd = sd(diff_p)))
      
      return(logL)
    }
    
    # optimizing
    optim.2 = optim(par = c(fix.m, fix.d),
                    fn = negloglik, method = "BFGS")
    if (dflag == 'original'){
      d = c(d, optim.2$par[2])
    }else{
      d = c(d, optim.2$par[2] + 1)
    }
    m = c(m, optim.2$par[1])
    ll.2 = loglik(optim.2$par)

    M = c(M, ll.1 + ll.2)
  }
  
  # tau and d at maximum log-likelihood
  M_df = data.frame(tau = tim[tau.idx], M = M, d = d, m = m)
  max.idx = which.max(M)
  max.tau = tim[tau.idx[max.idx]]
  max.d = d[max.idx]
  max.m = m[max.idx]
  
  return(list(M_df = M_df, res = res, tim = tim, tau.idx = tau.idx,
              tau = max.tau, d = max.d, m = max.m,
              idx = tau.idx[max.idx], logL = M[max.idx],
              dflag = dflag))
  
}

