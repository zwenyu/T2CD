require(arfima)
require(fracdiff)
require(splines)
require(ggplot2)
require(mvtnorm)
source('./method/bsWLS.R')
source('./method/t2cd_sigmoid.R')
source('./helperfunction/helperfn.R')

# multivariate implementation for T2CD-sigmoid
# optimize the likelihood for d and tau
# option to initialize tau at multiple indices
t2cd_sigmoid.mv = function(dat, t.max = 72, tau.range = c(10, 50), 
                           init.tau = c(15, 30, 45), deg = 3, C = 1000,
                           seqby = 1, resd.seqby = 5, use_scale = TRUE){
  # dat: input time series
  # t.max: cutoff for time considered
  # tau.range candidate change point range
  # init.tau: candidate taus to initialize learning
  # deg: degree for B-spline  
  # C: regularization coefficient
  # seqby, resd.seqby: interval between knots  
  # use_scale: if true, scale time series 
  
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
  N = ncol(res)
  p = nrow(res)
  
  # points in change range
  tau.idx = which(apply(tim >= tau.range[1] & tim <= tau.range[2], 2, all))
  tim_cp = matrix(tim[,tau.idx], nrow = nrow(tim))

  # optimize for spline component
  ll.1.mat = matrix(ncol = N, nrow = 0)
  for (k in 1:p){
    fit1 = refitWLS(tim[k,], res[k,], deg = deg)
    resd1 = res[k,] - fit1$fit.vals
    var.resd1 = fit1$var.resd
    ll.1.vec = dnorm(resd1, log = TRUE, sd = sqrt(var.resd1))  
    ll.1.mat = rbind(ll.1.mat, ll.1.vec)
  }
  
  # determine range of d to search and find initializing parameters
  init.param = matrix(NA, p, 4)
  init.d = rep(NA, p)
  for (k in 1:p){
    res_k = t2cd_sigmoid(list(res=matrix(res[k,], 1), tim=matrix(tim[k,], 1)), 
                         t.max, tau.range, init.tau, deg, C, 
                         seqby = seqby, resd.seqby = resd.seqby, use_scale = use_scale)
    init.param[k,] = res_k$par
    init.d[k] = res_k$d
  }

  dflag = 'original'
  x.2 = res_mean
  init.d.2 = init.d

  # loglikelihood, penalty to enforce tau within tau.range
  negloglik_partial_pen = function(param){
    m = param[1:p]  
    dfrac = param[p+1]
    if (dfrac<=-0.5 | dfrac>=1.5){
      return(1e+10)
    }
        
    # weights
    wt_cp = sigmoid(alpha0+alpha1*tim_cp) # 0 to 1
    wt = cbind(matrix(rep(wt_cp[,1], tau.idx[1]-1), p, byrow = F),
               wt_cp,
               matrix(rep(wt_cp[,ncol(wt_cp)], N-tau.idx[length(tau.idx)]), p, byrow = F))
    
    x.2m = x.2-m
    diff_p = t(diffseries_keepmean(t(wt*(x.2m)), dfrac))
    
    neglogL = -sum((1-wt)*ll.1.mat) + 0.5*log(2*pi)*sum(wt) + 0.5*sum(wt) +
      sum(0.5*rowSums(wt)*log(rowSums(wt*diff_p^2)/rowSums(wt))) - 
      p*C*sum(wt_cp[,ncol(wt_cp)] - wt_cp[,1])      
    
    return(neglogL)
  }
  
  # loglikelihood
  loglik = function(param){
    m = param[1:p]  
    dfrac = param[p+1]
    
    # weights
    wt_cp = sigmoid(alpha0+alpha1*tim_cp) # 0 to 1
    wt = cbind(matrix(rep(wt_cp[,1], tau.idx[1]-1), p, byrow = F),
               wt_cp,
               matrix(rep(wt_cp[,ncol(wt_cp)], N-tau.idx[length(tau.idx)]), p, byrow = F))
    
    x.2m = x.2-m
    diff_p = t(diffseries_keepmean(t(wt*(x.2m)), dfrac))
    
    logL = sum((1-wt)*ll.1.mat) - 0.5*log(2*pi)*sum(wt) - 0.5*sum(wt) -
      sum(0.5*rowSums(wt)*log(rowSums(wt*diff_p^2)/rowSums(wt)))
    
    return(logL)
  }
  
  alpha0 = init.param[,1]
  alpha1 = init.param[,2]
  m = init.param[,3]
  optim_params = optim(par = c(m, mean(init.d.2)),
                       fn = negloglik_partial_pen, method = "BFGS")
  
  opt_param = c(alpha0, alpha1, optim_params$par)
  opt_logL = loglik(optim_params$par)
  opt_d = opt_param[3*p+1]
  univ_d = init.param[,4]
  
  # weights
  wt_cp = sigmoid(alpha0+alpha1*tim_cp) # 0 to 1
  wt = cbind(matrix(rep(wt_cp[,1], tau.idx[1]-1), p, byrow = F),
             wt_cp,
             matrix(rep(wt_cp[,ncol(wt_cp)], N-tau.idx[length(tau.idx)]), p, byrow = F))
  opt_tau.idx = apply(wt, 1, function(x){return(which(x>=0.5, arr.ind = TRUE)[1]-1)})
  opt_tau = c()
  for (k in 1:p){
    opt_tau = c(opt_tau, tim[k,opt_tau.idx[k]])
  }
  
  return(list(res = res, tim = tim, tau.idx = tau.idx,
              d = opt_d, univ_d = univ_d, tau = opt_tau, param = opt_param, logL = opt_logL,
              dflag = dflag))
}

### helper functions

# sigmoid
sigmoid = function(a){return(1/(1+exp(-a)))}
softmax = function(a,b){
  return(exp(a)/(exp(a) + exp(b)))
}
