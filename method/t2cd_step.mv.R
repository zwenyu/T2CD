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
                        use_arf = TRUE, use_scale = TRUE){
  # dat: input time series
  # t.max: cutoff for time considered
  # tau.range candidate change point range
  # deg: degree for B-spline  
  # seqby, resd.seqby: interval between knots  
  # use_arf: if true, use arfima estimates from arfima package
  # use_scale: if true, scale time series
  
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
  
  # retain all data in regime 2
  reg2 = matrix(NA, nrow = p, ncol = N)
  init.d = init.m = tau = idx = rep(NA, p)
  for (k in 1:p){
    res_k = t2cd_step(list(res=res[k,], tim=tim[k,]), t.max, tau.range, deg,
                      seqby = seqby, resd.seqby = resd.seqby,
                      use_arf = use_arf, use_scale = use_scale) 
    reg2[k, (res_k$idx+1):N] = res[k, (res_k$idx+1):N]
    init.d[k] = res_k$d
    init.m[k] = res_k$m
    tau[k] = res_k$tau
    idx[k] = res_k$idx
  }
  
  # preprocessing
  if (mean(init.d) < 0.5){
    dflag = 'original'
    if (use_scale){
      reg2 = t(scale(t(reg2), center = F))
    }
    x.2 = reg2
    init.d.2 = init.d
  }else{
    dflag = 'fdiff'
    if (use_scale){
      reg2 = t(scale(t(reg2), center = F))
      x.2 = reg2[, -1] - res[, -ncol(reg2)]/attributes(reg2)$`scaled:scale`
    }else{
      x.2 = reg2[, -1] - res[, -ncol(reg2)] 
    }
    init.d.2 = init.d - 1
  }
  
  # helper functions for loglikelihood
  negloglik = function(param){
    m = param[1:p]
    dfrac = param[p+1]
    if (dfrac<=-0.5 | dfrac>=0.5){
      return(1e+10)
    }
    
    n.2 = rowSums(!is.na(x.2))
    diff_p = t(diffseries_keepmean(t(x.2-m), dfrac))
    
    diff_p[diff_p==0] = NA
    neglogL = sum(n.2*log(rowMeans(diff_p^2, na.rm = T)))
    
    return(neglogL)
  }
  
  # optimizing
  if (use_arf){
    if (dflag == 'original'){
      opt = optimize(fima.ll.auto, interval = c(-0.5, 0.5), x = reg2, 
                     maximum = TRUE, tol = .Machine$double.eps)
    }else{
      opt = optimize(fima.ll.auto, interval = c(0.5, 1.5), x = reg2, 
                     maximum = TRUE, tol = .Machine$double.eps)
    }
    d = opt$maximum
  }else{
    optim_params = optim(par = c(rowMeans(x.2, na.rm = T), mean(init.d.2)),
                         fn = negloglik, method = "BFGS")
    if (dflag == 'original'){
      d = optim_params$par[p + 1]
    }else{
      d = optim_params$par[p + 1] + 1  
    }
  }
  
  return(list(res = res, tim = tim, tau.idx = tau.idx,
              tau = tau, d = d, univ_d = init.d, idx = idx, dflag = dflag))
}
