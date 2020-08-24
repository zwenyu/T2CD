library(ecp)
source('./helperfunction/logL_fima.R')
source('./helperfunction/helperfn.R')

### ECP method

# univariate version
# segment by edivisive, d found by maximizing log-likelihood on regime 2
ecp_d = function(dat, t.max = 72, tau.range = c(10, 50), dflag = 'diff'){
  # dat: input time series
  # t.max: cutoff for time considered
  # tau.range candidate change point range
  # dflag: if true, take first difference of data
  
  # select data below t.max
  if (is.na(t.max)){
    t.max = min(apply(dat$tim, 1, max, na.rm = T))
  }
  
  tim.ind = !is.na(dat$tim) & dat$tim <= t.max
  t.maxidx = which(apply(tim.ind, 2, all) == T)
  
  res = matrix(dat$res[, t.maxidx], nrow = nrow(dat$res))
  tim = matrix(dat$tim[, t.maxidx], nrow = nrow(dat$tim))
  
  tim.ind = !is.na(tim) & tim <= tau.range[2] & tim >= tau.range[1]
  t1 = min(which(apply(tim.ind, 2, all) == T))
  t2 = max(which(apply(tim.ind, 2, all) == T))
  p = nrow(tim.ind)
  
  # segment
  k = 5
  if (dflag == 'original'){
    x = t(res)
    ecp = e.divisive(x, k = k)
    idx = which(ecp$order.found[3:(k+2)] >= t1 & ecp$order.found[3:(k+2)] <= t2)
    if (length(idx) > 0){
      tau.idx = ecp$order.found[3:(k+2)][idx[1]] - 1
    }else{
      tau.idx = ecp$order.found[3] - 1
    }
  }else{
    x = diff(t(res))
    ecp = e.divisive(x, k = k)
    idx = which(ecp$order.found[3:(k+2)] >= t1 & ecp$order.found[3:(k+2)] <= t2)
    if (length(idx) > 0){
      tau.idx = ecp$order.found[3:(k+2)][idx[1]]
    }else{
      tau.idx = ecp$order.found[3]
    }
  }
  tau = tim[,tau.idx]
  
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
  
  # original
  reg2 = matrix(res[,(t.max+1):ncol(res)], nrow = nrow(res))
  x.2 = t(scale(t(reg2), center = F))
  optim_params_orig = optim(par = c(rowMeans(x.2, na.rm = T), 0),
                            fn = negloglik, method = "BFGS")
  # diff
  reg2 = matrix(res[,(t.max):ncol(res)], nrow = nrow(res))
  reg2 = t(scale(t(reg2), center = F))
  x.2 = matrix(reg2[, -1] - reg2[, -ncol(reg2)], nrow = p)
  optim_params_diff = optim(par = c(rowMeans(x.2, na.rm = T), 0),
                             fn = negloglik, method = "BFGS")
  if (optim_params_orig$value < optim_params_diff$value){
    d = optim_params_orig$par[p + 1]
  }else{
    d = optim_params_diff$par[p + 1] + 1  
  }
  
  return(list(idx = tau.idx, tau = tau, d = d))
}

# multivariate version
# segment by edivisive, d found by maximizing log-likelihood on regime 2
ecp_d.mv = function(dat, t.max = 72, tau.range = c(10, 50), dflag = 'diff'){
  # dat: input time series
  # t.max: cutoff for time considered
  # tau.range candidate change point range
  # dflag: if true, take first difference of data  
  
  # select data below t.max
  if (is.na(t.max)){
    t.max = min(apply(dat$tim, 1, max, na.rm = T))
  }
  
  tim.ind = !is.na(dat$tim) & dat$tim <= t.max
  t.maxidx = which(apply(tim.ind, 2, all) == T)
  
  res = matrix(dat$res[, t.maxidx], nrow = nrow(dat$res))
  tim = matrix(dat$tim[, t.maxidx], nrow = nrow(dat$tim))
  
  tim.ind = !is.na(tim) & tim <= tau.range[2] & tim >= tau.range[1]
  N = ncol(tim.ind)
  p = nrow(tim.ind)
  
  reg2 = matrix(NA, nrow = p, ncol = N)
  init.d = tau = idx = rep(NA, p)
  for (k in 1:p){
    res_k = ecp_d(list(res=matrix(res[k,], 1), tim=matrix(tim[k,], 1)), 
                  t.max, tau.range, dflag)
    reg2[k, (res_k$idx+1):N] = res[k, (res_k$idx+1):N]
    init.d[k] = res_k$d
    tau[k] = res_k$tau
    idx[k] = res_k$idx
  }
  
  if (mean(init.d) < 0.5){
    dflag = 'original'
    x.2 = t(scale(t(reg2), center = F))
    init.d.2 = init.d
  }else{
    dflag = 'diff'
    init.d.2 = init.d - 1    
    reg2 = t(scale(t(reg2), center = F))
    x.2 = reg2[, -1] - res[, -ncol(reg2)]/attributes(reg2)$`scaled:scale`
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
  optim_params = optim(par = c(rowMeans(x.2, na.rm = T), mean(init.d.2)),
                       fn = negloglik, method = "BFGS")
  if (dflag == 'original'){
    d = optim_params$par[p + 1]
  }else{
    d = optim_params$par[p + 1] + 1  
  }
  
  return(list(tau = tau, d = d, univ_d = init.d))
}

### FixedTau and TrueTau method

# univariate version
# given tau, d found by maximizing log-likelihood on regime 2
estd = function(dat, tau, t.max = 72){
  # dat: input time series
  # t.max: cutoff for time considered
  # tau.range candidate change point range
  
  # select data below t.max
  if (is.na(t.max)){
    t.max = min(apply(dat$tim, 1, max, na.rm = T))
  }
  
  tim.ind = !is.na(dat$tim) & dat$tim <= t.max
  t.maxidx = which(apply(tim.ind, 2, all) == T)
  
  res = matrix(dat$res[, t.maxidx], nrow = nrow(dat$res))
  tim = matrix(dat$tim[, t.maxidx], nrow = nrow(dat$tim))
  
  tim.ind = !is.na(dat$tim) & dat$tim <= tau
  t.max = max(which(apply(tim.ind, 2, all) == T))
  p = nrow(tim.ind)
  
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
  
  # original
  reg2 = matrix(res[,(t.max+1):ncol(res)], nrow = nrow(res))
  x.2 = t(scale(t(reg2), center = F))
  optim_params_orig = optim(par = c(rowMeans(x.2, na.rm = T), 0),
                            fn = negloglik, method = "BFGS")
  # diff
  reg2 = matrix(res[,(t.max):ncol(res)], nrow = nrow(res))
  reg2 = t(scale(t(reg2), center = F))
  x.2 = matrix(reg2[, -1] - reg2[, -ncol(reg2)], nrow = p)
  optim_params_diff = optim(par = c(rowMeans(x.2, na.rm = T), 0),
                             fn = negloglik, method = "BFGS")
  if (optim_params_orig$value < optim_params_diff$value){
    d = optim_params_orig$par[p + 1]
  }else{
    d = optim_params_diff$par[p + 1] + 1  
  }
  
  return(list(d = d))
}

# multivariate version
# dimensions can have different tau
# given tau, d found by maximizing log-likelihood on regime 2
estd.mv = function(dat, tau, t.max = 72){
  # dat: input time series
  # t.max: cutoff for time considered
  # tau.range candidate change point range
  
  # selelct data below t.max
  if (is.na(t.max)){
    t.max = min(apply(dat$tim, 1, max, na.rm = T))
  }
  
  tim.ind = !is.na(dat$tim) & dat$tim <= t.max
  t.maxidx = which(apply(tim.ind, 2, all) == T)
  
  res = matrix(dat$res[, t.maxidx], nrow = nrow(dat$res))
  tim = matrix(dat$tim[, t.maxidx], nrow = nrow(dat$tim))
  N = ncol(tim.ind)
  p = nrow(tim.ind)
  
  reg2 = matrix(NA, nrow = p, ncol = N)
  init.d = rep(NA, p)
  for (k in 1:p){
    res_k = estd(list(res=matrix(res[k,], 1), tim=matrix(tim[k,], 1)), 
                 tau[k], t.max)
    tau_k.ind = tim[k,] <= tau[k]
    tau_k.idx = which.min(tau_k.ind)
    reg2[k, tau_k.idx:N] = res[k, tau_k.idx:N]
    init.d[k] = res_k$d
  }
  
  if (mean(init.d) < 0.5){
    dflag = 'original'
    x.2 = t(scale(t(reg2), center = F))
    init.d.2 = init.d
  }else{
    dflag = 'diff'
    reg2 = t(scale(t(reg2), center = F))
    x.2 = reg2[, -1] - res[, -ncol(reg2)]/attributes(reg2)$`scaled:scale`
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
  optim_params = optim(par = c(rowMeans(x.2, na.rm = T), mean(init.d.2)),
                       fn = negloglik, method = "BFGS")
  if (dflag == 'original'){
    d = optim_params$par[p + 1]
  }else{
    d = optim_params$par[p + 1] + 1  
  }
  
  return(list(d = d, univ_d = init.d))
}





