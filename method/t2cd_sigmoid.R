require(arfima)
require(fracdiff)
require(splines)
require(ggplot2)
require(mvtnorm)
source('./method/bsWLS.R')
source('./helperfunction/helperfn.R')

# T2CD-sigmoid method
# wrapper around search_dtau_sigmoid to fit with both ranges of d
# and pick the fit with highest likelihood
t2cd_sigmoid = function(dat, t.max = 72, tau.range = c(10, 50), 
                        init.tau = c(15, 30, 45), deg = 3, C = 1000,
                        seqby = 1, resd.seqby = 5, use_scale = TRUE){
  # dat: input time series
  # t.max: cutoff for time considered
  # tau.range candidate change point range
  # init.tau: candidate taus to initialize learning
  # deg: degree for B-spline  
  # C: regularization coefficient
  # segby, resd.seqby: interval between knots
  # use_scale: if true, scale time series  
  res1 = search_dtau_sigmoid(dat, t.max, tau.range, init.tau, deg, C, dflag = 'original',
                         seqby = seqby, resd.seqby = resd.seqby, use_scale = use_scale)
  res2 = search_dtau_sigmoid(dat, t.max, tau.range, init.tau, deg, C, dflag = 'fdiff',
                         seqby = seqby, resd.seqby = resd.seqby, use_scale = use_scale)
  if (res1$logL > res2$logL){
    return(res1)
  }else{
    return(res2)
  }
}

### helper and plotting functions

# sigmoid
sigmoid = function(a){return(1/(1+exp(-a)))}
softmax = function(a,b){
  return(exp(a)/(exp(a) + exp(b)))
}

# optimize the likelihood for d and tau
# option to initialize tau at multiple indices
search_dtau_sigmoid = function(dat, t.max = 72, tau.range = c(10, 50), 
                           init.tau = c(15, 30, 45), deg = 3, C = 1000, dflag = 'fdiff', 
                           seqby = 1, resd.seqby = 5, use_scale = T){
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
  p = nrow(res_mean)
  
  # points in change range
  tau.idx = which(apply(tim >= tau.range[1] & tim <= tau.range[2], 2, all))
  tim_cp = matrix(tim[,tau.idx], nrow = nrow(tim))
  
  tol = 0.01
  # determine range of d to search
  if (is.null(dflag)){
    tau_last = tau.idx[length(tau.idx)]
    dfrac_p = rep(0, p)
    for (k in 1:p){
      x.last = res_mean[k,(tau_last+1):N]
      arf = arfima::arfima(x.last)
      dfrac_p[k] = arf$modes[[1]]$dfrac
    }
    if ((0.5-mean(dfrac_p)) < tol){
      dflag = 'fdiff' # likelihood evaluated on first differences
    }else{
      dflag = 'original' # likelihood evaluated on original series
    }
  }

  # optimize for spline component
  ll.1.mat = matrix(ncol = N, nrow = 0)
  for (k in 1:p){
    fit1 = refitWLS(tim[k,], res_mean[k,], deg = deg, 
                    seqby = seqby, resd.seqby = resd.seqby)
    resd1 = res_mean[k,] - fit1$fit.vals
    var.resd1 = fit1$var.resd
    ll.1.vec = dnorm(resd1, log = TRUE, sd = sqrt(var.resd1))  
    ll.1.mat = rbind(ll.1.mat, ll.1.vec)
  }
  
  if (dflag == 'original'){
    x.2 = res_mean
    lastN = N
  }else{
    x.2 = t(diff(t(res_mean), 1))
    lastN = N-1
  }
  # loglikelihood, penalty to enforce tau within tau.range
  negloglik_pen = function(param){
    alpha0 = param[1]
    alpha1 = param[2]
    m = param[3]  
    dfrac = param[4]
    if (dfrac<=-0.5 | dfrac>=0.5){
      return(1e+10)
    }

    # weights
    wt_cp = sigmoid(alpha0+alpha1*tim_cp) # 0 to 1
    wt = cbind(matrix(rep(wt_cp[,1], tau.idx[1]-1), p, byrow = F),
               wt_cp,
               matrix(rep(wt_cp[,ncol(wt_cp)], N-tau.idx[length(tau.idx)]), p, byrow = F))

    if (dflag == 'original'){
      x.2m = x.2-m
    }else{
      x.2m = cbind(rep(0, p), x.2-m)
    }
    diff_p = t(diffseries_keepmean(t(wt*(x.2m)), dfrac))

    neglogL = -sum((1-wt)*ll.1.mat) + 0.5*log(2*pi)*sum(wt) + 0.5*sum(wt) +
      sum(0.5*rowSums(wt)*log(rowSums(wt*diff_p^2)/rowSums(wt))) - 
      p*C*sum(wt_cp[,ncol(wt_cp)] - wt_cp[,1])      
    
    return(neglogL)
  }
  
  # loglikelihood
  loglik = function(param){
    alpha0 = param[1]
    alpha1 = param[2]
    m = param[3]  
    dfrac = param[4]
    if (dfrac<=-0.5 | dfrac>=0.5){
      return(1e+10)
    }
    
    # weights
    wt_cp = sigmoid(alpha0+alpha1*tim_cp) # 0 to 1
    wt = cbind(matrix(rep(wt_cp[,1], tau.idx[1]-1), p, byrow = F),
               wt_cp,
               matrix(rep(wt_cp[,ncol(wt_cp)], N-tau.idx[length(tau.idx)]), p, byrow = F))
    
    if (dflag == 'original'){
      x.2m = x.2-m
    }else{
      x.2m = cbind(rep(0, p), x.2-m)
    }
    diff_p = t(diffseries_keepmean(t(wt*(x.2m)), dfrac))
    
    logL = sum((1-wt)*ll.1.mat) - 0.5*log(2*pi)*sum(wt) - 0.5*sum(wt) -
      sum(0.5*rowSums(wt)*log(rowSums(wt*diff_p^2)/rowSums(wt)))
    
    return(logL)
  }
  
  opt_logL = -Inf
  for (tau_i in init.tau){
    tau_i.idx = which.min(apply(tim <= tau_i, 2, all) == T)
    optim_i = optim(par = c(-tau_i, 1, mean(x.2[tau_i.idx:lastN]), 0),
                    fn = negloglik_pen, method = "BFGS")
    logL_i = loglik(optim_i$par)
    if (logL_i > opt_logL){
      opt_logL = logL_i
      optim_params = optim_i
    } 
  }
  
  opt_param = optim_params$par
  alpha0 = opt_param[1]
  alpha1 = opt_param[2]
  if (dflag == 'original'){
    opt_d = opt_param[4]  
  }else{
    opt_d = opt_param[4] + 1
  }
  
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
  
  # range of change location
  opt_taurange1.idx = apply(wt, 1, function(x){return(which(x>=0.1, arr.ind = TRUE)[1])})
  opt_taurange2.idx = apply(wt, 1, function(x){return(which(x>=0.9, arr.ind = TRUE)[1]-1)})
  opt_taurange1 = c()
  opt_taurange2 = c()
  for (k in 1:p){
    opt_taurange1 = c(opt_taurange1, tim[k,opt_taurange1.idx[k]])
    opt_taurange2 = c(opt_taurange2, tim[k,opt_taurange2.idx[k]])
  }
  opt_taurange1[is.na(opt_taurange1)] = tau.range[1]
  opt_taurange2[is.na(opt_taurange2)] = tau.range[2]
  
  return(list(res = res, tim = tim, tau.idx = tau.idx,
              d = opt_d, tau = opt_tau, 
              tau.range1 = opt_taurange1, tau.range2 = opt_taurange2,
              param = opt_param, logL = opt_logL, dflag = dflag))
}

# plot sequences and fitted lines
plot.t2cd_sigmoid = function(results, tau.range = c(10, 50), deg = 3, 
                         seqby = 1, resd.seqby = 5){
  res = results$res
  tim = results$tim
  tau.idx = results$tau.idx
  dflag = results$dflag
  res_mean = t(scale(t(res), center = F)) # scaling
  N = ncol(res_mean)
  p = nrow(res_mean)
  
  fit1 = var.resd1 = matrix(nrow = 0, ncol = N)
  for (k in 1:p){
    fitwls = refitWLS(tim[k,], res_mean[k,], deg = deg,
                      seqby = seqby, resd.seqby = resd.seqby)  
    fit1 = rbind(fit1, fitwls$fit.vals)
    var.resd1 = rbind(var.resd1, fitwls$var.resd)
  }
  
  # points in change range
  tim_cp = matrix(tim[,tau.idx], nrow = nrow(tim))
  
  # select optimal parameters
  opt_d = results$d
  opt_tau = results$tau
  opt_param = results$param
  
  # fitted values
  alpha0 = opt_param[1]
  alpha1 = opt_param[2]
  m = opt_param[3]  
  
  # weights
  wt_cp = sigmoid(alpha0+alpha1*tim_cp) # 0 to 1
  wt = cbind(matrix(rep(wt_cp[,1], tau.idx[1]-1), p, byrow = F),
             wt_cp,
             matrix(rep(wt_cp[,ncol(wt_cp)], N-tau.idx[length(tau.idx)]), p, byrow = F))
  
  # update variables if using original or first difference
  if (dflag == 'original'){
    x.2 = res_mean
    d = opt_d
    diff_p = t(diffseries_keepmean(t(wt*(x.2-m)), d))
  }else{
    x.2 = cbind(0, t(diff(t(res_mean), 1))-m)
    d = opt_d - 1
    diff_p = t(diffseries_keepmean(t(wt*x.2), d))
  }
  
  if (dflag == 'original'){
    mu.2 = wt*(x.2-m) - diff_p
    fit.vals = (mu.2 + wt*m)*attributes(res_mean)$'scaled:scale' + 
      (1-wt)*fit1*attributes(res_mean)$'scaled:scale'
  }else{
    diff.2 = wt*x.2 - diff_p
    mu.2 = matrix(nrow = p, ncol = 0)
    for (i in 1:N){
      mu.2 = cbind(mu.2, res_mean[,i] + diff.2[,i])
    }
    fit.vals = (wt*mu.2)*attributes(res_mean)$'scaled:scale' +
      (1-wt)*fit1*attributes(res_mean)$'scaled:scale'
  }
  
  # plotting
  plot(tim[1,], res[1,], ylim = c(min(c(res, fit.vals)), max(c(res, fit.vals))), type = 'l', 
       main = paste('Values fitted with d: ', round(opt_d,3), ' tau: ', round(mean(opt_tau),3)),
       xlab = 'Time (hour)', ylab = 'Resistance (ohm)')
  if (p > 1){
    for (k in 2:p){
      lines(tim[k,], res[k,])
    }
  }
  if (is.na(opt_tau[1])){
    lines(tim[1,], fit.vals[1,], col = "blue", lwd = 1)    
  }else{
    opt_tau.idx = which(tim[1,] == opt_tau[1])
    lines(tim[1,1:opt_tau.idx], fit.vals[1,1:opt_tau.idx], col = "blue", lwd = 1)  
    lines(tim[1,(opt_tau.idx):ncol(fit.vals)], fit.vals[1,(opt_tau.idx):ncol(fit.vals)], col = "green", lwd = 1)  
    abline(v = opt_tau[1], lty = 2, col = "red")
  }
  if (p > 1){
    for (k in 2:p){
      if (is.na(opt_tau[k])){
        lines(tim[k,], fit.vals[k,], col = "blue", lwd = 1)    
      }else{
        opt_tau.idx = which(tim[k,] == opt_tau[k])
        lines(tim[k,1:opt_tau.idx], fit.vals[k,1:opt_tau.idx], col = "blue", lwd = 1)  
        lines(tim[k,(opt_tau.idx):ncol(fit.vals)], fit.vals[k,(opt_tau.idx):ncol(fit.vals)], col = "green", lwd = 1)  
        abline(v = opt_tau[k], lty = 2, col = "red")
      }
    }
  }
  abline(v = tau.range, lty = 1, col = "red")
  
  return(list(fit.vals = fit.vals, 
              var.resd1 = var.resd1*attributes(res_mean)$'scaled:scale'^2))
}
