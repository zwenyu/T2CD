# fit first regime with model that accounts to heteroscedastic noise
# fit second regime with ARFIMA as in cpsearch2.R
# uses loglikelihood of original series

require(arfima)
require(fracdiff)
require(splines)
require(ggplot2)
require(mvtnorm)
source('./method/bsWLS.R')

# T2CD-step method
# wrapper around search_dtau_step to fit with both ranges of d
# and pick the fit with highest likelihood
t2cd_step = function(dat, t.max = 72, tau.range = c(10, 50), deg = 3,
                     t_resd = FALSE,
                     seqby = 1, resd.seqby = 5, 
                     use_arf = TRUE, use_scale = TRUE){
  # dat: input time series
  # t.max: cutoff for time considered
  # tau.range candidate change point range
  # deg: degree for B-spline  
  # segby, resd.seqby: interval between knots
  # use_arf: if true, use arfima estimates from arfima package
  # use_scale: if true, scale time series
  res1 = search_dtau_step(dat, t.max, tau.range, deg,
                          t_resd = t_resd,
                          dflag = 'original', seqby = seqby, resd.seqby = resd.seqby,
                          use_arf = use_arf, use_scale = use_scale)
  return(res1)
}

### helper and plotting functions

# grid search for d and tau
# return the likelihood for each combination
search_dtau_step = function(dat, t.max = 72, tau.range = c(10, 50), deg = 3,
                            t_resd = FALSE,
                            dflag = 'fdiff', seqby = 1, resd.seqby = 5,
                            use_arf = TRUE, use_scale = TRUE){
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
  if (tau.range[1] == tau.range[2]){
    diff = tim - tau.range[1]
    tau.idx = which(diff > 0)[1]
  }else{
    tau.idx = which(tim >= tau.range[1] & tim <= tau.range[2])    
  }
  M = c()
  d = c()
  m = c()
  t_df = c()
  t_scale = c()
  
  tol = 0.01
  if (is.null(dflag)){
    # determine range of d to search
    tau_last = tau.idx[length(tau.idx)]
    x.last = res_mean[(tau_last+1):N]
    arf = arfima::arfima(x.last)
    if ((0.5-arf$modes[[1]]$dfrac) < tol){
      dflag = 'fdiff' # likelihood evaluated on first differences
    }else{
      dflag = 'original' # likelihood evaluated on original series
    }
  }
  
  if (t_resd){
    # fit on post-50hr portion of sequence with normal errors for initialization
    res_init = search_dtau_step(dat, t.max, c(50, 50), deg,
                                t_resd = FALSE,
                                dflag = dflag, seqby = seqby, resd.seqby = resd.seqby,
                                use_arf = use_arf, use_scale = use_scale)  
    fit_init = plot.t2cd_step(res_init, tau.range = c(50, 50), deg = deg, 
                              use_arf = use_arf, use_scale = use_scale, return_plot = FALSE)
    r2 = res_mean[(res_init$idx+1):length(res_mean)] - fit_init$fit.vals2
    m_init = res_init$m
    d_init = res_init$d
    t_scale_init = sqrt(mean(r2^2))
    t_df_init = 3
  }  
  
  # iterate through each tau, return log-likelihood  
  for (j in 1:length(tau.idx)){
    tau_j = tau.idx[j]
    
    # optimize for polynomial component
    x.1 = res_mean[1:tau_j]
    t.1 = tim[1:tau_j]
    n.1 = length(x.1)
    
    fit = function(){
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
      if (t_resd){ # fit residual/t_scale using t-distribution with t_df
        negloglik = function(param){
          m = param[1]
          dfrac = param[2]
          t_df = param[4]
          t_scale = param[3]/sqrt(t_df/(t_df - 2))
          
          diff_p = c(diffseries_keepmean(matrix(x.2-m, ncol = 1), dfrac))
          logL = sum(log(gamma((t_df+1)/2)) - 0.5*log(t_df*pi) - log(gamma(t_df/2)) 
                     - 0.5*(t_df+1)*log(1+((diff_p/t_scale)^2)/t_df) - log(t_scale))
          neglogL = -logL
          
          return(neglogL)
        }
        loglik = function(param){
          m = param[1]
          dfrac = param[2]
          t_df = param[4]
          t_scale = param[3]/sqrt(t_df/(t_df - 2))
          
          diff_p = c(diffseries_keepmean(matrix(x.2-m, ncol = 1), dfrac))
          logL = sum(log(gamma((t_df+1)/2)) - 0.5*log(t_df*pi) - log(gamma(t_df/2)) 
                     - 0.5*(t_df+1)*log(1+((diff_p/t_scale)^2)/t_df) - log(t_scale)) 
          
          return(logL)
        }
                
      }else{ # Gaussian distribution
        negloglik = function(param){
          m = param[1]
          dfrac = param[2]
          
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
      }
      
      # optimizing
      if (use_arf){
        if (t_resd){
          stop('use_arf = TRUE option not implemented for t_resd = TRUE')
        }
          
        if (dflag == 'original'){
          arf = arfima::arfima(x.2)
          fit_d = arf$modes[[1]]$dfrac
        }else{
          arf = arfima::arfima(x.2)
          fit_d = arf$modes[[1]]$dfrac + 1
        }
        fit_m = arf$models[[1]]$muHat
        n.2 = length(x.2)
        ll.2 = arf$modes[[1]]$loglik - (n.2/2)*log(2*pi) - n.2/2
      }else{
        if (t_resd){
          eps = 1e-5
          optim.2 = optim(par = c(m_init, d_init, t_scale_init, t_df_init),
                          fn = negloglik, method = "L-BFGS-B",
                          lower = c(-Inf, -Inf, eps, 2+eps), 
                          upper = c(Inf, Inf, Inf, 300))
        }else{
          optim.2 = optim(par = c(mean(x.2), 0),
                          fn = negloglik, method = "BFGS")          
        }
        if (dflag == 'original'){
          fit_d = optim.2$par[2]
        }else{
          fit_d = optim.2$par[2] + 1
        }
        fit_m = optim.2$par[1]
        ll.2 = loglik(optim.2$par)[1]
        if (t_resd){
          fit_t_scale = optim.2$par[3]
          fit_t_df = optim.2$par[4]
        }else{
          fit_t_scale = -1
          fit_t_df = -1
        }
      }
      return(list(fit_M=ll.1 + ll.2, fit_d=fit_d, fit_m=fit_m,
                  fit_t_scale=fit_t_scale, fit_t_df = fit_t_df))
    }
    
    fit_res = tryCatch(fit(),
                    error = function(e) {return(NA)})
    if (any(is.na(fit_res))){
      M = c(M, -Inf)
      d = c(d, NA)
      m = c(m, NA)
      t_scale = c(m, NA)
      t_df = c(m, NA)
    }else{
      M = c(M, fit_res$fit_M)
      d = c(d, fit_res$fit_d)
      m = c(m, fit_res$fit_m)      
      t_scale = c(t_scale, fit_res$fit_t_scale)
      t_df = c(t_df, fit_res$fit_t_df)
    }
  }
  
  # tau and d at maximum log-likelihood
  M_df = data.frame(tau = tim[tau.idx], M = M, d = d, m = m)
  max.idx = which.max(M)
  max.tau = tim[tau.idx[max.idx]]
  max.d = d[max.idx]
  max.m = m[max.idx]
  max.t_df = t_df[max.idx]
  max.t_scale = t_scale[max.idx]/sqrt(max.t_df/(max.t_df - 2))
  
  return(list(M_df = M_df, res = res, tim = tim, tau.idx = tau.idx,
              tau = max.tau, d = max.d, m = max.m, t_scale = max.t_scale, t_df = max.t_df,
              idx = tau.idx[max.idx], logL = M[max.idx],
              dflag = dflag))
}

# plot sequences and fitted lines
plot.t2cd_step = function(results, tau.range = c(10, 50), deg = 3, 
                           use_arf = TRUE, use_scale = TRUE, return_plot = TRUE){
  res = results$res
  tim = results$tim
  tau.idx = results$tau.idx
  dflag = results$dflag
  if (use_scale){
    res_mean = scale(res, center = F) # scaling  
  }else{
    res_mean = matrix(res, ncol = 1)
  }
  N = length(res)
  
  # select optimal parameters
  opt_d = results$d
  opt_tau = results$tau
  opt_idx = results$idx
  
  ### fitted values for first regime
  fit1 = refitWLS(tim[1:opt_idx], res_mean[1:opt_idx], deg = deg)
  var.resd1.1 = fit1$var.resd
  mu = c(fit1$fit.vals, rep(NA, N-opt_idx))
  
  # fitted values for second regime
  if (use_arf){
    if (dflag == 'original'){
      x.2 = res_mean[(opt_idx+1):N]
      arf = arfima::arfima(x.2)
      fit = arf$modes[[1]]$muHat + arf$modes[[1]]$fitted
      mu[(opt_idx+1):N] = fit
    }else{
      x.2 = res_mean[(opt_idx):N]
      n.2 = length(x.2) - 1
      arf = arfima::arfima(diff(x.2, 1))
      diff.2 = arf$modes[[1]]$muHat + arf$modes[[1]]$fitted
      mu.2 = c()
      for (i in 1:n.2){
        mu.2 = c(mu.2, x.2[i] + diff.2[i])
      }
      mu[(opt_idx+1):N] = mu.2
    }
  }else{
    m = results$m
    if (dflag == 'original'){
      x.2 = res_mean[(opt_idx+1):N]
      diff_p = c(diffseries_keepmean(matrix(x.2-m, ncol = 1), opt_d))
    }else{
      x.2 = c(0, diff(res_mean[(opt_idx+1):N]-m, 1))
      diff_p = c(diffseries_keepmean(matrix(x.2, ncol = 1), opt_d - 1))
    }
    
    if (dflag == 'original'){
      mu.2 = (x.2 - m) - diff_p
      mu[(opt_idx+1):N] = mu.2 + m
    }else{
      diff.2 = x.2 - diff_p
      mu.2 = c()
      for (i in (opt_idx+1):N){
        mu.2 = c(mu.2, res_mean[i] + diff.2[i-opt_idx])
      }
      mu[(opt_idx+1):N] = mu.2
    }
  }
  
  if (use_scale){
    fit.vals = mu*attributes(res_mean)$'scaled:scale'  
    var.resd1 = fit1$var.resd*attributes(res_mean)$'scaled:scale'^2
  }else{
    fit.vals = mu
    var.resd1 = fit1$var.resd
  }
  
  # plotting
  if (return_plot){
    plot(tim, res, ylim = c(min(c(res, fit.vals)), max(c(res, fit.vals))), type = 'l', 
         main = paste('Values fitted with d: ', round(opt_d,3), ' tau: ', round(opt_tau,3)),
         xlab = 'Time (hour)', ylab = 'Resistance (ohm)')
    if (is.na(opt_tau)){
      lines(tim, fit.vals, col = "blue", lwd = 1)    
    }else{
      opt_tau.idx = which(tim == opt_tau)
      lines(tim[1:opt_idx], fit.vals[1:opt_idx], col = "blue", lwd = 1)  
      lines(tim[(opt_idx+1):N], fit.vals[(opt_idx+1):N], col = "green", lwd = 1)  
      abline(v = opt_tau, lty = 2, col = "red")
    }
    abline(v = tau.range, lty = 1, col = "red")    
  }
  
  return(list(fit.vals1 = mu[1:opt_idx], fit.vals2 = mu[(opt_idx+1):N],
              opt_idx = opt_idx, N = N, scaling = attributes(res_mean)$'scaled:scale',
              var.resd1 = fit1$var.resd))
}

# parametric bootstrap using outputs from t2cd_step and plot.t2cd_step
bootstrap_sample_step = function(results, plot_results, seed = 0, t_resd = FALSE){
  
  set.seed(seed)
  res = results$res
  tim = results$tim
  N = length(res)
  opt_idx = results$idx
  
  # regime 1
  fit.vals1 = plot_results$fit.vals1*plot_results$scaling
  var.resd1 = plot_results$var.resd1*plot_results$scaling^2
  noise1 = rnorm(opt_idx, 0, sqrt(var.resd1))
  
  # regime 2
  opt_d = results$d
  m = results$m
  res_mean = scale(res, center = F) # scaling
  diff_p = c(diffseries_keepmean(matrix(res_mean[(opt_idx+1):N]-m, ncol = 1), opt_d))
  if (t_resd){
    sim = sim.fi(N-opt_idx, opt_d, t_scale=results$t_scale, t_df=results$t_df, t_resd=TRUE)   
  }else{
    sd.resd2 = sqrt(mean(diff_p^2))
    sim = sim.fi(N-opt_idx, opt_d, sd.resd2, t_resd=FALSE)      
  }
  seq_fi = sim$s
  
  samp = c(fit.vals1 + noise1, (seq_fi + m)*plot_results$scaling)
  
  return(list(res=matrix(samp, nrow=1), tim=tim))  
}
