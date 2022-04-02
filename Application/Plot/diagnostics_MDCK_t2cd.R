#! /usr/bin/Rscript

# running univariate methods of MDCK data
# set working directory as T2CD

library(parallel)
source('./helperfunction/helperfn.R')
source('./method/t2cd_step.R')
source('./method/t2cd_sigmoid.R')

### load data
load("Cornell_MDCK_M_hominis_runs.RData")

### setup
dat.info = list(gel=gel, inf=inf, nor=nor, null=null, wou=wou)

# Calculate the number of cores
no_cores = detectCores() - 1
# Initiate cluster
cl = makeCluster(no_cores, type = 'FORK')
# each iteration of function output as a list entry
# each function output is a named list

# hypothesis testing by parametric bootstrap
boot_test = function(K, results, plot_results, t_resd=FALSE){
  
  test2_stat = c()
  for (k in 1:K){
    samp = bootstrap_sample_step(results, plot_results, seed = k, t_resd = t_resd)
    res_k = t2cd_step(samp, tau.range = c(50, 50), t_resd = t_resd, use_arf = F)
    res_mean = scale(res_k$res, center = F) 
    fit_step = plot.t2cd_step(res_k, tau.range = c(50, 50), use_arf = F, return_plot = FALSE)
    r2 = res_mean[(res_k$idx+1):length(res_mean)] - fit_step$fit.vals2
    test2_k = ks.test(r2/res_k$t_scale, "pt", res_k$t_df) # two-sided, exact
    test2_k_stat = test2_k$statistic
    test2_stat = c(test2_stat, test2_k_stat)
  }
  
  return(test2_stat)
}

# iterate through experiments, frequency, gel and inf/nor
ecisfreq = function(f){
  
  mat_step = matrix(NA, 0, 10)
  colnames(mat_step) = c('statistic1', 'p.value1', 'statistic2', 'p.value2', 
                         'expt', 'freq', 'gel', 'inf', 'm', 'time')

  for (x in 1:4){
    print(x)
    for (g in 0:1){
      for (i in 0:1){
        dat = extractdata(dat.com, dat.info, expt = x, freq = f, gel = g, inf = i, nor = (1-i))
        M = nrow(dat$res)
        for (m in 1:M){
          dat_m = list(res=dat$res[m,], tim=dat$tim[m,])
          
          # step method
          sink('aux')
          s = proc.time()
          res_step = t2cd_step(dat_m, tau.range = c(50, 50), t_resd = T, use_arf = F)
          ptime = proc.time() - s
          sink(NULL)
          res_mean = scale(res_step$res, center = F) 
          fit_step = plot.t2cd_step(res_step, tau.range = c(50, 50), use_arf = F, return_plot = FALSE)
          r1 = (res_mean[1:res_step$idx] - fit_step$fit.vals1)/sqrt(fit_step$var.resd1)
          test1 = shapiro.test(r1)
          r2 = res_mean[(res_step$idx+1):length(res_mean)] - fit_step$fit.vals2
          test2 = ks.test(r2/res_step$t_scale, "pt", res_step$t_df) # two-sided, exact
          test2_stat = test2$statistic
          
          test2_stat_boot = boot_test(100, res_step, fit_step, t_resd = T)
          test2_pval = mean(test2_stat_boot > test2_stat)
          
          mat_step = rbind(mat_step, c(test1$statistic, test1$p.value, test2_stat, test2_pval,
                                       res_step$tau, x, f, g, i, m, ptime[1]))
        }
      }
    }
  }
  
  return(list(mat_step = mat_step))
}

cptresults = parLapply(cl, 1, ecisfreq)

save.image('./Application/Univariate/diagnostics_tau50_tdist_mdck.RData')
stopCluster(cl)


### save to pdf
# library(e1071)
# 
# ecisfreq = function(f){
#   
#   for (x in 1:4){
#     print(x)
#     for (g in 0:1){
#       for (i in 0:1){
#         dat = extractdata(dat.com, dat.info, expt = x, freq = f, gel = g, inf = i, nor = (1-i))
#         M = nrow(dat$res)
#         for (m in 1:M){
#           dat_m = list(res=dat$res[m,], tim=dat$tim[m,])
#           
#           # step method
#           res_step = t2cd_step(dat_m, tau.range = c(50, 50), t_resd = T, use_arf = F)
#           res_mean = scale(res_step$res, center = F)
#           fit_step = plot.t2cd_step(res_step, tau.range = c(50, 50), use_arf = F, return_plot = FALSE)
#           r2 = res_mean[(res_step$idx+1):length(res_mean)] - fit_step$fit.vals2
#           t_df = 7
#           test2 = ks.test(r2, "pt", t_df) # two-sided, exact
#           var_r2 = var(r2)
#           ex_kurtosis = kurtosis(r2)
#           title = paste('Expt ', x, ', Gel ', g, ', Inf ', i, ', Series ', m, sep = '')
#           value = paste('Var = ', var_r2, ', Ex. Kurtosis = ', ex_kurtosis, sep = '')
#           test_value = paste('Stat = ', test2$statistic, ', P-value = ', test2$p.value, sep = '')
#           hist(r2, main = paste(title, value, test_value, sep='\n'))
#         }
#       }
#     }
#   }
#   
#   return(list(mat_step = mat_step))
# }
# 
# pdf(file='./Application/Univariate/diagnostics_tau50_tdist_mdck_residuals.pdf',
#     width=8, height=4) 
# res = ecisfreq(1)
# dev.off()
# 
