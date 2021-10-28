library(parallel)
source('./helperfunction/helperfn.R')
source('./method/t2cd_step.R')
source('./method/t2cd_sigmoid.R')

### MDCK
load("Cornell_MDCK_M_hominis_runs.RData")
dat.info = list(gel=gel, inf=inf, nor=nor, null=null, wou=wou)
num_boot = 500

# Calculate the number of cores
no_cores = 12
# Initiate cluster
cl = makeCluster(no_cores, type = 'FORK')
# each iteration of function output as a list entry
# each function output is a named list

# iterate through experiments, frequency, gel and inf/nor
boot_k = function(k, results, plot_results, methodname){
  
  if (methodname=='step'){
    samp = bootstrap_sample_step(results, plot_results, seed = k)
    res = t2cd_step(samp, use_arf = F)    
  }else if (methodname=='sigmoid'){
    samp = bootstrap_sample_sigmoid(results, plot_results, seed = k)
    res = t2cd_sigmoid(list(res=matrix(samp$res,1), tim=matrix(samp$tim,1)))
  }
  tau = res$tau
  d = res$d    

  return(list(tau=tau, d=d))
}

all_cptresults_step = list()
all_cptresults_sigmoid = list()
f = 1
count = 1
for (x in 1:4){
  print(x)
  for (g in 0:1){
    for (i in 0:1){
      dat = extractdata(dat.com, dat.info, expt = x, freq = f, gel = g, inf = i, nor = (1-i))
      M = nrow(dat$res)
      for (m in 1:M){
        dat_m = list(res=dat$res[m,], tim=dat$tim[m,])
        
        # first result is on original data, next num_boot on bootstrapped samples
        # step method
        res_step = t2cd_step(dat_m, use_arf = F)
        plot_results_step = plot.t2cd_step(res_step, use_arf = F, return_plot = FALSE)
        cptresults_step = parLapply(cl, 1:num_boot, boot_k, 
                                    res_step, plot_results_step, 'step')
        cptresults_step = c(list(list(tau=res_step$tau, d=res_step$d, expt=x, freq=f, gel=g, inf=i, m=m)), 
                            cptresults_step)
        
        # weighted method           
        res_sigmoid = t2cd_sigmoid(list(res=matrix(dat_m$res,1), tim=matrix(dat_m$tim,1)))
        plot_results_sigmoid = plot.t2cd_sigmoid(res_sigmoid, return_plot = FALSE)
        cptresults_sigmoid = parLapply(cl, 1:num_boot, boot_k,
                                    res_sigmoid, plot_results_sigmoid, 'sigmoid')
        cptresults_sigmoid = c(list(list(tau=res_sigmoid$tau, d=res_sigmoid$d, expt=x, freq=f, gel=g, inf=i, m=m)), 
                            cptresults_sigmoid)   
        
        # append to overall results list
        all_cptresults_step[[count]] = cptresults_step
        all_cptresults_sigmoid[[count]] = cptresults_sigmoid
        count = count + 1
      }
    }
  }
}

save.image('./Application/Univariate/mdck_boot.RData')
stopCluster(cl)

# saved image contains all_cptresults_step and all_cptresults_sigmoid
# which are lists of length = number of series in dataset
# each element is a list of length = 1 + num_boot
# where first result is fit on original data, and next num_boot are on bootstrapped samples

# for (res in all_cptresults_step){
#   x = res[[1]]$expt
#   g = res[[1]]$gel
#   i = res[[1]]$inf
#   print(paste('x',x,'g',g,'i',i))
#   print(res[[1]]$d) # d
#   boot = res[2:501]
#   df_boot = data.frame(matrix(unlist(boot), nrow=length(boot), byrow=TRUE))
#   colnames(df_boot) = c('tau', 'd')
#   print(quantile(df_boot$d, probs=c(0.05,0.95)))
# }
# 
# for (res in all_cptresults_sigmoid){
#   x = res[[1]]$expt
#   g = res[[1]]$gel
#   i = res[[1]]$inf
#   print(paste('x',x,'g',g,'i',i))
#   print(res[[1]]$d) # d
#   boot = res[2:501]
#   df_boot = data.frame(matrix(unlist(boot), nrow=length(boot), byrow=TRUE))
#   colnames(df_boot) = c('tau', 'd')
#   print(quantile(df_boot$d, probs=c(0.05,0.95)))
# }



