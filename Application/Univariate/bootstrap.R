library(parallel)
source('./helperfunction/helperfn.R')
source('./method/t2cd_step.R')

### MDCK
load("Cornell_MDCK_M_hominis_runs.RData")
dat.info = list(gel=gel, inf=inf, nor=nor, null=null, wou=wou)

num_boot = 500

# expt 1, freq 500, gel medium, infected
inf_dat = extractdata(dat.com, dat.info, expt = 1, freq = 1, gel = 1, inf = 1, nor = 0)
t.max = 72
tim.ind = !is.na(inf_dat$tim[1,]) & inf_dat$tim[1,] <= t.max
t.maxidx = which(tim.ind == T)
res = inf_dat$res[5,t.maxidx] # 2
tim = inf_dat$tim[5,t.maxidx]

# fit
df1 = data.frame(tim=tim, res=res)
results = t2cd_step(df1, use_arf = F)
plot_results = plot.t2cd_step(results, use_arf = F)

# Calculate the number of cores
no_cores = detectCores() - 1
# Initiate cluster
cl = makeCluster(no_cores, type = 'FORK')
# each iteration of function output as a list entry
# each function output is a named list

# iterate through experiments, frequency, gel and inf/nor
boot_k = function(k){
  
  samp = bootstrap_sample(results, plot_results, seed = k)
  
  # T2CD-step
  res_step = t2cd_step(samp, t.max = 70, use_arf = F)
  tau_step = res_step$tau
  d_step = res_step$d
  
  return(list(tau=tau_step, d=d_step))
}

cptresults = parLapply(cl, 1:num_boot, boot_k)

save.image('./Application/Univariate/mdck_boot.RData')
stopCluster(cl)

# retrieve results
# results$d # 1.118006
# results$tau # 39.25373
# # CI
# tau = c()
# d = c()
# for (i in 1:num_boot){
#   tau = c(tau, cptresults[[i]]$tau)
#   d = c(d, cptresults[[i]]$d)
# }
# quantile(d, probs=c(0.05,0.95)) # 0.9771801 1.1961741
# quantile(tau, probs=c(0.05,0.95)) # 23.75685 48.21105
