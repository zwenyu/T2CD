library(parallel)
source('./helperfunction/helperfn.R')
source('./method/t2cd_step.R')

num_boot = 500

d_list = seq(-0.45, 1.45, by = 0.2)
tau_list = seq(15, 45, by = 5)
taupercent_list = tau_list/70
truetau_list = tau_list

# small d: -0.45
# big d: 1.35
sim_ij = sim.simple(tau.percent=taupercent_list[3], d = -0.45, 
                    phip = 0, piq = 0,
                    hetero1 = 'TRUE', regime1 = 'gp',
                    seed = 1)

# fit
results = t2cd_step(sim_ij, use_arf = F)
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

save.image('./Simulations/Univariate/smalld_boot.RData')
stopCluster(cl)