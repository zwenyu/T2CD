# running bootstrap for multivariate methods of MDCK data
# slurm version
root <- "~/T2CD"

library(MASS)
require(arfima)
require(fracdiff)
require(splines)
require(ggplot2)
require(mvtnorm)

source(paste(root, "/method/bsWLS.R", sep = ""))
source(paste(root, "/method/t2cd_step.mv.R", sep = ""))
source(paste(root, "/method/t2cd_sigmoid.mv.R", sep = ""))
source(paste(root, "/helperfunction/helperfn.R", sep = ""))
source(paste(root, "/helperfunction/logL_fima.R", sep = ""))

boot_k_mv = function(k, results, plot_results, methodname){
  
  if (methodname=='step'){
    samp = bootstrap_sample_step_mv(results, plot_results, seed = k)
    res = t2cd_step_mv(samp, use_t = use_t, t_resd = T)
  }else if (methodname=='sigmoid'){
    samp = bootstrap_sample_sigmoid_mv(results, plot_results, seed = k)
    res = t2cd_sigmoid_mv(samp, use_t = use_t, t_resd = T)
  }
  tau = res$tau
  d = res$d
  
  return(list(tau=tau, d=d))
}

### MDCK
load(paste(root, "/Cornell_MDCK_M_hominis_runs.RData", sep = ""))
dat.info = list(gel=gel, inf=inf, nor=nor, null=null, wou=wou)

f = 1
numsim <- 500
for (x in 1:4){
  for (g in 0:1){
    for (i in 0:1){
      if (x == 1 & g == 0 & i == 0) {
        pars <- cbind(x, g, i)
      } else {
        pars <- rbind(pars, cbind(x, g, i))
      }
    }
  }
}

pars <- cbind(pars, 1:nrow(pars), 1)
colnames(pars) <- c("x", "g", "i", "par", "id")
sims <- cbind("id" = 1, "sim" = 1:numsim)
pars <- merge(pars, sims, by = "id")
pars <- pars[order(pars$sim, decreasing = FALSE), ]

par <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(par)

x <- pars[par, "x"]
g <- pars[par, "g"]
i <- pars[par, "i"]
sim <- pars[par, "sim"]
par <- pars[par, "par"]

dat = extractdata(dat.com, dat.info, expt = x, freq = f, gel = g, inf = i, nor = (1-i))

# step method
load(paste(root, "/Out/Data/mdck_step_estimates_mv_", par, ".RData", sep = ""))
cptresults_step <- boot_k.mv(k = sim, res_step, plot_results_step, 'step')
save(cptresults_step,
     file = paste(root, "/Out/Data/sim_mdck_step_estimates_mv_", par, "_", sim, ".RData", sep = ""))

# weighted method
load(paste(root, "/Out/Data/mdck_sigmoid_estimates_mv_", par, ".RData", sep = ""))
cptresults_sigmoid <- boot_k.mv(k = sim, res_sigmoid, plot_results_sigmoid, 'sigmoid')
save(cptresults_sigmoid,
     file = paste(root, "/Out/Data/sim_mdck_sigmoid_estimates_mv_",
                  par, "_", sim, ".RData", sep = ""))


