# running multivariate methods of MDCK data
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

### MDCK
load(paste(root, "/Cornell_MDCK_M_hominis_runs.RData", sep = ""))
dat.info = list(gel=gel, inf=inf, nor=nor, null=null, wou=wou)

f = 1
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


par <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

x <- pars[par, 1]
g <- pars[par, 2]
i <- pars[par, 3]

dat = extractdata(dat.com, dat.info, expt = x, freq = f, gel = g, inf = i, nor = (1-i))

# step method
res_step = t2cd_step_mv(dat, t_resd = T)
plot_results_step = plot.t2cd_step.mv(res_step, use_arf = F, return_plot = FALSE)
save(pars, par, res_step, plot_results_step, 
     file = paste(root, "/Out/Data/mdck_step_estimates_mv_", par, ".RData", sep = ""))


# weighted method           
res_sigmoid = t2cd_sigmoid_mv(dat, t_resd = T)
plot_results_sigmoid = plot.t2cd_sigmoid_mv(res_sigmoid, return_plot = FALSE)
save(pars, par, res_sigmoid, plot_results_sigmoid, 
     file = paste(root, "/Out/Data/mdck_sigmoid_estimates_mv_", par, ".RData", sep = ""))


