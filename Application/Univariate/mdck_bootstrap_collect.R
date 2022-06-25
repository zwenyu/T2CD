# collecting results for bootstrap for univariate methods of MDCK data
# slurm version
root <- "~/T2CD"

library(MASS)
require(arfima)
require(fracdiff)
require(splines)
require(ggplot2)
require(mvtnorm)

source(paste(root, "/helperfunction/helperfn.R", sep = ""))

### BSC
load(paste(root, "/Cornell_MDCK_M_hominis_runs.RData", sep = ""))
dat.info = list(gel=gel, inf=inf, nor=nor, null=null, wou=wou)

f = 1
numsim <- 200
for (x in 1:4){
  for (g in 0:1){
    for (i in 0:1){
      dat = extractdata(dat.com, dat.info, expt = x, freq = f, gel = g, inf = i, nor = (1-i))
      M = nrow(dat$res)
      if (x == 1 & g == 0 & i == 0) {
        pars <- cbind(x, g, i, 1:M)
      } else {
        pars <- rbind(pars, cbind(x, g, i, 1:M))
      }
    }
  }
}

pars <- cbind(pars, 1:nrow(pars), 1)
colnames(pars) <- c("x", "g", "i", "m", "par", "id")
sims <- cbind("id" = 1, "sim" = 1:numsim)
pars <- merge(pars, sims, by = "id")
pars <- pars[order(pars$sim, decreasing = FALSE), ]
pars <- data.frame(pars)
pars$d.step <- pars$d.sig <- pars$tau.step <- pars$tau.sig <- NA
pars.save <- pars

for (k in 1:nrow(pars.save)) {
  
  x <- pars.save[k, "x"]
  g <- pars.save[k, "g"]
  i <- pars.save[k, "i"]
  m <- pars.save[k, "m"]
  sim <- pars.save[k, "sim"]
  par <- pars.save[k, "par"]
  
  load(paste(root, "/Out/Data/sim_mdck_step_estimates_", par, "_", sim, ".RData", sep = ""))
  pars.save[k, c("tau.step", "d.step")] <- unlist(cptresults_step)
  
  load(paste(root, "/Out/Data/sim_mdck_sigmoid_estimates_", par, "_", sim, ".RData", sep = ""))
  pars.save[k, c("tau.sig", "d.sig")] <- unlist(cptresults_sigmoid)
}

save(pars.save, 
     file = paste(root, "/Out/Data/mdck_boot_collect.RData", sep = ""))


