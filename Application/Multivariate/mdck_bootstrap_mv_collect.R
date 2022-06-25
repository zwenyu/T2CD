# running bootstrap for multivariate methods of MDCK data
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
numsim <- 100
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
pars <- data.frame(pars)
pars$d.step <- pars$d.sig <- NA
pars.save <- pars

for (k in 1:nrow(pars.save)) {
  
  x <- pars.save[k, "x"]
  g <- pars.save[k, "g"]
  i <- pars.save[k, "i"]
  sim <- pars.save[k, "sim"]
  par <- pars.save[k, "par"]
  
  load(paste(root, "/Out/Data/sim_mdck_step_estimates_mv_", par, "_", sim, ".RData", sep = ""))
  pars.save[k, c("d.step")] <- unlist(cptresults_step)[1]
  
  load(paste(root, "/Out/Data/sim_mdck_sigmoid_estimates_mv_", par, "_", sim, ".RData", sep = ""))
  pars.save[k, c("d.sig")] <- unlist(cptresults_sigmoid)[1]
}

save(pars.save, 
     file = paste(root, "/Out/Data/mdck_boot_mv_collect.RData", sep = ""))


