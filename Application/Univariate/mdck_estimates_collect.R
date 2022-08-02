# collecting results for univariate methods of MDCK data
# slurm version
root <- "~/T2CD"

library(MASS)
require(arfima)
require(fracdiff)
require(splines)
require(ggplot2)
require(mvtnorm)

source(paste(root, "/helperfunction/helperfn.R", sep = ""))

### MDCK
load(paste(root, "/Cornell_MDCK_M_hominis_runs.RData", sep = ""))
dat.info = list(gel=gel, inf=inf, nor=nor, null=null, wou=wou)

f = 1
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

colnames(pars) <- c("x", "g", "i", "m")
pars <- data.frame(pars)
pars$d.step <- pars$tau.step <- pars$d.sig <- pars$tau.sig <- NA
pars.save <- pars

for (k in 1:nrow(pars.save)) {

  x <- pars.save[k, 1]
  g <- pars.save[k, 2]
  i <- pars.save[k, 3]
  m <- pars.save[k, 4]
  
  load(paste(root, "/Out/Data/mdck_step_estimates_", k, ".RData", sep = ""))
  pars.save[k, c("tau.step", "d.step")] <- c(res_step$tau, res_step$d)
  
  load(paste(root, "/Out/Data/mdck_sigmoid_estimates_", k, ".RData", sep = ""))
  pars.save[k, c("tau.sig", "d.sig")] <- c(res_sigmoid$tau, res_sigmoid$d)

}

save(pars.save, 
     file = paste(root, "/Out/Data/mdck_estimates_collect.RData", sep = ""))



