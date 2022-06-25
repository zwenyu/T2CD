# collecting results for multivariate methods of MDCK data
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
      if (x == 1 & g == 0 & i == 0) {
        pars <- cbind(x, g, i)
      } else {
        pars <- rbind(pars, cbind(x, g, i))
      }
    }
  }
}

colnames(pars) <- c("x", "g", "i")
pars <- data.frame(pars)
pars$d.step <- pars$d.sig <- NA
pars.save <- pars

for (k in 1:nrow(pars.save)) {

  x <- pars.save[k, 1]
  g <- pars.save[k, 2]
  i <- pars.save[k, 3]

  load(paste(root, "/Out/Data/mdck_step_estimates_mv_", k, ".RData", sep = ""))
  pars.save[k, c("d.step")] <- c(res_step$d)

  load(paste(root, "/Out/Data/mdck_sigmoid_estimates_mv_", k, ".RData", sep = ""))
  pars.save[k, c("d.sig")] <- c(res_sigmoid$d)

}

save(pars.save,
     file = paste(root, "/Out/Data/mdck_estimates_mv_collect.RData", sep = ""))



