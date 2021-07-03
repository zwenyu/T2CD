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

# iterate through experiments, frequency, gel and inf/nor
ecisfreq = function(f){
  
  mat_step = mat_sigmoid = matrix(NA, 0, 8)
  colnames(mat_step) = colnames(mat_sigmoid) = c('d', 'tau', 'expt', 'freq', 'gel', 'inf', 'm', 'time')

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
          res_step = t2cd_step(dat_m, use_arf = F)
          ptime = proc.time() - s
          sink(NULL)
          mat_step = rbind(mat_step, c(res_step$d, res_step$tau, x, f, g, i, m, ptime[1]))
          
          # weighted method        
          s = proc.time()
          res_sigmoid = t2cd_sigmoid(list(res=matrix(dat_m$res,1), tim=matrix(dat_m$tim,1)))
          ptime = proc.time() - s
          mat_sigmoid = rbind(mat_sigmoid, c(res_sigmoid$d, res_sigmoid$tau, x, f, g, i, m, ptime[1]))
        }
      }
    }
  }
  
  return(list(mat_step = mat_step, mat_sigmoid = mat_sigmoid))
}

cptresults = parLapply(cl, 1, ecisfreq)

save.image('./Application/Univariate/mdck.RData')
stopCluster(cl)

