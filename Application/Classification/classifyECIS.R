setwd('/Users/wenyuzhang/Desktop/Cornell/Change points/ECIS/')
source('./Applications2/Classification/feature.R')
source('./helperfunction/helperfn.R')
library(dplyr)
library(BSDA)

mdck_ftd_exact = extractFeat(cellline = 'mdck', freq = 1, method = 'mat_step')
bsc_ftd_exact = extractFeat(cellline = 'bsc', freq = 1, method = 'mat_step')
mdck_ftd_approx = extractFeat(cellline = 'mdck', freq = 1, method = 'mat_wt')
bsc_ftd_approx = extractFeat(cellline = 'bsc', freq = 1, method = 'mat_wt')
mdck_ftd_ecp = extractFeat(cellline = 'mdck', freq = 1, method = 'mat_ecp')
bsc_ftd_ecp = extractFeat(cellline = 'bsc', freq = 1, method = 'mat_ecp')
mdck_ftd_ecpfdiff = extractFeat(cellline = 'mdck', freq = 1, method = 'mat_ecp.fdiff')
bsc_ftd_ecpfdiff = extractFeat(cellline = 'bsc', freq = 1, method = 'mat_ecp.fdiff')

bsc_ftd_exact %>% group_by(expt, gel, inf) %>%
  summarise(mean = mean(d), sd = sd(d))
bsc_ftd_approx %>% group_by(expt, gel, inf) %>%
  summarise(mean = mean(d), sd = sd(d))
exact = mdck_ftd_exact
approx = mdck_ftd_approx
for (x in 1:4){
  for (g in 0:1){
    for (i in 0:1){
      exact_d = exact[exact$expt==x & exact$gel==g & exact$inf==i, 'd']
      approx_d = approx[approx$expt==x & approx$gel==g & approx$inf==i, 'd']
      test = SIGN.test(exact_d - approx_d)
      if (test$p.value < 0.05){
        print(paste(x, g, i))
      }
    }
  }
}

mdck_acc_exact = classify(mdck_ftd_exact)
bsc_acc_exact = classify(bsc_ftd_exact)
mdck_acc_approx = classify(mdck_ftd_approx)
bsc_acc_approx = classify(bsc_ftd_approx)
mdck_acc_ecp = classify(mdck_ftd_ecp)
bsc_acc_ecp = classify(bsc_ftd_ecp)
mdck_acc_ecpfdiff = classify(mdck_ftd_ecpfdiff)
bsc_acc_ecpfdiff = classify(bsc_ftd_ecpfdiff)

# summarizing accuracy in terms of mean and variance
acc_summary = function(mat_list){
  res = lapply(mat_list, function(mat){
    mat_mean = colMeans(mat)
    mat_sd = apply(mat, 2, sd)
    return(rbind(mat_mean, mat_sd))
  })
  return(res)
}

mdck_exact = acc_summary(mdck_acc_exact)
bsc_exact = acc_summary(bsc_acc_exact)
mdck_approx = acc_summary(mdck_acc_approx)
bsc_approx = acc_summary(bsc_acc_approx)
mdck_ecp = acc_summary(mdck_acc_ecp)
bsc_ecp = acc_summary(bsc_acc_ecp)
mdck_ecpfdiff = acc_summary(mdck_acc_ecpfdiff)
bsc_ecpfdiff = acc_summary(bsc_acc_ecpfdiff)

mdck_exact
bsc_exact
mdck_approx
bsc_approx
mdck_ecp
bsc_ecp
mdck_ecpfdiff
bsc_ecpfdiff



