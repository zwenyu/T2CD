# obtain classification results for infection status
# set working directory as T2CD

require(dplyr)
source('./Application/Classification/feature.R')
source('./helperfunction/helperfn.R')

# extract features
mdck_ftd_step = extractFeat(cellline = 'mdck', freq = 1, method = 'mat_step')
bsc_ftd_step = extractFeat(cellline = 'bsc', freq = 1, method = 'mat_step')
mdck_ftd_sigmoid = extractFeat(cellline = 'mdck', freq = 1, method = 'mat_sigmoid')
bsc_ftd_sigmoid = extractFeat(cellline = 'bsc', freq = 1, method = 'mat_sigmoid')

mdck_ftd_step %>% group_by(expt, gel, inf) %>%
  summarise(mean = mean(d), sd = sd(d))
mdck_ftd_sigmoid %>% group_by(expt, gel, inf) %>%
  summarise(mean = mean(d), sd = sd(d))

bsc_ftd_step %>% group_by(expt, gel, inf) %>%
  summarise(mean = mean(d), sd = sd(d))
bsc_ftd_sigmoid %>% group_by(expt, gel, inf) %>%
  summarise(mean = mean(d), sd = sd(d))

# perform classfication by LDA
mdck_acc_step = classify(mdck_ftd_step)
bsc_acc_step = classify(bsc_ftd_step)
mdck_acc_sigmoid = classify(mdck_ftd_sigmoid)
bsc_acc_sigmoid = classify(bsc_ftd_sigmoid)

# summarizing accuracy in terms of mean and standard deviation
acc_summary = function(mat_list){
  res = lapply(mat_list, function(mat){
    mat_mean = colMeans(mat)
    mat_sd = apply(mat, 2, sd)
    return(rbind(mat_mean, mat_sd))
  })
  return(res)
}

mdck_step = acc_summary(mdck_acc_step)
bsc_step = acc_summary(bsc_acc_step)
mdck_sigmoid = acc_summary(mdck_acc_sigmoid)
bsc_sigmoid = acc_summary(bsc_acc_sigmoid)

# classification performance for MDCK
mdck_step
mdck_sigmoid
# classification performance for BSC
bsc_step
bsc_sigmoid





