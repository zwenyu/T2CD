# classification helper functions for classifyECIS.R

library(klaR)
library(pracma)

# extract feature from one sequence
# 1. resh - resistance at selected time
# 2. maxres - the maximum resistance across the time series 
# 3. endofrun - average resistance over last 5 time points in time series 
extractFeat1 = function(res, tim, cellline){
  # take simple moving average with window length 5
  res_sma = movavg(res, 5, type = 's')
  tim_sma = movavg(tim, 5, type = 's')
  # 1. resh - resistance at 17 hours for mdck, 2 hours for bsc
  if (cellline == 'mdck'){
    resh = res_sma[which.min(abs(tim_sma - 17))]  
  }else if (cellline == 'bsc'){
    resh = res_sma[which.min(abs(tim_sma - 2))]
  }
  # 2. maxres - the maximum resistance across the time series 
  maxres = max(res_sma)
  # 3. endofrun - average resistance over last 5 time points in time series 
  endofrun = res_sma[length(res_sma)]
  return(c(resh, maxres, endofrun))
} 

# extract features from all data, and combine with d and tau
extractFeat = function(method, cellline = c('mdck', 'bsc'), freq = 1){
  if (cellline == 'mdck'){
    load("Cornell_MDCK_M_hominis_runs.RData")
  }else if (cellline == 'bsc'){
    load("Cornell_BSC_M_hominis_runs.RData")
  }
  dat.info = list(gel=gel, inf=inf, nor=nor, null=null, wou=wou)
  t.max = 72
  
  # extract features
  featmat = matrix(NA, 0, 8)
  colnames(featmat) = c('resh', 'maxres', 'endofrun', 'expt', 'freq', 'gel', 'inf', 'm')
  for (x in 1:4){
    for (g in 0:1){
      for (i in 0:1){
        dat = extractdata(dat.com, dat.info, expt = x, freq = freq, gel = g, inf = i, nor = (1-i))
        M = nrow(dat$res)
        for (m in 1:M){
          tim.ind = !is.na(dat$tim[m,]) & dat$tim[m,] <= t.max
          t.maxidx = which(tim.ind == T)
          res = dat$res[m,t.maxidx]
          tim = dat$tim[m,t.maxidx]
          feat1 = extractFeat1(res, tim, cellline)
          featmat = rbind(featmat, c(feat1, x, freq, g, i, m))
        }
      }
    }
  }
  
  # combine with d and tau
  if (cellline == 'mdck'){
    load("./Application/Univariate/mdck.RData")
  }else if (cellline == 'bsc'){
    load("./Application/Univariate/bsc.RData")
  }
  mat_ftd = merge(featmat, cptresults[[freq]][[method]])
  return(mat_ftd)
}

# classification accuracy with a given expt as test set
classify1 = function(mat_ftd, testexpt,
                     classifier = c('lda', 'qda', 'rda')){
  # set rda parameters based on classifier picked
  if (classifier == 'lda'){
    gamma = 0
    lambda = 1
  }else if (classifier == 'qda'){
    gamma = 0
    lambda = 0
  }else{
    # restrict to cases where rda is between lda and qda
    gamma = 0
    lambda = NA
  }
  
  # train-test split
  train = data.frame(mat_ftd[mat_ftd$expt != testexpt, ])
  test = data.frame(mat_ftd[mat_ftd$expt == testexpt, ])
  
  # features according to whether tau and d are included
  # using original features
  set.seed(0)
  model_orig = rda(inf ~ resh + maxres + endofrun, data = train,
                   gamma = gamma, lambda = lambda)
  pred_orig = as.numeric(predict(model_orig, test)$class) - 1
  acc_orig = 1 - mean(abs(pred_orig - test$inf))
  lambda_orig = model_orig$regularization['lambda']
  # using only tau
  set.seed(0)
  model_t = rda(inf ~ tau, data = train,
                gamma = gamma, lambda = lambda)
  pred_t = as.numeric(predict(model_t, test)$class) - 1
  acc_t = 1 - mean(abs(pred_t - test$inf))
  lambda_t = model_t$regularization['lambda']
  # using only d
  set.seed(0)
  model_d = rda(inf ~ d, data = train,
                gamma = gamma, lambda = lambda)
  pred_d = as.numeric(predict(model_d, test)$class) - 1
  acc_d = 1 - mean(abs(pred_d - test$inf))
  lambda_d = model_d$regularization['lambda']
  # using d and tau and one other original feature
  #
  model_td1 = rda(inf ~ d + tau + resh, data = train,
                  gamma = gamma, lambda = lambda)
  pred_td1_train = as.numeric(predict(model_td1, train)$class) - 1
  acc_td1_train = 1 - mean(abs(pred_td1_train - train$inf))
  pred_td1 = as.numeric(predict(model_td1, test)$class) - 1
  acc_td1 = 1 - mean(abs(pred_td1 - test$inf))
  #
  model_td2 = rda(inf ~ d + tau + maxres, data = train,
                  gamma = gamma, lambda = lambda)
  pred_td2_train = as.numeric(predict(model_td2, train)$class) - 1
  acc_td2_train = 1 - mean(abs(pred_td2_train - train$inf))
  pred_td2 = as.numeric(predict(model_td2, test)$class) - 1
  acc_td2 = 1 - mean(abs(pred_td2 - test$inf))
  #
  model_td3 = rda(inf ~ d + tau + endofrun, data = train,
                  gamma = gamma, lambda = lambda)
  pred_td3_train = as.numeric(predict(model_td3, train)$class) - 1
  acc_td3_train = 1 - mean(abs(pred_td3_train - train$inf))
  pred_td3 = as.numeric(predict(model_td3, test)$class) - 1
  acc_td3 = 1 - mean(abs(pred_td3 - test$inf))
  idx.train = which.max(c(acc_td1_train, acc_td2_train, acc_td3_train))
  acc_tdw = c(acc_td1, acc_td2, acc_td3)[idx.train]
  lambda_tdw = c(model_td1$regularization['lambda'], model_td2$regularization['lambda'],
                 model_td3$regularization['lambda'])[idx.train]
  # using all features
  set.seed(0)
  model_all = rda(inf ~ resh + maxres + endofrun + tau + d, data = train,
                  gamma = gamma, lambda = lambda)
  pred_all = as.numeric(predict(model_all, test)$class) - 1
  acc_all = 1 - mean(abs(pred_all - test$inf))
  lambda_all = model_all$regularization['lambda']
  
  return(c(acc_orig = acc_orig, acc_t = acc_t, acc_d = acc_d,
           acc_tdw = acc_tdw, acc_all = acc_all, 
           lambda_orig, lambda_t, lambda_d, lambda_tdw, lambda_all))
}

# classification with each expt taking turns being the test set
classify = function(mat_ftd, classifier){
  lda_mat = qda_mat = rda_mat = matrix(NA, 0, 10)
  for (x in 1:4){
    lda1 = classify1(mat_ftd, x, classifier)
    lda_mat = rbind(lda_mat, lda1)
  }
  
  return(list(lda_mat = lda_mat))
}


