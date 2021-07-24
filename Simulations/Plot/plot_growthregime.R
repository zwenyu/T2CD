# set working directory as T2CD

require(dplyr)
require(ggplot2)
require(reshape2)
require(RColorBrewer)
load('./Simulations/Univariate/simulate_growth_gp.RData')

# simulation settings
simnum = 100
tau_list = seq(15, 45, by = 5)

### plotting rmse
rmse.full = fit1[[1]]$rmse.full
rmse.part = fit1[[1]]$rmse.part
for (i in 2:simnum){
  rmse.full = rbind(rmse.full, fit1[[i]]$rmse.full)  
  rmse.part = rbind(rmse.part, fit1[[i]]$rmse.part)  
}

fit1_df = data.frame(matrix(nrow = length(tau_list)*2, ncol = 5))
names(fit1_df) = c('tau', 'median', 'lower', 'upper', 'method')
fit1_df$tau = rep(tau_list, 2)
fit1_df$median = c(apply(rmse.full, 2, median), apply(rmse.part, 2, median))
fit1_df$lower = c(apply(rmse.full, 2, quantile, probs = 0.25), 
                  apply(rmse.part, 2, quantile, probs = 0.25))
fit1_df$upper = c(apply(rmse.full, 2, quantile, probs = 0.75), 
                  apply(rmse.part, 2, quantile, probs = 0.75))
fit1_df$method = rep(c('fit on whole series', 'fit on regime 1'), each = length(tau_list))

cairo_ps(file = './Simulations/Univariate/Figures/fit1_gp.eps', width = 10, height = 3, pointsize = 12)
p = ggplot(fit1_df, aes(x=tau, y=median, group=method, shape=method)) + 
  geom_line(aes(color=method)) +
  geom_point(aes(color=method), size = 4) + 
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=method), alpha = 0.15)
p + labs(x=expression(paste('True ', tau)), y = expression(paste("RMSE on regime 1")))+
  theme_bw(base_size = 20)
dev.off()

### plotting difference in spline coefficients
valscoeff.diff = fit1[[1]]$valscoeff.diff
valscoeff.diff_df = melt(valscoeff.diff)
colnames(valscoeff.diff_df) = c('d', 'tau', 'diff')
cairo_ps(file = './Simulations/Univariate/Figures/fit1_gp_valscoeff.eps', width = 5, height = 3, pointsize = 12)
p = ggplot(valscoeff.diff_df, aes(x=factor(tau), y=diff)) + 
  geom_boxplot() +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1))
p + labs(x=expression(paste('True ', tau)), y = expression(paste("Diff in coefficient")))+
  theme_bw(base_size = 20)
dev.off()

resdcoeff.diff = fit1[[1]]$resdcoeff.diff
resdcoeff.diff_df = melt(resdcoeff.diff)
colnames(resdcoeff.diff_df) = c('d', 'tau', 'diff')
cairo_ps(file = './Simulations/Univariate/Figures/fit1_gp_resdcoeff.eps', width = 5, height = 3, pointsize = 12)
p = ggplot(resdcoeff.diff_df, aes(x=factor(tau), y=diff)) + 
  geom_boxplot() +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1))
p + labs(x=expression(paste('True ', tau)), y = expression(paste("Diff in coefficient")))+
  theme_bw(base_size = 20)
dev.off()

### plotting differences in lambda hyperparameters
valslambda.diff = fit1[[1]]$valslambda.diff
valslambda.diff_df = melt(valslambda.diff)
colnames(valslambda.diff_df) = c('d', 'tau', 'diff')
cairo_ps(file = './Simulations/Univariate/Figures/fit1_gp_valslambda.eps', width = 5, height = 3, pointsize = 12)
p = ggplot(valslambda.diff_df, aes(x=factor(tau), y=diff)) + 
  geom_boxplot() +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1))
p + labs(x=expression(paste('True ', tau)), y = expression(paste("Diff in hyperparam")))+
  theme_bw(base_size = 20)
dev.off()

resdlambda.diff = fit1[[1]]$resdlambda.diff
resdlambda.diff_df = melt(resdlambda.diff)
colnames(resdlambda.diff_df) = c('d', 'tau', 'diff')
cairo_ps(file = './Simulations/Univariate/Figures/fit1_gp_resdlambda.eps', width = 5, height = 3, pointsize = 12)
p = ggplot(resdlambda.diff_df, aes(x=factor(tau), y=diff)) + 
  geom_boxplot() +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1))
p + labs(x=expression(paste('True ', tau)), y = expression(paste("Diff in hyperparam")))+
  theme_bw(base_size = 20)
dev.off()
