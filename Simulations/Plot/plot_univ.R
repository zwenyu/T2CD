# set working directory as T2CD

require(dplyr)
require(ggplot2)
require(reshape2)
require(RColorBrewer)
load('./Simulations/Univariate/simulate_gpTRUE00.RData')

# names for methods
methodnames = c('step', 'sigmoid', 'ecp', 'ecp.diff', 'seg', 'truetau')
methodnames_ext = c('T2CD-step', 'T2CD-sigmoid', 'ECP', 'ECP.diff', 'FixedTau', 'TrueTau')
tau_comnames = paste('tau', methodnames[-c(5,6)], sep = '_')
d_comnames = paste('d', methodnames, sep = '_')
# color map
colmap = brewer.pal(n = length(methodnames_ext), name = "Dark2")
names(colmap) = methodnames_ext
# shape map
shapemap = 0:(length(methodnames_ext)-1)
names(shapemap) = methodnames_ext

# simulation settings
d_list = seq(-0.45, 1.45, by = 0.2)
tau_list = seq(15, 45, by = 5)


### plotting tau estimates of T2CD methods across all tau and d values tested

# combining results
tau_com = vector("list", length(tau_comnames))
names(tau_com) = tau_comnames

for (i in 1:simnum){
  for (m in tau_comnames){
    tau_com[[m]] = rbind(tau_com[[m]], 
                         t(t(cptresults[[i]][[m]])))
  }
}
for (m in tau_comnames){
  colnames(tau_com[[m]]) = tau_list
}

gr = expand.grid(d_list, tau_list)
slopeline_tau = data.frame(gr, slope = rep(5, nrow(gr)), intercept = rep(10, nrow(gr)))
colnames(slopeline_tau) = c('d', 'tau', 'slope', 'intercept')

slopeline_tau_stationary = slopeline_tau[slopeline_tau$d<0.5,]
slopeline_tau_nonstationary = slopeline_tau[slopeline_tau$d>0.5,]

# tau
tau_step = melt(tau_com$tau_step)
tau_sigmoid = melt(tau_com$tau_sigmoid)
colnames(tau_step) = colnames(tau_sigmoid) = c('d', 'tau', 'est')
tau_step$Method = 'T2CD-step'
tau_sigmoid$Method = 'T2CD-sigmoid'

tau_stationary = rbind(tau_step[tau_step$d<0.5,], tau_sigmoid[tau_sigmoid$d<0.5,])
tau_nonstationary = rbind(tau_step[tau_step$d>0.5,], tau_sigmoid[tau_sigmoid$d>0.5,])
d_list_stationary = d_list[d_list<0.5]
d_label_stationary = paste('True d = ', d_list_stationary, sep = '')
names(d_label_stationary) = d_list_stationary
d_list_nonstationary = d_list[d_list>0.5]
d_label_nonstationary = paste('True d = ', d_list_nonstationary, sep = '')
names(d_label_nonstationary) = d_list_nonstationary
labeller_stationary = labeller(.cols=d_label_stationary)
labeller_nonstationary = labeller(.cols=d_label_nonstationary)

# plotting
cairo_ps(file = './Simulations/Univariate/Figures/univgp_stationary_tau.eps', width = 9, height = 5, pointsize = 12)
p = ggplot(tau_stationary, aes(x = factor(tau), y = est, fill = factor(d))) + 
  geom_boxplot() +
  facet_grid(Method~factor(d),
             labeller = labeller_stationary) + 
  geom_abline(data = slopeline_tau_stationary, aes(slope = slope, intercept = intercept),
              linetype = 'dashed')
p + labs(x = expression(paste("True ", tau)), y = expression(paste("Estimated ", tau)), fill = 'True d') + 
  scale_fill_grey() + scale_color_grey() + 
  theme_bw(base_size = 20) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        strip.text.x = element_text(size = 10))
dev.off()

cairo_ps(file = './Simulations/Univariate/Figures/univgp_nonstationary_tau.eps', width = 9, height = 5, pointsize = 12)
p = ggplot(tau_nonstationary, aes(x = factor(tau), y = est, fill = factor(d))) + 
  geom_boxplot() +
  facet_grid(Method~factor(d),
             labeller = labeller_nonstationary) + 
  geom_abline(data = slopeline_tau_nonstationary, aes(slope = slope, intercept = intercept),
              linetype = 'dashed')
p + labs(x = expression(paste("True ", tau)), y = expression(paste("Estimated ", tau)), fill = 'True d') + 
  scale_fill_grey() + scale_color_grey() + 
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        strip.text.x = element_text(size = 10))
dev.off()


### plotting comparison of T2CD methods and FixedTau

# combining results
d_error = vector("list", length(d_comnames))
names(d_error) = d_comnames

for (i in 1:simnum){
  for (m in d_comnames){
    d_error[[m]] = cbind(d_error[[m]],
                         abs(cptresults[[i]][[m]] - d_list))
  }
}

d_step_error = melt(d_error$d_step)
d_sigmoid_error = melt(d_error$d_sigmoid)
d_seg_error = melt(d_error$d_seg)
colnames(d_step_error) = colnames(d_sigmoid_error) = colnames(d_seg_error) =
  c('d', 'tau', 'error')

d_step_avg = d_step_error %>% 
  group_by(tau) %>% dplyr::summarize(Mean = mean(error))
d_sigmoid_avg = d_sigmoid_error %>% 
  group_by(tau) %>% dplyr::summarize(Mean = mean(error))
d_seg_avg = d_seg_error %>% 
  group_by(tau) %>% dplyr::summarize(Mean = mean(error))

d_step_avg = d_step_error %>% 
  group_by(tau) %>% 
  do(data.frame(t(quantile(.$error, probs = c(0.25, 0.5, 0.75)))))
d_sigmoid_avg = d_sigmoid_error %>% 
  group_by(tau) %>% 
  do(data.frame(t(quantile(.$error, probs = c(0.25, 0.5, 0.75)))))
d_seg_avg = d_seg_error %>% 
  group_by(tau) %>% 
  do(data.frame(t(quantile(.$error, probs = c(0.25, 0.5, 0.75)))))
colnames(d_step_avg) = colnames(d_sigmoid_avg) = colnames(d_seg_avg) =
  c('tau', 'lower', 'median', 'upper')
d_step_avg$method = 'T2CD-step'
d_sigmoid_avg$method = 'T2CD-sigmoid'
d_seg_avg$method = 'FixedTau'

d_methods_error = rbind(d_step_avg, d_sigmoid_avg, d_seg_avg)

cairo_ps(file = './Simulations/Univariate/Figures/univgp_d_error.eps', width = 6, height = 3, pointsize = 12)
ggplot(d_methods_error, aes(x = tau, y = median, shape = method)) +
  geom_line(aes(color=method)) +
  geom_point(aes(color=method), size = 4) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=method), alpha = 0.15) +
  labs(x = expression(paste('True ', tau)),
       y = 'Error in d estimate') +
  theme_bw(base_size = 20) +
  scale_color_manual(values=colmap) +
  scale_shape_manual(values=shapemap) +
  scale_fill_manual(values=colmap)
dev.off()

