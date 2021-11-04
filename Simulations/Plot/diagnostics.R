source('./helperfunction/helperfn.R')
source('./method/t2cd_sigmoid.R')

require(ggplot2)
require(reshape2)
require(RColorBrewer)

simnum = 100
tau_list = seq(15, 45, by = 5)
taupercent_list = tau_list/70
truetau_list = tau_list

# small d (stationary)
sim_ij = sim.simple(tau.percent=taupercent_list[3], d = -0.45, 
                    phip = 0, piq = 0,
                    hetero1 = 'TRUE', regime1 = 'gp',
                    seed = 1)

# fit
res_sigmoid = t2cd_sigmoid(sim_ij)
fit1 = plot.t2cd_sigmoid(res_sigmoid)
setEPS()
cairo_ps(file = './Simulations/Univariate/Figures/fitsigmoid_smalld.eps', width = 6, height = 3, pointsize = 12)
ggplot() +
  geom_line(data = data.frame(tim = res_sigmoid$tim[1,], res = res_sigmoid$res[1,]), 
            aes(x = tim, y = res)) +
  geom_vline(xintercept = res_sigmoid$tau, linetype = 'dashed') +
  labs(y = "Resistance (ohm)", x = 'Time (hour)') +
  theme_bw(base_size = 20)
dev.off()

cairo_ps(file = './Simulations/Univariate/Figures/wt_smalld.eps', width = 6, height = 3, pointsize = 12)
ggplot(data = data.frame(tim = res_sigmoid$tim[1,], wt = fit1$wt[1,]),
       aes(x = tim, y = wt)) +
  geom_line() + 
  geom_vline(xintercept = res_sigmoid$tau, linetype = 'dashed') + 
  labs(y = expression(paste('w(t,', hat(alpha), ')')), x = 'Time (hour)') +
  theme_bw(base_size = 20)
dev.off()

### large d (nonstationary)
sim_ij = sim.simple(tau.percent=taupercent_list[3], d = 1.35, 
                     phip = 0, piq = 0,
                     hetero1 = 'TRUE', regime1 = 'gp',
                     seed = 1)

# fit
res_sigmoid = t2cd_sigmoid(sim_ij)
fit1 = plot.t2cd_sigmoid(res_sigmoid)
setEPS()
cairo_ps(file = './Simulations/Univariate/Figures/fitsigmoid_larged.eps', width = 6, height = 3, pointsize = 12)
ggplot() +
  geom_line(data = data.frame(tim = res_sigmoid$tim[1,], res = res_sigmoid$res[1,]), 
            aes(x = tim, y = res)) +
  geom_vline(xintercept = res_sigmoid$tau, linetype = 'dashed') +
  labs(y = "Resistance (ohm)", x = 'Time (hour)') +
  theme_bw(base_size = 20)
dev.off()

cairo_ps(file = './Simulations/Univariate/Figures/wt_larged.eps', width = 6, height = 3, pointsize = 12)
ggplot(data = data.frame(tim = res_sigmoid$tim[1,], wt = fit1$wt[1,]),
       aes(x = tim, y = wt)) +
  geom_line() + 
  geom_vline(xintercept = res_sigmoid$tau, linetype = 'dashed') + 
  labs(y = expression(paste('w(t,', hat(alpha), ')')), x = 'Time (hour)') +
  theme_bw(base_size = 20)
dev.off()


