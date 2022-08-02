source('./helperfunction/helperfn.R')
source('./method/t2cd_step.R')

require(ggplot2)
require(reshape2)
require(RColorBrewer)
library(FinTS)

### MDCK
load("Cornell_MDCK_M_hominis_runs.RData")
dat.info = list(gel=gel, inf=inf, nor=nor, null=null, wou=wou)

# expt 1, freq 500, gel medium, infected
inf_dat = extractdata(dat.com, dat.info, expt = 1, freq = 1, gel = 1, inf = 1, nor = 0)
t.max = 72
tim.ind = !is.na(inf_dat$tim[1,]) & inf_dat$tim[1,] <= t.max
t.maxidx = which(tim.ind == T)
res = inf_dat$res[2,t.maxidx]
res_mean = scale(res, center = F) 
tim = inf_dat$tim[2,t.maxidx]

# fit
df1 = data.frame(tim=tim, res=res)
res_step = t2cd_step(df1, use_arf = F)
fit1 = plot.t2cd_step(res_step, use_arf = F)
setEPS()
cairo_ps(file = './Application/Plot/Figures/mdck_inf_res.eps', width = 6, height = 3, pointsize = 12)
ggplot() +
  geom_line(data = data.frame(tim = res_step$tim, res = res_step$res), 
            aes(x = tim, y = res)) +
  labs(y = "Resistance (ohm)", x = 'Time (hour)') +
  theme_bw(base_size = 20)
dev.off()

# residuals
r1 = (res_mean[1:res_step$idx] - fit1$fit.vals1)/sqrt(fit1$var.resd1)
Box.test(r1, lag = 20, type = "Ljung-Box")
cairo_ps(file = './Application/Plot/Figures/mdck_inf_residual1_sq.eps', width = 4, height = 3, pointsize = 12)
plot(tim[1:res_step$idx], r1, xlab = 'Time (hour)', ylab = 'Standardized Residuals')
dev.off()
cairo_ps(file = './Application/Plot/Figures/mdck_inf_qqplot1_sq.eps', width = 4, height = 3, pointsize = 12)
qqnorm(r1)
qqline(r1)
dev.off()
cairo_ps(file = './Application/Plot/Figures/mdck_inf_acf1_sq.eps', width = 4, height = 3, pointsize = 12)
Acf(r1, main = 'ACF', lag.max=50, ci.col='grey')
dev.off()
r2 = res_mean[(res_step$idx+1):length(res)] - fit1$fit.vals2
cairo_ps(file = './Application/Plot/Figures/mdck_inf_residual2_sq.eps', width = 4, height = 3, pointsize = 12)
plot(tim[(res_step$idx+1):length(res)], r2, xlab = 'Time (hour)', ylab = 'Standardized Residuals')
dev.off()
cairo_ps(file = './Application/Plot/Figures/mdck_inf_qqplot2_sq.eps', width = 4, height = 3, pointsize = 12)
qqnorm(r2)
qqline(r2)
dev.off()
cairo_ps(file = './Application/Plot/Figures/mdck_inf_acf2_sq.eps', width = 4, height = 3, pointsize = 12)
Acf(r2, main = 'ACF', lag.max=50, ci.col='grey')
dev.off()