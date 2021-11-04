# set working directory as T2CD

require(ggplot2)
require(reshape2)
load('./Simulations/Multivariate/simulate_mvTRUE3.RData')

### plot comparing T2CD estimate of d from univariate and multivariate methods

# combining results
d_list = seq(-0.45, 1.45, by = 0.2)
poolnames = c('mult', 'univ', 'mean')
d_pool = vector("list", length(poolnames))
names(d_pool) = poolnames

for (i in 1:simnum){
  m = 'mult'
  d_pool[[m]] = cbind(d_pool[[m]], 
                      abs(cptresults_comp[[i]]$d_step - d_list))
  
  m = 'univ'
  d_pool[[m]] = cbind(d_pool[[m]], 
                      rowMeans(abs(cptresults_comp[[i]]$d_univ_step - d_list)))
}

d_pool_t2cd = melt(d_pool$mult)[,c(1,3)]
d_pool_univ = melt(d_pool$univ)[,c(1,3)]
colnames(d_pool_t2cd) = colnames(d_pool_univ) =
  c('d', 'error')
d_pool_t2cd$method = 'multivariate'
d_pool_univ$method = 'univariate'

d_pool_methods = rbind(d_pool_t2cd, d_pool_univ)

cairo_ps(file = './Simulations/Multivariate/Figures/mv_d_error.eps', width = 7, height = 3, pointsize = 12)
ggplot(d_pool_methods, aes(x=factor(d), y=error, fill = method)) + 
  geom_boxplot() +
  ylim(0, 0.3) +
  labs(x = 'True d',
       y = 'Error in d estimate') +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1))
dev.off()
