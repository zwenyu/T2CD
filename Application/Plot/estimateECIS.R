# plot tau and d estimates 
# set working directory as T2CD

# plot for MDCK cell line
load("./Application/Univariate/mdck.RData")
load("./Application/Multivariate/mdck.mv.RData")
# switch to plot for BSC cell line
# load("./Application/Univariate/bsc.RData")
# load("./Application/Multivariate/bsc.mv.RData")

plot_univ = function(method, freq = 1, tau.range = c(10, 50), mvline = TRUE){
  
  d.range = c(0.5, 1.5)
  mat_f = cptresults[[freq]][[method]]
  
  # plot of d by group
  mat00 = mat_f[mat_f[,'gel']==0 & mat_f[,'inf']==0, 'd' ]
  mat10 = mat_f[mat_f[,'gel']==1 & mat_f[,'inf']==0, 'd' ]
  mat01 = mat_f[mat_f[,'gel']==0 & mat_f[,'inf']==1, 'd' ]
  mat11 = mat_f[mat_f[,'gel']==1 & mat_f[,'inf']==1, 'd' ]
  mat = data.frame(d = c(mat00, mat10, mat01, mat11), 
                   group = c(rep('BSA-Nor',length(mat00)), rep('Gel-Nor',length(mat10)), 
                             rep('BSA-Inf',length(mat01)), rep('Gel-Inf',length(mat11))))
  plot_d = ggplot(mat, aes(x=group, y=d)) + geom_boxplot() +
    theme_bw(base_size = 14)
  
  # plot of tau by group
  mat00t = mat_f[mat_f[,'gel']==0 & mat_f[,'inf']==0, 'tau' ]
  mat10t = mat_f[mat_f[,'gel']==1 & mat_f[,'inf']==0, 'tau' ]
  mat01t = mat_f[mat_f[,'gel']==0 & mat_f[,'inf']==1, 'tau' ]
  mat11t = mat_f[mat_f[,'gel']==1 & mat_f[,'inf']==1, 'tau' ]
  matt = data.frame(tau = c(mat00t, mat10t, mat01t, mat11t), 
                    group = c(rep('BSA-Nor',length(mat00t)), rep('Gel-Nor',length(mat10t)), 
                              rep('BSA-Inf',length(mat01t)), rep('Gel-Inf',length(mat11t))))
  plot_tau = ggplot(matt, aes(x=group, y=tau)) + geom_boxplot() +
    theme_bw(base_size = 14)
  
  # scatter of d vs tau, shape by expt, colored by group
  mat_dt = data.frame(mat_f[,c('d', 'tau', 'expt','gel','inf')])
  mat_dt$group = paste(ifelse(mat_dt$gel==1, 'Gel', 'BSA'), 
                       ifelse(mat_dt$inf==1, 'Inf', 'Nor'), sep = '-')
  mat_dt$expt = factor(paste('Expt', mat_dt$expt))
  mat_dt$Infection = factor(ifelse(mat_dt$inf==1, 'Yes', 'No'))
  mat_dt$serum = factor(ifelse(mat_dt$gel==1, 'Gel', 'BSA'))
  plot_scatter = ggplot(mat_dt, aes(x=tau, y=d, shape=expt, color=group)) + geom_point() +
    scale_color_manual(values = c('#FF66CC', '#0099FF', '#CC0033', '#0000CC')) +
    theme_bw(base_size = 14)
  
  # scatter, separated by experiment and serum type
  if (mvline){
    mat_f.mv = cptresults.mv[[freq]][[method]]
    mat_dt.mv = data.frame(mat_f.mv[,c('d', 'tau', 'expt','gel','inf')])
    mat_dt.mv$group = paste(ifelse(mat_dt.mv$gel==1, 'Gel', 'BSA'), 
                         ifelse(mat_dt.mv$inf==1, 'Inf', 'Nor'), sep = '-')
    mat_dt.mv$expt = factor(paste('Expt', mat_dt.mv$expt))
    mat_dt.mv$Infection = factor(ifelse(mat_dt.mv$inf==1, 'Yes', 'No'))
    mat_dt.mv$serum = factor(ifelse(mat_dt.mv$gel==1, 'Gel', 'BSA'))
    plot_grid = ggplot(mat_dt, aes(x=tau, y=d, shape=Infection, color=Infection)) + geom_point() +
      geom_hline(data = mat_dt.mv, aes(yintercept=d, linetype=Infection, color=Infection), alpha = 0.5) +
      facet_grid(serum~expt) + scale_color_manual(values = c("#3333CC", '#FF3333')) +
      scale_linetype_manual(values = c("solid", "longdash")) + 
      scale_shape_manual(values = c(15, 19)) + 
      theme_bw(base_size = 14) + labs(x=expression(tau))
  }else{
    plot_grid = ggplot(mat_dt, aes(x=tau, y=d, shape=Infection, color=Infection)) + geom_point() +
      facet_grid(serum~expt) + scale_color_manual(values = c("#3333CC", '#FF3333')) +
      scale_shape_manual(values = c(15, 19)) + 
      theme_bw(base_size = 14) + labs(x=expression(tau))
  }

  return(list(plot_d = plot_d, plot_tau = plot_tau, 
              plot_scatter = plot_scatter, plot_grid = plot_grid))
}


# plotting tau and d estimates for T2CD_step
plot_step = plot_univ(method = 'mat_step', mvline=TRUE)
plot_step$plot_grid

# plotting tau and d estimates for T2CD_sigmoid
plot_sigmoid = plot_univ(method = 'mat_sigmoid', mvline=TRUE)
plot_sigmoid$plot_grid


