T2CD
--------------------------------------------------

Code implementing T2CD (Trend-to-Confluence Detector) procedure in paper Modeling Nonlinear Trend Followed by Long-Memory Equilibrium with Unknown Change Point (https://arxiv.org/abs/2007.09417).

Implementation
--------------------------------------------------

To reproduce results on ECIS data presented in the paper, go to the Application folder. To implement univariate methods, run `Application/Univariate/MDCK_t2cd.R` and `Application/Univariate/BSC_t2cd.R`. To implement multivariate methods, run `Application/Multivariate/MDCK_t2cd.mv.R` and `Application/Multivariate/BSC_t2cd.mv.R`. To get classification results, run `Application/Classification/classifyECIS.R`.

Example
--------------------------------------------------

```r
# simulate data where first regime is a Gaussian process
simdata = sim.simple(Tend=70, N=400, tau.percent=0.25, d=0.15, regime1='gp', sig1=2, sig2=0.5, seed=1)

# T2CD-step
res_t2cd_step = t2cd_step(simdata, t.max = 70, use_arf = F) 
res_t2cd_step$tau
res_t2cd_step$d
fit_t2cd_step = plot.t2cd_step(res_t2cd_step, use_arf = F)

# T2CD-sigmoid
res_t2cd_sigmoid = t2cd_sigmoid(list(res=matrix(simdata$res,1), tim=matrix(simdata$tim,1)), t.max = 70)
res_t2cd_sigmoid$tau
res_t2cd_sigmoid$d
fit_t2cd_sigmoid = plot.t2cd_sigmoid(res_t2cd_sigmoid)
```