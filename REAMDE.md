T2CD
--------------------------------------------------

Code implementing T2CD (Trend-to-Confluence Detector) procedure in paper Modeling Nonlinear Trend Followed by Long-Memory Equilibrium with Unknown Change Point (https://arxiv.org/abs/2007.09417).

Code Description
--------------------------------------------------

The repository contains code implementing both univariate and multivariate versions of the T2CD-step and T2CD-sigmoid methods. T2CD (Trend-to-Confluence Detector) detects the change point between the trend regime and the confluence regime in a series. Measurements during the trend period are independent deviations from a smooth nonlinear function of time, and measurements during the equilibrium period are characterized by a simple long memory model. T2CD simultaneously estimates the parameters of the trend and equilibrium processes and locate the change point between the two. The repository also contains experiments on both simulations and real-world ECIS data.

ECIS Data Description
--------------------------------------------------

The data consists of four experiments for two cell types, Madin-Darby Canine Kidney (MDCK) cells and epithelial cells of African green monkey kidney origin (BSC-1 cells), using ECIS measurements. Each experiment corresponds to ECIS measurements on cells on a single tray of 96 wells obtained over the course of at least 72 hours. Of the 96 wells, 16 are left empty, 32 contain uncontaminated cells and 48 contain cells contaminated by M. hominis, a species of mycoplasma. Half of the wells were prepared using bovine serum albumin (BSA), and half were prepared using gel.

Data is available upon request.

Implementation
--------------------------------------------------

To reproduce results on simulated data presented in the paper, go to the Simulations folder. To implement univariate methods, run `Simulations/Univariate/simulate_gp.R`. To implement multivariate methods, run `Simulations/Multivariate/simulate.mv.R`.

To reproduce results on ECIS data presented in the paper, go to the Application folder. To implement univariate methods, run `Application/Univariate/MDCK_t2cd.R` for MDCK cell line and `Application/Univariate/BSC_t2cd.R` for BSC cell line. To implement multivariate methods, run `Application/Multivariate/MDCK_t2cd.mv.R` for MDCK cell line and `Application/Multivariate/BSC_t2cd.mv.R` for BSC cell line. To get classification results, run `Application/Classification/classifyECIS.R`.

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