Miscellaneous (mostly) R Code
====================

This is a place for miscellaneous R and other code I've put together for clients, co-workers or myself for learning and demonstration purposes. The attempt is made to put together some well-commented and/or conceptually clear code from scratch, though most functionality is readily available in any number of well-developed R packages.  Typically examples are provided using such packages for comparison of results.


Model Fitting
-------------
Code related to fitting of various models. 

[standard regression](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/standardlm.R), 
[penalized regression](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/penalizedML.R), 
one factor random effects [(R)](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/onefactorRE.R) 
[(Julia)](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/onefactorRE.jl) 
[(Matlab)](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/onefactorRE.m), 
two factor random effects [(R)](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/twofactorRE.R) 
[(Julia)](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/twofactorRE.jl) 
[(Matlab)](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/twofactorRE.m), 
[cubic spline](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/cubicsplines.R), 
[hurdle poisson](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/hurdle.R), 
[hurdle negbin](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/hurdle.R), 
[zero-inflated poisson](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/poiszeroinfl.R), 
[zero-inflated negbin](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/NBzeroinfl.R), 
[Cox survival](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/survivalCox.R),
[confirmatory factor analysis](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/cfa_ml.R),
[EM mixture univariate](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/EM%20Examples/EM%20Mixture.R),
[EM mixture multivariate](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/EM%20Examples/EM%20Mixture%20MV.R),
[EM probit](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/EM%20Examples/EM%20algorithm%20for%20probit%20example.R),
[EM pca](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/EM%20Examples/EM%20for%20pca.R),
[EM probabilistic pca](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/EM%20Examples/EM%20algorithm%20for%20ppca.R),
[EM state space model](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/EM%20Examples/EM%20for%20state%20space%20unobserved%20components.R),
[Gaussian Process noisy](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/gp%20Examples/gaussianprocessNoisy.R),
[Gaussian Process noise-free](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/gp%20Examples/gaussianprocessNoiseFree.R)...

### Bayesian (mostly with Stan/rstan)
[BEST t-test](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/rstant_testBEST.R),
[linear regression](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/rstan_linregwithprior.R)
(Compare with [BUGS version](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/bugs_linreg.R), [JAGS](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/jags_linreg.R)),
[mixed model](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/rstan_MixedModelSleepstudy.R), 
[mixed model with correlated random effects](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/rstan_MixedModelSleepstudy_withREcorrelation.R), 
[beta regression](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/rstanBetaRegression.R),
mixed model with beta response [(Stan)](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/rstan_MixedModelBetaRegression.R) [(JAGS)](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/jags_MixedModelBetaRegression.R),
[mixture model](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/rstan_MixtureModel.R)
...

SC and TR
---------
Code specific to short courses and technical reports I put together from time to time.

[Introduction to R](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/SC%20and%20TR/coursecode.r),
[Generalized Additive Models](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/SC%20and%20TR/GAMS.R),
[Machine Learning](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/SC%20and%20TR/MLcode.R),
...

Other
-----
Random shenanigans.

[ggplot2 theme](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/Other/ggtheme.R),
FizzBuzz test [(R)](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/Other/fizzbuzz.R) [(julia)](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/Other/fizzbuzz.jl) [(Python)](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/Other/fizzbuzz.py),
Scrape xkcd [(R)](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/Other/xkcdscrape.R) [(Python)](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/Other/xkcdscrape.py), 
[Shakespearean Insulter](https://github.com/mclark--/Miscellaneous-R-Code/blob/master/Other/shakespeareanInsulter.R), 
...



