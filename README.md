Miscellaneous (mostly) R Code
====================

This is a place for miscellaneous R and other code I've put together for clients, co-workers or myself for learning and demonstration purposes. The attempt is made to put together some well-commented and/or conceptually clear code from scratch, though most functionality is readily available in any number of well-developed R packages.  Typically examples are provided using such packages for comparison of results.

More recently, I usually just create a document of some kind rather than a standard `*.R` file, so you might check out the [docs](https://github.com/m-clark/docs) repo as well.


Model Fitting
-------------

Code related to fitting of various models. 

[standard regression](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/standardlm.R), 
[penalized regression](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/penalizedML.R), 
[gradient descent regression](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/lm_gradientdescent.R) [(online)](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/stochastic_gradientdescent.R), 
one factor random effects [(R)](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/onefactorRE.R) 
[(Julia)](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/onefactorRE.jl) 
[(Matlab)](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/onefactorRE.m), 
two factor random effects [(R)](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/twofactorRE.R) 
[(Julia)](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/twofactorRE.jl) 
[(Matlab)](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/twofactorRE.m), 
[mixed model via ML](https://m-clark.github.io/docs/mixedModels/mixedModelML.html)
[cubic spline](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/cubicsplines.R), 
[hurdle poisson](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/hurdle.R), 
[hurdle negbin](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/hurdle.R), 
[zero-inflated poisson](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/poiszeroinfl.R), 
[zero-inflated negbin](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/NBzeroinfl.R), 
[Cox survival](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/survivalCox.R),
[confirmatory factor analysis](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/cfa_ml.R),
[EM mixture univariate](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/EM%20Examples/EM%20Mixture.R),
[EM mixture multivariate](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/EM%20Examples/EM%20Mixture%20MV.R),
[EM probit](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/EM%20Examples/EM%20algorithm%20for%20probit%20example.R),
[EM pca](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/EM%20Examples/EM%20for%20pca.R),
[EM probabilistic pca](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/EM%20Examples/EM%20algorithm%20for%20ppca.R),
[EM state space model](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/EM%20Examples/EM%20for%20state%20space%20unobserved%20components.R),
[Gaussian Process noisy](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/gp%20Examples/gaussianprocessNoisy.R),
[Gaussian Process noise-free](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/gp%20Examples/gaussianprocessNoiseFree.R), 
[reproducing kernel hilbert space regression](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/RKHSReg/RKHSReg.md), 
[stochastic volatility](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/stochasticVolatility.R),
[bivariate probit](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/bivariateProbit.R),
[quantile regression](http://htmlpreview.github.io/?https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/quantileRegression.html),
[ordinal regression](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/ordinal_regression.R),
[naive bayes](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/naivebayes.R),
[extreme learning machine](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/elm.R),
[Chinese restaurant process, Indian buffet process](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/crp.R), 
[One-line models (an exercise)](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/one_line_models.R),
...


### Bayesian (mostly with Stan/rstan)

[BEST t-test](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/rstant_testBEST.R),
[linear regression](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/rstan_linregwithprior.R)
(Compare with [BUGS version](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/bugs_linreg.R), [JAGS](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/jags_linreg.R)),
[mixed model](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/rstan_MixedModelSleepstudy.R), 
[mixed model with correlated random effects](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/rstan_MixedModelSleepstudy_withREcorrelation.R), 
[beta regression](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/rstanBetaRegression.R),
mixed model with beta response [(Stan)](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/rstan_MixedModelBetaRegression.R) [(JAGS)](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/jags_MixedModelBetaRegression.R),
[mixture model](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/rstan_MixtureModel.R),
[topic model](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/topicModelgibbs.R),
[multilevel mediation](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/rstan_multilevelMediation.R), 
[variational bayes regression](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/variationalBayesRegression.Rmd), 
[gaussian process](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting//gp%20Examples/gaussianProcessStan.Rmd),
[horseshoe prior](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/horseshoe/README.md),
[item response theory](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/StanBugsJags/IRT_models), ...


SC and TR
---------

This part of the repository is deprecated, but used to be a section of 'short courses' and 'technical reports'.  See the [Workshops](https://github.com/m-clark/Workshops) or [docs](https://github.com/m-clark/docs) repositories instead, or go to the [workshops](http://m-clark.github.io/workshops/) and [documents](http://m-clark.github.io/documents/) sections of the website where you can see finished products ...


Other
-----

Random shenanigans.

FizzBuzz test [(R)](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/Other/fizzbuzz.R) [(julia)](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/Other/fizzbuzz.jl) [(Python)](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/Other/fizzbuzz.py),
Reverse a string recursively [(R)](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/Other/Programming_Shenanigans/stringReverseRecursively.R) [(Python)](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/Other/Programming_Shenanigans/stringReverseRecursively.py),
Recursive Word Wrap [(R)](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/Other/Programming_Shenanigans/wordWrap.R) [(Python)](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/Other/Programming_Shenanigans/wordWrap.py),
[calculate compound interest recursively](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/Other/Programming_Shenanigans/compound.R),
[get US Congress roll call data](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/Other/getRollCall.R),
Scrape xkcd [(R)](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/Other/xkcdscrape.R) [(Python)](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/Other/xkcdscrape.py), 
[Shakespearean Insulter](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/Other/shakespeareanInsulter.R), 
[ggplot2 theme](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/Other/ggtheme.R),
[spurious correlation with ratios](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/Other/spuriousCorrelationwithRatios.R),
[R matrix speedups](https://github.com/m-clark/Miscellaneous-R-Code/blob/master/Other/Programming_Shenanigans/matrixOperations.md), ...



