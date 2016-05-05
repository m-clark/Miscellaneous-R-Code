# Comparing performance of growth curve models in an SEM framework vs standard mixed models.

# 5 within cluster sample sizes
# 5 cluster sample sizes
# 5 correlation sizes
# balanced vs. not  # maybe in the future

# ints and Times
dataGen <- function(nclusters, nwithin, corr=0, balanced=T) {
  # setup
  nclus = nclusters                                                       # number of groups
  clus = factor(rep(1:nclus, each=nwithin))                               # cluster variable
  n = length(clus)                                                        # total n
                        
  # parameters                      
  sigma = 1                                                               # residual sd
  psi = matrix(c(1,corr,corr,1), 2, 2)                                    # re covar
  gamma_ = MASS::mvrnorm(nclus, mu=c(0,0), Sigma=psi, empirical=TRUE)     # random effects
  e = rnorm(n, mean=0, sd=sigma)                                          # residual error
  intercept = 3                                                           # fixed effects
  b1 = .75                            
                              
  # data                            
  x = rep(1:nwithin-1, times=nclus)                                       # covariate
  y = intercept+gamma_[clus,1] + (b1+gamma_[clus,2])*x  + e               # see model 1
  d = data.frame(time=x, y, clus=clus)
}

runNLME <- function(data) {
  nlmemod = lme(y ~ time, data=data, random =  ~time|clus, 
                control=list(maxIter=1000, msMaxIter=1000, msMaxEval=1000, returnObject=T), 
                weights=varIdent(form = ~1|time), method='ML')
  varRes = coef(nlmemod$modelStruct$varStruct, unconstrained =FALSE,allCoef=T)*nlmemod$sigma
  varRE  = as.numeric(VarCorr(nlmemod)[1:2,1])
  corRE = as.numeric(VarCorr(nlmemod)['time','Corr'])
  fixed = fixef(nlmemod)
  list(varRes=varRes, varRE=varRE, corRE=corRE, fixed=fixed)
}

runGC <- function(data) {
  ntime = unique(data$time)
  data$time = factor(data$time)
  dataWide = tidyr::spread(data, time, y)
  colnames(dataWide)[-1] = paste0('y', colnames(dataWide)[-1])

  IModel = paste0('I =~ 1*y0 ', paste0('+ 1*', colnames(dataWide)[-c(1:2)], collapse=''), '\n')
  SModel = paste0('S =~ 0*y0 ', paste0('+ ', 1:dplyr::last(ntime), '*', colnames(dataWide)[-c(1:2)], collapse=''), '\n')
  CenterY = paste0('y0', paste0(' + ', colnames(dataWide)[-c(1:2)], collapse=''), ' ~ 0*1')
  LVmodel = paste0(IModel, SModel, CenterY)

  suppressWarnings({semres = growth(LVmodel, data=dataWide)})
  varRes = sqrt(coef(semres)[ntime+1])
  varRE = coef(semres)[c('I~~I','S~~S')]
  covRE = coef(semres)[c('I~~S')]
  corRE = covRE/prod(sqrt(varRE))
  fixed = coef(semres)[c('I~1', 'S~1')]
  list(varRes=varRes, varRE=varRE, corRE=corRE, fixed=fixed)
}


# Data setup and generation
nclusters = c(10, 25, 50)
withinSizes = c(5, 10, 25)
corrs = seq(-.5,.5,.25)
grid = expand.grid(nclusters, withinSizes, corrs); colnames(grid) = c('nClusters', 'nWithinCluster', 'corRE')



library(parallel)
set.seed(1234)
clus = makeCluster(7)
clusterEvalQ(clus, library(nlme))
clusterEvalQ(clus, library(lavaan))
clusterExport(clus, c('runNLME', 'runGC', 'dataGen'))

dataList = parApply(clus, grid, 1, function(x) replicate(500, dataGen(nclusters=x[1], nwithin=x[2], corr=x[3]), simplify=F))

# because nlme keeps having issues in parallel (possibly only due to when I included nclus=5)
mixedResults = vector('list', length(dataList))

for (i in 1:length(dataList)) {
  mixedResults[[i]] = parLapply(clus, dataList[[i]], runNLME)
}

save(mixedResults, file='../data/growthvsMixedResults.RData')

growthResults = sapply(dataList, function(dat) parLapply(clus, dat, runGC), simplify=F)


# mixedResults = sapply(dataList, function(dat) parLapply(clus, dat, runNLME), simplify=F)

save(growthResults, mixedResults, file='../data/growthvsMixed_ModelResults.RData')

# summarize results growth
# fixed effects
growthFE0 = parSapply(clus, growthResults, function(x) sapply(x, function(res) res$fixed), simplify=F)
growthFE = lapply(growthFE0, function(res) c(rowMeans(res), c(apply(res, 1, quantile, p=c(.025,.975)))))
growthFE = do.call('rbind', growthFE); colnames(growthFE) = c('Int','Time', 'LL_Int', 'UL_Int', 'LL_Time', 'UL_Time')
growthFE = data.frame(grid, round(growthFE, 2))

# residual variance
growthvarRes0 = parSapply(clus, growthResults, function(x) sapply(x, function(res) res$varRes), simplify=F)

# random effects variance
growthvarRE0 = parSapply(clus, growthResults, function(x) sapply(x, function(res) res$varRE), simplify=F)
growthvarRE = lapply(growthvarRE0, function(res) c(rowMeans(res), c(apply(res, 1, quantile, p=c(.025,.975)))))
growthvarRE = do.call('rbind', growthvarRE); colnames(growthvarRE) = c('Int','Time', 'LL_Int', 'UL_Int', 'LL_Time', 'UL_Time')
growthvarRE = data.frame(grid, round(growthvarRE, 2))

# cor random effects
growthcorRE0 = parSapply(clus, growthResults, function(x) sapply(x, function(res) res$corRE), simplify=F)
growthcorRE = lapply(growthcorRE0, function(res) c(mean(res, na.rm=T), quantile(res, p=c(.025,.975), na.rm=T)))
growthcorRE = do.call('rbind', growthcorRE); colnames(growthcorRE) = c('corRE_est', 'LL_corRE', 'UL_corRE')
growthcorRE = data.frame(grid, round(growthcorRE, 2))

# summarize results mixed
# fixed effects
mixedFE0 = parSapply(clus, mixedResults, function(x) sapply(x, function(res) res$fixed), simplify=F)
mixedFE = lapply(mixedFE0, function(res) c(rowMeans(res), c(apply(res, 1, quantile, p=c(.025,.975)))))
mixedFE = do.call('rbind', mixedFE); colnames(mixedFE) = c('Int','Time', 'LL_Int', 'UL_Int', 'LL_Time', 'UL_Time')
mixedFE = data.frame(grid, round(mixedFE, 2))

# residual variance
mixedvarRes0 = parSapply(clus, mixedResults, function(x) sapply(x, function(res) res$varRes), simplify=F)

# random effects variance
mixedvarRE0 = parSapply(clus, mixedResults, function(x) sapply(x, function(res) res$varRE), simplify=F)
mixedvarRE = lapply(mixedvarRE0, function(res) c(rowMeans(res), c(apply(res, 1, quantile, p=c(.025,.975)))))
mixedvarRE = do.call('rbind', mixedvarRE); colnames(mixedvarRE) = c('Int','Time', 'LL_Int', 'UL_Int', 'LL_Time', 'UL_Time')
mixedvarRE = data.frame(grid, round(mixedvarRE, 2))

# cor random effects
mixedcorRE0 = parSapply(clus, mixedResults, function(x) sapply(x, function(res) res$corRE), simplify=F)
mixedcorRE = lapply(mixedcorRE0, function(res) c(mean(res), quantile(res, p=c(.025,.975))))
mixedcorRE = do.call('rbind', mixedcorRE); colnames(mixedcorRE) = c('corRE_est', 'LL_corRE', 'UL_corRE')
mixedcorRE = data.frame(grid, round(mixedcorRE, 2))

stopCluster(clus)

# Comparison of results ---------------------------------------------------
library(dplyr)

### Fixed effects
fixedEffects = left_join(growthFE, mixedFE, by=c('nClusters', 'nWithinCluster', 'corRE')) %>% 
  arrange(nClusters, nWithinCluster, corRE) %>% 
  mutate(widthGrowth_Int = UL_Int.x-LL_Int.x,
         widthMixed_Int =  UL_Int.y-LL_Int.y,
         mixedWidthMinusgrowthWidth_Int = widthMixed_Int-widthGrowth_Int,
         widthGrowth_Time = UL_Time.x-LL_Time.x,
         widthMixed_Time =  UL_Time.y-LL_Time.y,
         mixedWidthMinusgrowthWidth_Time = widthMixed_Time-widthGrowth_Time)

### random effects
randomEffects = left_join(growthvarRE, mixedvarRE, by=c('nClusters', 'nWithinCluster', 'corRE')) %>% 
  arrange(nClusters, nWithinCluster, corRE) %>% 
  mutate(widthGrowth_Int = UL_Int.x-LL_Int.x,
         widthMixed_Int =  UL_Int.y-LL_Int.y,
         mixedWidthMinusgrowthWidth_Int = widthMixed_Int-widthGrowth_Int,
         widthGrowth_Time = UL_Time.x-LL_Time.x,
         widthMixed_Time =  UL_Time.y-LL_Time.y,
         mixedWidthMinusgrowthWidth_Time = widthMixed_Time-widthGrowth_Time)


### cor of random effects
randomEffectsCor = left_join(growthcorRE, mixedcorRE, by=c('nClusters', 'nWithinCluster', 'corRE')) %>% 
  arrange(nClusters, nWithinCluster, corRE) %>% 
  mutate(widthGrowth = UL_corRE.x-LL_corRE.x,
         widthMixed =  UL_corRE.y-LL_corRE.y,
         mixedWidthMinusgrowthWidth = widthMixed-widthGrowth)

feEst = fixedEffects %>% select(one_of('nClusters', 'nWithinCluster', 'corRE'), Int.x, Time.x, Int.y, Time.y)
reEst = randomEffects %>% select(one_of('nClusters', 'nWithinCluster', 'corRE'), Int.x, Time.x, Int.y, Time.y)
corEst = randomEffectsCor %>% select(one_of('nClusters', 'nWithinCluster', 'corRE'),  corRE_est.x,  corRE_est.y)

biasFE = feEst %>% 
  arrange(nClusters, nWithinCluster, corRE) %>% 
  mutate(biasLGC_Int = Int.x-3,
         biasLGC_Time = Time.x-.75,
         biasMM_Int = Int.x-3,
         biasMM_Time = Time.y-.75)

biasRE = reEst %>% 
  arrange(nClusters, nWithinCluster, corRE) %>% 
  mutate(biasLGC_Int = Int.x-1,
         biasLGC_Time = Time.x-1,
         biasMM_Int = Int.x-1,
         biasMM_Time = Time.y-1)

biasREcor = corEst %>% 
  arrange(nClusters, nWithinCluster, corRE) %>% 
  mutate(biasLGC_corRE = corRE_est.x-corRE,
         biasMM_corRE = corRE_est.y-corRE)

# biasREcor %>% 
#   select(contains('bias')) %>% 
#   d3heatmap::d3heatmap(Rowv=F, Colv=F)
#          
# largely matched except for small samples where nlme has smaller intervals
save(fixedEffects, randomEffects, randomEffectsCor, 
     biasFE, biasRE, biasREcor,
     file='../data/growthvsMixed_EstResults.RData')

save(fixedEffects, randomEffects, randomEffectsCor, 
     biasFE, biasRE, biasREcor,
     file='SC and TR/mixedModels/growth_vs_mixed_files/growthvsMixed_EstResults.RData')

