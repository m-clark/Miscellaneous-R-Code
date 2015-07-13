# -*- coding: utf-8 -*-
"""
Created on Sat Jun 20 17:25:18 2015

@author: MC
"""
import numpy as np
import pandas as pd
import scipy as sp
import statsmodels.api as sm
import statsmodels.formula.api as smf

np.random.seed(1234)
N = 1000
nGroups = 40
nPerGroup = N//nGroups
x = np.random.normal(size=N)
group = np.repeat(np.arange(0, nGroups), nPerGroup)
ranEff = np.random.normal(scale=.5, size=nGroups)
coefs = np.array([2,.2]).reshape(2,1)
randInts = coefs[0] + ranEff[group]

X = np.array([randInts, x]).T

y = randInts + x*coefs[1] +  np.random.normal(scale=.75, size=N)
y.shape

df = pd.DataFrame({'x':x, 'group':group, 'y':y})

### Mode
mod0 = smf.mixedlm('y ~ x', data=df, groups=df['group']) 
mod = mod0.fit() 
print(mod.summary())  # intercept RE is in sd 

reModel = mod.random_effects['Intercept'][:]
sm.qqplot(reModel, sp.stats.norm, line='45', fit=True, scale=.25)
