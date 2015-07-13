# -*- coding: utf-8 -*-
"""
Created on Sat Jun 20 09:08:20 2015

@author: MC
"""

import numpy as np
import pandas as pd
import scipy as sp
from sklearn.preprocessing import scale
import statsmodels.api as sm
import statsmodels.formula.api as smf
# import seaborn as sea
import matplotlib.pyplot as plt

# Regression example using np arrays, np matrices and statsmodels; serves as a 
# quick reminder for me of some basic statistical modeling in python; 
# Unfortunately numpy doesn't deal with vectors as matrices, assumes row
# entry when creating arrays, and elementwise rather than matrix operations, 
# which combined means regularly transposing and reshaping, using built in 


### Data setup
np.random.seed(1234)

N = 100
x = np.random.normal(size=N)
X = np.array((np.ones(N), (x))).T
# compare to R: X = cbind(1, rnorm(N))

# true value of coefficients
beta = np.array((5, .5)).reshape(2,1)

# outcome of interest
y = np.dot(X,beta) + np.random.normal(scale=.5, size=N).reshape(N,1)

# example scale
scale(y)

### Inspect
X.shape
X[:5, 1]
beta.shape

plt.plot(x, y, 'ro')


### Regression models

## numpy linalg
np.linalg.lstsq(X, y)[0]

## Normal equations numpy array
# note that Python will finally get a matrix multiplier (@ symbol) in 3.5
xx = np.linalg.inv(np.dot(X.T, X))
xy = np.dot(X.T, y)

# alternative
xx = np.linalg.inv(X.T.dot(X))
xy = X.T.dot(y)

np.dot(xx, xy)

## Normal equations numpy matrix
# is only 2d but perhaps more intuitive and multiplication is matrix multiplication by default
Xmat = np.matrix(X) 
np.linalg.inv(Xmat.T * Xmat) *  np.dot(Xmat.T, y)

## statsmodels 
res0 = sm.OLS(y, X)
res = sm.OLS(y, X).fit()

res.summary()
res.aic
res.fittedvalues
res.rsquared_adj

## using scipy linregress for simple regression
slope, intercept, r_value, p_value, std_err = sp.stats.linregress(x, y[:,0])
[intercept, slope]

## pandas data frames with formula approach
df = pd.DataFrame(X)
df['y'] = y
df.columns = ['Intercept', 'x', 'y']
df.head()

res = smf.ols('y ~ x',  data=df).fit()
res.summary()

plt.plot(x, y, 'ro',
         x, res.fittedvalues)

