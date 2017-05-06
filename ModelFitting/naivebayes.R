# Naive bayes demo for binary data


# Initialization ----------------------------------------------------------

# generate some data

set.seed(123)
x = matrix(sample(0:1, 50, replace = T), ncol=5)
xf = data.frame(lapply(data.frame(x), factor))
y = sample(0:1, 10, prob=c(.25, .75), replace=T)


# use e1071 for comparison

library(e1071)
m = naiveBayes(xf, y)
m



# Base R approach  -----------------------------------------------------------

lapply(xf, function(var) t(prop.table(table(' '=var, y), margin=2)))

m

