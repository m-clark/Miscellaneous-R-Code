# a function that will take an input n and generate n x y and z values. The
# correlation of interest is between x and y, when either or both are a ratio
# based on z.

spurcorr = function(n){
  x = rpois(n, 10)
  y = rpois(n, 10)
  z = runif(n, 50, 100)
  
  cbind(N=n, origCorr=cor(x,y), xzyzCorr=cor(x/z, y/z), xzyCorr=cor(x/z,y), xyzCorr=cor(x,y/z))
}

# replicate at various sample sizes
out = replicate(1000, sapply(c(10, 100, 500, 1000, 2500), spurcorr))

origcorrs = apply(out, 3, function(mat) mat[2,]) ;rownames(origcorrs) = paste(c(10, 100, 500, 1000, 2500))
spurrcorrsBoth = apply(out, 3, function(mat) mat[3,]);rownames(spurrcorrsBoth) = paste(c(10, 100, 500, 1000, 2500))
spurrcorrsXonly = apply(out, 3, function(mat) mat[4,]);rownames(spurrcorrsXonly) = paste(c(10, 100, 500, 1000, 2500))
spurrcorrsYonly = apply(out, 3, function(mat) mat[5,]);rownames(spurrcorrsYonly) = paste(c(10, 100, 500, 1000, 2500))

data.frame(origZero=rowMeans(origcorrs), Xratio=rowMeans(spurrcorrsXonly), 
           Yratio=rowMeans(spurrcorrsYonly), Both=rowMeans(spurrcorrsBoth))


