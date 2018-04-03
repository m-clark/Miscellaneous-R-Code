#' Nelder Mead algorithm
#'
#' @param f  function to optimize, must return a scalar score and operate over a
#'   numpy array of the same dimensions as x_start
#' @param x_start initial position
#' @param step look-around radius in initial step
#' @param no_improve_thr See no_improv_break
#' @param no_improv_break break after no_improv_break iterations with an
#'   improvement lower than no_improv_thr
#' @param max_iter always break after this number of iterations. Set it to 0 to
#'   loop indefinitely.
#' @param alpha parameters of the algorithm (see Wikipedia page for reference)
#' @param gamma parameters of the algorithm (see Wikipedia page for reference)
#' @param rho parameters of the algorithm (see Wikipedia page for reference)
#' @param sigma parameters of the algorithm (see Wikipedia page for reference)
#' @param verbose Print current iteration?
#'
#' @details This is based on the pure Python implementation by François Chollet
#'   found at https://github.com/fchollet/nelder-mead. To be honest, this is
#'   just an academic exercise.  I'm not sure how much one would use NM in
#'   practice. In my limited experience BFGS would be faster, more accurate, and
#'   less sensitive to starting values for the types of problems I've played
#'   around with. Others who actually spend their time researching such things
#'   seem to agree.
#'
#'   There were two issues on
#'   (GitHub)[https://github.com/fchollet/nelder-mead/issues/2] regarding this,
#'   and I've implemented them with notes.  The current code is not very R-like,
#'   as the goal was to keep more similar to the original Python, which used a
#'   list approach (e.g. I'd never felt a need to use breaks and nexts in R in
#'   my entire decade-plus use of it).  I'll provide a cleaner version at some
#'   point.
#'
#' @return best parameter array, best score
#' @export
#'
#' @examples
nelder_mead = function(f, 
                       x_start,
                       step=0.1, 
                       no_improve_thr=1e-12,
                       no_improv_break=10, 
                       max_iter=0,
                       alpha=1, 
                       gamma=2, 
                       rho=0.5, 
                       sigma=0.5, 
                       verbose=FALSE) {
  # init
  dim = length(x_start)
  prev_best = f(x_start)
  no_improv = 0
  res = list(list(x_start=x_start, prev_best=prev_best))
  
  
  for (i in 1:dim) {
    x = x_start
    x[i] = x[i] + step
    score = f(x)
    res = append(res, list(list(x_start=x, prev_best=score)))
  }
  
  # simplex iter
  iters = 0
  while (TRUE) {
    # order
    idx = sapply(res, `[[`, 2)
    res = res[order(idx)]   # ascending order
    best = res[[1]][[2]]
    
    # break after max_iter
    if (max_iter > 0 & iters >= max_iter) return(res[[1]])
    iters = iters + 1
    
    # break after no_improv_break iterations with no improvement
    if (verbose) message(paste('...best so far:', best))
    
    if (best < (prev_best - no_improve_thr)) {
      no_improv = 0
      prev_best = best
    } else {
      no_improv = no_improv + 1
    }
    
    if (no_improv >= no_improv_break) return(res[[1]])
    
    # centroid
    x0 = rep(0, dim)
    for (tup in 1:(length(res)-1)){
      for (i in 1:dim){
        x0[i] = x0[i] + res[[tup]][[1]][i] / (length(res)-1)
      }
    }
    
   # reflection
   xr = x0 + alpha*(x0 - res[[length(res)]][[1]])
   rscore = f(xr)
   if (res[[1]][[2]] <= rscore & 
       rscore < res[[length(res)-1]][[2]]) {
     res[[length(res)]] = list(xr, rscore)
     next
   }
     
   # expansion
   if (rscore < res[[1]][[2]]) {
     # xe = x0 + gamma*(x0 - res[[length(res)]][[1]])   # issue with this
     xe = x0 + gamma*(xr - x0)   
     escore = f(xe)
     if (escore < rscore) {
       res[[length(res)]] = list(xe, escore)
       next
     } else {
       res[[length(res)]] = list(xr, rscore)
       next
     }
   }
   
   # contraction
   # xc = x0 + rho*(x0 - res[[length(res)]][[1]])  # issue with wiki consistency for rho values (and optim)
   xc = x0 + rho*(res[[length(res)]][[1]] - x0)
   cscore = f(xc)
   if (cscore < res[[length(res)]][[2]]) {
     res[[length(res)]] = list(xc, cscore)
     next
   }
   
   # reduction
   x1 = res[[1]][[1]]
   nres = list()
   for (tup in res) {
     redx = x1 + sigma*(tup[[1]] - x1)
     score = f(redx)
     nres = append(nres, list(list(redx, score)))
   }
   res = nres
  }
}




# Example -----------------------------------------------------------------


f = function(x) {
  sin(x[1]) *cos(x[2]) * (1 / (abs(x[3]) + 1))
}
# debug(nelder_mead) 
nelder_mead(f, 
            c(0,0,0), 
            max_iter = 1000, 
            no_improve_thr = 1e-12)
optimx::optimx(c(0,0,0), f, 
      method = "Nelder-Mead",
      control = list(alpha=1, 
                     gamma=2,
                     beta=0.5,     #rho
                     maxit=1000,
                     reltol=1e-12)  
      )



# A Regression Model ------------------------------------------------------

set.seed(8675309)
N = 500
npreds = 5
X = cbind(1, matrix(rnorm(N*npreds), ncol=npreds))
beta = runif(ncol(X), -1, 1)
y = X %*% beta + rnorm(nrow(X))


f = function(b) {
  crossprod(y - X %*% b)[,1]  # if using optimx need scalar
}


lm.fit(X, y)$coef

nelder_mead(f, 
            runif(ncol(X)), 
            max_iter = 2000,
            no_improve_thr = 1e-12)

# run model via least squares loss
opt_out = optimx::optimx(
  runif(ncol(X)),
  fn = f,  # model function
  method = 'Nelder-Mead',
  control = list(alpha=1, 
                 gamma=2,
                 beta=0.5,     #rho
                 maxit=2000,
                 reltol = 1e-12)
)
opt_out

