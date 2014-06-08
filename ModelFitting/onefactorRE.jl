############################################################################
### see onefactorRE.R for details and onefactorRE.m for matlab code. The ###
### commented part of the function reflects the R code.                  ###
############################################################################



#####################
### Main function ###
#####################
function sfran_loglike(par::Vector)
    d, ni = size(y)
    mu = par[1]
    sigma2_mu = par[2]
    sigma2 = par[3]

    Sigmai = sigma2*eye(ni) + sigma2_mu*ones(ni,ni)
    l = -(ni*d)/2*log(2*pi) - d/2*log(det(Sigmai))

    for i in 1:d
      yi = y[i,:]'
      l = l - .5*(yi-mu)'* (Sigmai\(yi-mu))
    end

#    l = zeros(10)

#    for i in 1:d
#        yi = y[i,:]'
#        l[i,:] =  .5* (yi-mu)' * (Sigmai\(yi-mu))
#    end
#    l =  -(ni*d)/2*log(2*pi) - d/2*log(det(Sigmai)) - sum(l)

    l = -l[1]  # having to do this line hurts
    return l
end


###################
### Data set up ###
###################
y = [22.6 20.5 20.8
     22.6 21.2 20.5
     17.3 16.2 16.6
     21.4 23.7 23.2
     20.9 22.2 22.6
     14.5 10.5 12.3
     20.8 19.1 21.3
     17.4 18.6 18.6
     25.1 24.8 24.9
     14.9 16.3 16.6]


################################
### Starting values and test ###
################################
mu0 = mean(y)
sigma2_mu0 = var(mean(y,2))
sigma20 = mean(var(y, 2))
theta0 = [mu0, sigma2_mu0, sigma20]

### test
sfran_loglike(theta0)


###########
### Run ###
###########
using Optim
res = optimize(sfran_loglike, theta0, method=:l_bfgs)


