############################################################################
### see twofactorRE.R for details and twofactorRE.m for matlab code.     ###
############################################################################



function sfran2_loglike(par::Vector)
    n = length(y)
    mu = par[1]
    
    sigma2_alpha = exp(par[2])
    sigma2_gamma = exp(par[3])
    sigma2 = exp(par[4])

    Sigma = sigma2*eye(n) + sigma2_alpha*(Xalpha * Xalpha') + sigma2_gamma * (Xgamma * Xgamma')

    l = -n/2*log(2*pi) - sum(log(diag(chol(Sigma)))) - .5*(y-mu)' * (Sigma\(y-mu))

    l = -l[1]
    return l
end


##################
### Data setup ###
##################
y = [1.39,1.29,1.12,1.16,1.52,1.62,1.88,1.87,1.24,1.18,
     .95,.96,.82,.92,1.18,1.20,1.47,1.41,1.57,1.65]

# See the R file for a conceptual data representation

Xalpha = kron(eye(5),  ones(4,1))

Xgamma = kron(eye(10), ones(2,1))


################################
### Starting values and test ###
################################
yhat = mean(reshape(y, 4, 5), 1)

theta0 = [mean(y), log(var(yhat)), log(var(y)/3), log(var(y)/3)]

sfran2_loglike(theta0)


###########
### Run ###
###########
using Optim
res = optimize(sfran2_loglike, theta0, method = :l_bfgs)
exp([-2.92,-3.44, -6.079])


