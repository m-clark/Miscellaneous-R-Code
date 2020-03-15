% matlab from Statistical Modeling and Computation (2014 p 311).  See the 
% associated twofactorRE.R file for details.

function one_factor_re_loglike(mu, sigma2_mu, sigma2, y)
	[d ni] = size(y);
	Sigmai = sigma2*eye(ni) + sigma2_mu*ones(ni,ni);
	l = -(ni*d) / 2*log(2*pi) - d / 2*log(det(Sigmai));
	for i=1:d
	  yi = y(i, :)';
	  l = l - .5*(yi - mu)' * (Sigmai\(yi - mu));
	end
end


y = [22.6 20.5 20.8;
     22.6 21.2 20.5;
     17.3 16.2 16.6;
     21.4 23.7 23.2;
     20.9 22.2 22.6;
     14.5 10.5 12.3;
     20.8 19.1 21.3;
     17.4 18.6 18.6;
     25.1 24.8 24.9;
     14.9 16.3 16.6];


f = @(theta) -one_factor_re_loglike(theta(1), theta(2), theta(3), y);
ybar = mean(y, 2);
theta0 = [mean(ybar) var(ybar) mean(var(y, 0, 2))];
thetahat = fminsearch(f, theta0);
