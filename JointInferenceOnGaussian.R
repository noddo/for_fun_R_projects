#connjugate prior:
# tau (precision = 1/sigma^2) ~ Gamma(v0/2,s20*v0/2) --> mu ~ N(mu0,tau) --> X1, X2,...Xn ~iid N(mu,tau)
# Estimate (point estimates and credible intervals) mu and tau from gaussian using priors
# Posterior predictive distribution for a new observation.

# mu ~ N(mu_0,n_0*tau) mu_0 is the mean of the prior sample and n_0 is the sample size of the prior observations.
# 1/sigma = tau ~ Gamma(v0/2,v0*Sigma0/2) - v0 is the prior sample size for the variance, sigma0 is the sample variance of prior obs.

#postierior distributions:
    # tau|data ~ gamma(vn/2, vn*sigma_n/2); vn = v0+n, sigma_n = 1/vn * [v0*sigma0 + (n-1)*s2 + (ybar - mu0)^2*k0*n/kn]
    # mu|sigma, data ~ N( (k0*mu0 + n*ybar)/kn ,1/tau*kn) ; kn = k0 + n

#metaparameters:
# sigma:
    v0 = 1
    sigma0 = 0.01
#mu:
    k0 = 1
    mu0 = 1.9

#data (midge data)
    x = c(1.64, 1.70, 1.72, 1.74, 1.82, 1.82, 1.82, 1.90, 2.08)

    n = length(x)
    xbar = mean(x)
    s2 = var(x)

#posterior params
#mu|data:
    kn = k0 + n
    mun = (k0*mu0 + n*xbar)/kn


# sigma|data:
    vn = v0 + n
    sigma_n = 1/vn * (v0*sigma0 + (n-1)*s2 + (xbar - mu0)^2*k0*n/kn)

#Sample posterior distributions
sigma = 1/rgamma(1000, mun/2 , sigma_n*mun/2)
mu = rnorm(1000,mun,sqrt(sigma_n/kn))

#credible intervals:
paste('mu 95% credible interval: [', quantile(mu,c(.025, .975)),']')
paste('sigma 95% credible interval: [', quantile(sigma,c(.025, .975)),']')


#posterior distributions
post.x = rnorm(1000,mu,sqrt(sigma))

sum(ifelse(post.x <= 1.6,1.0,0.0))/1000

