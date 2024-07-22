#connjugate prior:

# Estimate mu and sigma^2 from gaussian using priors
# Posterior predictive distribution for a new observation.

# mu ~ N(mu_0,sigma_0)
# 1/sigma = tau ~ Gamma(alpha,beta)


# sigma|data ~ gamma()
# mu|sigma, data ~ N( [tau*n* xbar + mu_0/sigma_0] / (n*tau + 1/sigma_0),(n*tau + 1/sigma_0))
