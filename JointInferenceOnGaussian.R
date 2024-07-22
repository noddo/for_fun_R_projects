#connjugate prior:

# Estimate mu and sigma^2 from gaussian using priors
# Posterior predictive distribution for a new observation.

# mu ~ N(mu_0,n_0*tau)
# 1/sigma = tau ~ Gamma(alpha,beta)


# sigma|data ~ gamma(alpha + n/2 , beta + 1/2 (n-1)s^2 + [n*n_0/2(n + n_0)] * (xbar - mu_0)^2)
# mu|sigma, data ~ N( [tau*n * xbar + tau*n_0 * mu_0] / (n*tau + n_0*tau),(n*tau + n_0*tau))

#metaparameters:
mu_0 = 1.9
n_0 = 1
alpha = .01
beta = .01

#data (midge data)
x = c(1.64, 1.70, 1.72, 1.74, 1.82, 1.82, 1.82, 1.90, 2.08)

n = length(x)
xbar = mean(x)
s2 = (n-1)*var(x)

#plot prior distributions
tau = rgamma(100, shape = alpha, scale = beta)
mu_0 = rnorm(100,mu_0,n_0*1/tau)

#posterior distributions