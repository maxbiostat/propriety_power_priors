true.mu  <- -100
true.sigma <- .1
N0 <- 2000
y0 <- rnorm(N0, mean = true.mu, sd = true.sigma)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)


mu_0 <- -99
kappa0 <- 200
alpha_0 <- 5
beta_0 <- 5

###
get_c_a0_gaussian <- function(y0, n0, alpha0, beta0, m0, k0, a_0){
  nstar <- a_0 * n0
  ybar <- mean(y0)
  s <- mean( (y0-ybar)^2 )
  kappa_n <- k0 + nstar  
  alpha_n <- alpha0 + nstar/2
  beta_n <- beta0 + .5 * (nstar * s +  (k0 * nstar * (ybar - m0)^2 )/kappa_n ) 
  ans <- lgamma(alpha_n)-lgamma(alpha0)
  ans <- ans + alpha0 * log(beta0) - alpha_n * log(beta_n)
  ans <- ans + .5 *( log(k0) - log(kappa_n) )-nstar/2 * log(2*pi)
  return(ans)
}
#
c_a0 <- function(x) {
  get_c_a0_gaussian(
    y0 = y0,
    n0 = N0,
    alpha0 = alpha_0,
    beta0 = beta_0,
    m0 = mu_0,
    k0 = kappa0,
    a_0 = x
  )
}
c_a0 <- Vectorize(c_a0)

################

ca0.data <- read.csv("../../data/Gaussian_logCA0.csv")


###########################################

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

gp.model <- stan_model(file = "gp_2nd_deriv.stan")
expose_stan_functions(gp.model)

# 
# library(brms)
# m0 <- brm(lc_a0 ~ gp(a0), data = ca0.data)

K <- 10
maxA <- max(ca0.data$a0) + 2
forpred <- seq(0, maxA, length.out = K)

# K <- nrow(ca0.data)
# forpred <- ca0.data$a0

gp.data <- list(
  N = nrow(ca0.data),
  x = ca0.data$a0,
  y = ca0.data$lc_a0,
  K = K,
  m = rep(1, K),
  s = forpred,
  m = rep(1, K),
  epsilon = 1E-8,
  N_pred = length(forpred),
  x_pred = forpred
)

# M <- full_cov_mat(x = gp.data$x, s = gp.data$s, alpha = 1, rho = .5, epsilon = gp.data$epsilon)
# fields::image.plot(M)

opt <- optimizing(gp.model, data = gp.data)
opt

plot(gp.data$x_pred, opt$par[grep("f_pred", names(opt$par))])

f_means <-  opt$par[grep("f\\[", names(opt$par))][1:gp.data$N]
sigma.hat <- opt$par["sigma"]
samps <- matrix(NA, ncol = gp.data$N, nrow = 1E6)
for(i in 1:gp.data$N) samps[, i] = rnorm(1E6, mean = f_means[i] , sd = sigma.hat)
plot(gp.data$x, colMeans(samps), type = "l")
points(lc_a0 ~ a0, ca0.data, add = TRUE, pch = 16, col = 2)

chain <- sampling(gp.model, data = gp.data, control = list(adapt_delta = .95, max_treedepth = 10))
chain
stan_trace(chain,pars = c("rho", "alpha", "sigma"))
pairs(chain, pars = c("rho", "alpha", "sigma"))
# 
preds <- extract(chain, 'f')$f
# 
plot(gp.data$x, colMeans(preds)[1:11], type = "l", xlab = expression(a[0]), ylab = expression(log(c(a[0]))))
points(lc_a0 ~ a0, ca0.data, add = TRUE, pch = 16, cex = 2)
