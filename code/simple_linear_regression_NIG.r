N0 <- 1000
true.betas <- c(-1, 1, .5, -.5) 
P <- 4
X0 <- matrix(NA, nrow = N0, ncol = P)
for(i in 1:P) X0[, i] <- rnorm(N0)
sy <- 2
y0 <- rnorm(N0, mean = X0%*%true.betas, sd = sy)

summary(lm(y0 ~ -1 + X0))


library(rstan)
rstan_options(auto_write = TRUE)

rstan:::rstudio_stanc("stan/simple_linear_regression_NIG_prior.stan")

compiled.model.prior <- stan_model("stan/simple_linear_regression_NIG_prior.stan")

as <- .5
bs <- 2
vb <- 1.5
lm.data <- list(
  N_0 = N0,
  X_0 = X0,
  y_0 = y0,
  mu_beta = rep(0, P),
  lambda_0 = solve(vb * diag(P)),
  alpha0 = as,
  beta0 = bs,
  a_0 = .15
)

prior.lm <- sampling(compiled.model.prior, data = lm.data)
prior.lm
log.c_a0 <- bridgesampling::bridge_sampler(prior.lm)

get_mal_NIG_regression <- function(y0, X0, n0, mu0, lambda0, alpha0, beta0, a_0){
  P <- ncol(X0)
  if(length(mu0) != P) stop("mu0 is not the same dimension as X")
  Xstar <- sqrt(a_0) * X0
  ystar <- sqrt(a_0) * y0
  lambda_n <- t(Xstar)%*%Xstar + lambda0   
  mu_n <- solve(lambda_n) %*% (lambda0%*%mu0 + t(Xstar)%*%ystar)
  alpha_n <-  as + (n0*a_0)/2
  beta_n <- beta0 + .5 * ( t(ystar)%*%ystar + t(mu0)%*%lambda0%*%mu0 - t(mu_n)%*%lambda_n%*%mu_n )
  det0 <- det(lambda0)
  det1 <- det(lambda_n)   
  ans <- -(a_0 * n0/2) * log(2*pi) + .5 * (log(det0) - log(det1)) + (alpha0 * log(beta0) - alpha_n * log(beta_n) ) + (lgamma(alpha_n)-lgamma(alpha0))  
  return(ans)
}

invlambda0 <- solve(lm.data$lambda_0)
get_mal_NIG_regression(
  y0 = lm.data$y_0,
  X0 = lm.data$X_0,
  n0 = lm.data$N_0,
  mu0 = lm.data$mu_beta,
  lambda0 = invlambda0,
  alpha0 = lm.data$alpha0,
  beta0 = lm.data$beta0,
  a_0 = lm.data$a_0
)
log.c_a0

anaughts <- seq(0, 2, length.out = 100)
c_a0 <- function(x) {
  get_mal_NIG_regression(
    y0 = lm.data$y_0,
    X0 = lm.data$X_0,
    n0 = lm.data$N_0,
    mu0 = lm.data$mu_beta,
    lambda0 = invlambda0,
    alpha0 = lm.data$alpha0,
    beta0 = lm.data$beta0,
    a_0 = x
  )
}

MaLs <- sapply(anaughts, c_a0)
plot(anaughts, MaLs, xlab = expression(a[0]), main = "Normal inverse Gamma regression",
     ylab = expression(log(c(a[0]))), type = "l")
############
### Unnormalised posterior

### will draw from the same process for simplicity
N <- 100
X <- matrix(NA, nrow = N0, ncol = P)
for(i in 1:P) X[, i] <- rnorm(N)
y <-  rnorm(N, mean = X%*%true.betas, sd = sy)

compiled.model.unnorm.posterior <- stan_model("stan/simple_linear_regression_NIG_posterior_unnormalised.stan")

delta <- 1
nu <- 1
lm.data.forposterior <- list(
  N_0 = N0,
  X_0 = X0,
  y_0 = y0,
  mu_beta = rep(0, P),
  lambda_0 = solve(vb * diag(P)),
  alpha0 = as,
  beta0 = bs,
  delta = delta,
  nu = nu,
  X = X0,
  N = N0,
  y = y0
)

unnorm.posterior.lm <- sampling(compiled.model.unnorm.posterior, data = lm.data.forposterior)
unnorm.posterior.lm
pairs(unnorm.posterior.lm, pars = c("a_0", "beta"))

## Now the normalised version

compiled.model.norm.posterior <- stan_model("stan/simple_linear_regression_NIG_posterior_normalised.stan")

norm.posterior.lm <- sampling(compiled.model.norm.posterior, data = lm.data.forposterior)
norm.posterior.lm
pairs(norm.posterior.lm, pars = c("a_0", "beta"))

##### Analysis of p(a_0 | data)
a0.unnorm <- extract(unnorm.posterior.lm, 'a_0')$a_0
a0.norm <- extract(norm.posterior.lm, 'a_0')$a_0

a0.dt <- data.frame(a0 = c(a0.unnorm, a0.norm),
                    normalisation = c( rep("unnormalised", length(a0.unnorm)),
                                       rep("normalised", length(a0.norm)) )
)
library(ggplot2)
ggplot(a0.dt, aes(x = a0, fill = normalisation, colour = normalisation)) +
  geom_density() +
  stat_function(fun = function(x) dbeta(x, delta, nu),
                geom = "line", colour = "black", linetype = "longdash") +
  ggtitle("Regression") +
  facet_grid(normalisation~., scales = "free") +
  scale_y_continuous("Density", expand = c(0, 0)) +
  scale_x_continuous(expression(a[0])) +
theme_bw()