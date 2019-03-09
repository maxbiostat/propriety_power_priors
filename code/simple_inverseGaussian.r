true.lambda <- 2.4
true.mu <- .05
N0 <- 100
y0 <- statmod::rinvgauss(N0, mean = true.mu, shape = true.lambda)

library(rstan)
rstan_options(auto_write = TRUE)

options(mc.cores = 4)
rstan:::rstudio_stanc("stan/simple_inverse_gaussian.stan")

#### Sampling from the "prior"
compiled.model.prior <- stan_model("stan/simple_inverse_gaussian.stan")

a.l <- .2
b.l <- .2
a.m <- 1
b.m <- 2
ig.data <- list(
  N0 = N0,
  y0 = y0,
  alpha_l = a.l,
  beta_l = b.l,
  alpha_m = a.m,
  beta_m = b.m,
  a_0 = .2
)

prior.invgauss <- sampling(compiled.model.prior, data = ig.data)
prior.invgauss
pairs(prior.invgauss, pars = c("mu", "lambda"))

log.c_a0 <- bridgesampling::bridge_sampler(prior.invgauss)
log.c_a0$logml

get_c_a0_invgauss_brute_force <- function(x){
  ig.data <- list(
    N0 = N0,
    y0 = y0,
    alpha_l = a.l,
    beta_l = b.l,
    alpha_m = a.m,
    beta_m = b.m,
    a_0 = x
  )
  prior.invgauss <- sampling(compiled.model.prior, data = ig.data)
  return(bridgesampling::bridge_sampler(prior.invgauss)$logml)
}

a0s <- c(.1, .2, .5, .75, .9)
mls <- sapply(a0s, get_c_a0_invgauss_brute_force)

plot(mls ~ a0s)

get_c_a0_invgauss <- function(y0, n0, alpha.l, beta.l, alpha.m, beta.m, a_0){
  logP <- sum(log(y0))
  S <- sum(y0)
  Sprime <- sum(1/y0)
  A <- (a_0*n0)/2 + alpha.l
  lconst <- -0.5 *n0*a_0*log(2*pi) - (3*a_0/2)*logP + lgamma(A)
  ###
  require(Rmpfr)
  igk <- function(mu){
    Delta <- -(a_0*N0)/mu +(a_0*S)/(2*mu^2) + (a_0*Sprime)/2  + beta.l
    ldens <- mpfr(-A * log(Delta) + (alpha.m-1)*log(mu) -beta.m*mu, 200)
    return(exp(ldens))  
  }
  logKInt <- log( integrateR(igk, 0.000001, 500, ord = 20)$value )
  ###
  ans <- lconst + logKInt
  return(ans)
}
system.time(
  exact.log.c_a0 <- get_c_a0_invgauss(
    y0 = ig.data$y0,
    n0 = ig.data$N0,
    alpha.l = ig.data$alpha_l,
    beta.l = ig.data$beta_l,
    alpha.m = ig.data$alpha_m,
    beta.m = ig.data$beta_m,
    a_0 = ig.data$a_0
  )
)
exact.log.c_a0
log.c_a0

anaughts <- seq(0, 2, length.out = 100)
c_a0 <- function(x) {
  get_c_a0_invgauss(
    y0 = ig.data$y0,
    n0 = ig.data$N0,
    alpha.l = ig.data$alpha_l,
    beta.l = ig.data$beta_l,
    alpha.m = ig.data$alpha_m,
    beta.m = ig.data$beta_m,
    a_0 = x
  )
} 
MaLs <- unlist(parallel::mclapply(anaughts, c_a0, mc.cores = 4))
MaLs <- unlist(lapply(MaLs, Rmpfr::asNumeric))
plot(anaughts, MaLs, xlab = expression(a[0]), main = "Inverse Gaussian",
     ylab = expression(log(c(a[0]))), type = "l")
