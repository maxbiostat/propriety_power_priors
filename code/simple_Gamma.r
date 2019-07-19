true.alpha <- 2.4
true.beta <- .5
N0 <- 500
y0 <- rgamma(N0, shape = true.alpha, rate = true.beta)

library(rstan)
rstan_options(auto_write = TRUE)

options(mc.cores = 4)
rstan:::rstudio_stanc("stan/simple_gamma_prior.stan")

#### Sampling from the "prior"
compiled.model.prior <- stan_model("stan/simple_gamma_prior.stan")

etaA <- .2
nuA <- .2
etaB <- 1
nuB <- 2
ga.data <- list(
  N0 = N0,
  y0 = y0,
  eta_a = etaA,
  nu_a = nuA,
  eta_b = etaB,
  nu_b = nuB,
  a_0 = 9
)

prior.gamma <- sampling(compiled.model.prior, data = ga.data)
prior.gamma
pairs(prior.gamma, pars = c("alpha", "beta"))

log.c_a0 <- bridgesampling::bridge_sampler(prior.gamma)

get_c_a0_gamma <- function(y0, n0, eta.a, nu.a, eta.b, nu.b, a_0){
  logP <- sum(log(y0))
  S <- sum(y0)
  Delta <- nu.a - a_0*logP
  gk <- function(alpha){
    term1 <- lgamma(a_0*n0*alpha + eta.b) - a_0*n0*lgamma(alpha)
    term2 <- (eta.a-1)*log(alpha)-Delta*alpha
    term3 <- (a_0*n0*alpha + eta.b)*log(a_0*S + nu.b)
    require(Rmpfr)
    ldens <- mpfr(term1 + term2 -term3, 200) 
    return(exp(ldens))
  }
  # gk <- Vectorize(gk)
  logKInt <- log(Rmpfr::integrateR(gk, Rmpfr::mpfr(0.00000001, 200), 500, abs.tol = 1e-10)$value )
  # gk <- function(alpha, args){
  #   mock <- args[1] + args[2]
  #   term1 <- lgamma(a_0*n0*alpha + eta.b) - a_0*n0*lgamma(alpha)
  #   term2 <- (eta.a-1)*log(alpha)-Delta*alpha
  #   term3 <- (a_0*n0*alpha + eta.b)*log(a_0*S + nu.b)
  #   ldens <- term1 + term2 -term3
  #   return(mock * ldens)
  # }
  # gk <- Vectorize(gk)
  # require(reticulate)
  # lint <- import("lintegrate", convert = FALSE)
  # intt <- lint$lcquad(py_func(gk), a=r_to_py(0), b=r_to_py(10),
  #                     args=c(1, 0), epsabs = r_to_py(1e-20))
  # logKInt <- reticulate::py_to_r(intt[0])
  # ans <- -a_0*logP + logKInt
  ans <- -a_0*logP + logKInt 
  return(ans)
}
system.time(
  exact.log.c_a0 <- get_c_a0_gamma(
    y0 = ga.data$y0,
    n0 = ga.data$N0,
    eta.a = ga.data$eta_a,
    nu.a = ga.data$nu_a,
    eta.b = ga.data$eta_b,
    nu.b = ga.data$nu_b,
    a_0 = ga.data$a_0
  )
)

Rmpfr::asNumeric(exact.log.c_a0)
log.c_a0


# anaughts <- seq(0, 2, length.out = 100)
# c_a0 <- function(x) {
#   get_c_a0_gamma(
#     y0 = ga.data$y0,
#     n0 = ga.data$N0,
#     eta.a = ga.data$eta_a,
#     nu.a = ga.data$nu_a,
#     eta.b = ga.data$eta_b,
#     nu.b = ga.data$nu_b,
#     a_0 = x
#   )
# } 
# MaLs <- unlist(parallel::mclapply(anaughts, c_a0, mc.cores = 4))
# numMals <- unlist(lapply(MaLs, Rmpfr::asNumeric))
# plot(anaughts, numMals, xlab = expression(a[0]), main = "Gamma",
#      ylab = expression(log(c(a[0]))), type = "l")

#####################
get_estimates_a0 <- function(x){
  ga.data <- list(
    N0 = N0,
    y0 = y0,
    eta_a = etaA,
    nu_a = nuA,
    eta_b = etaB,
    nu_b = nuB,
    a_0 = x
  )
  res <- sampling(compiled.model.prior, data = ga.data, refresh = 0)
  return(res)
}
K <- 25
maxA <- 2
a0s <- seq(0.05, maxA, length.out = K)
# a0s <- pracma::logseq(0.01, maxA, n = K)
# a0s <- c(seq(0.05, 1, length.out = K-3), seq(1.2, maxA, length.out = 3))
system.time(
  results <- lapply(a0s, get_estimates_a0)  
)

### Sensitivity to a_0: parameter estimates
mean_pars <- data.frame(
  do.call(rbind,
          lapply(results, function(r) summary(r)$summary[1:2, "mean"])        
  )
) 
mean_pars$a0 <- a0s
round(mean_pars, 3)
#
sd_pars <- data.frame(
  do.call(rbind,
          lapply(results, function(r) summary(r)$summary[1:2, "sd"])        
  )
) 
sd_pars$a0 <- a0s
round(sd_pars, 4)


library(ggplot2)

ggplot(data = mean_pars, aes(x = a0, y = alpha)) +
  geom_line() +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  scale_y_continuous(expression(Mean(alpha)), expand = c(0, 0),
                     limits = c(true.alpha - .1 * true.alpha,  true.alpha + .1 * true.alpha ) ) +
  theme_bw(base_size = 20) +
  geom_hline(yintercept = true.alpha, linetype = "dashed") +
  NULL

ggplot(data = mean_pars, aes(x = a0, y = beta)) +
  geom_line() +
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  scale_y_continuous(expression(Mean(beta)), expand = c(0, 0),
                     limits = c(true.beta - .1 * true.beta,  true.beta + .1 * true.beta ) ) +
  theme_bw(base_size = 20) +
  geom_hline(yintercept = true.beta, linetype = "dashed") +
  NULL

ggplot(data = sd_pars, aes(x = a0, y = alpha)) +
  geom_line()+
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  scale_y_continuous(expression(SD(alpha)), expand = c(0, 0)) +
  theme_bw(base_size = 20) +
  NULL

ggplot(data = sd_pars, aes(x = a0, y = beta)) +
  geom_line()+
  scale_x_continuous(expression(a[0]), expand = c(0, 0)) +
  scale_y_continuous(expression(SD(beta)), expand = c(0, 0)) +
  theme_bw(base_size = 20) +
  NULL