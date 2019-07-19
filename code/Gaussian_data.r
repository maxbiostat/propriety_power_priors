
set.seed(666)


true.mu <- -100
true.sigma <- .01
N_0 <- 2000

y_0 <- rnorm(N_0, mean = true.mu, sd = true.sigma)

mu_0 <- -99
kappa_0 <- 200
alpha_0 <- 5
beta_0 <- 5

N <- 50
y <- rnorm(N, mean = true.mu, sd = true.sigma)

delta <- 1
nu <- 1