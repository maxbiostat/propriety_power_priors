set.seed(666)

true.lambda <- 2
N_0 <- 200
y_0 <- rpois(N_0, lambda = true.lambda)

N <- 100
y <- rpois(N, lambda = true.lambda)

alpha_0 <- 2
beta_0 <- 2

nu <- 1
eta <- 1
