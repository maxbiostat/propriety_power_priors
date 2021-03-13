set.seed(6666)

N_0 <- 1000
true.alpha <- 1.2
true.betas <- c(-1, 1, .5, -.5) 
P <- 4
X_0 <- matrix(NA, nrow = N_0, ncol = P)
for(i in 1:P) X_0[, i] <- rnorm(N_0)
mu_0 <- true.alpha  + X_0%*%true.betas
y_0 <- rbinom(N_0, size = 1, prob = arm::invlogit(mu_0)  )

### will draw from the same process for simplicity
N <- 100
X <- matrix(NA, nrow = N, ncol = P)
for(i in 1:P) X[, i] <- rnorm(N)
mu <- true.alpha  + X%*%true.betas
y <- rbinom(N, size = 1, prob = arm::invlogit(mu))

eta <- 1
nu <- 1

