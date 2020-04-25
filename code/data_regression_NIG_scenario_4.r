set.seed(666)

N_0 <- 1000
P <- 100
true.betas <- sample(c(-1, -1/2, 1/2, 1), P, replace = TRUE)
X_0 <- matrix(NA, nrow = N_0, ncol = P)
for(i in 1:P) X_0[, i] <- rnorm(N_0)
sy <- 2
y_0 <- rnorm(N_0, mean = X_0%*%true.betas, sd = sy)

as <- .5
bs <- 2
vb <- 1.5

### will draw from the same process for simplicity
N <- 100
X <- matrix(NA, nrow = N, ncol = P)
for(i in 1:P) X[, i] <- rnorm(N)
y <- rnorm(N, mean = X%*%true.betas, sd = sy)


nu <- 1
eta <- 1
