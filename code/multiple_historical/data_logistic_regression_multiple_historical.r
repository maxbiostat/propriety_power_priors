## Setup
set.seed(6666)
verbose <- FALSE
D <- 5
N_0 <- 1000
true.alpha <- 1.2
true.betas <- c(-1, 10, .5, -5) 
P <- 4
dataSets <- vector(D, mode = "list")

for (d in 1:D){
  ### will draw from the same process for simplicity
  X_0 <- matrix(NA, nrow = N_0, ncol = P)
  for(i in 1:P) X_0[, i] <- rnorm(N_0)
  mu_0 <- true.alpha  + X_0%*%true.betas
  y0 <- rbinom(N_0, size = 1, prob = arm::invlogit(mu_0)  )
  dataSets[[d]] <- data.frame(X_0, y0 = y0)
}
rm(X_0)
rm(mu_0)
rm(y0)

if(verbose) lapply(dataSets, function(x) summary(glm(y0 ~., data = x, family = "binomial")))

## Current data
N <- 100
X <- matrix(NA, nrow = N, ncol = P)
for(i in 1:P) X[, i] <- rnorm(N)
mu <- true.alpha  + X%*%true.betas
y <- rbinom(N, size = 1, prob = arm::invlogit(mu))

eta <- rep(1, D)
nu <- rep(1, D)
