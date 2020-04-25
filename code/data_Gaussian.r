set.seed(666)

true.mu <- -.1
true.sigma <- .001
N_0 <- 50

y_0 <- rnorm(N_0, mean = true.mu, sd = true.sigma)

mu_0 <- 0
kappa_0 <- 5
alpha_0 <- 1
beta_0 <- 1

copy <- FALSE

if(copy){
  N <- N_0
  y <- y_0
}else{
  N <- 200
  y <- rnorm(N, mean = true.mu, sd = true.sigma)  
}

nu <- 1
eta <- 1
