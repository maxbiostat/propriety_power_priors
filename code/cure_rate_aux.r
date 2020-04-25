# rpois_trunc <- function(n, lambda, lb = 0){
#   thres <- ppois(lb,  lambda = lambda)
#   us <- runif(n, min = thres, max = 1)
#   samples <- qpois(us, lambda = lambda)
#   return(samples)
# }
generate_one_data_point <- function(theta, alpha, lambda){
  # N <- rpois_trunc(1, lambda = theta) ## generate from a zero-truncated Poisson
  N <- rpois(1, lambda = theta) ## generate from a zero-truncated Poisson
  if(N > 0){
    ## convert between parametrisations
    a <- alpha
    b <- exp(-lambda/alpha)
    Z <- rweibull(N, shape = a, scale = b)
  }else{
    Z <- Inf
  }
  return(min(Z))
}
#
r_cure_rate <- function(n, theta, alpha, lambda){
  ys <- sapply(1:n, function(i) generate_one_data_point(theta = theta, alpha = alpha, lambda = lambda))
  return(ys)
}
#
generate_n_data_points <- function(n, beta, alpha, lambda, X, tcs){
  theta <-  exp(X%*%beta)
  samples <- rep(NA, n)
  for(i in 1:n){
    samples[i] <- r_cure_rate(n = 1, theta = theta[i], alpha = alpha, lambda = lambda)
  }
  deltas <- as.numeric(samples < tcs)
  samples_cens <- ifelse(deltas == 1, samples, tcs)
  return(
    list(
      y = samples,
      y_censored = samples_cens,
      delta_censored = deltas
    )
  )
}