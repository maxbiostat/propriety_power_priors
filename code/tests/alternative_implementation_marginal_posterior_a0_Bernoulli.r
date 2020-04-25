get_l_a0_bernoulli <- function(y0, n0, cc, dd, a_0){
  ans <- lbeta(a_0 * y0 + cc, a_0 *(n0 -y0) + dd)-lbeta(cc, dd)
  return(ans)
}
posterior_a0_Bernoulli_2 <- function(a_0, y0, n0, y, n, cc, dd, delta, nu, log = FALSE){
  term1 <- get_l_a0_bernoulli(y0 = y0, n0 = n0, cc = cc, dd = dd, a_0 = a_0)
  term2 <- dbeta(a_0, shape1 = delta, shape2 = nu, log = TRUE)
  term3 <- lbeta(a_0 * y0 + y + cc - 1, a_0 * (n0 - y0) + dd + (n-y) - 1)
  ans <- -term1 + term2 + term3
  if(!log) ans <- exp(ans)
  return(ans)
}
post_a0_2 <- function(x) {
  posterior_a0_Bernoulli_2(a_0 = x, y0 = y_0, n0 = N_0,
                           y = y, n = N, cc = cc, dd = dd, nu = nu, delta = delta)
}
post_a0_2  <- Vectorize(post_a0_2) 
post_a0_2  <- Vectorize(post_a0_2) 
curve(post_a0_2)

K2 <- integrate(post_a0_2, 0, 1)$value
norm_post_a0_2 <- function(x) post_a0_2(x)/K2
norm_post_a0_2  <- Vectorize(norm_post_a0_2)
curve(norm_post_a0_2)