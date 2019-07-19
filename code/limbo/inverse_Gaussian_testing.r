invgaus_kernel_forint <- function(u, a, b, c, d, s, r){
  (a/u^2 - b/u + c)^-d * u^(s-1) * exp(-r * u)
}
invgaus_kernel_forint <- Vectorize(invgaus_kernel_forint)

invgaus_kernel_forint_2 <- function(u, a, b, c, d, s, r){
  (a/(c * u^2) - b/(c* u) + 1)^-d * u^(s-1) * exp(-r * u)
}
invgaus_kernel_forint_2 <- Vectorize(invgaus_kernel_forint_2)
########################
true.lambda <- 2.4
true.mu <- 10
a.l <- .2
b.l <- .2
a.m <- 1
b.m <- 2
N0 <- 2000
y0 <- statmod::rinvgauss(N0, mean = true.mu, shape = true.lambda)
a0 <- .3
########################
s <- sum(y0)
sp <- sum(1/y0)
########################
u <- 10
a <- a0*s/2
b <- a0*N0
c <- a0*sp/2 + b.l
d <- (a0*N0 + 2*a.l)/2
s <- a.m
r <- b.m   

integrate(function(x) invgaus_kernel_forint(x, a, b, c, d, s, r), 0, Inf)  
integrate(function(x) invgaus_kernel_forint_2(x, a, b, c, d, s, r), 0, Inf)  
-d * log(c) + log(integrate(function(x) invgaus_kernel_forint_2(x, a, b, c, d, s, r), 0, Inf)$value)
