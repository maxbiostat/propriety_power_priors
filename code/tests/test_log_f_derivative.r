f <- function(x) (x + 3)^2
f <- Vectorize(f)
lf <- function(x) log(f(x))
lf <- Vectorize(lf)

lf_prime_num <- function(x) numDeriv::grad(lf, x)

lf_prime_analytical <- function(x) (2*(x+3))/f(x)

curve(lf_prime_num(x), 0, 10)
curve(lf_prime_analytical(x), 0 , 10, lwd = 2, lty = 2, col = 2, add = TRUE)
