set.seed(666)
true.lambda <- 2
true.mu <- .1
N0 <- 200
y0 <- statmod::rinvgauss(N0, mean = true.mu, shape = true.lambda)

a.l <- 1
b.l <- 1
a.m <- 2
b.m <- 2

N <- 100
y <-  statmod::rinvgauss(N, mean = true.mu, shape = true.lambda)

delta <- 1
nu <- 1