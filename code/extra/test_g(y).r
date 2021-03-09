## let L(D|theta) be a likelihood > 0 for all thera.
## let Y = [L(D|theta)]^2
## then the prior pi(theta) induces a measure on the space Y

### Binomial

n <- 10
x <- 3

ap <- 2
bp <- 2

p.samples <- rbeta(n = 1e6, shape1 = ap, shape2 = bp)

X <- dbinom(x = x, size = n, prob = p.samples)

hist(X, main = "Binomial")

Y <- X^2

hist(Y, main = "Binomial")

###############

#### Normal, fixed variance

x <- -1
s <- .1


mu.samples <- rnorm(n = 1E6)

X <- dnorm(x = x, mean = mu.samples, sd = s)

hist(X, main = "Gaussian")

Y <- X^2

hist(Y, main = "Gaussian")

###############

#### Beta with fixed shape1

psi <- .5
s1 <- 2

s2.samples <- rgamma(n = 1E6, shape = 2, rate = 1)

X <- dbeta(x = psi, shape1 = s1, shape2 = s2.samples)

hist(X, main = "Beta")

Y <- X^2

hist(Y, main = "Beta")

