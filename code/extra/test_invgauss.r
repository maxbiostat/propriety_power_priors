true.lambda <- 2.4
true.mu <- .000000005
N0 <- 2
y0 <- statmod::rinvgauss(N0, mean = true.mu, shape = true.lambda)

a_0 <- .9

log(prod(statmod::dinvgauss(y0, mean = true.mu, shape = true.lambda)^a_0) )
a_0* sum(statmod::dinvgauss(y0, mean = true.mu, shape = true.lambda, log = TRUE))
sum(a_0 * statmod::dinvgauss(y0, mean = true.mu, shape = true.lambda, log = TRUE))

logP <- sum(log(y0))
S <- sum(y0)
Sprime <- sum(1/y0)

a_0 * N0 * .5 * ( log(true.lambda) - log(2*pi) ) -3/2 * a_0 * logP  + (a_0*N0*true.lambda)/true.mu -(a_0*S*true.lambda)/(2*true.mu^2)-(a_0*Sprime*true.lambda)/2
a_0 * N0 * .5 * ( log(true.lambda) - log(2*pi) ) -3/2 * a_0 * logP  + true.lambda * ( (a_0*N0)/true.mu -(a_0*S)/(2*true.mu^2)-(a_0*Sprime)/2   ) 
a_0 * N0 * .5 * ( log(true.lambda) - log(2*pi) ) -3/2 * a_0 * logP  - ( -(a_0*N0)/true.mu +(a_0*S)/(2*true.mu^2) + (a_0*Sprime)/2 ) * true.lambda 

-(a_0*N0)/true.mu +(a_0*S)/(2*true.mu^2) + (a_0*Sprime)/2
