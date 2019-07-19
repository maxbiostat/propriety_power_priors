# plot(anaughts, numDeriv::grad(c_a0, anaughts), type = "l")
# post.dens <- function(p, n, y, cc, dd, a_0){
#   exp( 
#     a_0 * dbinom(x = y, size = n, prob = p, log = TRUE) + dbeta(p, cc, dd, log = TRUE)
#   )
# }
# post.dens <- Vectorize(post.dens)
# 
# expec_ll <- function(y, N, cc, dd, a_0){
#   integrate(
#     function(p) post.dens(p, n = N, y = y, cc = cc, dd = dd, a_0 = a_0) * dbinom(x = y, size = N, prob = p, log = TRUE),
#     0, 1
#   )$value
# }
# 
# e_c0 <- function(x) exp(c_a0(x))
# expec_ll(y = y, N = N, cc = cc, dd = dd, a_0 = .5)
# numDeriv::grad(e_c0, .5)
