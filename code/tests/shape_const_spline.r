smoothWithDer <- function(t, y, d, m = 3,
                          Hstar = c(3, 0.2, 0.1)^2, sigma2eta = 1.0^2) {
  
  ## define the SSM matrices, depending on 'delta_k' or on 'd_k'
  Tfun <- function(delta) {
    mat <-  matrix(0, nrow = m, ncol = m)
    for (i in 0:(m-1)) {
      mat[col(mat) == row(mat) + i] <- delta^i / gamma(i + 1)
    }
    mat
  }
  Qfun <- function(delta) {
    im <- (m - 1):0
    x <- delta^im / gamma(im + 1)
    mat <- outer(X = x, Y = x, FUN = "*")
    im2 <- outer(im, im, FUN = "+")
    sigma2eta * mat * delta / (im2 + 1) 
  }
  Zfun <-  function(d) {
    Z <- matrix(0.0, nrow = 1, ncol = m)
    Z[1, d + 1] <- 1.0
    Z
  }
  Hfun <- function(d) ifelse(d >= 0, Hstar[d + 1], 0.0)
  Rfun <- function() diag(x = 1.0, nrow = m)
  
  ## define arrays by stacking the SSM matrices. We need one more
  ## 'delta' at the end of the series
  n <- length(t)
  delta <-  diff(t)
  delta <- c(delta, mean(delta))
  
  Ta <- Qa <- array(0.0, dim = c(m, m, n))
  Za <- array(0.0, dim = c(1, m, n))
  Ha <- array(0.0, dim = c(1, 1, n))
  Ra <-  array(0.0, dim = c(m, m, n))
  
  for (k in 1:n) {
    Ta[ , , k] <- Tfun(delta[k])
    Qa[ , , k] <- Qfun(delta[k])
    Za[ , , k] <- Zfun(d[k])
    Ha[ , , k] <- Hfun(d[k])
    Ra[ , , k] <- Rfun()
  }
  
  require(KFAS)
  ## define the SSM and perform Kalman Filtering and smoothing
  mod <- SSModel(y ~ SSMcustom(Z = Za, T = Ta, R = Ra, Q = Qa, n = n,
                               P1 = matrix(0, nrow = m, ncol = m),
                               P1inf = diag(1.0, nrow = m), 
                               state_names = paste0("d", 0:(m-1))) - 1)
  out <- KFS(mod, smoothing = "state")
  list(t = t, filtered = out$att, smoothed = out$alphahat)
  
}
Hst <- c(10, 1, 1)^2
siget <- (.01)^2

a0s <- c(d15$a0, d15$a0)
res1 <- smoothWithDer(t = a0s,
                      y = c(d15$lc_a0, d15$deriv_lc),
                      d = c(rep(0, nrow(d15)), c(-1, rep(1, nrow(d15)-1)) ), m = 5,
                      Hstar = Hst, sigma2eta = siget)

curve(l_a0, 0, 10)
curve(l_a0_p, 0, 10, col = 2, add =  TRUE)
matplot(res1$t, res1$smoothed[, 1:2], pch = 16, cex = 0.7, ylab = "", xlab = "", add = TRUE)

curve(l_a0)
curve(l_a0_p, col = 2, add =  TRUE)
matplot(res1$t, res1$smoothed[, 1:2],
        pch = 16, cex = 0.7, ylab = "", xlab = "", add = TRUE)
