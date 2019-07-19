source("Bernoulli_data.r")
###
bb.data <- list(
  N0 = N_0,
  y0 = y_0,
  c = cc,
  d = dd,
  a_0 = .4
)
get_c_a0_bernoulli <- function(y0, n0, cc, dd, a_0){
  ans <- lbeta(a_0 * y0 + cc, a_0 *(n0 -y0) + dd)
  return(ans)
}

c_a0 <- function(x) {
  get_c_a0_bernoulli(
    y0 = bb.data$y0,
    n0 = bb.data$N0,
    cc = bb.data$c,
    dd = bb.data$d,
    a_0 = x
  )
} 

constant_data <- read.csv("../data/Bernoulli_logCA0.csv")

library(mgcv)
fit.gam <- gam(lc_a0 ~ s(a0), data = constant_data)

K <- 1E4
pred_a0s <- seq(0, max(constant_data$a0), length.out = K)

a0_grid <- data.frame(a0 = pred_a0s, lc_pred = predict(fit.gam, newdata = data.frame(a0 = pred_a0s)))

curve(c_a0, 0, 5, xlab = expression(a[0]), ylab = expression(log(c(a[0]))), lwd = 3, lty = 2)
points(a0_grid$a0, a0_grid$lc_pred)

get_approx_lc <- function(x, grid){
  y <- grid[, 1]
  i <- which.min(abs(x-y))
  if(i != 1){
    x1 <- grid[i, 1]
    x2 <- grid[i + 1, 1]
    y1 <- grid[i, 2]
    y2 <- grid[i + 1, 2]
    ans <- y1 + (y2-y1) * (x-x1)/(x2-x1)
  }else{
    ans <- grid[i, 2]
  }
  return(ans)
}
app_f <- function(x) get_approx_lc(x, grid = a0_grid)
app_f <- Vectorize(app_f)

curve(c_a0, 0, 5, xlab = expression(a[0]), ylab = expression(log(c(a[0]))), lwd = 3, lty = 2)
curve(app_f, 0, 5, lwd = 3, col = 2, add = TRUE)
lines(a0_grid$a0, a0_grid$lc_pred, lwd = 2, lty = 2, col = 3)
