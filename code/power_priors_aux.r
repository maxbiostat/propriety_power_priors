get_band <- function(fit, xpred){
  inipred <- predict(fit, newdata = data.frame(a0 = xpred),  se.fit = TRUE)
  ans <- data.frame(
    mean = inipred$fit,
    lwr =  inipred$fit - (2 * inipred$se.fit),
    upr =  inipred$fit + (2 * inipred$se.fit)
  )
  return(ans)
}