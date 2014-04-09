# ------- Functions for making quantile predictions one step ahead----------

# AR1 quantile regression
ar1_qreg <- function(x, pval){
  n <- length(x)
  fit <- rq(x[-1]~abs(x[-n]), tau=pval)
  
  sum(coefficients(fit)*c(1,abs(x[n])))
}

# HAR quantile regression
har_qreg <- function(x, pval, ff=abs){
  n <- length(x)
  nn <- n - 19
  y <- x[-(1:20)]
  xx <- ff(x)
  x1 <- xx[-c(1:19)]
  x2 <- mavg(xx[-(1:15)], 5)
  x3 <- mavg(xx, 20)
  fit <- rq(y~x1[-nn]+x2[-nn]+x3[-nn], tau=pval)
  
  sum(coefficients(fit)*c(1, x1[nn], x2[nn], x3[nn]))
}

# Historical standard deviation
stde <- function(x, pval){
  qnorm(pval, mean=mean(x), sd=sd(x))
}

# Simple historical simulation
hisi <- function(x, pval){
  quantile(x, pval)
}
# ------------------------------------------------------------------------------
