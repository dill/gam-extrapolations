# how do extrapolations look in 2D?

library(mvtnorm)
library(mgcv)
#library(rgl)

set.seed(123)

# simulate the data
gsim <- function(x1, x2, scale){

  f <- dmvnorm(cbind(x1, x2), c(0, -1), matrix(c(0.5, 0.2, 0.2, 0.5), 2, 2)) +
        dmvnorm(cbind(x1, x2), c(0.5, 0.5), matrix(c(0.1, 0, 0, 0.1), 2, 2))

  #f  <-  5*plogis(2*x1)
  #f  <-  2*x1

  y <- f+rnorm(length(x1), 0, scale)

  data <- data.frame(y=y, x1=x1, x2=x2, f=f)
  return(data)
}

# generate some data
n <- 200
dat <- gsim(runif(n, -1, 1), runif(n, -1, 1), 2)


# SNR
cat("SNR=", with(dat, cor(y, f)), "\n")



# predict and plot function
pred_and_plot <- function(model){

  # prediction grid
  x1g <- seq(-10, 10, len=100)
  x2g <- seq(-10, 10, len=100)
  xgrid <- expand.grid(x1=x1g, x2=x2g)

  zlims <- range(dat$f) + c(-2, 2)
  preds <- predict(model, xgrid, se.fit=TRUE)
  pred.mat <- matrix(preds$fit, 100, 100)
  pred.mat_p <- matrix(preds$fit, 100, 100)+2*matrix(preds$se.fit,100,100)
  pred.mat_m <- matrix(preds$fit, 100, 100)-2*matrix(preds$se.fit,100,100)

  persp3d(z=pred.mat, x=x1g, y=x2g, zlim=zlims, col="green")
  persp3d(z=pred.mat_m, x=x1g, y=x2g, zlim=zlims, col="blue", add=TRUE)
  persp3d(z=pred.mat_p, x=x1g, y=x2g, zlim=zlims, col="red", add=TRUE)

}


# tprs fit
b_tprs <- gam(y~s(x1, x2), data=dat,  method="REML")
#pred_and_plot(b_tprs)

### Duchon fit (s=0.5)
b_duchon <- gam(y~s(x1, x2, bs="ds", m=c(1,0.5)), data=dat, method="REML")
#pred_and_plot(b_duchon)

### use Eric's trick -- extra zero-weight data
#xtra <- expand.grid(x1=c(-10, 10), x2=c(-10, 10))
#xtra$y <- 1e6
#dat2  <-  rbind(dat[,c("x1", "x2", "y")],
#                xtra)
#b_weights <- gam(y~s(x1, x2, bs="ds", m=c(1, 0.5)), data=dat2, method="REML",
#                 weights=c(rep(1, nrow(dat)), rep(0, nrow(xtra)) ))
#pred_and_plot(b_weights)
#
#
#
### Simon's trick -- add in the x1 and x2 terms again
#b_ex <- gam(y~s(x1, x2, bs="ds", m=c(1,0.5), k=30) + x1 + x2, data=dat, method="REML")
#pred_and_plot(b_ex)
#
