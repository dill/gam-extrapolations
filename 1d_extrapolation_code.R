sim_post_curves= function(model, pdata, n=25){
  lp <- predict(model, newdata = pdata, type = "lpmatrix")
  coefs <- coef(model)
  vc <- vcov(model)
  sim_coef = mvrnorm(n, coefs, vc)
  fits = lp %*% t(sim_coef)
  return(fits)
}


library(mgcv)
library(MASS)
library(scales)
library(microbenchmark)

ff <- function(x2, func=1) {
  if(func==1){
    0.2 * x2^11 * (10 * (1 - x2))^6 + 10 * (10 *
                                            x2)^3 * (1 - x2)^10
  } else if(func==2){
  20/(1+exp(0.5*20)*exp(-20*x2))
  }else if (func==3){
    10*sin(x*1.5*pi)
  }else if (func==4){
    20*x^0.8
  }
}
n <- 200
n_psuedo = 10
k = 50
func_type = 4
sigma = 1
#maximum and minimum ranges to extrapolate over
min_exp = -1
max_exp = 2
x <- runif(n,0,1)
y <- rnorm(n,ff(x,func = func_type),sigma)
dat_psuedo = data.frame(x=  c(seq(min_exp,0,length=n_psuedo),
                               x, seq(1, max_exp,length=n_psuedo)),
                        y= c(rep(0,times=n_psuedo),y,rep(0,times=n_psuedo)),
                        w = c(rep(0,n_psuedo),rep(1,n),rep(0,n_psuedo)))

pd <- data.frame(x=seq(min_exp,max_exp,length=200))

b1 <- gam(y~s(x, m=2, k=k), method="REML")
b2 <- gam(y~s(x, m=1, k=k)+x, method="REML")
b3 <- gam(y~s(x, m=1, k=k*2), weights=w,
          data=dat_psuedo, method="REML")
b4 <- gam(y~s(x, m=2, k=k*2), weights=w,
          data=dat_psuedo, method="REML")
b5 <- gam(y~s(x, m=1, k=k*2), method="REML", knots= list(x=seq(min_exp, max_exp, length=k*2)))
b6 <- gam(y~s(x, m=2, k=k*2), method="REML", knots= list(x=seq(min_exp, max_exp, length=k*2)))


#fv1 <- predict(b1,pd,se=TRUE)
#fv2 <- predict(b2,pd,se=TRUE)
#fv3 <- predict(b3,pd,se=TRUE)
#fv4 <- predict(b4,pd,se=TRUE)
#fv5 <- predict(b5,pd,se=TRUE)
#fv6 <- predict(b6,pd,se=TRUE)


par(mfrow=c(2,3))
n_sims = 50

# set same ylims across plots
if(func_type != 4){
  ylim <- c(-125,125)
}else{
  ylim <- c(-50, 50)
}


plotter <- function(model, pd, n_sims){

  title <- model$title
  model <- model$model

  fv1 <- predict(model, pd, se=TRUE)
  ll <- fv1$fit - fv1$se.fit*2
  ul <- fv1$fit + fv1$se.fit*2

  post_sim <- sim_post_curves(model, pd, n_sims)

  plot(pd$x,fv1$fit,ylim=ylim,type="l",xlab="x",ylab="y",
       main = title,lwd=2, col="blue")
  lines(pd$x,ll,lty=2,lwd=2, col="blue")
  lines(pd$x,ul,lty=2,lwd=2, col="blue")
  matlines(pd$x, post_sim, col=alpha("black",0.25),lty=1)
  points(x,y,pch=19,cex=.25, col="red")

}

models <- list(
  list(model=b1, title="gam(y~s(x,m=2,k=k))"),
  list(model=b2, title="gam(y~s(x,m=2,k=k*2),\nweights=w, data=dat_psuedo)"),
  list(model=b3, title="gam(y~s(x,m=2,k=k*2),\nknots= seq(min_exp,max_exp,length=k*2)"),
  list(model=b4, title="gam(y~s(x,m=1,k=k)+x)"),
  list(model=b5, title="gam(y~s(x,m=1,k=k*2),\nweights=w,data=dat_psuedo)"),
  list(model=b6, title="gam(y~s(x,m=1,k=k*2),\nknots= seq(min_exp,max_exp,length=k*2)"))

ll <- lapply(models, plotter, pd=pd, n_sims=n_sims)

#testing computation time for each extrapolation solution (assuming m=1):
linear_trend_time = microbenchmark(gam(y~s(x,m=1,k=k)+x,method="REML"),times = 25,unit = "ms")
pseudodata_time  = microbenchmark(gam(y~s(x,m=1,k=k*2),weights=w,
                                      data=dat_psuedo,method="REML"),
                                  times = 25,unit = "ms")
ext_knots_time  = microbenchmark(gam(y~s(x,m=1,k=k*2),method="REML",
                                     knots= list(x=seq(min_exp,max_exp,length=k*2))),
                                  times = 25,unit = "ms")


