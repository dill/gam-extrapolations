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
k = 50
func_type = 4
sigma = 1
#maximum and minimum ranges to extrapolate over
min_exp = -1 
max_exp = 2
x <- runif(n,0,1)
y <- rnorm(n,ff(x,func = func_type),sigma)
dat_psuedo = data.frame(x=  c(seq(min_exp,0,length=50),
                               x, seq(1, max_exp,length=50)),
                        y= c(rep(0,times=50),y,rep(0,times=50)),
                        w = c(rep(0,50),rep(1,n),rep(0,50)))

pd <- data.frame(x=seq(min_exp,max_exp,length=200))

b1 <- gam(y~s(x,m=2,k=k),method="REML")
b2 <- gam(y~s(x,m=1,k=k)+x,method="REML")
b3 <- gam(y~s(x,m=1,k=k*2),weights=w,
          data=dat_psuedo,method="REML")
b4 <- gam(y~s(x,m=2,k=k*2),weights=w, 
          data=dat_psuedo,method="REML")

fv1 <- predict(b1,pd,se=TRUE)
fv2 <- predict(b2,pd,se=TRUE)
fv3 <- predict(b3,pd,se=TRUE)
fv4 <- predict(b4,pd,se=TRUE)

par(mfrow=c(2,2))
n_sims = 50
 ## black order 2 penalty...
 ll <- fv1$fit-fv1$se.fit*2;ul <- fv1$fit+fv1$se.fit*2
 post_sim = sim_post_curves(b1,pd,n_sims)
 plot(pd$x,fv1$fit,ylim=range(c(ul,ll)),type="l",xlab="x",ylab="y",
      main  ="gam(y~s(x,m=2,k=k))",lwd=2)
 lines(pd$x,ll,lty=2,lwd=2);lines(pd$x,ul,lty=2,lwd=2)
 matlines(pd$x, post_sim, col=alpha("black",0.25),lty=1)
 points(x,y,pch=19,cex=.4)
 
 ## blue order 2 penalty + pseudodata...
 ll <- fv4$fit-fv4$se.fit*2;ul <- fv4$fit+fv4$se.fit*2
 post_sim = sim_post_curves(b4,pd,n_sims)
 plot(pd$x,fv4$fit,type="l",ylim=range(c(ul,ll)),xlab="x",ylab="y",
      main  ="gam(y~s(x,m=2,k=k*2),\nweights=w, data=dat_psuedo)",
      col=4,lwd=2)
 lines(pd$x,ll,lty=2,col=4,lwd=2);lines(pd$x,ul,lty=2,col=4,lwd=2)
 matlines(pd$x, post_sim, col=alpha(4,0.25),lty=1)
 points(x,y,pch=19,cex=.4)

## red order 1 penalty + linear effect...
ll <- fv2$fit-fv2$se.fit*2;ul <- fv2$fit+fv2$se.fit*2
post_sim = sim_post_curves(b2,pd,n_sims)
plot(pd$x,fv2$fit,ylim=range(c(ul,ll)),type="l",xlab="x",ylab="y",
     main  ="gam(y~s(x,m=1,k=k)+x)",col=2,lwd=2)
lines(pd$x,ll,lty=2,col=2,lwd=2);lines(pd$x,ul,lty=2,
                                       col=2,lwd=2)
matlines(pd$x, post_sim, col=alpha(2,0.25),lty=1)
points(x,y,pch=19,cex=.4)

## green order 1 penalty + pseudodata...
ll <- fv3$fit-fv3$se.fit*2;ul <- fv3$fit+fv3$se.fit*2
post_sim = sim_post_curves(b3,pd,n_sims)
plot(pd$x,fv3$fit,type="l",ylim=range(c(ul,ll)),xlab="x",ylab="y",
     main  ="gam(y~s(x,m=1,k=k*2),\nweights=w,data=dat_psuedo)",
     col=3,lwd=2)
lines(pd$x,ll,lty=2,col=3,lwd=2);lines(pd$x,ul,lty=2,col=3,lwd=2)
matlines(pd$x, post_sim, col=alpha(3,0.25),lty=1)
points(x,y,pch=19,cex=.4)

