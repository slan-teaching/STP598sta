#########################################################################################
# Adapted from code by: Hedibert Freitas Lopes
#########################################################################################
##-- Notation --##
#           y(t) = x(t) + v(t);    v(t) ~ iid N(0,V)                     
#           x(t) = x(t-1) + w(t);  w(t) ~ iid N(0,W)                        
#  priors:  x(0) ~ N(m0,C0);  V ~ IG(a,b);  W ~ IG(c,d)
#    FFBS:  x(t|t) ~ N(m,C);  x(t|n) ~ N(mm,CC);  x(t|t+1) ~ N(a,R)  
##-- 
ffbs = function(y,V,W,m0,C0){
  n  = length(y);  a  = rep(0,n);  R  = rep(0,n)
  m  = rep(0,n);   C  = rep(0,n);  B  = rep(0,n-1)     
  H  = rep(0,n-1); mm = rep(0,n);  CC = rep(0,n)
  x  = rep(0,n); llike = 0.0
  for (t in 1:n){
    if(t==1){a[1] = m0; R[1] = C0 + W
    }else{ a[t] = m[t-1]; R[t] = C[t-1] + W }
    f      = a[t]
    Q      = R[t] + V
    A      = R[t]/Q
    m[t]   = a[t]+A*(y[t]-f)
    C[t]   = R[t]-Q*A**2
    B[t-1] = C[t-1]/R[t]
    H[t-1] = C[t-1]-R[t]*B[t-1]**2
    llike  = llike + dnorm(y[t],f,sqrt(Q),log=TRUE) }
  mm[n] = m[n]; CC[n] = C[n]
  x[n]  = rnorm(1,m[n],sqrt(C[n]))
  for (t in (n-1):1){
    mm[t] = m[t] + C[t]/R[t+1]*(mm[t+1]-a[t+1])
    CC[t] = C[t] - (C[t]^2)/(R[t+1]^2)*(R[t+1]-CC[t+1])
    x[t]  = rnorm(1,m[t]+B[t]*(x[t+1]-a[t+1]),sqrt(H[t]))  }
  return(list(x=x,m=m,C=C,mm=mm,CC=CC,llike=llike))   }
# Simulate states and data
set.seed(1); W = 0.5; V = 1.0
n  = 100; m0 = 0.0; C0 = 10.0; x0 = 0
w  = rnorm(n,0,sqrt(W))
v  = rnorm(n,0,sqrt(V))
x  = y = rep(0,n)
x[1] = x0   + w[1]
y[1] = x[1] + v[1]
for (t in 2:n){
  x[t] = x[t-1] + w[t]
  y[t] = x[t] + v[t]   }
# actual smoother (for plotting)
ks = Ksmooth0(num=n, y, A=1, m0, C0, Phi=1, cQ=sqrt(W), cR=sqrt(V))
xsmooth = as.vector(ks$xs)
#
run = ffbs(y,V,W,m0,C0)
m   = run$m; C = run$C; mm = run$mm
CC  = run$CC; L1 = m-2*C; U1  = m+2*C
L2  = mm-2*CC; U2 = mm+2*CC
N   = 50
Vs  = seq(0.1,2,length=N)
Ws  = seq(0.1,2,length=N)
likes = matrix(0,N,N)
for (i in 1:N){
  for (j in 1:N){
    V   = Vs[i]
    W   = Ws[j]
    run = ffbs(y,V,W,m0,C0)    
    likes[i,j] = run$llike  }  }  
# Hyperparameters
a = 0.01; b = 0.01; c = 0.01; d = 0.01
# MCMC step
set.seed(90210)
burn  = 10;  M = 1000
niter = burn + M
V1    = V;  W1 = W
draws = NULL
all_draws = NULL
for (iter in 1:niter){
  run   = ffbs(y,V1,W1,m0,C0)
  x     = run$x
  V1    = 1/rgamma(1,a+n/2,b+sum((y-x)^2)/2)
  W1    = 1/rgamma(1,c+(n-1)/2,d+sum(diff(x)^2)/2)
  draws = rbind(draws,c(V1,W1,x))    }
all_draws = draws[,1:2]
q025  = function(x){quantile(x,0.025)}
q975  = function(x){quantile(x,0.975)}
draws = draws[(burn+1):(niter),]
xs    = draws[,3:(n+2)]
lx    = apply(xs,2,q025)
mx    = apply(xs,2,mean)
ux    = apply(xs,2,q975)
##  plot of the data
par(mfrow=c(2,2), mgp=c(1.6,.6,0), mar=c(3,3.2,1,1))
ts.plot(ts(x), ts(y), ylab='', col=c(1,8), lwd=2)
points(y)
legend(0, 11, legend=c("x(t)","y(t)"), lty=1, col=c(1,8), lwd=2, bty="n", pch=c(-1,1))
contour(Vs, Ws, exp(likes), xlab=expression(sigma[v]^2), ylab=expression(sigma[w]^2), 
        drawlabels=FALSE, ylim=c(0,1.2))
points(draws[,1:2], pch=16, col=rgb(.9,0,0,0.3), cex=.7)
hist(draws[,1], ylab="Density",main="", xlab=expression(sigma[v]^2))
abline(v=mean(draws[,1]), col=3, lwd=3)
hist(draws[,2],main="", ylab="Density", xlab=expression(sigma[w]^2))
abline(v=mean(draws[,2]), col=3, lwd=3)
## plot states
par(mgp=c(1.6,.6,0), mar=c(2,1,.5,0)+.5)
plot(ts(mx), ylab='', type='n', ylim=c(min(y),max(y)))
grid(lty=2); points(y)
lines(xsmooth, lwd=4, col=rgb(1,0,1,alpha=.4))
lines(mx, col= 4) 
xx=c(1:100, 100:1)
yy=c(lx, rev(ux))
polygon(xx, yy, border=NA, col= gray(.6,alpha=.2))
lines(y, col=gray(.4))
legend('topleft', c('true smoother', 'data', 'posterior mean', '95% of draws'), lty=1, 
       lwd=c(3,1,1,10), pch=c(-1,1,-1,-1), col=c(6, gray(.4), 4, gray(.6, alpha=.5)), 
       bg='white' )  