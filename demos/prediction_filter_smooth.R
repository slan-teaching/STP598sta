library(astsa)
# generate data 
set.seed(1)  
num = 50
w = rnorm(num+1,0,1)
v = rnorm(num,0,1)
                               
mu = cumsum(w)  # states:  mu[0], mu[1], . . ., mu[50] 
y = mu[-1] + v  # obs:  y[1], . . ., y[50]

# filter and smooth (Ksmooth0 does both)
mu0 = 0; sigma0 = 1;  phi = 1; cQ = 1; cR = 1   
ks = Ksmooth0(num, y, 1, mu0, sigma0, phi, cQ, cR)   

# pictures 
par(mfrow=c(3,1))
Time = 1:num

plot(Time, mu[-1], main="Prediction", ylim=c(-5,10))      
  lines(ks$xp)
  lines(ks$xp+2*sqrt(ks$Pp), lty="dashed", col="blue")
  lines(ks$xp-2*sqrt(ks$Pp), lty="dashed", col="blue")

plot(Time, mu[-1], main="Filter", ylim=c(-5,10))
  lines(ks$xf)
  lines(ks$xf+2*sqrt(ks$Pf), lty="dashed", col="blue")
  lines(ks$xf-2*sqrt(ks$Pf), lty="dashed", col="blue")

plot(Time, mu[-1],  main="Smoother", ylim=c(-5,10))
  lines(ks$xs)
  lines(ks$xs+2*sqrt(ks$Ps), lty="dashed", col="blue")
  lines(ks$xs-2*sqrt(ks$Ps), lty="dashed", col="blue") 

mu[1]; ks$x0n; sqrt(ks$P0n)   # initial value info

# In case you can't see the differences in the figures...
# ... either get new glasses or ... 
# ... plot them on the same graph (not shown)
dev.new()
plot(Time, mu[-1], type='n')
abline(v=Time, lty=3, col=8)
abline(h=-1:5, lty=3, col=8)
lines(ks$xp, col=4, lwd=5)
lines(ks$xf, col=3, lwd=5) 
lines(ks$xs, col=2, lwd=5)
points(Time, mu[-1], pch=19, cex=1.5)
names = c("predictor","filter","smoother")
legend("bottomright", names, col=4:2, lwd=5, lty=1, bg="white")