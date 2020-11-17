library(depmixS4)
##  note - as of updata 1.4-2, you can get standard errors ##
model <- depmix(EQcount ~1, nstates=2, data=data.frame(EQcount), family=poisson('identity'), respstart=c(15,25))
set.seed(90210)
summary(fm <- fit(model))   # estimation results  
standardError(fm)           # with standard errors

##-- Get + Display Parameters --##
para.mle = as.vector(getpars(fm))[3:8]  
( mtrans = matrix(para.mle[1:4], byrow=TRUE, nrow=2) )
( lams   = para.mle[5:6] )  
( pi1    = mtrans[2,1]/(2 - mtrans[1,1] - mtrans[2,2]) )  
( pi2    = 1 - pi1 )

#-- Graphics --##
par(mfrow=c(3,1))
# data and states
tsplot(EQcount, main="", ylab='EQcount', type='h', col=gray(.7))
text(EQcount, col=6*posterior(fm)[,1]-2, labels=posterior(fm)[,1], cex=.9)
# prob of state 2
tsplot(ts(posterior(fm)[,3], start=1900), ylab = expression(hat(pi)[~2]*'(t|n)'));  abline(h=.5, lty=2)
# histogram
hist(EQcount, breaks=30, prob=TRUE, main="")
xvals = seq(1,45)
u1 = pi1*dpois(xvals, lams[1])  
u2 = pi2*dpois(xvals, lams[2])
lines(xvals, u1, col=4)   
lines(xvals, u2, col=2)