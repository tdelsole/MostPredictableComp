rm(list=ls())

lplotfile     = TRUE

ndim          = 5
delta         = 0.1
l0.real       = 1
ntot          = 200
fvec          = rep(1,ndim)
dt            = 0.2
lag.max       = 40

set.seed(8)

ftitle.meta = bquote(Delta==.(delta))


dir.Rlib  = NULL
source(paste(dir.Rlib,'pdf.eps.R',sep=''))
source(paste(dir.Rlib,'gev.R',sep=''))
source(paste(dir.Rlib,'gev.complex.R',sep=''))
source(paste(dir.Rlib,'lyap.R',sep=''))


lambda        = l0.real + delta * (1:ndim-1)
qmat          = fvec %*% t(fvec)

#######################################
######## COMPUTE MAXIMUM APT
#######################################
emat = array(NA,dim=c(ndim,ndim))
for (j in 1:ndim) for (i in 1:ndim) emat[i,j] = 1/(lambda[i] + lambda[j])

apt.exact = as.numeric(rep(NA,ndim))
q.exact   = array(NA,dim=c(ndim,ndim))
for (j in 2:ndim) {
	apt.list        = gev((emat * emat)[1:j,1:j],emat[1:j,1:j])
	q.exact[1:j,j]  = apt.list$q[,1]/apt.list$q[1,1]
	apt.exact[j]    = 2 * apt.list$lambda[1]
}

q.pascal = array(NA,dim=c(ndim,ndim))
for (j in 2:ndim) q.pascal[1:j,j] = (-1)^(1:j-1) * choose(j-1,1:j-1)

#######################################
######## COMPUTE EQUIVALENT AR MODEL
#######################################
lmat               = diag(-lambda)
amat               = diag(exp(-lambda*dt))

gamma              = array(NA,dim=c(ndim,ndim))
for (j in 1:ndim) for (i in 1:ndim) gamma[i,j] = - ( (exp(-(lambda[i] + lambda[j]) * dt) - 1 ) * fvec[i] * fvec[j])/(lambda[i] + lambda[j])
gamma.eigen        = eigen(gamma)
gamma.eigen$values = ifelse(abs(gamma.eigen$values)<1.e-12,0,gamma.eigen$values)
gamma.sqrt         = gamma.eigen$vector %*% diag(sqrt(gamma.eigen$values))

#######################################
######## COMPUTE AUTOCORRELATION FUNCTION
#######################################
lmat.eigen    = eigen(lmat)
covm.inf      = lyap(lmat.eigen$vectors, lmat.eigen$values, qmat)
covm.yx       = array(NA,dim=c(ndim,ndim,lag.max+1))
corr.tau      = array(NA,dim=c(ndim,ndim,lag.max+1))
for (lag in 0:lag.max) covm.yx[,,lag+1] = diag(exp(-lambda*dt*lag)) %*% covm.inf
for (j in 1:ndim) for (i in 1:ndim) corr.tau[i,j,] = covm.yx[i,j,]/sqrt(covm.yx[i,i,1]*covm.yx[j,j,1])

covm.exact    = array(NA,dim=c(lag.max+1,ndim))
covm.pascal   = array(NA,dim=c(lag.max+1,ndim))
for (n in 1:ndim) for (lag in 0:lag.max) covm.exact [lag+1,n] = q.exact [1:n,n] %*% covm.yx[1:n,1:n,lag+1] %*% q.exact [1:n,n] 
for (n in 1:ndim) for (lag in 0:lag.max) covm.pascal[lag+1,n] = q.pascal[1:n,n] %*% covm.yx[1:n,1:n,lag+1] %*% q.pascal[1:n,n] 

#######################################
######## COMPUTE PREDICTABILITY
#######################################
covm.tau = array(NA,dim=c(ndim,ndim,lag.max+1))
covm.tau[,,1] = covm.inf
for (lag in 1:lag.max) covm.tau[,,lag+1] = diag(exp(-lambda*dt)) %*% covm.tau[,,lag] %*% diag(exp(-lambda*dt))
# for (lag in 0:lag.max) covm.tau[,,lag+1] = diag(exp(-lambda*dt*lag)) %*% covm.inf %*% diag(exp(-lambda*dt*lag))


p.tau.eigen  = array(NA,dim=c(ndim,lag.max+1))
p.tau.exact  = array(NA,dim=c(ndim,lag.max+1))
p.tau.pascal = array(NA,dim=c(ndim,lag.max+1))
for (n in 1:ndim) for (lag in 0:lag.max) p.tau.eigen [n,lag+1] = covm.tau[n,n,lag+1]
for (n in 1:ndim) for (lag in 0:lag.max) p.tau.exact [n,lag+1] = q.exact [1:n,n] %*% covm.tau[1:n,1:n,lag+1] %*% q.exact [1:n,n]
for (n in 1:ndim) for (lag in 0:lag.max) p.tau.pascal[n,lag+1] = q.pascal[1:n,n] %*% covm.tau[1:n,1:n,lag+1] %*% q.pascal[1:n,n]
p.tau.eigen  = p.tau.eigen /p.tau.eigen [,1]
p.tau.exact  = p.tau.exact /p.tau.exact [,1]
p.tau.pascal = p.tau.pascal/p.tau.pascal[,1]

### ESTIMATE APT
apt.brute.eigen  = as.numeric(rep(NA,ndim))
apt.brute.exact  = as.numeric(rep(NA,ndim))
apt.brute.pascal = as.numeric(rep(NA,ndim))
for (n in 1:ndim) apt.brute.eigen [n] = (p.tau.eigen [n,1] + 2 * sum(p.tau.eigen [n,1:lag.max+1])) * dt
for (n in 1:ndim) apt.brute.exact [n] = (p.tau.exact [n,1] + 2 * sum(p.tau.exact [n,1:lag.max+1])) * dt
for (n in 1:ndim) apt.brute.pascal[n] = (p.tau.pascal[n,1] + 2 * sum(p.tau.pascal[n,1:lag.max+1])) * dt

#######################################
######## GENERATE TIME SERIES
#######################################
yvec     = array(NA,dim=c(ntot,ndim))
yvec[1,] = 0
for (nt in 2:ntot) yvec[nt,] = amat %*% yvec[nt-1,] + gamma.sqrt %*% rnorm(ndim)

y.max.exact    = array(NA,dim=c(ntot,ndim))
for (n in 1:ndim) y.max.exact[,n] = yvec[,1:n,drop=FALSE] %*% q.exact[1:n,n]

y.max.pascal   = array(NA,dim=c(ntot,ndim))
for (n in 1:ndim) y.max.pascal[,n] = yvec[,1:n,drop=FALSE] %*% q.pascal[1:n,n]

#######################################
######## PLOTS
#######################################
fout = './figures/NDim.TimeSeries'
if (lplotfile) pdf.eps(fout,'pdf',height = 10,width=8.5)
par(mfcol=c(3,1),mar=c(4.5,4.5,3,0.2))
cex.say = 1.7
par(cex.lab=cex.say,cex.axis=cex.say,cex.main=cex.say)

###### PLOT TIME SERIES
y.plot = array(NA,dim=c(ntot,ndim+ndim))
y.plot[,1:ndim] = yvec

offset = 0.7
offset.one.time = 1
ynorm  = as.numeric(rep(NA,ndim))
for (n in 1:ndim) {
	ynorm[n] = max(abs(y.max.exact[,n]))
	y.plot[,n+ndim] = y.max.exact[,n]/ynorm[n] - (n-1) * offset - offset.one.time
}

yrange = range(y.plot,na.rm=TRUE)
xrange = c(1,ntot+20)
plot(1,1,type='n',xlim=xrange,ylim=yrange,xlab='time',ylab='Y')
for (n in 1:ndim) lines(y.plot[,n     ],lwd=1,col='grey50')
for (n in 2:ndim) lines(y.plot[,n+ndim],lwd=2,col=n-1)
for (n in 2:ndim) text(ntot,y.plot[ntot,n+ndim],paste('x',round(1/ynorm[n])),pos=4,col=n-1)
text(ntot,y.plot[ntot,1],'modes',pos=4,col='black')
title(main='Time Series Realizations From a Stochastic Model with Rank-1 Forcing',line=2.0)
title(main=ftitle.meta,line=0.5)
mtext('(a)',side=3,adj=0,cex=1.3,line=0.5)

### PLOT ACF
yrange = c(0,1)
x.plot = dt*(0:lag.max)
xrange = range(x.plot)
ylab.say = bquote(P[tau])
plot(1,1,type='n',xlim=xrange, ylim=yrange,xlab='time lag',ylab=ylab.say)

for (n in 1:ndim) lines(x.plot,p.tau.eigen[n,],lwd=2,col='grey50')
for (n in 1:ndim) lines(x.plot,p.tau.exact[n,],lwd=2,col=n-1)
title(main='Predictability',line=0.5)
legend.say = paste('APT optimal for ',2:ndim,'D',sep='')
legend.say = c('individual eigenmodes',legend.say)
legend('topright',legend=legend.say,col=c('grey50',1:ndim),lwd=2,cex=1.3)
mtext('(b)',side=3,adj=0,cex=1.3,line=0.5)

ylab.say = 'ACF'
plot(1,1,type='n',xlim=xrange, ylim=yrange,xlab='time lag',ylab=ylab.say)
for (i in 1:ndim) lines(x.plot,corr.tau[i,i,],lwd=2,col='grey50')
for (n in 2:ndim) lines(x.plot,covm.exact[,n]/covm.exact[1,n],lwd=2,col=n-1)
# for (n in 2:ndim) lines(x.plot,covm.pascal[,n]/covm.pascal[1,n],lwd=2,col=n-1)
title(main='Autocorrelation Function',line=0.5)
legend.say = paste('APT optimal for ',2:ndim,'D',sep='')
legend.say = c('individual eigenmodes',legend.say)
legend('topright',legend=legend.say,col=c('grey50',1:ndim),lwd=2,cex=1.3)
mtext('(c)',side=3,adj=0,cex=1.3,line=0.5)

## PLOT SAMPLE ACFS
# acf.eig1  = acf(yvec[,1],lag.max=lag.max,plot=FALSE)
# acf.max1  = acf(y.max.exact[,ndim],lag.max=lag.max,plot=FALSE)
# lines(x.plot,acf.eig1$acf,col='red' ,lwd=2,lty='dashed')
# lines(x.plot,acf.max1$acf,col='blue',lwd=2,lty='dashed')

if (lplotfile) dev.off()
