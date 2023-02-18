rm(list=ls())

lplotfile      = TRUE
neof.max       = 50
area.name      = 'PACIFIC30S30N'
area.name      = 'NASST'
area.name      = 'NPSST'
ntrun          = 5
npoly          = 3
time.window    = 'seasonal'
time.window    = 'annual'
time.window    = 'monthly'

iyst           = 1950
iynd           = 2021

dir.rlib       = NULL
dir.ersst      = '/Users/delsole/data/ersstv5/'

nyrs           = iynd - iyst + 1

ftitle.meta = paste(ntrun,' EOFs; ERSSTv5; ',area.name,'; ',iyst,'-',iynd,'; ',time.window,sep='')


library(base)
library(fields)
library(ncdf4)
library(maps)
library(rworldmap)
library(RNetCDF)
library(lubridate)

source(paste(dir.rlib,'plot_latlon_v4.R',sep=''))
source(paste(dir.rlib,'plot_latlon_contour.R',sep=''))
source(paste(dir.rlib,'eof.latlon.R',sep=''))
source(paste(dir.rlib,'plot_eof.R',sep=''))
source(paste(dir.rlib,'index.climate.v2.R',sep=''))
source(paste(dir.rlib,'pdf.eps.R',sep=''))
source(paste(dir.rlib,'lyap.R',sep=''))
source(paste(dir.rlib,'gev.R',sep=''))
source(paste(dir.rlib,'gev.complex.R',sep=''))


##################################
######### READ DATA
##################################
filename   = paste(dir.ersst,'sst.mnmean.nc',sep='')
ncin       = nc_open(filename)
lon        = ncvar_get(ncin,'lon')
lat        = ncvar_get(ncin,'lat')
sst.all    = ncvar_get(ncin,'sst')
time.coord = ncvar_get(ncin,'time')
time.units = ncatt_get(ncin,'time','units')
nc_close(ncin)

time.sst   = utcal.nc(time.units$value,time.coord,'c')
year.sst   = year(time.sst)
month.sst  = month(time.sst)

nlon       = length(lon)
nlat       = length(lat)

ntime.get  = year.sst >= iyst & year.sst <= iynd
sst.all    = sst.all[,,ntime.get]

ntot       = sum(ntime.get)
nyrs       = ntot / 12
year.say   = seq(from=iyst,by=1/12,length.out=ntot)

plot_latlon_v4(lon,lat,sst.all[,,100])
ocean = !is.na(sst.all[,,1])

area.say = NA
if (area.name == 'NASST') area.say = 'North Atlantic SST'
if (area.name == 'PACIFIC30S30N') area.say = 'Tropical Pacific'
if (area.name == 'NPSST') area.say = 'North Pacific'
if (is.na(area.say)) stop('give your domain a name')

##################################
######### DEFINE DOMAIN
##################################
lgood = index.climate.v2(lon,lat,area.name=area.name,ocean=ocean)$view
plot_latlon_v4(lon,lat,lgood)
dim(sst.all) = c(nlon*nlat,ntot)
sst.all[!lgood,] = NA

##################################
######### REMOVE ANNUAL CYCLE
##################################
dim(sst.all) = c(nlon*nlat*12,nyrs)
sst.mean     = rowMeans(sst.all)
sst.all      = sst.all - sst.mean
dim(sst.all) = c(nlon,nlat,ntot)
plot_latlon_v4(lon,lat,sst.all[,,100])

##################################
######### REMOVE POLYNOMIAL
##################################
poly.time    = poly(year.say,npoly)
dim(sst.all) = c(nlon*nlat,ntot)
sst.all      = sst.all - (sst.all %*% poly.time) %*% t(poly.time)

##################################
######### CONVERT TO SEASONAL OR ANNUAL MEANS
##################################
if (time.window == 'monthly') {
	nchunk = 1
	tunit.say = 'months'
} else if (time.window == 'seasonal') {
	nchunk = 3
	tunit.say = 'seasons'
} else if (time.window == 'annual') {
	nchunk = 12
	tunit.say = 'years'
} else stop('do not understand time.window')

if (nchunk != 1) {
	dum          = sst.all
	dim(dum)     = c(nlon*nlat,nchunk,ntot/nchunk)
	sst.all      = array(NA,dim=c(nlon*nlat,ntot/nchunk))
	for (n in 1:(ntot/nchunk)) sst.all[,n] = rowMeans(dum[,,n])
	dim(sst.all) = c(nlon,nlat,ntot/nchunk)
	ntot         = ntot/nchunk
	year.say     = seq(from=iyst,by=nchunk/12,length.out=ntot)
	rm(dum)
}


##################################
######### COMPUTE EOFS
##################################
eof.list = eof.latlon(lon,lat,sst.all)

nplot.pic = 1
fout = paste('./figures/rank1.',area.name,'.eof',nplot.pic,sep='')
if (lplotfile) pdf.eps(fout,'pdf',height=8)
par(mfcol=c(2,1),mar=c(6,5,4,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
plot_latlon_contour(lon,lat,eof.list$eof[, nplot.pic],shrinkdomain=TRUE)
ftitle.top = paste('EOF',nplot.pic,sep='')
title(ftitle.top,line=2.0)
title(ftitle.meta,line=0.5)
par(mar=c(5,5,3,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
plot(year.say,eof.list$pc[, nplot.pic],type='l',xlab='year',ylab='')
title('PC time series',line=0.5)
if (lplotfile) dev.off()


##################################
######### COMPUTE LIM
##################################
npic   = 2:ntot
y.all  = eof.list$pc[npic  ,1:ntrun]
x.all  = eof.list$pc[npic-1,1:ntrun]
lim.lm = lm(y.all~x.all)
lmat   = lim.lm$coef[-1,]
qmat   = ( t(residuals(lim.lm)) %*% residuals(lim.lm) ) / lim.lm$df.residual

##################################
######### COMPUTE WHITENED Q
##################################
lmat.eigen   = eigen(lmat)
smat.inv     = solve(lmat.eigen$vectors)

q.tilde      = smat.inv %*% qmat %*% t(Conj(smat.inv))

q.eigen = eigen(q.tilde,symmetric=TRUE)
print(q.eigen$values)

fout = paste('./figures/rank1.',area.name,'.E',ntrun,'.Qeigen',sep='')
if (lplotfile) pdf.eps(fout,'pdf',height=6)
par(mfcol=c(1,1),mar=c(5,5,4,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
plot(q.eigen$values,type='l',xlab='eigenmode',ylab='eigenvalue')
points(q.eigen$values,pch=19)
title("Eigenvalues of Transformed Noise Covariance Matrix",line=2.0)
title(ftitle.meta,line=0.5)
if (lplotfile) dev.off()

### SOMETIMES THE AR-EIGENVALUE IS NEGATIVE, SO LOG(EVALUE) = LOG(-EVALUE) + i * PI.
# amat       = lmat.eigen$vectors %*% diag(eval.lim) %*% smat.inv
if (any(Re(lmat.eigen$values) < 0)) stop('violently oscillatory mode exists')


##################################
######### COMPUTE APT
##################################
eval.lim   = log (lmat.eigen$values)
covm       = lyap(lmat.eigen$vectors,eval.lim,qmat)
phim       = lyap(lmat.eigen$vectors,eval.lim,qmat,e.square=TRUE)

if (any(abs(Im(covm)) > 1.e-13)) stop('imaginary values in covm') else covm = Re(covm)
if (any(abs(Im(phim)) > 1.e-13)) stop('imaginary values in phim') else phim = Re(phim)

apt.eigen        = gev(phim,covm)
q.max            = apt.eigen$q[,1]/apt.eigen$q[1,1]
apt.eigen$lambda = 2 * apt.eigen$lambda

### check
apt.max.check    = 2 * (Conj(q.max) %*% phim %*% q.max) / (Conj(q.max) %*% covm %*% q.max)

apt.mode         = -1/Re(eval.lim[1])
apt.mode.check2  = 2 * (Conj(smat.inv[1,]) %*% phim %*% smat.inv[1,]) / (Conj(smat.inv[1,]) %*% covm %*% smat.inv[1,])

##################################
######### COMPUTE RANK1 UPPER BOUND
##################################
lambda = eval.lim
emat = array(NA,dim=c(ntrun,ntrun))
for ( j in 1:ntrun) for ( i in 1:ntrun) emat[i,j] = -1/(lambda[i] + Conj(lambda[j]))
apt.upper.bound = gev.complex(2 * emat * emat,emat)$values

### check
q.rank1    = array(1,dim=c(ntrun,ntrun))
covm.rank1 = lyap(lmat.eigen$vectors,eval.lim,q.rank1)
phim.rank1 = lyap(lmat.eigen$vectors,eval.lim,q.rank1,e.square=TRUE)

if (any(abs(Im(covm.rank1)) > 1.e-13)) stop('imaginary values in covm.rank1') else covm.rank1 = Re(covm.rank1)
if (any(abs(Im(phim.rank1)) > 1.e-13)) stop('imaginary values in phim.rank1') else phim.rank1 = Re(phim.rank1)

apt.rank1  = gev(2*phim.rank1,covm.rank1)$lambda

## UNPROVEN UPPER BOUND
apt.max.unproven = -sum(4/(lambda))


##################################
######### PASCAL MODE
### since Q is not rank1, technically there is not a Pascal model
### have to choose leading eigenvector of Q.  
### PASCAL MODE SHOULD NOT BE PLOTTED!
##################################
if (Im(eval.lim)[1] != 0) stop('first eigenmode is not pure real')
cmat         = array(NA,dim=c(ntrun,ntrun))
cmat[1,]     = c(1,rep(0,ntrun-1))
for (np in 0:(ntrun-2)) cmat[np+2,] = (eval.lim/Re(eval.lim)[1])^np
q.moment     = solve(cmat,c(1,rep(0,ntrun-1)))

fvec         = q.eigen$vectors[,1]
q.pascal     = t(Conj(smat.inv)) %*% diag(1/Conj(fvec)) %*% q.moment
q.pascal     = q.pascal/q.pascal[1]
if (any(abs(Im(q.pascal)) > 1.e-5)) stop('q.pascal has imaginary values')
q.pascal     = as.numeric(Re(q.pascal))

apt.pascal   = 2*as.numeric( ( q.pascal %*% phim %*% q.pascal ) / ( q.pascal %*% covm %*% q.pascal ) )

q.pascal.check = complex(ntrun)
for (j in 1:ntrun) q.pascal.check[j] = 1/prod(eval.lim[j] - eval.lim[-j])
q.pascal.check = q.pascal.check/q.pascal.check[1]
print(all.equal(q.moment,q.pascal.check))


##################################
######### LOADING VECTORS
##################################
p.max          = covm %*% q.max
p.pascal       = covm %*% q.pascal
p.max.space    = eof.list$eof [,1:ntrun] %*% p.max
p.pascal.space = eof.list$eof [,1:ntrun] %*% p.pascal
q.max.space    = eof.list$eofi[,1:ntrun] %*% q.max
q.pascal.space = eof.list$eofi[,1:ntrun] %*% q.pascal

if (sum(p.pascal      ,na.rm=TRUE) < 0) p.pascal       = -p.pascal
if (sum(p.max.space   ,na.rm=TRUE) < 0) p.max.space    = -p.max.space
if (sum(p.pascal.space,na.rm=TRUE) < 0) p.pascal.space = -p.pascal.space
if (sum(q.max.space   ,na.rm=TRUE) < 0) q.max.space    = -q.max.space
if (sum(q.pascal.space,na.rm=TRUE) < 0) q.pascal.space = -q.pascal.space

##################################
######### P TAU
##################################
tau.max     = 48
ntau        = 201
tau         = seq(from=0,to=tau.max,length.out=ntau)
ptau.exact  = as.numeric(rep(NA,ntau))
ptau.mode   = as.numeric(rep(NA,ntau))
ptau.pascal = as.numeric(rep(NA,ntau))
ptau.max    = as.numeric(rep(NA,ntau))
for (qtype in 1:4) for (nt in 1:ntau) {
	numerator   = 0
	denominator = 0
	if (qtype == 1) {
		q = q.max 
	} else if (qtype == 2) {
		q = smat.inv[1,]
	} else if (qtype == 3) {
		q = q.pascal
	} else if (qtype == 4) {
		q = NA
	}
	
	etau = array(NA,dim=c(ntrun,ntrun))
	for (j in 1:ntrun) for (i in 1:ntrun) etau[i,j] = -exp(tau[nt]*(lambda[i]+Conj(lambda[j])))/(lambda[i]+Conj(lambda[j]))
	emat = array(NA,dim=c(ntrun,ntrun))
	for (j in 1:ntrun) for (i in 1:ntrun) emat[i,j] = -1/(lambda[i]+Conj(lambda[j]))
	
	if (qtype != 4) {
		shq         = Conj(t(lmat.eigen$vectors)) %*% q
		denominator = Conj(t(shq)) %*% ( q.tilde * emat ) %*% shq		
		numerator   = Conj(t(shq)) %*% ( q.tilde * etau ) %*% shq	    
		dum = Re(numerator)/Re(denominator)		
	} else {
		gev.list = gev.complex(etau * q.tilde,emat * q.tilde)
		ptau.max[nt] = gev.list$values[1]
	}
	
	    
	if (qtype == 1) {
		ptau.exact[nt] = dum 
	} else if (qtype == 2) {
		ptau.mode[nt] = dum
	} else if (qtype == 3) {
		ptau.pascal[nt] = dum
	}
}



##### INTEGRAL CHECK
apt.exact.check  = (1 + 2 * sum(ptau.exact[-1]) ) * (tau[2]-tau[1])
apt.mode.check   = (1 + 2 * sum(ptau.mode[-1])  ) * (tau[2]-tau[1])
apt.pascal.check = (1 + 2 * sum(ptau.pascal[-1])) * (tau[2]-tau[1])

apt.archive = c(apt.eigen$lambda[1],apt.exact.check)
apt.archive = rbind(apt.archive,c(apt.mode,apt.mode.check))
colnames(apt.archive) = c('exact','finite integral')
rownames(apt.archive) = c('optimal','l.d.mode')
print(apt.archive)

ptau.mode.check  = exp(2*tau*Re(lambda[1]))

fout = paste('./figures/rank1.',area.name,'.E',ntrun,'.ptau',sep='')
if (lplotfile) pdf.eps(fout,'pdf',height=8)
par(mfcol=c(1,1),mar=c(6,5,4,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
yrange = c(0,1)
xrange = c(0,tau.max)
plot(1,1,type='n',xlim=xrange,ylim=yrange,xlab='lead time (months)',ylab=bquote(P[tau]))
lines(tau,ptau.exact,col='black',lwd=2)
lines(tau,ptau.mode ,col='red',lwd=2)
# lines(tau,ptau.pascal,col='green',lwd=2)
ftitle.top = bquote(P[tau]~'of the Most Predictable Component and Least Damped mode')
title(ftitle.top,line=2.0)
title(ftitle.meta,line=0.5)
legend('topright',legend=c('Most Predictable Component','Least Damped Mode'),lwd=2,col=c('black','red'),cex=1.2)
if (lplotfile) dev.off()



fout = paste('./figures/rank1.',area.name,'.E',ntrun,'.q',sep='')
if (lplotfile) pdf.eps(fout,'pdf',height=8)
par(mfcol=c(2,1),mar=c(6,5,4,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
dum = plot_latlon_contour(lon,lat,q.max.space,shrinkdomain=TRUE)
ftitle.top = paste('Projection Vector of Most Predictable Component',sep='')
title(ftitle.top,line=2.0)
title(ftitle.meta,line=0.5)
dum = plot_latlon_contour(lon,lat,q.pascal.space,shrinkdomain=TRUE)
ftitle.top = paste('Projection Vector of the Pascal Mode',sep='')
title(ftitle.top,line=0.5)
if (lplotfile) dev.off()

fout = paste('./figures/rank1.',area.name,'.E',ntrun,'.p',sep='')
if (lplotfile) pdf.eps(fout,'pdf',height=8)
par(mfcol=c(2,1),mar=c(6,5,4,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
dum = plot_latlon_contour(lon,lat,p.max.space,shrinkdomain=TRUE)
ftitle.top = paste('Loading Vector of Most Predictable Component',sep='')
title(ftitle.top,line=2.0)
title(ftitle.meta,line=0.5)
dum = plot_latlon_contour(lon,lat,p.pascal.space,shrinkdomain=TRUE)
ftitle.top = paste('Loading Vector of the Pascal Mode',sep='')
title(ftitle.top,line=0.5)
if (lplotfile) dev.off()

fout = paste('./figures/rank1.',area.name,'.E',ntrun,'.APT',sep='')
if (lplotfile) pdf.eps(fout,'pdf',height=6)
par(mfcol=c(1,1),mar=c(5,5,4,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
ylab.say = paste('APT (',tunit.say,')',sep='')
plot(apt.eigen$lambda,type='l',xlab='eigenmode',ylab=ylab.say)
points(apt.eigen$lambda,pch=19)
title("Maximized APTs",line=2.0)
title(ftitle.meta,line=0.5)
if (lplotfile) dev.off()

fout = paste('./figures/rank1.',area.name,'.E',ntrun,'.L.eval',sep='')
if (lplotfile) pdf.eps(fout,'pdf',height=8)
par(mfcol=c(2,1),mar=c(5,5,4,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
ylab.say = bquote('real eigenvalue ('*.(tunit.say)^{-1}*')')
plot(Re(eval.lim),type='l',xlab='eigenmode',ylab=ylab.say)
points(Re(eval.lim),pch=19)
title("Dynamical Eigenvalues (Real)",line=2.0)
title(ftitle.meta,line=0.5)
xlab.say = bquote('real ('*.(tunit.say)^{-1}*')')
ylab.say = bquote('imaginary ('*.(tunit.say)^{-1}*')')
plot(Re(eval.lim),Im(eval.lim),xlab=xlab.say,ylab=ylab.say,pch=19)
title("Dynamical Eigenvalues (Complex)",line=0.5)
abline(h=0,col='grey')
abline(v=0,col='grey')
if (lplotfile) dev.off()


fout = paste('./figures/rank1.',area.name,'.E',ntrun,'.qvec',sep='')
if (lplotfile) pdf.eps(fout,'pdf',height=8)
par(mfcol=c(1,1),mar=c(5,5,4,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
ylab.say = 'q'
yrange = range(q.max,q.pascal)
plot(q.max,type='l',xlab='vector element',ylab=ylab.say,ylim=yrange,lwd=2)
points(q.max,pch=19)
abline(h=0,col='grey')
lines (q.pascal,col='red',lwd=2)
points(q.pascal,col='red',pch=19)
title("Projection Vector for Max APT",line=2.0)
title(ftitle.meta,line=0.5)
legend('topright',legend=c('exact','pascal'),col=c('black','red'),lwd=2)
if (lplotfile) dev.off()

fout = paste('./figures/rank1.',area.name,'.E',ntrun,'.L.evec',sep='')
if (lplotfile) pdf.eps(fout,'pdf',height=8)
par(mfcol=c(2,1),mar=c(6,5,4,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
L.evec = eof.list$eof [,1:ntrun] %*% Re(lmat.eigen$vectors[1:ntrun,1])
if (sum(L.evec,na.rm=TRUE) < 0) L.evec = - L.evec
dum = plot_latlon_contour(lon,lat,L.evec,shrinkdomain=TRUE)
ftitle.top = paste('Least Damped Eigenmode (real)',sep='')
title(ftitle.top,line=2.0)
title(ftitle.meta,line=0.5)
dum = plot_latlon_contour(lon,lat,p.pascal.space,shrinkdomain=TRUE)
ftitle.top = paste('Loading Vector of the Pascal Mode',sep='')
title(ftitle.top,line=0.5)
if (lplotfile) dev.off()


#################################
######## COMPUTE PREDICTABILTY VERSUS EOF TRUNCATION
#################################
ntrun.max = 20
npic      = 2:ntot
apt.rank1.coalesce = as.numeric(rep(NA,ntrun.max))
apt.rank1.actual   = as.numeric(rep(NA,ntrun.max))
apt.actual.actual  = as.numeric(rep(NA,ntrun.max))
apt.mode           = as.numeric(rep(NA,ntrun.max))
for (nt in 2:ntrun.max) {
	y.all  = eof.list$pc[npic  ,1:nt]
	x.all  = eof.list$pc[npic-1,1:nt]
	lim.lm = lm(y.all~x.all)
	lmat   = lim.lm$coef[-1,]
	qmat   = ( t(residuals(lim.lm)) %*% residuals(lim.lm) ) / lim.lm$df.residual
	
	lmat.eigen   = eigen(lmat)
	smat.inv     = solve(lmat.eigen$vectors)
	q.tilde      = smat.inv %*% qmat %*% t(Conj(smat.inv))
	
	eval.lim   = log (lmat.eigen$values)
	covm       = lyap(lmat.eigen$vectors,eval.lim,qmat)
	phim       = lyap(lmat.eigen$vectors,eval.lim,qmat,e.square=TRUE)
	
	if (any(abs(Im(covm)) > 1.e-13)) stop('imaginary values in covm') else covm = Re(covm)
	if (any(abs(Im(phim)) > 1.e-13)) stop('imaginary values in phim') else phim = Re(phim)
	
	lambda     = eval.lim
	
	if (abs(Im(eval.lim[1])) < 1.e-13) {
		### RANK1 + LAMBDA-COALSENCE
		k                  = 0:nt
		char.poly          = choose(nt,k)^2 * factorial(k) * (-1)^k
		char.roots         = polyroot(rev(char.poly))
		if (any(abs(Im(char.roots)) > 1.e-5)) stop('complex roots detected') else char.roots = Re(char.roots)
		apt.rank1.coalesce[nt] = -char.roots[nt]/Re(eval.lim[1])
		
		### RANK1 + ACTUAL EIGENVALUES
		emat  = array(NA,dim=c(nt,nt))
		emat2 = array(NA,dim=c(nt,nt))
		for (j in 1:nt) for (i in 1:nt) emat [i,j] = -1/(lambda[i] + Conj(lambda[j]))
		for (j in 1:nt) for (i in 1:nt) emat2[i,j] =  1/(lambda[i] + Conj(lambda[j]))^2
		gev.list = gev.complex(2 * emat2, emat)
		apt.rank1.actual[nt] = gev.list$values[1]
		
		### ACTUAL Q + ACTUAL EIGENVALUES
		gev.list = gev.complex(2 * (q.tilde * emat2), q.tilde * emat)
		apt.actual.actual[nt] = gev.list$values[1]
		
		## check
		gev.list = gev.complex(2*phim,covm)
	
		### MODE
		apt.mode[nt] = -1/Re(eval.lim[1])
	
	} else print('first dynamical eigenvalue is complex')
	
}

fout = paste('./figures/rank1.',area.name,'.E',ntrun,'.APT.vs.D',sep='')
if (lplotfile) pdf.eps(fout,'pdf',height=5)
par(mfcol=c(1,1),mar=c(5,5,3,1))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
yrange = range(apt.rank1.coalesce,apt.rank1.actual,apt.actual.actual,apt.mode,na.rm=TRUE)
if (area.name == 'NASST') yrange = c(0,40) else yrange = c(0,100)
xrange = c(1,ntrun.max)
plot(1,1,type='n',xlim=xrange,ylim=yrange,xlab='dimension D',ylab='Maximum APT (months)')
lines(apt.rank1.coalesce,col='red',lwd=2)
lines(apt.rank1.actual,col='darkgreen',lwd=2)
lines(apt.actual.actual,col='black',lwd=2)
lines(apt.mode,col='blue',lwd=2)
legend('topright',legend=c(
  'Rank1 Q; Coalesce E.values',
  'Rank1 Q; Actual E.values',
  'Actual Q; Actual E.values',
  'least damped mode'),
  lwd=2,col=c('red','darkgreen','black','blue'))
title(main=paste('Maximum APTs for LIMs for',area.say),line=0.7)
if (lplotfile) dev.off()


#################################
######## PLOT APT AND P-TAU ON SAME FIGURE
#################################
fout = paste('./figures/rank1.',area.name,'.E',ntrun,'.APT.Ptau',sep='')
if (lplotfile) pdf.eps(fout,'pdf',height=4)
par(mfcol=c(1,2),mar=c(4.5,4.5,3,0.5))
cex.all = 1.1
par(cex.lab=cex.all,cex.axis=cex.all,cex.main=cex.all)
yrange = range(apt.rank1.coalesce,apt.rank1.actual,apt.actual.actual,apt.mode,na.rm=TRUE)
if (area.name == 'NASST') yrange = c(0,40) else yrange = c(0,100)
xrange = c(1,ntrun.max)
plot(1,1,type='n',xlim=xrange,ylim=yrange,xlab='dimension D',ylab='Maximum APT (months)')
lines(apt.rank1.coalesce,col='gold3',lwd=2)
lines(apt.rank1.actual,col='darkgreen',lwd=2)
lines(apt.actual.actual,col='black',lwd=2)
lines(apt.mode,col='red',lwd=2)
legend('topright',legend=c(
  'Rank1 Q; Coalescent lambda',
  'Rank1 Q; Actual lambda',
  'Actual Q; Actual lambda',
  'least damped mode'),
  lwd=2,col=c('gold3','darkgreen','black','red'),cex=0.8,bg='white')
title(main=paste('Maximum APT'),line=2.0)
title(main=paste('LIMs of',area.say),line=0.7)
mtext('(a)',side=3,adj=0,cex=1.0,line=0.5)

yrange = c(0,1)
xrange = c(0,tau.max)
par(cex.lab=cex.all,cex.axis=cex.all,cex.main=cex.all)
plot(1,1,type='n',xlim=xrange,ylim=yrange,xlab='lead time (months)',ylab=bquote(P[tau]))
lines(tau,ptau.exact,col='black',lwd=2)
lines(tau,ptau.mode ,col='red',lwd=2)
lines(tau,ptau.max  ,col='blue',lwd=2,lty='dashed')
# lines(tau,ptau.pascal,col='green',lwd=2)
ftitle.top = bquote(P[tau]~'of the Most Predictable Component and Least Damped mode')
ftitle.top = bquote(P[tau]~'based on'~.(ntrun)~'EOFs')
title(ftitle.top ,line=0.7)
# title(ftitle.meta,line=0.5)
legend('topright',legend=c('APT optimal','finite-time optimal','Least Damped Mode'),lwd=2,col=c('black','blue','red'),cex=0.8)
mtext('(b)',side=3,adj=0,cex=1.0,line=0.5)
if (lplotfile) dev.off()

