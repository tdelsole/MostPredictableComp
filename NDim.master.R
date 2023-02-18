rm(list=ls())

lplotfile     = TRUE

ndelta        = 150
delta.max     = 3
l0.real       = 1
l.imag        = 0.01
real.imag.vec = c('rr','rc','rrc','rcc')
real.imag.vec = c('rr','rrr','rrrr','rrrrr')

dir.Rlib  = NULL
source(paste(dir.Rlib,'pdf.eps.R',sep=''))
source(paste(dir.Rlib,'gev.R',sep=''))
source(paste(dir.Rlib,'gev.complex.R',sep=''))



ncases        = length(real.imag.vec)
delta.all     = seq(from=0,to=delta.max,length.out=ndelta+1)[-1]

ndim.max      = sum(ifelse(strsplit(real.imag.vec[ncases],split='')[[1]] == 'r',1,2))

q.pascal      = array(NA,dim=c(ndim.max,ndim.max))
for (nd in 1:ndim.max) q.pascal[1:nd,nd] = (-1)^(1:nd-1) * choose(nd-1,1:nd-1)

apt.pascal    = array(NA,dim=c(ndim.max,ndelta,ncases))
apt.max       = array(NA,dim=c(ndelta,ncases))
apt.mode      = array(NA,dim=c(ndelta,ncases))
apt.upper     = array(NA,dim=c(ndelta,ncases))
apt.pascal.it = array(NA,dim=c(ndelta,ncases))
qvec.max      = array(NA,dim=c(ndim.max,ndelta,ncases))
lambda.at.max = array(NA,dim=c(ndim.max,ncases))
for (nc in 1:ncases) for (nd in 1:ndelta) {
	real.imag.split = strsplit(real.imag.vec[nc],split='')[[1]]
	if (real.imag.split[1] != 'r') stop('first eigenvalue assumed to be real')
	delta      = delta.all[nd]
	lambda     = l0.real
	for (ctype in real.imag.split[-1]) {
		nl = length(lambda)
		if (ctype == 'r') {
			lambda = c(lambda,Re(lambda[nl]) + delta)
		} else {
			lambda = c(lambda,Re(lambda[nl]) + delta + 1i * l.imag, Re(lambda[nl])+ delta - 1i * l.imag)
		}
	}
	
	
	
	ndim = length(lambda)
	if (ndim != nc + 1) stop('miscounted dimension')
	if (nd == 1) lambda.at.max[1:ndim,nc] = lambda

	emat = array(NA,dim=c(ndim,ndim))
	for (j in 1:ndim) for (i in 1:ndim) emat[i,j] = 1/(lambda[i] + Conj(lambda[j]))
	amat = emat * emat
	bmat = emat
	
	apt.list  = gev.complex(amat,bmat)

	### CHECK EIGENVECTOR ###
	lhs = amat %*% apt.list$vectors
	rhs = bmat %*% apt.list$vectors * rep(apt.list$values,each=ndim)
	# if (!all.equal(lhs,rhs)) stop("inconsistent GEV solution")
	
	apt.max[nd,nc]         = 2 * apt.list$values[1]
	qvec.max[1:ndim,nd,nc] = Re(apt.list$vectors[,1]/apt.list$vectors[1,1])
	
	### PASCAL MODE
	for (n in 2:ndim) apt.pascal[n,nd,nc] = 2*
	    (q.pascal[1:n,n] %*% amat[1:n,1:n] %*% q.pascal[1:n,n] ) / 
	    (q.pascal[1:n,n] %*% bmat[1:n,1:n] %*% q.pascal[1:n,n])

	apt.mode[nd,nc] = sum(1/lambda)
	apt.pascal.it[nd,nc] = 2*(q.pascal[1:ndim,ndim] %*% amat %*% q.pascal[1:ndim,ndim])/(q.pascal[1:ndim,ndim] %*% bmat %*% q.pascal[1:ndim,ndim])
	
	apt.upper[nd,nc] = 2 * sum(1/lambda)
}

################################
########## SOLVE APT ANALYTICALLY FOR 2D
################################
apt.2D = array(NA,dim=c(ndelta,2))
for (nd in 1:ndelta) {
	lambda = c(l0.real,l0.real + delta.all[nd])
	i = 1
	j = 2
	beta           = (4*lambda[i]*lambda[j] + (lambda[i] + lambda[j])^2)/(4*lambda[i]*lambda[j]*abs(lambda[i]+lambda[j]))
	bcoef          = -2*beta
	ccoef          =    beta/abs(lambda[i]+lambda[j])
	apt.2D[nd,] = -bcoef + c(1,-1) * sqrt(bcoef^2-4*ccoef)
}


################################
########## PLOTS
################################
fout = './figures/NDim.master'
if (lplotfile) pdf.eps(fout,'pdf',height = 5,width=8.5)
par(mfcol=c(1,2),mar=c(4,4.5,3,0.2))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
yrange = c(0.5,max(apt.max)+0.5)
xrange = c(0,delta.max)
ysay   = bquote('maximum APT '*lambda[0])
ysay   = 'maximum APT'
plot(1,1,type='n',xlab=bquote(Delta),ylab=ysay,xlim=xrange,ylim=yrange)
abline(h=1,col='grey50',lwd=2,lty='dashed')
mtext('(a)',side=3,adj=0,cex=1.3,line=0.5)
for (nc in 1:ncases) lines(delta.all,apt.max[,nc],type='l',col=nc,lwd=2)
for (nc in 1:ncases) points(delta.all[1],apt.max[1,nc],pch=19,cex=0.7,col=nc)
for (nc in 1:ncases) text(delta.all[1],apt.max[1,nc],round(apt.max[1,nc],1),pos=3,col=nc)
delta.choose = 0.75
nd.pic       = max(which(delta.all < delta.choose))
for (nc in 1:ncases) text(delta.all[nd.pic],apt.max[nd.pic,nc],paste('D=',nc+1,sep=''),pos=3,col=nc)
title(main='Maximum Average\nPredictability Time',line=0.5)

nd.pic = seq(from=1,to=ndelta,by=10)
points(delta.all[nd.pic],apt.2D[nd.pic,1],pch=4,col='black',cex=1.3,lwd=2)


ymax   = max(abs(qvec.max[,1,]),na.rm=TRUE)
yoff   = 7
yrange = c(-2,yoff*ncases+ymax-1)
xrange = c(1,ndim.max)
ylab.say = bquote(hat(q))
par(mar=c(4,3,3,0.2))
plot(1,1,type='n',xlab='vector element',ylab=ylab.say,xlim=xrange,ylim=yrange,yaxt='n')
mtext(ylab.say,side=2,cex=1.5)
for (nc in 1:ncases) lines(qvec.max[,1,nc]+yoff*(nc-1),col=nc,lwd=2)
for (nc in 1:ncases) for (nd in 1:(nc+1)) text(nd,qvec.max[nd,1,nc]+yoff*(nc-1),round(qvec.max[nd,1,nc],1),pos=ifelse(nd %% 2,3,1),col=nc)
mtext('(b)',side=3,adj=0,cex=1.3,line=0.5)
title.top = bquote('Optimal Projection Vector')
title.bot = bquote(Delta==.(delta.all[1]))
title(main=title.top,line=2.0)
title(main=title.bot,line=0.7)

if (lplotfile) dev.off()