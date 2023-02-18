rm(list=ls())

lplotfile     = TRUE

ndelta        = 150
# ndelta        = 30
delta.max     = 3
l0.real       = 1
l.imag        = 0.01
power         = c(1,0.5,1.5)

ndim          = 3
npower        = length(power)

dir.Rlib  = NULL
source(paste(dir.Rlib,'pdf.eps.R',sep=''))
source(paste(dir.Rlib,'gev.R',sep=''))
source(paste(dir.Rlib,'gev.complex.R',sep=''))


delta.all     = seq(from=0,to=delta.max,length.out=ndelta+1)[-1]

apt.max       = array(NA,dim=c(ndelta,npower))
apt.pascal    = array(NA,dim=c(ndelta,npower))
apt.nifty     = array(NA,dim=c(ndelta,npower))
q.exact       = array(NA,dim=c(ndim,ndelta,npower))
q.pascal      = array(NA,dim=c(ndim,ndelta,npower))
apt.asymp     = array(NA,dim=c(npower))

for (np in 1:npower) for (nd in 1:ndelta) {
	delta   = delta.all[nd]
	lambda  = l0.real + delta * (1:ndim-1)^power[np]

	# SETUP MATRICES
	emat = array(NA,dim=c(ndim,ndim))
	for (j in 1:ndim) for (i in 1:ndim) emat[i,j] = 1/(lambda[i] + Conj(lambda[j]))
	amat = emat * emat
	bmat = emat
	
	# EXACT MAXIMUM
	apt.list         = gev.complex(amat,bmat)
	apt.max [ nd,np] = 2 * apt.list$values[1]
	q.exact [,nd,np] = Re(apt.list$vectors[,1]/apt.list$vectors[1,1])
	
	# GENERALIZED PASCAL MODE
	fmat               = array(NA,dim=c(ndim,ndim))
	fmat[1,]           = c(1,rep(0,ndim-1))
	for (i in 2:ndim)  fmat[i,] = lambda^(i-2)
	rhs                = c(1,rep(0,ndim-1))
	q.pascal  [,nd,np] = solve(fmat,rhs)
	apt.pascal[ nd,np] = 2*(q.pascal[,nd,np] %*% amat %*% q.pascal[,nd,np] ) / (q.pascal[,nd,np] %*% bmat %*% q.pascal[,nd,np])
	
	### NIFTY CALCULATION
	apt.nifty[ nd,np] = 0
	for (j in 1:ndim) for (i in j:ndim) apt.nifty[nd,np] = apt.nifty[nd,np] + 2/(lambda[i] + lambda[j])
	
	# if (nd == ndelta & np == 2) stop()
	
	# GENERALIZED GENERALIZED PASCAL MODE
	if (nd == 1) {
		# fmat               = array(NA,dim=c(ndim,ndim))
		# for (i in 1:ndim-1) fmat[i+1,] = lambda^i
		# fmat.inv           = solve(fmat)
		# for (i in 1:ndim)  fmat.inv[,i] = fmat.inv[,i]/fmat.inv[1,i]
		gmat1              = array(NA,dim=c(ndim,ndim))
		gmat2              = array(NA,dim=c(ndim,ndim))
		# for (j in 1:ndim-1) for (i in 1:ndim-1) gmat1[i+1,j+1] = factorial(i+j+1)  
		# for (j in 1:ndim-1) for (i in 1:ndim-1) gmat2[i+1,j+1] = factorial(i+j  ) 
		for (j in 1:ndim-1) for (i in 1:ndim-1) gmat1[i+1,j+1] = choose(i+j,i) * (i+j+1) 
		for (j in 1:ndim-1) for (i in 1:ndim-1) gmat2[i+1,j+1] = choose(i+j,i)
		gmat.apt           = gev(gmat1,gmat2)
		apt.asymp[np]      = gmat.apt$lambda[1]		
	}
	
}

fout = './figures/NDim.GenPascal'
if (lplotfile) pdf.eps(fout,'pdf',height =6,width=8.5)
par(mfcol=c(1,1),mar=c(4,4.5,1,0.2))
par(cex.lab=1.3,cex.axis=1.3,cex.main=1.3)
yrange = range(apt.max,apt.pascal,apt.asymp[1])
xrange = c(0,delta.max)
plot(1,1,type='n',xlim=xrange,ylim=yrange,xlab=bquote(Delta),ylab='APT')
for (np in 1:npower) lines(delta.all,apt.max   [,np],col=np,lwd=2)
for (np in 1:npower) lines(delta.all,apt.pascal[,np],col=np,lwd=2,lty='dashed')
legend.say = c(paste('power= ',power,'; exact',sep=''),paste('power= ',power,'; Generalized Pascal',sep=''))
legend.say = c(paste('APT optimal; p= ',power,sep=''),paste('Pascal; p= ',power,sep=''))
col.say    = c(1:npower,1:npower)
lty.say    = c(rep('solid',npower),rep('dashed',npower))
legend('topright',legend=legend.say,lwd=2,col=col.say,lty=lty.say,ncol=1,cex=1.2)

points(0,apt.asymp[1],pch=19,cex=1.3)
points(0,2*ndim-1,pch=19,cex=1.3)

# mtext('(a)',side=3,adj=0,cex=1.3,line=0.5)
# yrange = range(q.exact,q.pascal)
# plot(1,1,type='n',xlim=xrange,ylim=yrange,xlab=bquote(Delta),ylab='q')
# for (nd in 2:ndim) for (np in 1:npower) lines(delta.all,q.exact [nd,,np],col=np,lwd=2)
# for (nd in 2:ndim) for (np in 1:npower) points(0,q.pascal[nd,1,np],col=np,lwd=2,pch=19)
# mtext('(b)',side=3,adj=0,cex=1.3,line=0.5)
# abline(h=0,col='grey',lwd=2)
if (lplotfile) dev.off()

gmat.say = array(paste(gmat1,'-a',gmat2,' & '),dim=c(ndim,ndim))