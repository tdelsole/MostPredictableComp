gev.complex = function(amat,bmat) {
### SOLVES GENERALIZED EIGENVALUE PROBLEM Ax = lambda Bx
	b.chol     = cholesky.complex(bmat)
	xmat       = forwardsolve.complex(b.chol,amat)
	xmat       = forwardsolve.complex(b.chol,t(Conj(xmat)))
	xmat.eigen = eigen(xmat,symmetric=TRUE)
	xmat.eigen$vectors  = backwardsolve.complex(t(Conj(b.chol)),xmat.eigen$vectors)
	return(xmat.eigen)
}



cholesky.complex     = function(bmat) {
	nrow.b = dim(bmat)[1]
	for (j in 1:nrow.b) {
		n1 = j:nrow.b
		n2 = 1:(j-1)
		if (j > 1) bmat[n1,j] = bmat[n1,j,drop=FALSE] - bmat[n1,n2,drop=FALSE] %*% t(Conj(bmat[j,n2,drop=FALSE]))
		bmat[n1,j] = bmat[n1,j]/sqrt(bmat[j,j])
	}
	for (j in 1:(nrow.b-1)) for (i in (j+1):nrow.b) bmat[j,i] = 0
	return(bmat)
}

forwardsolve.complex = function(lower,rhs) {
	nrow.low = dim(lower)[1]
	ncol.rhs = length(rhs)/nrow.low
	dim(rhs) = c(nrow.low,ncol.rhs)
	xmat     = rhs/lower[1,1]
	for (nc in 1:ncol.rhs) for (m in 2:nrow.low) xmat[m,nc] = (rhs[m,nc] - sum( lower[m,1:(m-1)] * xmat[1:(m-1),nc] ) ) / lower[m,m]
	return(xmat)
}

backwardsolve.complex = function(upper,rhs) {
	nrow.upp = dim(upper)[1]
	ncol.rhs = length(rhs)/nrow.upp
	dim(rhs) = c(nrow.upp,ncol.rhs)
	xmat     = rhs/upper[nrow.upp,nrow.upp]
	for (nc in 1:ncol.rhs) for (m in (nrow.upp-1):1) xmat[m,nc] = (rhs[m,nc] - sum( upper[m,(m+1):nrow.upp] * xmat[(m+1):nrow.upp,nc])) / upper[m,m]
	return(xmat)
}
