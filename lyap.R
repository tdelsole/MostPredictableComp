lyap = function(e.vectors,e.values,qq,lstop.unstable=TRUE,e.square=FALSE) {
		
	jdim = length(e.values)
	
	if ( !all(Re(e.values)<0) ) {
		if (lstop.unstable) stop("system is unstable") else dum1=array(NA,dim=c(jdim,jdim))
	} else {
		### SOLVE S^-1 Q S^-H #####################
		dum1 = solve(e.vectors,qq)
		dum2 = solve(e.vectors,Conj(t(dum1)))
		dum1 = Conj(t(dum2))
	
		### CALCULATE SCHUR PRODUCT ###############
#		for ( j in 1:jdim) for (i in 1:jdim) dum2[i,j]=-dum1[i,j]/(e.values[i]+Conj(e.values[j]))

		x = rep(e.values,jdim)
		y = rep(Conj(e.values),each=jdim)
		dum2 = -dum1 /(x + y)
		
		if (e.square) dum2 = -dum2/(x+y)
	
		### INVERT E.MODE TRANSFORMATION ##########
		dum1 = e.vectors %*% dum2 %*% Conj(t(e.vectors))
	}
	
	dum1
	
}