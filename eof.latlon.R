eof.latlon <- function(lon, lat, data.array, neof=30, idate=NULL, dname=NULL, alpha=0.05, SubtractMean=TRUE) {

# COMPUTE PRINCIPAL COMPONENTS OF A DATA ARRAY.  

# INPUT:
#	LON[NLON]: VECTOR SPECIFYING LONGITUDES
#	LAT[NLAT]: VECTOR SPECIFYING LATITUDES
#	DATA.ARRAY[NLON,NLAT,NTOT] OR DATA.ARRAY[NLON*NLAT,NTOT]:  
#	  THE DATA ARRAY IN [SPACE1,SPACE2,TIME] OR [SPACE,TIME] FORMAT 
#	NEOF = NUMBER OF SPATIAL EOFS (WITH MASK) TO BE OUTPUTTED.  DEFAULT = 30
#   IDATE = A VECTOR INDICATING THE STARTING TIME AND INCREMENT (E.G., IDATE = C(1, "JAN", 1979, "1MO") )
#   DNAME = NAME OF DATA (E.G., HADSST, ERSSTV, NCEP/NCAR REAN TEMP, ETC)
#   ALPHA = SIZE OF THE CONFIDENCE INTERVAL OF EXPLAINED VARIANCES: (1-ALPHA)*100%

# OUTPUT LIST:
#	$EOF[NLON*NLAT,NEOF]: THE FIRST NEOF SCALED EOFS 
#	$PC[NTOT,NEOF]: THE PCS, NORMALIZED TO UNIT VARIANCE 
#	$FEXPVAR: FRACTION OF EXPLAINED VARIANCE FOR EACH EOF (ALL OF THEM).
#	$FEXPVAR.CI: CONFIDENCE INTERVALS FOR FRACTION OF EXPLAINED VARIANCE.
#	$EOFI[NLON*NLAT,NEOF]: PSEUDO INVERSE OF EOF (I.E.,T(EOFI) %*% EOF = I)
#	$NEOF: MINIMUM OF (NEOF IN ARGUMENT LIST, RANK OF DATA.ARRAY)
#	$BAD[NLON*NLAT]: LOGICAL ARRAY INDICATING POINTS DROPPED FROM THE ANALYSIS
#	$WEIGHT[NLON*NLAT]: THE WEIGHTS USED TO COMPUTE THE EOF
#	$LON: THE LONGITUDES OF THE SPATIAL FIELD
#	$LAT: THE LATITUDES OF THE SPATIAL FIELD
#	$SVAL: THE SINGULAR VALUES OF THE SCALED ANOMALIES

# DEFINE PARAMETERS
nlon = length(lon)
nlat = length(lat)
ntot = length(data.array)/(nlon*nlat)
if ( length(data.array) %% nlon*nlat != 0 ) stop("data.array dimension not integral multiple of nlon*nlat")

# DEFINE WEIGHTING
weight =rep(sqrt(cos(lat*pi/180)),each=nlon)

# IDENTIFY MISSING OR 0-WEIGHTED DATA
dim(data.array) = c(nlon*nlat,ntot)
lbad = is.na(rowSums(data.array)) | weight == 0 

# REMOVE CLIMATOLOGY
if (SubtractMean) {
	clim       = rowMeans(data.array)
	data.array = data.array - clim	
} else {
	clim       = rep(NA,nlon*nlat)
}

# NORMALIZE/WEIGHT THE VARIABLES
data.array = data.array * weight

# COMPRESS DATA BY ELIMINATING MISSING GRID POINTS
data.array = data.array[!lbad,] 
ndef = sum(!lbad)
                
# COMPUTE SVD
mmin = min(ndef,ntot,neof)   
data.svd = svd(data.array,nu=mmin,nv=mmin)

# COMPUTE FEXPVAR
fexpvar = data.svd$d^2/sum(data.svd$d^2)
    
# COMPUTE CONFIDENCE INTERVAL OF FEXPVAR
stderr     = qnorm(alpha/2,lower.tail=FALSE) * sqrt(2/ntot)
stderr     = sqrt(2/ntot)
fexpvar.ci = cbind ( fexpvar * (1 - stderr) , fexpvar * (1 + stderr) )
    
# COMPUTE PCS
pc = sqrt(ntot-1)*data.svd$v[,(1:mmin)]
    
# FILL IN MISSING POINTS IN EOF 
eof  = array(NA,dim=c(nlon*nlat,mmin))
eofi = array(NA,dim=c(nlon*nlat,mmin))        
for ( n in 1:mmin ) eof [!lbad,n] = data.svd$u[,n]/weight[!lbad]*data.svd$d[n]/sqrt(ntot-1)
for ( n in 1:mmin ) eofi[!lbad,n] = data.svd$u[,n]*weight[!lbad]/data.svd$d[n]*sqrt(ntot-1)
    
list(eof=eof,pc=pc,fexpvar=fexpvar,fexpvar.ci=fexpvar.ci,eofi=eofi,neof=mmin,lbad=lbad,weight=weight,lon=lon,lat=lat,
   sval=data.svd$d,idate=idate,dname=dname,clim=clim)
}

