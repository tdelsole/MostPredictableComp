plot_eof <- function(data.eof,npic,data.lon,data.lat) {
	
library(maps)
library(fields)
source("/Users/delsole/R/statbook/eof/bluered.R")

expvar = round(100*data.eof$sval[npic]^2/sum(data.eof$sval^2))
rg <- range(data.eof$eof[,npic],na.rm=TRUE)
ntot = length(data.eof$pc[,1])

n = 10
lowerlim <- rg[1]
upperlim <- rg[2]
maximum <- signif(max(abs(c(lowerlim,upperlim))),2)
breaks <- seq(-maximum,maximum,(maximum-(-maximum))/n)


# breaks = seq(from=-1,to=1,by=0.2)
# n = length(breaks)-1

colours = bluered(seq(0,n-1,by=1), "linear", yellow =TRUE)


layout(matrix(1:2,ncol=1,nrow=2),heights=c(4,4))
par(mar=c(4,5,4,1))

dum = array(data.eof$eof[,npic],dim=c(length(data.lon),length(data.lat)))

title = paste("EOF-",npic," of ",data.eof$dname," (",expvar,"%)",sep="")
image.plot(data.lon,data.lat,dum,col=colours,xlab="lon",ylab="lat",breaks=breaks,main=title)
contour(data.lon,data.lat,dum,add=T,xlab=lon,ylab=lat,labcex=0.9)

if ( data.lon[1] < 0 ) map('world',interior=FALSE,add=T) else map('world2',interior=FALSE,add=T)

if ( data.eof$idate[4] == "1mo") incr = 1/12
if ( data.eof$idate[4] == "1yr") incr = 1
if ( data.eof$idate[4] == "6mo") incr = 1/2


par(mar=c(6,5,3,1))
year = seq(from=as.numeric(data.eof$idate[3]),by=incr,length.out=ntot)
# plot(year[year>=1979 & year < 2001],data.eof$pc[year>=1979 & year < 2001,npic],type="l",xlab="Year",ylab="PC")
plot(year,data.eof$pc[,npic],type="l",xlab="Year",ylab="PC")
}


