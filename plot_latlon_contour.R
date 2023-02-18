plot_latlon_contour = function(lon,lat,field,nbreaks=10,plot.legend=TRUE,title.say=NULL,
    lmaskzero=TRUE,shrinkdomain=FALSE,plot.continents=TRUE,force.pm=FALSE,
    cbar.list=NULL,suppress.axislab=FALSE,rotate.land=FALSE,lcontour=FALSE,contour.labels=FALSE) {
	
## GENERATES SPATIAL PLOTS ON A GLOBE WITH SUPERIMPOSED CONTINENTAL MAP
## LON[NLON]: LONGITUDE COORDINATES
## LAT[NLAT]: LATITUDE COORDINATES
## NBREAKS: NUMBER OF BREAK POINTS FOR COLOR BAR (DEFAULT=10)
## PLOT.LEGEND: =TRUE TO PLOT THE LEGEND, =FALSE TO SUPPRESS LEGEND
## TITLE.SAY: CHARACTER STRING FOR THE TITLE
## LMASKZERO: =TRUE TO SET ZERO TO 'WHITE', OTHERWISE SET TO 'GREY'
## SHRINKDOMAIN: =TRUE TO RE-ADJUST LAT/LON TO SHOW ONLY DEFINED POINTS
## PLOT.CONTINENTS: = TRUE TO PLOT CONTINENTS, OTHERWISE CONTINENTS ARE NOT PLOTTED
## FORCE.PM: FORCE COLOR BAR TO RANGE FROM NEGATIVE TO POSITIVE VALUES

nlon = length(lon)
nlat = length(lat)
if ( length(field) %% (nlon*nlat) != 0) stop('field not dimensioned correctly')


all.na        = all(is.na(field))
if (!all.na) {
	if (all(field[!is.na(field)] == 0)) {
		all.constant = TRUE
	} else {
		all.constant  = (var(as.numeric(field),na.rm=TRUE)/max(abs(field),na.rm=TRUE) <= 1.e-15 ) 
	}
} else all.constant = TRUE
# all.constant  = (var(as.numeric(field),na.rm=TRUE)/max(abs(field),na.rm=TRUE) <= 1.e-15 ) & !all.na
nbreaks.half  = floor(nbreaks/2)

#####################################################################
## IMPORT COLOR BAR PARAMETERS
#####################################################################
if (!is.null(cbar.list)) {
	plot.cols   = cbar.list$plot.cols
	plot.breaks = cbar.list$plot.breaks
}

#####################################################################
## IF FIELD IS ALL NAS, MAKE UP DATA SO THAT A BLANK PLOT IS FORMED
#####################################################################
if (all.na) {
	dim(field) = c(length(lon),length(lat))
	field[1,1] = 0
	plot.breaks = c(-1,0,1)
	plot.cols   = c("white","white")	

} else if (is.null(cbar.list)) {
	#####################################################################
	## IF FIELD IS A CONSTANT, THEN MAKE UP A RANGE OF DATA
	#####################################################################
	if ( all.constant) {
		constant.value = mean(as.numeric(field),na.rm=TRUE)
		if (constant.value == 0 ) {
			plot.breaks = seq(from=0,to=1,length.out=nbreaks.half) 
		} else { 
			plot.breaks = seq(from = 0, to = constant.value,length.out = nbreaks.half)
			plot.breaks[nbreaks.half] = plot.breaks[nbreaks.half]*1.001
		}
	} else {
		plot.breaks = pretty(abs(field),n=floor(nbreaks.half))
	}
	
	#####################################################################
	## GENERATE SYMMETRIC BREAKS ABOUT 0, OR INSERT ZERO IN THE BREAKS
	#####################################################################
	l.pm   = min(field,na.rm=TRUE)*max(field,na.rm=TRUE) < 0
	if (l.pm | force.pm) {
		plot.breaks = c(-rev(plot.breaks),plot.breaks)
		plot.breaks = plot.breaks[ plot.breaks != 0 ] ; # get rid of zero break points
	} else {
		plot.breaks = c(0,plot.breaks)
	}
	plot.breaks = sort(union(plot.breaks,plot.breaks)); #get rid of redundant breaks

	#####################################################################
	## GENERATE COLORS
	#####################################################################
	if (!all.na) {
		nb          = length(plot.breaks)
		plot.cols   = colorRampPalette(c('blue','deepskyblue','white','orange','red'),space='rgb')(nb-1)
	}

	#####################################################################
	## SET ZERO INTERVAL TO WHITE IF DESIRED
	#####################################################################
	nb                 = length(plot.breaks)
	break.include.zero = plot.breaks[1:(nb-1)] <=0 & plot.breaks[2:nb] >= 0
	if (lmaskzero) {
		plot.cols[break.include.zero] = 'white'
	} else {
		plot.cols[break.include.zero] = 'grey'
	}
	
}

#####################################################################
## DEFINE DUMMY ARRAY FOR PLOTTING
#####################################################################
plot.var = array(field,dim=c(length(lon),length(lat)))

#####################################################################
## REVERSE LATITUDES IF THEY GO FROM NORTH TO SOUTH
#####################################################################
if ( lat[2] - lat[1] < 0 ) {
	plot.var = plot.var[,length(lat):1]
	lat.say  = lat[length(lat):1]
	print('latitudes reversed')
} else {
	lat.say  = lat
}

#####################################################################
## ROTATE TO AVOID CUTTING ACROSS LAND MASSES
#####################################################################
if (rotate.land & lon[1] >= 0 ) {
	nlon     = length(lon)
	nst      = max(which(lon <= 180))
	npic     = c((nst+1):nlon,1:nst)
	plot.var = plot.var[npic,]
	lon      = c(lon[(nst+1):nlon] - 360,lon[1:nst])
}


#####################################################################
## SHRINK DOMAIN IF DATA IS MISSING 
#####################################################################
if (suppress.axislab) {xaxt='n';yaxt='n'} else {xaxt='s';yaxt='s'}
if (shrinkdomain) {
	good  = which(!is.na(plot.var),arr.ind=TRUE)
	xpic  = range(good[,'row'])
	ypic  = range(good[,'col'])
	xlim = lon[xpic]; ylim = lat.say[ypic]
	# image(lon,lat.say,plot.var,col=plot.cols,breaks=plot.breaks,xlab="",ylab="",xlim=xlim,ylim=ylim,xaxt=xaxt,yaxt=yaxt)
	
	# ###changed by Xiaoqin
	image(lon,lat.say,array(1,dim=dim(plot.var)),col=rgb(0.9,0.9,0.9),xlab="",ylab="",axes=FALSE)
	par(new=TRUE)
	# image(lon,lat.say,plot.var,col=plot.cols,breaks=plot.breaks,xlab="",ylab="",xlim=xlim,ylim=ylim,axes=FALSE)
    contour(lon,lat.say,plot.var,col=plot.cols,levels = plot.breaks,
       nlevels=length(plot.breaks),labcex=0.65,xlab='',ylab='',xlim=xlim,ylim=ylim,
       xaxs='i',yaxs='i',xaxt=xaxt,yaxt=yaxt,frame.plot=TRUE,drawlabels=contour.labels)
    par(new=TRUE)
	# .filled.contour(lon,lat.say,plot.var,col=plot.cols,levels = plot.breaks,nlevels=length(plot.breaks),xlab='',ylab='')
	.filled.contour(lon,lat.say,plot.var,col=plot.cols,levels = plot.breaks)
	par(new=TRUE)
	if(lcontour) contour(lon,lat.say,plot.var,col='black',levels = plot.breaks,
	  nlevels=length(plot.breaks),labcex=0.65,xlab='',ylab='',xlim=xlim,ylim=ylim,xaxt=xaxt,yaxt=yaxt,drawlabels=contour.labels)
	par(new=TRUE)
	
} else {
	# image(lon,lat.say,plot.var,col=plot.cols,breaks=plot.breaks,xlab="",ylab="",xaxt=xaxt,yaxt=yaxt)

	# ###changed by Xiaoqin
	image(lon,lat.say,array(1,dim=dim(plot.var)),col=rgb(0.9,0.9,0.9),xlab="",ylab="",axes=FALSE)
	par(new=TRUE)
	# image(lon,lat.say,plot.var,col=plot.cols,breaks=plot.breaks,xlab="",ylab="",axes=FALSE)
	
    contour(lon,lat.say,plot.var,col=plot.cols,levels = plot.breaks,nlevels=length(plot.breaks),labcex=0.65,xlab='',ylab='',xlim=xlim,ylim=ylim,xaxt=xaxt,yaxt=yaxt)
    par(new=TRUE)
	# .filled.contour(lon,lat.say,plot.var,col=plot.cols,levels = plot.breaks,nlevels=length(plot.breaks),xlab='',ylab='')
	.filled.contour(lon,lat.say,plot.var,col=plot.cols,levels = plot.breaks)
	par(new=TRUE)
	contour(lon,lat.say,plot.var,col=plot.cols,levels = plot.breaks,nlevels=length(plot.breaks),labcex=0.65,xlab='',ylab='',xlim=xlim,ylim=ylim,xaxt=xaxt,yaxt=yaxt)
	par(new=TRUE)
}


#####################################################################
## ADD COLOR BAR 
#####################################################################
if (plot.legend) addMapLegend(cutVector=plot.breaks,colourVector=plot.cols,legendLabels='all',
   legendMar=5,horizontal=TRUE,labelFontSize=1.5,catMethod='pretty')

#####################################################################
## ADD TITLE
#####################################################################
if (!is.null(title.say)) title(main=as.character(title.say))

#####################################################################
## ADD CONTINENTS
#####################################################################
if (plot.continents) if ( lon[1] < 0 ) map('world',interior=FALSE,add=T) else map('world2',interior=FALSE,add=T)

#####################################################################
## OUTPUT PLOTTING PARAMETERS
#####################################################################
list(plot.breaks=plot.breaks,plot.cols=plot.cols,xlim=xlim,ylim=ylim)


}