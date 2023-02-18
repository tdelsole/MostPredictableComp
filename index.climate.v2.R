index.climate.v2 = function(lon,lat,field=NULL,area.name=NULL,ocean=NULL,lplot=FALSE) {
## THIS FUNCTION DEFINES A MASK FOR SELECTING GEOGRAPHIC DOMAINS
## AND COMPUTES AREA AVERAGE WITHIN THAT DOMAIN
## INPUT:
##  LON[NLON]: LONGITUDE VALUES OF THE FIELD
##  LAT[NLAT]: LATITUDE VALUES OF THE FIELD
##  FIELD[NLON,NLAT,NSTEP]: FIELD TO AVERAGE (MISSING POINTS REMAIN MISSING)
##  AREA.NAME: A CHARACTER STRING DEFINING THE GEOGRAPHIC DOMAIN TO SELECT
##  OCEAN[NLON,NLAT]: A LOGICAL MASK = TRUE AT OCEAN POINTS; IF DEFINED, FALSE POINTS WILL BE OMITTED
##  LPLOT: LOGICAL INDICATING WHETHER LOGICAL MASK SHOULD BE PLOTTED
## OUTPUT:
##  INDEX[NSTEP]: INDEX CALCULATED AT EACH TIME STEP OF THE FIELD
##  VIEW[NLON,NLAT]: LOGICAL ARRAY.  
##     =TRUE IF THE GRID POINT IS TO BE USED FOR ANALYSIS
##     =FALSE IF THE GRID POINT IS TO BE IGNORED FOR ANALYSIS
## 
## COMMENTS: 
##  1) DEPENDENCIES: PLOT_LATLON_V4 IF LPLOT=TRUE

master.list = NULL
master.list = rbind(master.list,list(area.name='NASST'         ,lon.mn=285,lon.mx=352,lat.mn=  0,lat.mx=60))
master.list = rbind(master.list,list(area.name='NINO34'        ,lon.mn=190,lon.mx=240,lat.mn= -5,lat.mx= 5))
master.list = rbind(master.list,list(area.name='SST30S60N'     ,lon.mn=  0,lon.mx=360,lat.mn=-30,lat.mx=60))
master.list = rbind(master.list,list(area.name='SST40S60N'     ,lon.mn=  0,lon.mx=360,lat.mn=-40,lat.mx=60))
master.list = rbind(master.list,list(area.name='NPSST'         ,lon.mn=100,lon.mx=265,lat.mn=  0,lat.mx=60))
master.list = rbind(master.list,list(area.name='GLO'           ,lon.mn=  0,lon.mx=360,lat.mn=-90,lat.mx=90))
master.list = rbind(master.list,list(area.name='OCEAN'         ,lon.mn=  0,lon.mx=360,lat.mn=-90,lat.mx=90))
master.list = rbind(master.list,list(area.name='VIRGINIA'      ,lon.mn=276,lon.mx=284,lat.mn= 36,lat.mx=40))
master.list = rbind(master.list,list(area.name='NAMERICA'      ,lon.mn=190,lon.mx=310,lat.mn= 15,lat.mx=70))
master.list = rbind(master.list,list(area.name='INDOPACIFIC'   ,lon.mn= 25,lon.mx=300,lat.mn=-90,lat.mx=90))
master.list = rbind(master.list,list(area.name='PACIFIC'       ,lon.mn=125,lon.mx=300,lat.mn=-90,lat.mx=90))
master.list = rbind(master.list,list(area.name='NOTPACIFIC'    ,lon.mn= 25,lon.mx=300,lat.mn=-90,lat.mx=90))
master.list = rbind(master.list,list(area.name='TEXAS'         ,lon.mn=254,lon.mx=266,lat.mn= 26,lat.mx=36))
master.list = rbind(master.list,list(area.name='NOANTARCTIC'   ,lon.mn=  0,lon.mx=360,lat.mn=-60,lat.mx=90))
master.list = rbind(master.list,list(area.name='PACIFIC30S60N' ,lon.mn=125,lon.mx=300,lat.mn=-30,lat.mx=60))
master.list = rbind(master.list,list(area.name='PACIFIC30S30N' ,lon.mn=125,lon.mx=300,lat.mn=-30,lat.mx=30))
master.list = rbind(master.list,list(area.name='ATLANTIC30S60N',lon.mn=275,lon.mx=352,lat.mn=-30,lat.mx=60))
master.list = rbind(master.list,list(area.name='INDOPAC30S60N' ,lon.mn= 25,lon.mx=300,lat.mn=-30,lat.mx=60))
master.list = rbind(master.list,list(area.name='CONUS'         ,lon.mn=190,lon.mx=310,lat.mn= 25,lat.mx=50))

npic = which(area.name == as.character(master.list[,"area.name"]))
if ( is.null(area.name) | length(npic) == 0) {print(master.list); stop('select area.name from the above list')}

nlon = length(lon); nlat = length(lat)
lonp = lon; lonp[lon<0] = lon[lon<0] + 360

lon.pic = ( lonp >= master.list[npic,'lon.mn'] & lonp <= master.list[npic,'lon.mx'])
lat.pic = ( lat  >= master.list[npic,'lat.mn'] & lat  <= master.list[npic,'lat.mx'])

### DEFINE MASTER MASK
view = array(FALSE,dim=c(nlon,nlat))
view[lon.pic,lat.pic] = TRUE

if (!is.null(ocean)) land = !ocean

if (!is.null(ocean)) {
## REMOVE MEDITERRANEAN, BALTIC, CASPIAN SEAS FROM OCEAN FIELD
	dim(ocean) = c(nlon,nlat)
	euro.lon = lon >= 0  & lon <= 60
	euro.lat = lat >= 10 & lat <= 50
	ocean[euro.lon,euro.lat] = FALSE 
	
	caspian.lon = lon >= 250 & lon <= 290
	caspian.lat = lat >= 40  & lat <= 70
	ocean[caspian.lon,caspian.lat] = FALSE
	
	baltic.lon  = lon >= 10 & lon <= 30
	baltic.lat  = lat >= 50 & lat <= 70
	ocean[baltic.lon,baltic.lat] = FALSE
	
	## INCLUDE ONLY OCEAN, DEPENDING ON DOMAIN
	if ( area.name == 'NINO34'    | 
	     area.name == 'NASST'     | 
	     area.name == 'SST30S60N' | 
	     area.name == 'SST40S60N' |
	     area.name == 'NPSST'     | 
	     area.name == 'OCEAN'     |
	     area.name == 'PACIFIC'   |
	     area.name == 'INDOPACIFIC'|
	     area.name == 'PACIFIC30S60N' |
	     area.name == 'ATLANTIC30S60N' |
	     area.name == 'INDOPAC30S60N' ) view = view & ocean
}

if ( area.name == 'INDOPACIFIC' | area.name == 'NOTPACIFIC' | area.name == 'PACIFIC' | area.name == 'PACIFIC30S60N' |
     area.name == 'INDOPAC30S60N') {
	atlsw.lon = lon >= 275 & lon <= 350
	atlsw.lat = lat >=  10 & lat <= 90
	view[atlsw.lon,atlsw.lat] = FALSE
	
	gulfmex.lon = lon >= 260 & lon <= 350
	gulfmex.lat = lat >=  18 & lat <= 90
	view[gulfmex.lon,gulfmex.lat] = FALSE
}

if ( area.name == 'NOTPACIFIC' ) view = ocean & !view

if ( area.name == 'ATLANTIC30S60N' ) {
	peru.lon = lon > 265 & lon < 300
	peru.lat = lat > -40 & lat < 10
	view[peru.lon,peru.lat] = FALSE
}

if ( area.name == 'NAMERICA' | area.name == 'CONUS') {
	cuba.lon = lon >= 275 & lon <= 310
	cuba.lat = lat >= 15  & lat <= 25
	view[cuba.lon,cuba.lat] = FALSE
	
	bahama.lon = lon >= 280 & lon <= 300
	bahama.lat = lat >=  22 & lat <=  30
	view[bahama.lon,bahama.lat] = FALSE
	
	baha.lon = lon >= 240 & lon <= 248
	baha.lat = lat >=  22 & lat <=  30
	view[baha.lon,baha.lat] = FALSE
		
	baha.lon = lon >= 240 & lon <= 251
	baha.lat = lat >=  22 & lat <=  27
	view[baha.lon,baha.lat] = FALSE
	
	greenland.lon = lon >= 300 & lon <= 360
	greenland.lat = lat >= 60  & lat <= 90
	view[greenland.lon,greenland.lat] = FALSE
	
	newfoundland.lon = lon >= 300 & lon <= 360
	newfoundland.lat = lat >=  45 & lat <=  52.5
	view[newfoundland.lon,newfoundland.lat] = FALSE
	
	hawaii.lon = lon >= 180 & lon <= 210
	hawaii.lat = lat >= 15 & lat <= 25
	view[hawaii.lon,hawaii.lat] = FALSE
	
	if (!is.null(ocean)) view = view & land
}

if (area.name == 'TEXAS') {
	gulf.lon = lon >= 262 & lon <= 266
	gulf.lat = lat >= 26  & lat <= 30
	view[gulf.lon,gulf.lat] = FALSE
}


## RESHAPE FIELD, DEFINE 'VIEW' TO INCLUDE ONLY NON-MISSING GRID POINTS OF FIELD
if ( !is.null(field)) {
	if ( length(field) %% nlon*nlat != 0 ) stop('dimension of field inconsistent with grid')
	nstep      = length(field)/nlon/nlat
	dim(field) = c(nlon*nlat,nstep)
	view       = view & !is.na(field[,1]) # mask out based on 'view' and missing data in 'field'
} 

weight  = rep(cos(lat*pi/180),each=nlon)
# weight  = rep(1,nlon*nlat)
weight  = weight / sum(weight[view])
weight[!view] = NA

if (!is.null(field)) index = colSums(field[view,,drop=FALSE]*weight[view]) else index = NULL

if ( lplot ) plot_latlon_v4(lon,lat,view)

list(index=index,view=view,weight=weight,metadata=master.list[npic,])

}