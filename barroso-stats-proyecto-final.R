# ######################################### #
# k-means clustering of Mexican earthquakes #
# ######################################### #

# Made by: Roberto Barroso Fernández
# Instituto de Geociencias UNAM
# 
# Proyecto final: Herramientas Estadísticas de las Geociencias

# ######################################### #
#            Setting up the code            #
# ######################################### #

# Reemplazar path por el preferido
setwd('/Users/roberto/Documents/stats-final-project')

library(geodata)
library(ggplot2)
library(sf)
library(dplyr)
library(scales)

######### ######### ######### ######### ######### ######### ######### ######### 

# ######################################### #
#            Defining functions             #
# ######################################### #

plims <- function(x,p=0.1) {
	# limits of plot axis
	# x - numeric
	# p - percentage (from 0 to 1) of padding
	maxval <- max(x)
	minval <- min(x)
	span   <- maxval-minval
	lim_hi <- maxval + 0.5 * p * span
	lim_lo <- minval - 0.5 * p * span
	return (c(lim_lo,lim_hi))
}

magfreqs <- function(m,mbin=0.1) {
	# Generates data.frame of number of events against magnitude
	# Input:
	# m - list of magnitudes
	# mbin - binning width (difference)
	# Output:
	# df$magnitude - 
	# df$amount - 
	mags = seq(round(min(m),digits=1), round(max(m),digits=1), mbin)
	amt = 1:length(mags)
	for (i in 1:length(mags)){
		cutmag = mags[i]
		amt[i] = length(m[m>cutmag])
	}
	df <- data.frame(magnitude=mags,amount=amt)
	return (df)
}

kagan <- function(x) {
	# Kagan relation for maximum aftershock affectation area
	# x - array (can be single-element)
	# returns influence radius in kilometers
	kgn <- 20 * 10 ^ ((x - 6)/2)
	return (kgn)
}

sph_to_rct <- function(lat,lon,depth,r=6371) {
	# function to extract Cartesian coordinates from earthquake
	# epicenter and depth
	# Input (all arrays)
	# lat - latitude, in degrees
	# lon - longitude, in degrees
	# depth - in kilometers
	# r - Earth's mean radius, in kilometers
	# returns data.frame with x, y, z coordinates in kilometers
	rho   = r - depth
	theta = lon * pi / 180  # longitude in radians
	phi   = (90-lat) * pi / 180  # 90º-latitude in radians
	x <- rho * sin(phi) * cos(theta)
	y <- rho * sin(phi) * sin(theta)
	z <- rho * cos(phi)
	df <- data.frame(X=x,Y=y,Z=z)
	return (df)
}

timediff <- function(date1,date2) {
	# converts time difference to days
	# 1 day = 24 hrs = 1440 mins = 86400 secs
	delta <- date2 - date1
	delta_val <- as.numeric(delta) 
	delta_unit <- attr(delta,'unit') 
	if (delta_unit == 'secs') {
		delta_val <- delta_val / 86400
	} else if (delta_unit == 'mins') {
		delta_val <- delta_val / 1440 
	} else if (delta_unit == 'hours') {
		delta_val <- delta_val / 24
	} 
	return (delta_val)
}


######### ######### ######### ######### ######### ######### ######### ######### 

# Definin Gardner-Knopoff (1974) windows

# If declustering this way, it is assumed that events that aren't
# clustered are completely independent, which is probably not the case
# for Mexico (or at least, the Mexican Subduction Zone, MSZ). 
# Declustering based on Kagan radii (Avila-Barrientos et al., 2015), 
# then some other models, would be better, but due to time constraints,
# I am using Gardner-Knopoff. Hope this helps :)

# magnitudes
gk_m <- c(2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0) 
# radii in kilometers
gk_l <- c(19.5,22.5,26,30,35,40,47,54,61,70,81,94)
# time intervals in days
gk_t <- c(6,11.5,22,42,83,155,290,510,790,915,960,985)

gk <- data.frame(m=gk_m,l=gk_l,t=gk_t)

######### ######### ######### ######### ######### ######### ######### ######### 

# ######################################### #
#           Processing the catalog          #
# ######################################### #

# Loading the catalog
sismos <- read.csv("catalogo_sismos.csv",skip=4)[-(304906:304912),] %>%
          filter(Magnitud != 'no calculable')
sismos$utc <- paste(sismos$Fecha.UTC,sismos$Hora.UTC) %>% 
              strptime('%Y-%m-%d %H:%M:%S') %>%
              as.POSIXct(tz = 'UTC')
sismos$Magnitud <- as.numeric(sismos$Magnitud)
sismo_coords <- sph_to_rct(
	sismos$Latitud, sismos$Longitud, sismos$Profundidad
)
sismos$X <- sismo_coords$X
sismos$Y <- sismo_coords$Y
sismos$Z <- sismo_coords$Z 

# relative time to first event in catalog, in days
time <- 1:length(sismos$X)
for (i in 1:length(sismos$X)) {
	date1 = sismos$utc[1]
	date2 = sismos$utc[i]
	time[i] <- timediff(date1,date2) 
}
sismos$time <- time

# First step: taking into account completeness of catalog
mc = 3.1
mc_bool = sismos$Magnitud >= mc
sismos_mc <- sismos[mc_bool,]   # complete for older data
sismos_mc_len = length(sismos_mc[,1])
sismos_mc$idx <- 1:sismos_mc_len

# Second step: removing aftershocks
# sismos_mc_hi$kagan <- kagan(sismos_mc_hi$Magnitud)

# ######################################### #
#                Declustering               #
# ######################################### #

# Using the Gardner-Knopoff algorithm

gk_rad <- spline(x=gk$m, y=gk$l, method='natural', xout=sismos_mc$Magnitud)$y
gk_tdy <- spline(x=gk$m, y=gk$t, method='natural', xout=sismos_mc$Magnitud)$y

# 0 if not an aftershock, 1 if yes
aftershocks <- rep(0,sismos_mc_len)
for (ii in 1:sismos_mc_len) {
	if (ii %% 5000 == 0){
		print(ii)
	}
	if (ii == sismos_mc_len){
		sprintf('DONE!')
		break
	}
	mag <- sismos_mc$Magnitud[ii]
	rad <- gk_rad[ii]  # radius
	tdy <- gk_tdy[ii]  # time
	# if the event isn't already marked as an aftershock
	if (aftershocks[ii] == 0) {
		jj <- ii + 1
		# it has to fall within the acceptable time window for
		# declustering
		while ((sismos_mc$time[jj] - sismos_mc$time[ii]) < tdy) {
			# could be more elegant
			x1 = sismos_mc$X[ii] 
			y1 = sismos_mc$Y[ii] 
			z1 = sismos_mc$Z[ii]
			x2 = sismos_mc$X[jj] 
			y2 = sismos_mc$Y[jj] 
			z2 = sismos_mc$Z[jj]
			# distance between events
			rad_diff = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)
			# it has to fall withing the acceptable radius for
			# it to be marked as an aftershock  
			if (rad_diff < rad) {
				aftershocks[jj] = 1
			}
			if (jj == sismos_mc_len){
				break
			}
			jj <- jj + 1
		}
	}
}

sismos_mc$aftershocks <- aftershocks

sismos_clust <- sismos_mc[sismos_mc$aftershocks==1,]
sismos_bckgd <- sismos_mc[sismos_mc$aftershocks==0,]

######### ######### ######### ######### ######### ######### ######### #########

# Now, the statistical analysis part

# ######################################### #
#             k-means clustering            #
# ######################################### #

# Splitting catalog into two parts, one for more superficial events,
# the other for deeper events
# threshold set at 40 km, after Zuñiga et al. (2017)

# Choosing X, Y, Z, magnitude

# Magnitude is a function of fault size anyway...

cutoff_depth = 40

bool_depth = sismos_bckgd$Profundidad<=cutoff_depth

shallow <- sismos_bckgd[bool_depth,
                        c('Longitud','Latitud','Magnitud','X','Y','Z')]
deep    <- sismos_bckgd[!bool_depth,
                        c('Longitud','Latitud','Magnitud','X','Y','Z')]

shallow_kmeans <- kmeans(
	shallow[,c('Magnitud','X','Y','Z')],10,nstart=20
)
plot(shallow[,c('Longitud','Latitud')],col=shallow_kmeans$cluster)

######### ######### ######### ######### ######### ######### ######### #########

bckgnd_countries = c('US','GT','BZ','SV','HN','NI','CR','PA')
bckgnd <- do.call('rbind',lapply(bckgnd_countries,
                  function(x) gadm(country=x,level=0,resolution=2,
                  path=getwd()))) %>%
          st_as_sf(crs = st_crs('EPSG:6365'))
mexico <- gadm(country='MX',level=1,resolution=2,path=getwd()) %>%
          st_as_sf(crs = st_crs('EPSG:6365'))
ev_pts <- st_as_sf(sismos_bckgd, coords=c('Longitud','Latitud'), 
                   crs=st_crs('EPSG:6365'))

######### ######### ######### ######### ######### ######### ######### #########




# No ejecutar esto a menos que esté bien arreglado
magmap1 <- ggplot() +
           geom_sf(data=bckgnd,fill='lightgrey',color='grey') +
           geom_sf(data=mexico,fill='white',color='grey') +
           geom_sf(data=ev_pts,alpha=0.2,
                   aes(color=Profundidad)) +
           coord_sf(xlim=plims(sismos$Longitud,p=-0.1),
                    ylim=plims(sismos$Latitud,p=-0.1)) 
magmap1
# No ejecutar esto a menos que esté bien arreglado
