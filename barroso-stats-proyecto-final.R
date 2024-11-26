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

######### ######### ######### ######### ######### ######### ######### ######### 

######### ######### ######### ######### ######### ######### ######### ######### 

# ######################################### #
#           Processing our catalog          #
# ######################################### #

# Loading the catalog
sismos <- read.csv("catalogo_sismos.csv",skip=4)[-(304906:304912),] %>%
          filter(Magnitud != 'no calculable')
sismos$utc <- paste(sismos$Fecha.UTC,sismos$Hora.UTC) %>% 
              strptime('%Y-%m-%d %H:%M:%S') %>%
              as.POSIXct(tz = 'UTC')
sismos$Magnitud <- as.numeric(sismos$Magnitud)

# First step: taking into account completeness of catalog
mc_bool = sismos$Magnitud >= 3.1
sismos_mc_hi <- sismos[mc_bool,]   # complete for older data



######### ######### ######### ######### ######### ######### ######### #########

bckgnd_countries = c('US','GT','BZ','SV','HN','NI','CR','PA')
bckgnd <- do.call('rbind',lapply(bckgnd_countries,
                  function(x) gadm(country=x,level=0,resolution=2,
                  path=getwd()))) %>%
          st_as_sf(crs = st_crs('EPSG:6365'))
mexico <- gadm(country='MX',level=1,resolution=2,path=getwd()) %>%
          st_as_sf(crs = st_crs('EPSG:6365'))
ev_pts <- st_as_sf(sismos_mc_hi, coords=c('Longitud','Latitud'), 
                   crs=st_crs('EPSG:6365'))

######### ######### ######### ######### ######### ######### ######### #########




# DECLUSTERING: Removal of aftershock from catalog



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
