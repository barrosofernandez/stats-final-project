# ######################################### #
# k-means clustering of Mexican earthquakes #
# ######################################### #

# Made by: Roberto Barroso Fernández
# Instituto de Geociencias UNAM
#
# Proyecto final: Herramientas Estadísticas de las Geociencias

# WARNING WARNING WARNING WARNING WARNING

# Only for plotting

# ######################################### #
#            Setting up the code            #
# ######################################### #

# Reemplazar path por el preferido
setwd("/Users/roberto/Documents/stats-final-project")
# setwd("/Users/robby/Documents/stats-final-project")

library(geodata)
library(ggplot2)
library(sf)
library(dplyr)
library(scales)
library(cluster)  # silhouette
library(ggpubr)
library(dbscan)

plims <- function(x, p = 0.1) {
    # limits of plot axis
    # x - numeric
    # p - percentage (from 0 to 1) of padding
    maxval <- max(x)
    minval <- min(x)
    span   <- maxval - minval
    lim_hi <- maxval + 0.5 * p * span
    lim_lo <- minval - 0.5 * p * span
    return(c(lim_lo, lim_hi))
}

magfreqs <- function(m, mbin = 0.1) {
	# Generates data.frame of number of events against magnitude
	# Input:
	# m - list of magnitudes
	# mbin - binning width (difference)
	# Output:
	# df$magnitude -
	# df$amount -
	mags = seq(round(min(m), digits = 1), round(max(m), digits = 1), mbin)
	amt = 1:length(mags)
	for (i in 1:length(mags)){
		cutmag = mags[i]
		amt[i] = length(m[m>cutmag])
	}
	df <- data.frame(magnitude = mags, amount = amt)
	return(df)
}


sismos <- read.csv("catalogo_sismos.csv", skip = 4)[-(304906:304912), ] %>%
    filter(Magnitud != "no calculable")
sismos$utc <- paste(sismos$Fecha.UTC,sismos$Hora.UTC) %>%
    strptime("%Y-%m-%d %H:%M:%S") %>%
    as.POSIXct(tz = "UTC")
sismos$Magnitud <- as.numeric(sismos$Magnitud)

ord <- order(sismos$Magnitud, decreasing = FALSE)
sismos <- sismos[ord, ]

mfqs <- magfreqs(sismos$Magnitud)

bckgnd_countries <- c("US", "GT", "BZ", "SV", "HN", "NI", "CR", "PA")
bckgnd <- do.call("rbind", lapply(bckgnd_countries,
                  function(x) {gadm(country = x, level = 0,
				  resolution = 2, path = getwd())})) %>%
    st_as_sf(crs = st_crs("EPSG:6365"))
mexico <- gadm(country = "MX", level = 1, resolution = 2, path = getwd()) %>%
    st_as_sf(crs = st_crs("EPSG:6365"))
ev_pts <- st_as_sf(sismos, coords = c("Longitud", "Latitud"),
                   crs = st_crs("EPSG:6365"))

evmap1 <- ggplot() +
          geom_sf(data=bckgnd,fill='lightgrey',color='grey') +
          geom_sf(data=mexico,fill='white',color='grey') +
          geom_sf(data=ev_pts,alpha=0.8,aes(color=Profundidad,size=Magnitud)) +
          scale_size_binned("Magnitud", breaks = c(2, 4, 6, 8),
          labels = c("2", "4", "6", "8"), range = c(0.1, 10),
          transform = "exp") + 
          scale_color_viridis_c(option = "H") + 
          coord_sf(xlim=plims(sismos$Longitud,p=-0.1),
                   ylim=plims(sismos$Latitud,p=-0.1)) +
		  labs(title = "Catálogo original",
		       subtitle = paste("n=", length(sismos$Longitud), sep = "")) +
		  theme_classic()

magplot1 <- ggplot(mfqs,aes(x=magnitude,y=amount)) +
    geom_line() +
    geom_point() +
    scale_y_continuous(trans="log10",
                       breaks=trans_breaks("log10",function(x) 10^x),
                       labels=trans_format("log10",math_format(10^.x))) +
    geom_vline(xintercept = 3.1, alpha = 0.5, linewidth = 2, color = "red") +
    theme_classic() +
    xlab("Magnitud") +
    ylab("No. de eventos 1988-2024") +
    labs(title = "Catálogo original",
		       subtitle = paste("n=", length(sismos$Longitud), sep = ""))

ggexport(evmap1, filename = "ppt-eventos.png", 
         width = 3250, height = 2750, res = 500, pointsize = 14)

ggexport(magplot1, filename = "ppt-magplot1.png",
         width = 3250, height = 2750, res = 500, pointsize = 14)