# ######################################### #
# k-means clustering of Mexican earthquakes #
# ######################################### #

# Made by: Roberto Barroso Fernández
# Instituto de Geociencias UNAM
#
# Proyecto final: Herramientas Estadísticas de las Geociencias

# WARNING WARNING WARNING WARNING WARNING

# This is an abridged version of the code. You may only run this code
# once you have the following files: "background_seismicity.csv" and
# "cluster_tests_per_k.csv". Otherwise, it just won't work, and you'd
# have to deal with the not-so-"goodly" programmed full code, 
# "barroso-stats-proyecto-final.R"

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

######### ######### ######### ######### ######### ######### ######### #########

# ######################################### #
#            Defining functions             #
# ######################################### #

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

rct_to_sph <- function(X, Y, Z,r=6371) {
	# function to extract geographical coordinates from Cartesian
	# inputs of this function are outputs of sph_to_rct
	# and viceversa
	rho <- sqrt(X^2 + Y^2 + Z^2)
	theta <- atan(Y / X)
	phi <- acos(Z / rho)

	lat <- 90 - (phi * 180 / pi)
	lon <- theta * 180 / pi - 180  # ONLY FOR MEXICO !!!!!
	depth <- r - rho

	df <- data.frame(lat=lat,lon=lon,depth=depth)
	return(df)
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
	return(delta_val)
}

b_value <- function(mags,amount) {
	M <- head(mags, -1)
	F <- head(amount, -1)
	F_log <- log10(F)
	model <- lm(F_log~M)
	smmry <- summary(model)
	aval <- smmry$coefficients[1, 1]
	asd <- smmry$coefficients[1, 2]
	bval <- -smmry$coefficients[2, 1]
	bsd <- smmry$coefficients[2, 2]
	results <- list(a.value = aval, a.sd = asd, b.value = bval, b.sd = bsd)
	return(results)
}

######### ######### ######### ######### ######### ######### ######### #########

# Loading the catalog
sismos <- read.csv("catalogo_sismos.csv", skip = 4)[-(304906:304912), ] %>%
    filter(Magnitud != "no calculable")
sismos$utc <- paste(sismos$Fecha.UTC,sismos$Hora.UTC) %>%
    strptime("%Y-%m-%d %H:%M:%S") %>%
    as.POSIXct(tz = "UTC")
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
    time[i] <- timediff(date1, date2)
}
sismos$time <- time

# ######################################### #
#             k-means clustering            #
# ######################################### #

# Actually computing our clusters of interest

clust_tests <- read.csv("cluster_tests_per_k.csv")
sismos_bckgd <- read.csv("background_seismicity.csv")

cutoff_depth <- 40

bool_depth <- sismos_bckgd$Profundidad <= cutoff_depth

shallow <- sismos_bckgd[bool_depth, c("Longitud", "Latitud", "Magnitud",
                                      "Profundidad","X", "Y", "Z")]
deep <- sismos_bckgd[!bool_depth, c("Longitud", "Latitud", "Magnitud",
                                    "Profundidad","X", "Y", "Z")]

# I calculated stuff using the base-10 logarithm, but I should've
# used the natural logarithm instead. Luckily, due to the nature
# of logarithms, this is a quite easy fix

clust_tests$gap <- clust_tests$gap * log(10)
clust_tests$sk <- clust_tests$sk * log(10)

# Shallow only
elbow_plot <- ggplot(clust_tests, aes(x = k, y = elbow)) +
    geom_vline(xintercept = 10, linewidth = 2, alpha = 0.2, color = 'blue') +
    geom_line(linewidth = 0.3) +
    geom_point(size = 0.5) +
    ggtitle("Prueba de codo") +
    xlab("N. de grupos") +
    ylab("TWCSS") +
	theme_classic() +
	ylim(0,2e+09)
elbow_plot

gap_plot <- ggplot(clust_tests, aes(x = k, y = gap)) +
    geom_vline(xintercept = 10, linewidth = 2, alpha = 0.2, color = 'blue') +
    geom_line(linewidth = 0.3) +
    geom_point(size = 0.5) +
    geom_errorbar(aes(ymin = gap - sk, ymax = gap + sk), 
                  width = .5, ) +
    ggtitle('Prueba de hueco') +
    xlab('N. de grupos') +
    ylab('Estadístico gap') +
	theme_classic()

sil_plot <- ggplot(clust_tests, aes(x = k, y = silhouette)) +
    geom_vline(xintercept = 10, linewidth = 2, alpha = 0.2, color = 'blue') +
    geom_line(linewidth = 0.3) +
    geom_point(size = 0.5) +
    ggtitle('Prueba de silueta') +
    xlab('N. de grupos') +
    ylab('Anchura de silueta promedio') +
	theme_classic()

ggarrange(elbow_plot, gap_plot, sil_plot, ncol = 3, nrow = 1,
          labels = "AUTO", font.label = list(size = 12, face = "bold",
                                             color = "red")) %>%
    ggexport(filename = 'elbow_gap_sil.png', width = 3250, height = 1250,
             res = 500, pointsize = 10)

# El k más pequeño para el cual gap(k) >= gap(k+1) - sk(k+1) es k=9
# Sin embargo, k=10 también cumple con la condición
# Para siluetas, el valor máximo es para k=2, donde los otros picos se dan
# para k=6 y k=10
# Por ende, se escoge k=10

clust_10 <- kmeans(shallow[, c("X", "Y", "Z")], 10, iter.max = 20, nstart = 25)

# Statistics of each cluster because it is necessary
clust_x <- clust_10$centers[, 1]
clust_y <- clust_10$centers[, 2]
clust_z <- clust_10$centers[, 3]
center_coords <- rct_to_sph(clust_x, clust_y, clust_z)



######### ######### ######### ######### ######### ######### ######### #########

# ######################################### #
#             Obtaining b-values            #
# ######################################### #

n_events <- clust_10$size
centroid_lat <- round(center_coords$lat, 3)
centroid_lon <- round(center_coords$lon, 3)
mag_min  <- 1:10
mag_max  <- 1:10
mag_mean <- 1:10
a_val <- 1:10
a_sd <- 1:10
b_val <- 1:10
b_sd <- 1:10
for (ii in 1:10){
	clust_ref <- shallow[clust_10$cluster == ii, ]
	mfqs <- magfreqs(clust_ref$Magnitud)
	m_clust <-mfqs$magnitude
	f_clust <- mfqs$amount / 36  # normalize to year span of data
	
	mag_min[ii] <- min(clust_ref$Magnitud)
	mag_max[ii] <- max(clust_ref$Magnitud)
	mag_mean[ii] <- round(mean(clust_ref$Magnitud), 3)
	# avoiding the non-linear parts of the curve 
	magbol <- m_clust < 5.2 & m_clust > 3.3 
	gutricht <- b_value(m_clust[magbol],f_clust[magbol]) 
	a_val[ii] <- round(gutricht$a.value, 3)
	a_sd[ii]  <- round(gutricht$a.sd, 3)
	b_val[ii] <- round(gutricht$b.value, 3)
	b_sd[ii]  <- round(gutricht$b.sd, 3)
}

clust_base_df <- data.frame(N = 1:10, Eventos = n_events, Latitud = centroid_lat,
                            Longitud = centroid_lon, MagMin = mag_min,
							MagMax = mag_max, MagMean = mag_mean, 
							a_val = a_val, a_sd = a_sd,
							b_val = b_val, b_sd = b_sd)

write.csv(clust_base_df, file = "clust_base_df.csv")

######### ######### ######### ######### ######### ######### ######### #########

bckgnd_countries <- c("US", "GT", "BZ", "SV", "HN", "NI", "CR", "PA")
bckgnd <- do.call("rbind", lapply(bckgnd_countries,
                  function(x) {gadm(country = x, level = 0,
				  resolution = 2, path = getwd())})) %>%
    st_as_sf(crs = st_crs("EPSG:6365"))
mexico <- gadm(country = "MX", level = 1, resolution = 2, path = getwd()) %>%
    st_as_sf(crs = st_crs("EPSG:6365"))
ev_pts <- st_as_sf(sismos, coords = c("Longitud", "Latitud"),
                   crs = st_crs("EPSG:6365"))
ev_pts_bckgd <- st_as_sf(sismos_bckgd, coords = c("Longitud", "Latitud"),
                   crs = st_crs("EPSG:6365"))
ev_pts_shallow <- st_as_sf(shallow, coords = c("Longitud", "Latitud"),
                   crs = st_crs("EPSG:6365"))
ev_pts_deep <- st_as_sf(deep, coords = c("Longitud", "Latitud"),
                   crs = st_crs("EPSG:6365"))

######### ######### ######### ######### ######### ######### ######### #########


# No ejecutar esto a menos que esté bien arreglado
evmap1 <- ggplot() +
          geom_sf(data=bckgnd,fill='lightgrey',color='grey') +
          geom_sf(data=mexico,fill='white',color='grey') +
          geom_sf(data=ev_pts,alpha=0.5,pch=16,size=0.05) +
          coord_sf(xlim=plims(sismos_bckgd$Longitud,p=-0.1),
                   ylim=plims(sismos_bckgd$Latitud,p=-0.1)) +
		  labs(title = "Catálogo original",
		       subtitle = paste("n=", length(sismos$X), sep = "")) +
		  theme_classic()

evmap2 <- ggplot() +
          geom_sf(data=bckgnd,fill='lightgrey',color='grey') +
          geom_sf(data=mexico,fill='white',color='grey') +
          geom_sf(data=ev_pts_bckgd,alpha=0.5,pch=16,size=0.05) +
          coord_sf(xlim=plims(sismos_bckgd$Longitud,p=-0.1),
                   ylim=plims(sismos_bckgd$Latitud,p=-0.1)) +
		  labs(title = "Sismicidad de fondo",
		       subtitle = paste("n=", length(sismos_bckgd$X), sep = "")) +
		  theme_classic()

evmap3 <- ggplot() +
          geom_sf(data=bckgnd,fill='lightgrey',color='grey') +
          geom_sf(data=mexico,fill='white',color='grey') +
          geom_sf(data=ev_pts_shallow,alpha=0.5,pch=16,size=0.05) +
          coord_sf(xlim=plims(sismos_bckgd$Longitud,p=-0.1),
                   ylim=plims(sismos_bckgd$Latitud,p=-0.1)) +
		  labs(title = "Sismicidad superficial",
		       subtitle = paste("n=", length(shallow$X), sep = "")) +
		  theme_classic()

evmap4 <- ggplot() +
          geom_sf(data=bckgnd,fill='lightgrey',color='grey') +
          geom_sf(data=mexico,fill='white',color='grey') +
          geom_sf(data=ev_pts_deep,alpha=0.5,pch=16,size=0.05) +
          coord_sf(xlim=plims(sismos_bckgd$Longitud,p=-0.1),
                   ylim=plims(sismos_bckgd$Latitud,p=-0.1)) +
		  labs(title = "Sismicidad profunda",
		       subtitle = paste("n=", length(deep$X), sep = "")) +
		  theme_classic()

ggarrange(evmap1, evmap2, evmap3, evmap4, ncol = 2, nrow = 2,
          labels = "AUTO", font.label = list(size = 12, face = "bold",
                                             color = "red")) %>%
    ggexport(filename = 'evmaps_full.png', width = 3250, height = 2750,
             res = 500, pointsize = 10)

# magmap1 <- ggplot() +
#            geom_sf(data=bckgnd,fill='lightgrey',color='grey') +
#            geom_sf(data=mexico,fill='white',color='grey') +
#            geom_sf(data=ev_pts,alpha=0.2,
#                    aes(color=Profundidad)) +
#            coord_sf(xlim=plims(sismos_bckgd$Longitud,p=-0.1),
#                     ylim=plims(sismos_bckgd$Latitud,p=-0.1))
			
# magmap1

clust_labels <- c("SubGO", "SubCG", "SubJL", "SubCM", "SubOC",
                  "SubGR", "SMOr", "BCS", "FVTM", "BCN")

clustmap1 <- ggplot() +
           geom_sf(data=bckgnd,fill='lightgrey',color='grey') +
           geom_sf(data=mexico,fill='white',color='grey') +
           geom_sf(data=ev_pts_shallow,pch=16,size=1,
                   aes(color=as.factor(clust_10$cluster))) +
		   scale_color_brewer(palette = "Paired") +
           coord_sf(xlim=plims(sismos_bckgd$Longitud,p=-0.1),
                    ylim=plims(sismos_bckgd$Latitud,p=-0.1)) +
		   labs(title = "Agrupación de eventos", color = "Grupo") +
		   annotate("text", x = clust_base_df$Longitud, 
		            y = clust_base_df$Latitud, label = clust_labels, 
					size=3) +
		   theme_classic() +
		   xlab(NULL) + ylab(NULL)

ggexport(clustmap1, filename = 'clustmap1.png', width = 3250, height = 2750,
             res = 500, pointsize = 10)


shallow$cluster <- clust_10$cluster

subcm <- shallow[clust_10$cluster == 4 & shallow$Magnitud < 5 & shallow$Magnitud > 3.3, ]
magscm <- head(magfreqs(subcm$Magnitud), -1)

clustmags <- ggplot(magscm, aes(magnitude, amount)) +
    geom_point() +
	scale_y_log10() +
	geom_smooth(method = lm) +
	theme_classic()
clustmags
			
# clustmap1
# No ejecutar esto a menos que esté bien arreglado
