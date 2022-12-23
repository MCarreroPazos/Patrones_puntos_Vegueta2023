# Script para reproducir los análisis del trabajo:
## Análisis de patrones espaciales de puntos para el estudio de tendencias 
## localización en distribuciones de yacimientos arqueológicos
# Autor: Miguel Carrero Pazos
# Email: miguel.carrero.pazos@gmail.com

# El siguiente código permite reproducir los análisis de la figura 3

# Definir el espacio de trabajo
setwd("~/Desktop/Patrones_puntos_Vegueta2023/datos") #MacOSX

# Cargar los paquetes de trabajo
library("raster")
library("rgdal")
library("readr")
library("spatstat")
library("maptools")

# Cargar los yacimientos, área de estudio y variable de distancia a líneas de cuenca hidrográfica
tumulos <- read.table(file ="datos_figura_3/sites.csv", header=TRUE, sep=";")
dist_cuencas <- raster("datos_figura_3/wsheddists.tif")
area_estudio <- readOGR(dsn ="datos_figura_3/studyarea.shp", layer="studyarea")

## Datos procedentes de Carrero-Pazos et al. 2019, JAS 
## https://github.com/MCarreroPazos/MontePenide)

# Preparar los datos para la función rhohat (Spatstat)
## Transformar el polígono de área en ventana de análisis
area <- as(area_estudio,"owin")
## Crear el patrón espacial de puntos (tumulos)
tumulos_sppp <- ppp(x=tumulos$x, y=tumulos$y, window=area)
## Convertir la variable raster a imagen (requisito para rhohat)
dist_cuencas_im <- as.im(as(dist_cuencas,"SpatialGridDataFrame"))

# Calcular la función rhohat de los túmulos y la covariable
rh <- rhohat(tumulos_sppp, dist_cuencas_im)
# Calcular la intensidad de túmulos basándose en la estimación de ρ (distancia a cuencas hidrográficas)
prh <- predict(rh)

# Crear la figura 3
png(file="figuras/Figura 3.png",
    width=13, height=6, units="in",res=300)
par(mar=c(5,5,5,2))
par(mfrow=c(1,3))
plot(dist_cuencas, main = "A", cex.main = 2, font.main = 2, axes = FALSE)
plot(tumulos_sppp, add = TRUE, cex = .8, pch = 20, col="black")
legend("topleft", legend="Distancia a ejes de\ncuenca hidrográfica y\ntúmulos megalíticos", cex = 1.5, bty='n', text.font = 2)
plot(rh, main="B", cex.main = 2, xlab="Distancia en metros", ylab="", ylim=c(0,2.2e-06), xlim=c(0,6000), legend=FALSE, cex.axis = 2, font.main = 2)
legend("topleft", legend="rhohat", cex = 1.5, bty='n', text.font = 2)
plot(prh*10000, main="")
title(main = "C", cex.main = 2, cex.axis = 2, font.main = 2)
plot(tumulos_sppp, add = TRUE, cex = .8, pch = 21, col="white", bg="black")
legend("topleft", legend="Intensidad de túmulos basada\n en la estimación de la covariable", cex = 1.5, bty='n', text.font = 2)
par(mfrow=c(1,1))
dev.off()
