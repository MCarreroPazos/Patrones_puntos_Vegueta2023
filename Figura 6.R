# Script para reproducir los análisis del trabajo:
## Análisis de patrones espaciales de puntos para el estudio de tendencias 
## localización en distribuciones de yacimientos arqueológicos
# Autor: Miguel Carrero Pazos
# Email: miguel.carrero.pazos@gmail.com

# El siguiente código permite reproducir los análisis de la figura 6
# Definir el espacio de trabajo
setwd("~/Desktop/Patrones_puntos_Vegueta2023/datos") #MacOSX

# Cargar los paquetes de R
spatpack<-c("raster","spatstat","rgdal","maptools", "sp", "sf")
lapply(spatpack, require, character.only=TRUE)

# Parte 1. Estudio del patrón monovariante
# Cargar datos de base y definir el área de estudio
lcps <- raster("grid/densidad_lcps_grid_1000.tiff")
area_estudio <-readOGR(dsn="shp/area_estudio.shp", layer="area_estudio")
area_w <- as.owin((area_estudio))
sites <-readOGR(dsn="shp/megalitos.shp", layer="megalitos")

## Cortar el raster de tránsito con la máscara del área de estudio
lcps_crop <- crop(lcps, extent(area_estudio))
lcps.m <- mask(lcps_crop, area_estudio)

# Crear el patrón espacial de puntos
sppp <- ppp(x=sites$UMTX, y=sites$UMTY, window=area_w)

# Convertir el raster a imagen para rhohat
lcps_im <- as.im(as(lcps.m,"SpatialGridDataFrame"))

# Parte 1. Estudiar el patrón monovariante (rhohat)
transito.rh <- rhohat(sppp, lcps_im, confidence=0.95)
## Aplicar las pruebas Z1 y Z2 de Berman para la dependencia de la covariable seleccionada
Z1 <- berman.test(sppp, lcps_im)
Z2 <- berman.test(sppp, lcps_im, which = "Z2")
print(Z1)
print(Z2)

# Visualizar el gráfico rhohat
par(mfrow=c(1,1))
plot(transito.rh, main="", xlab="Densidad", ylab="", legend=FALSE, cex.axis=0.7)
legend("topleft", legend="Tránsito", cex=0.7, bty='n', text.font=2)

# Parte 2. Simulaciones poisson no homogéneas
# Cargar raster de tránsito y definir la nueva área de estudio (zonas altas de la sierra)
lcps <- raster("grid/densidad_lcps_grid_1000.tiff")
area_estudio <-readOGR(dsn="shp/area_red.shp", layer="area_red")
area_red <- as.owin((area_estudio))

## Cortar el raster de tránsito con la máscara del área de estudio
lcps_crop <- crop(lcps, extent(area_estudio))
lcps.m <- mask(lcps_crop, area_estudio)
lcps_im_red <- as.im(as(lcps.m,"SpatialGridDataFrame"))

# Generar un patrón de puntos Poisson no homogéneo influenciado por el tránsito
sppp_ipoiss <- rpoispp(lambda = predict(rhohat(sppp,
                                               lcps_im_red, 
                                               method="ratio",
                                               confidence=0.95, 
                                               eps=50)), nsim= 99)

# Ejecutar la simulación (función de correlación par)
sims <- 99
rango <- round((sims + 1) / 100 * 2.5, 0) # intervalo de confianza del 95%
sim_PCF_prim_orden <- envelope(sppp,
                               W = area_red,
                               pcfinhom,
                               simulate = sppp_ipoiss,
                               nsim = sims,
                               rank = rango)

# Guardar el resultado de la función
save(sim_PCF_prim_orden, file = "~/Desktop/Patrones_puntos/datos/Rdata/sim_PCF_prim_orden_99simulaciones.RData")

# Cargar el resultado de la función con 999 simulaciones
load(file="~/Desktop/Patrones_puntos/datos/Rdata/sim_PCF_prim_orden999simulaciones.RData")

# Crear la figura 6
png(file="figuras/Figura 6.png",
    width=11, height=6, units="in",res=300)
l <- layout(matrix(c(1, 2, 3, 
                     4, 5, 6),
                   nrow = 2,
                   ncol = 3,
                   byrow = TRUE))
layout.show(l)
par(mar=c(4,2,2,2))
plot(transito.rh, main="A. Rhohat tránsito", xlab="Densidad de tránsito", ylab="", legend=FALSE, cex.axis=0.7)
plot(Z1, main="B. Test de Berman", xlab="Densidad de tránsito", ylab="Probabilidad", legend=FALSE, cex.axis=0.7)
plot(lcps_im, main="C. Monumentos megalíticos", cex.main = 1.5)
plot(sppp, add=T, pch = 21, col = "white", bg="black")
plot(lcps_im, main="D. Simulación no. 23", cex.main = 1.5)
plot(sppp_ipoiss$`Simulation 23`, add=T, pch = 20, col = "white")
plot(lcps_im, main="E. Simulación no. 38", cex.main = 1.5)
plot(sppp_ipoiss$`Simulation 38`, add=T, pch = 20, col = "white")
par(mar=c(3,3,1,3))
plot(sim_PCF_prim_orden,ylim=c(0,50), xlim=c(0, 2000),legend=FALSE,xlab="",ylab="",bty="l",main="F. Función de correlación par")
title(ylab=expression(italic("g(r)")), line=2,cex.lab=1)
title(xlab=expression(italic("Distancia entre puntos (metros)")), line=2,cex.lab=1)
legend(500,35, legend=c("Patrón de puntos observado","Estadístico aleatorio esperado"),lty=c(1,2), col=c("black","red"),bty="n",cex=1)
legend(555,27,legend=c("Simulación de Monte Carlo"),pch=0,fill=c("grey"),col=NA,border=NA,bty="n",cex=1)
dev.off()
