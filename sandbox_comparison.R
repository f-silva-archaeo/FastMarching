
## Comparison of movecost, fastmarching and fasthiking functions
library(raster)


# movecost ----------------------------------------------------------------
library(movecost)
volc <- raster::raster(system.file("external/maungawhau.grd", package="gdistance"))
# volc[20,20:40] <- 1000
values(volc) <- 100
data(volc.loc)
data(destin.loc)
system.time(res.mc <- movecost(dtm=volc, origin=volc.loc, destin=destin.loc, breaks=0.05))



# fasthiking ------------------------------------------------------------
#devtools::install_github('f-silva-archaeo/fastmaRching')
library(fastmaRching)
seeds <- volc.loc
seeds@data <- data.frame(Id=0, incept=0, off.path=F)
system.time(res.fh <- spFastHike(volc, seeds))



# fastmarching ------------------------------------------------------------
slope <- terrain(volc, unit='degrees')
dd <- hiking.speed(slope)
seeds@data <- data.frame(Id=0, incept=0, speed=1)
res.fm <- spFastMarch(dd, seeds)



# Comparison --------------------------------------------------------------
par(mfrow=c(1,3), mar=c(1,1,1,3))
plot(res.mc$accumulated.cost.raster, axes=F); contour(res.mc$accumulated.cost.raster, add=TRUE)
plot(res.fh$arrival.time, axes=F); contour(res.fh$arrival.time, add=TRUE)
plot(res.fm$arrival.time, axes=F); contour(res.fm$arrival.time, add=TRUE)


### TODO
# Implement different processes and inception times competing for space as in Mod Fast Marching
