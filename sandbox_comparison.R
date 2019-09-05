
## Comparison of movecost, fastmarching and fasthiking functions
library(raster)


# movecost ----------------------------------------------------------------
library(movecost)
volc <- raster::raster(system.file("external/maungawhau.grd", package="gdistance"))
values(volc) <- 100
#volc[20,20:40] <- 1000
# data(volc.loc)
# data(destin.loc)
# res.mc <- movecost(dtm=volc, origin=volc.loc, destin=destin.loc, breaks=0.05)
#


# fasthiking ------------------------------------------------------------
#devtools::install_github('f-silva-archaeo/fastmaRching')
library(fastmaRching)


# Voronoi diagram demo ----------------------------------------------------
seeds <- destin.loc
seeds@data <- data.frame(Id=1:length(seeds), incept=rep(0,length(seeds)), off.path=rep(F,length(seeds)))
res.fh <- spFastHike(volc, seeds)

par(mfrow=c(1,2), mar=c(1,1,1,3))
plot(res.fh$arrival.time, axes=F); contour(res.fh$arrival.time, add=TRUE)
plot(res.fh$process, axes=F); contour(res.fh$arrival.time, add=TRUE)



# First Arrival demo ------------------------------------------------------
seeds <- destin.loc[c(4,8),]
seeds@data <- data.frame(Id=1:2, incept=rep(0,2), off.path=rep(F,2), v0=c(2,NA))
res.fh <- spFastHike(volc, seeds)

par(mfrow=c(1,2), mar=c(1,1,1,3))
plot(res.fh$arrival.time, axes=F); contour(res.fh$arrival.time, add=TRUE)
plot(res.fh$process, axes=F); contour(res.fh$arrival.time, add=TRUE)



ttt <- spSpath(res.fh$arrival.time, destin.loc[1])   ## TODO make this work

