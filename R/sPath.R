#' Shortest Path on a gridded domain
#'
#' TODO description
#' @param cdist
#' @param drop
#' @param spatial.res (Optional) Spatial resolution of the grid in metres.
#' @export
#' @examples
#'
gridSpath <- function(cdist, drop, spatial.res=1) {
  if (class(cdist)=='fastmaRching') {
    cdist <- cdist$cost.distance
    spatial.res <- cdist$spatial.res
  }

  # prep cost surface
  dd <- cdist
  dd <- rbind(rep(NA,NCOL(dd)),dd,rep(NA,NCOL(dd)))
  dd <- cbind(rep(NA,NROW(dd)),dd,rep(NA,NROW(dd)))

  # init droplet
  dist <- spatial.res*rbind(c(sqrt(2*spatial.res^2),spatial.res,sqrt(2*spatial.res^2)),c(spatial.res,NA,spatial.res),c(sqrt(2*spatial.res^2),spatial.res,sqrt(2*spatial.res^2)))
  drop <- drop + c(1,1)
  sPath <- c(NA,NA)

  sink <- which(dd==min(dd, na.rm=T), arr.ind = T)

  while(sum((drop-sink)^2) > 0) {
    sPath <- rbind(sPath,drop)

    neigh <- dd[(drop[1]-1):(drop[1]+1),(drop[2]-1):(drop[2]+1)]
    aux <- neigh*dist
    ind <- which( aux==min(aux, na.rm=T), arr.ind = T)[1,] ## TODO one path was chosen
    drop <- ind-2+drop
  }
  sPath <- rbind(sPath,drop)

  sPath <- sPath[-1,]
  sPath <- sPath - 1
  rownames(sPath) <- c()

  # image(fm$cost.distance)
  # points(sPath[,1]/10, sPath[,2]/10)
  return(sPath)
}


#' Shortest Path on a spatial domain
#'
#' TODO description
#' @param cdist
#' @param drop
#' @param spatial.res (Optional) Spatial resolution of the raster in metres.
#' @export
#' @examples
#'
spSpath <- function(cdist, drop) {
  if (class(cdist)=='fastmaRching') { cdist <- cdist$cost.distance }

  cdist.mat <- raster::as.matrix(cdist)
  drop.rp <- sp::spTransform(drop, raster::crs(cdist))

  # Check if seeds are inside domain ----------------------------------------
  test <- raster::extract(cdist, drop.rp); length(which(test==0 | is.na(test)))>0
  if (length(which(test==0 | is.na(test)))>0) { stop('Droplet(s) not inside valid domain. Please check and rerun.') }

  drop.mat <- raster::as.matrix(raster::rasterize(drop.rp, cdist, field=1, background=NA))
  aux <- t(cbind(drop.mat[which(!is.na(drop.mat))], which(!is.na(drop.mat), arr.ind=TRUE))); drop.mat <- aux[-1,]

  sPath <- gridSpath(cdist.mat, drop.mat)

  mm <- matrix(NA,nrow=NROW(cdist.mat), ncol=NCOL(cdist.mat))

  for (i in 1:NROW(sPath)) {
    mm[sPath[i,1], sPath[i,2]] <- i
  }
  rr <- raster(mm, template=cdist)
  pp <- rasterToPoints(rr)
  pp <- pp[sort(pp[,3], index.return=T)$ix,]
  pp <- pp[,-3]

  ll <- Line(pp)
  S1 <- Lines(list(ll), ID = "line 1")
  Sl <- SpatialLines(list(S1), proj4string = crs(cdist))

  return(Sl)
}
