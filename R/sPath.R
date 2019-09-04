#' Shortest Path on a gridded domain
#'
#' This function returns the shortest-path from \emph{droplet} location to
#' source of dispesal surface.
#' @param surface Grid (matrix) of arrival times or \emph{fastmaRching} object.
#' @param droplet Location of droplet
#' @param spatial.res (Optional) Spatial resolution of the grid in metre.
#' See example below.
#' @export
#' @examples
gridSpath <- function(surface, droplet, spatial.res) {
  if (class(surface)=='fastmaRching') {
    surface <- surface$cost.distance
    spatial.res <- surface$spatial.res
  }

  # prep cost surface
  dd <- surface
  dd <- rbind(rep(NA,NCOL(dd)),dd,rep(NA,NCOL(dd)))
  dd <- cbind(rep(NA,NROW(dd)),dd,rep(NA,NROW(dd)))

  # init droplet
  dist <- spatial.res*rbind(c(sqrt(2*spatial.res^2),spatial.res,sqrt(2*spatial.res^2)),c(spatial.res,NA,spatial.res),c(sqrt(2*spatial.res^2),spatial.res,sqrt(2*spatial.res^2)))
  droplet <- droplet + c(1,1)
  sPath <- c(NA,NA)

  sink <- which(dd==min(dd, na.rm=T), arr.ind = T)

  while(sum((droplet-sink)^2) > 0) {
    sPath <- rbind(sPath,droplet)

    neigh <- dd[(droplet[1]-1):(droplet[1]+1),(droplet[2]-1):(droplet[2]+1)]
    aux <- neigh*dist
    ind <- which( aux==min(aux, na.rm=T), arr.ind = T)[1,] ## TODO one path was chosen
    droplet <- ind-2+droplet
  }
  sPath <- rbind(sPath,droplet)

  sPath <- sPath[-1,]
  sPath <- sPath - 1
  rownames(sPath) <- c()

  # image(fm$cost.distance)
  # points(sPath[,1]/10, sPath[,2]/10)
  return(sPath)
}



#' Shortest Path on a spatial domain
#'
#' This function returns the shortest-path from \emph{droplet} location to
#' source of dispesal surface.
#' @param surface A \code{\link[raster]{raster}} object with arrival times
#' or \emph{fastmaRching} object.
#' @param droplet A \code{\link[sp]{SpatialPointsDataFrame}} object with
#' droplet location.
#' @export
#' @examples
spSpath <- function(surface, droplet) {
  if (class(surface)=='fastmaRching') { surface <- surface$cost.distance }

  surface.mat <- raster::as.matrix(surface)
  droplet.rp <- sp::spTransform(droplet, raster::crs(surface))

  # Check if seeds are inside domain ----------------------------------------
  test <- raster::extract(surface, droplet.rp); length(which(test==0 | is.na(test)))>0
  if (length(which(test==0 | is.na(test)))>0) { stop('dropletlet(s) not inside valid domain. Please check and rerun.') }

  droplet.mat <- raster::as.matrix(raster::rasterize(droplet.rp, surface, field=1, background=NA))
  aux <- t(cbind(droplet.mat[which(!is.na(droplet.mat))], which(!is.na(droplet.mat), arr.ind=TRUE))); droplet.mat <- aux[-1,]

  sPath <- gridSpath(surface.mat, droplet.mat)

  mm <- matrix(NA,nrow=NROW(surface.mat), ncol=NCOL(surface.mat))

  for (i in 1:NROW(sPath)) {
    mm[sPath[i,1], sPath[i,2]] <- i
  }
  rr <- raster(mm, template=surface)
  pp <- rasterToPoints(rr)
  pp <- pp[sort(pp[,3], index.return=T)$ix,]
  pp <- pp[,-3]

  ll <- Line(pp)
  S1 <- Lines(list(ll), ID = "line 1")
  Sl <- SpatialLines(list(S1), proj4string = crs(surface))

  return(Sl)
}
