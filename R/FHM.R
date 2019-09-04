#' Hiking Speed Function(s)
#'
#' This function outputs the hiking speed in km/hr dependant on slope for a
#' variety of functions (see below).
#' @param slope Slope values in either degrees, radians or unitless. Can be array,
#' matrix or raster object.
#' @param units (Optional) Either 'degree', 'radian' or 'nounit'. Default is 'degree'.
#' @param fun (Optional) Hiking function. Either 'Tobler', 'Naismith', 'Langmuir1' or
#' 'Langmuir2'. Default is 'Tobler'.
#' @param v0 (Optional) Speed at flat surface. If given it will normalise hiking
#' function outputs to the value given.
#' @param off.path (Optional) Boolean for Tobler hiking function. Default is FALSE.
#' @param horseback (Optional) Boolean for Tobler hiking function.Default is FALSE.
#' @param cutoff (Optional) If a value is given any slopes above the cutoff will output speed zero.
#' @references Sethian, J.A. (1996), A fast marching level set method for
#' monotonically advancing fronts, \emph{Proc. Natl. Acad. Sci.} 93 (4),
#' 1591-1595.
#' @references Tobler, Waldo (1993), Three presentations on geographical analysis and modeling:
#' Non-isotropic geographic modeling speculations on the geometry of geography global spatial
#' analysis, \emph{Technical report. National center for geographic information and analysis} 93 (1).
#' @references Naismith, W. W. (1892), Excursions: Cruach Ardran, Stobinian, and Ben More,
#' \emph{Scottish Mountaineering Club Journal} 2 (3): 136.
#' @references Langmuir, Eric (1984), \emph{ Mountaincraft and Leadership: Official Handbook of the
#' Mountain Leader Training Boards of Great Britain and Northern Ireland}. Edinburgh Scotland: Britain
#' & Scottish Sports Council. (Langmuir1)
#' @references Langmuir, Eric (2013). \emph{ Mountaincraft and Leadership: A Handbook for Mountaineers
#' and Hillwalking Leaders in the British Isles (Fourth ed.)}. Mountain Training England; Mountain Training Scotland.
#' @export
#' @examples
#' hiking.speed(20, fun='Tobler')
#' hiking.speed(20, fun='Tobler', horseback=TRUE)
#' hiking.speed(-10, fun='Langmuir1')
#' plot(seq(-20,20,1), hiking.speed(seq(-20,20,1), fun='Langmuir1'))
#'
#' dem <- raster::raster(system.file("external/test.grd", package="raster"))
#' slope <- raster::terrain(dem)
#' speed1 <- hiking.speed(slope)
#' speed2 <- hiking.speed(slope, fun='Naismith')

hiking.speed <- function(slope, units='degree', fun='Tobler', v0=NULL, off.path=F, horseback=F, cutoff) {
  if (class(slope)=='RasterLayer') {
    #require(raster)
    r <- slope
    slope <- raster::as.matrix(r)
    rast <- TRUE
  } else { rast <- FALSE }

  slp <- switch( units,
                 degree = tan(slope/180*pi),
                 radian = tan(slope),
                 nounit = slope )
  if (is.null(slp)) { stop('Units not recognised. Please use either degree, radian or nounit (for unitless).') }


  if (fun=='Tobler') {
    if (!is.null(v0)) { factor <- v0/exp(-3.5 * abs(0.05)) } else { factor <- 6 }
    if (off.path) { factor <- factor * 3/5 }
    if (horseback) { factor <- factor * 5/4 }
    speed <- factor * exp(-3.5 * abs(slp + 0.05))
  } else

    if (fun=='ToblerMod') {
      if (!is.null(v0)) { factor <- v0/exp(-5.3 * abs(0.03)) } else { factor <- 4.8 }
      speed <- factor * exp(-5.3 * abs(0.7*slp + 0.03))
    } else

    if (fun=='Naismith') {
      if (missing(v0)) { v0 <- 4.83 }

      dx <- 1
      dy <- slp
      dh <-  sqrt(dy^2 + dx^2)
      dt <- dx / v0 + dy / 0.61
      speed <- dh / dt

      speed[slope < 0] <- NA
    } else

      if (fun=='Langmuir1') {
        if (is.null(v0)) { v0 <- 4.83 }

        dx <- matrix(1,nrow=NROW(slp), ncol=NCOL(slp))
        dy <- slp
        dh <-  sqrt(dy^2 + dx^2)
        dt <- dx / v0
        ind <- which(slp > 0); if (length(ind)>0) { dt[ind] <- dt[ind] + abs(dy[ind]) / 0.61 }
        ind <- which(slp <= -0.08748866 & slp > -0.2125566); if (length(ind)>0) { dt[ind] <- dt[ind] - 10/60 * abs(dy[ind]) / 0.3 }
        ind <- which(slp <= -0.2125566); if (length(ind)>0) { dt[ind] <- dt[ind] + 10/60 * abs(dy[ind]) / 0.3 }
        speed <- dh / dt
      }  else

        if (fun=='Langmuir2') {
          if (is.null(v0)) { v0 <- 4 }

          dx <- matrix(1,nrow=NROW(slp), ncol=NCOL(slp))
          dy <- slp
          dh <-  sqrt(dy^2 + dx^2)
          dt <- dx / v0
          ind <- which(slp > 0); if (length(ind)>0) { dt[ind] <- dt[ind] + abs(dy[ind]) / 0.45 }
          ind <- which(slp <= -0.08748866 & slp > -0.2125566); if (length(ind)>0) { dt[ind] <- dt[ind] - 10/60 * abs(dy[ind]) / 0.3 }
          ind <- which(slp <= -0.2125566); if (length(ind)>0) { dt[ind] <- dt[ind] + 10/60 * abs(dy[ind]) / 0.3 }
          speed <- dh / dt
        } else { stop('Function not recognised. Please use Tobler, Naismith, Langmuir1 or Langmuir2 only.')}

  if (!missing(cutoff)) { speed[slope > cutoff] <- 0 }

  if (rast) {
    raster::values(r) <- speed
    speed <- r
  }
  return(speed)
}



#' Gridded Modified Fast Hiking Method
#' @noRd
ModFastHiking <- function(domain, seeds, spatial.res, fun) {

  # Seed cleanup ------------------------------------------------------------
  if ("v0" %in% names(seeds)) { v0 <- seeds$v0 } else { v0 <- rep(NULL,NROW(seeds)) }
  if ("off.path" %in% names(seeds)) { off.path <- seeds$off.path } else { off.path <- rep(F,NROW(seeds)) }
  if ("horseback" %in% names(seeds)) { horseback <- seeds$horseback } else { horseback <- rep(F,NROW(seeds)) }
  seeds <- t(seeds[,1:3])
  seeds[3,] <- seeds[3,]
  incept <- matrix(c(1:(ncol(seeds)+1),c(seeds[3,],Inf)), nrow=2, byrow=TRUE) # incept times for seeds


  # Intialization of Fast Hiking grids ------------------------------------
  Map <- matrix(NA, nrow=NROW(domain)+4, ncol=NCOL(domain)+4)
  Map[3:(NROW(domain)+2), 3:(NCOL(domain)+2)] <- domain
  Arrival <- matrix(Inf, nrow=NROW(Map), ncol=NCOL(Map))
  Frozen <- Map==0
  TT <- Arrival
  dist <- rbind( c(sqrt(8), sqrt(5), 2, sqrt(5), sqrt(8)),
                 c(sqrt(5), sqrt(2), 1, sqrt(2), sqrt(5)),
                 c(2,       1,       0, 1      , 2),
                 c(sqrt(5), sqrt(2), 1, sqrt(2), sqrt(5)),
                 c(sqrt(8), sqrt(5), 2, sqrt(5), sqrt(8)) )*spatial.res

  dist <- rbind( c(Inf      , sqrt(5), Inf , sqrt(5), Inf),
                 c(sqrt(5), sqrt(2), 1 , sqrt(2), sqrt(5)),
                 c(Inf      ,       1, 0 , 1      , Inf),
                 c(sqrt(5), sqrt(2), 1 , sqrt(2), sqrt(5)),
                 c(Inf      , sqrt(5), Inf , sqrt(5), Inf))*spatial.res

  # Incepts seeds to initiate processes at incept time ----------------------
  ind <- seeds[1:2] + c(2,2)
  Frozen[ind[1], ind[2]] <- 1

  # Calculate arrival time for all neighbours (knight's move)
  Patch <- Map[(ind[1]-2):(ind[1]+2), (ind[2]-2):(ind[2]+2)]
  slope <- (Patch - Patch[3,3])/dist
  speed <- hiking.speed(slope, 'nounit', fun, v0, off.path, horseback)*1000 # in metres per hour
  time <- dist/speed
  time[Frozen[(ind[1]-2):(ind[1]+2), (ind[2]-2):(ind[2]+2)]==1] <- NA
  Arrival[(ind[1]-2):(ind[1]+2), (ind[2]-2):(ind[2]+2)] <- time

  TT <- Arrival
  TT[sort(TT, index.return=T, na.last=T)$ix] <- seq(1,length(Arrival))
  TT[!is.finite(Arrival)] <- NA

  # Loop through unitl boundary has covered entire domain
  F1 <- length(which(Frozen == 1))
  pb <- utils::txtProgressBar(max=length(which(Frozen==0)), style=3)
  while(sum(Frozen==0, na.rm=T) > 0) {

    ## Pick fastest pixel and make it current
    ind <- which(TT==min(TT, na.rm=TRUE), arr.ind = TRUE)
    TT[ind] <- NA

    Patch <- Map[(ind[1]-2):(ind[1]+2), (ind[2]-2):(ind[2]+2)]
    slope <- (Patch - Patch[3,3])/dist
    speed <- hiking.speed(slope, 'nounit', fun, v0, off.path, horseback)*1000 # in metres per hour
    time <- dist/speed + Arrival[ind[1], ind[2]]

    Arrival2 <- matrix(NA, nrow=NROW(Map), ncol=NCOL(Map))
    time[Frozen[(ind[1]-2):(ind[1]+2), (ind[2]-2):(ind[2]+2)]==1] <- NA
    Arrival2[(ind[1]-2):(ind[1]+2), (ind[2]-2):(ind[2]+2)] <- time

    ## Any smaller than Arrival, then this is better route
    Arrival[which(Arrival2 < Arrival)] <- Arrival2[which(Arrival2 < Arrival)]

    Frozen[ind[1], ind[2]] <- 1
    TT <- Arrival
    TT[!is.finite(Arrival)] <- NA
    TT[Frozen==1] <- NA
    TT[sort(TT, index.return=T, na.last=T)$ix] <- seq(1,length(Arrival))
    TT[!is.finite(Arrival)] <- NA
    TT[Frozen==1] <- NA
    ## TODO clean up above

    utils::setTxtProgressBar(pb, length(which(Frozen==1))-F1)
  }

  # cleanup
  Arrival[!is.finite(Arrival)] <- 0
  Arrival <- Arrival[3:(NROW(Arrival)-2), 3:(NCOL(Arrival)-2)]

  # Output ------------------------------------------------------------------
  out <- list(domain = domain, seeds = seeds, spatial.res = spatial.res, arrival.time = Arrival)
  class(out) <- 'fastmaRching'

  return(out)
}




#' Fast Hiking Method on a gridded domain
#'
#' This function runs the Fast Hiking Method on a gridded domain.
#' Output arrival time is in hours.
#' @param domain Grid (matrix) of chosen dimension with diffusivity values
#'  for every grid cell. Values above 1 will boost diffusivity, below 1 will
#'  inhibit it. Values of 0 should mark cells that block diffusion.
#' @param seeds A (4 x n) array containing the x-coordinate, y-coordinate,
#' incept time and  \emph{hiking.speed} parameters for each of the n seeds.
#' @param spatial.res (Optional) Spatial resolution of the grid in metre.
#' See example below.
#' @param fun (Optional) Hiking function to use, See \link{hiking.speed} for details.
#' @export
#' @examples
#' # Single process
#' grid <- matrix(1,10,10)
#' seed <- data.frame(row=5, col=5, incept=0)
#' fm <- gridFastHike(grid, seed)
#' image(fm$arrival.time)
#'
#' # Two processes with same incept time
#' seeds <- data.frame(row=c(7,2), col=c(7,2), incept=c(0,0))
#' fm2 <- gridFastHike(grid, seeds)
#' par(mfrow=c(1,3))
#' image(fm2$process, main='process')
#' image(fm2$arrival.time, main='arrival time')
#' image(fm2$cost.distance, main='distance')
#'
#' # Same as before but changing spatial.res parameter
#' fm3 <- gridFastHike(grid, seeds, spatial.res=10)
#'
#' # Same as before but with a barrier in middle
#' grid[5,2:9] <- 0
#' fm4 <- gridFastHike(grid, seeds, spatial.res=10)
#' par(mfrow=c(1,3))
#' image(fm4$process, main='process')
#' image(fm4$arrival.time, main='arrival time')
#' image(fm4$cost.distance, main='distance')
#'
#' # Same as before but with different incept times and speeds
#' seeds <- data.frame(row=c(7,2), col=c(7,2), incept=c(0,1), v0=c(2,5))
#' fm5 <- gridFastHike(grid, seeds, spatial.res=10)
#' par(mfrow=c(1,3))
#' image(fm5$process, main='process')
#' image(fm5$arrival.time, main='arrival time')
#' image(fm5$cost.distance, main='distance')
gridFastHike <- function(domain, seeds, spatial.res=1, fun='Tobler') {
  compiler::enableJIT(3)
  gFH <- compiler::cmpfun(ModFastHiking)
  return(gFH(domain, seeds, spatial.res, fun))
}



#' Fast Hiking Method on a spatial domain
#'
#' This function runs the Fast Hiking Method from \emph{sp} and
#' \emph{raster} objects and outputs results in the same formats, making
#' it more convenient for (geo)spatial analyses and simulation. Output
#' arrival time is in hours.
#' @param dem A \code{\link[raster]{raster}} object of chosen dimension
#' and resolution with elevation data.
#' @param seeds A \code{\link[sp]{SpatialPointsDataFrame}} object containing
#' parameter values for each of the n seeds in its data.frame,
#'  in colums named exactly \emph{incept} (for incept time), \emph{off.path},
#'  \emph{horseback} (see Tobler's function). This object will be automatically
#'   transformed to the projection of \emph{domain}.
#' @param spatial.res (Optional) Spatial resolution of the raster in metres.
#' Defaults to that of the raster used.
#' @param fun (Optional) Hiking function to use, See \link{hiking.speed} for details.
#' @export
#' @examples
#' library(raster); library(sp); library(rgdal)
#' domain <- raster(system.file("external/test.grd", package="raster")) # sample raster
#' coords <- cbind(c(179000,181200), c(330000, 333000)) # coordinates for seeds
#' seed.df <- data.frame(incept=c(0,1)) # parameters for each seed
#' seeds <- SpatialPointsDataFrame(coords, seed.df, proj4string=crs(domain))
#'
#' fm <- spFastHike(domain, seeds)
#  plot(domain)
#  contour(fm$arrival.time, levels=seq(1,5,1), add=T)
spFastHike <- function(dem, seeds, spatial.res, fun='Tobler') {
  # Convert Raster to Matrix ------------------------------------------------
  domain.grid <- raster::as.matrix(dem)
  domain.grid[is.na(domain.grid)] <- 0
  if (missing(spatial.res)) { spatial.res <- mean(raster::res(dem)) } # in meters


  # Convert SpatialPointsDataFrame to seeds array ---------------------------
  seeds.rp <- sp::spTransform(seeds, raster::crs(dem))
  seeds.matrix <- raster::as.matrix(raster::rasterize(seeds.rp, dem, background=NA)$ID)
  aux <- t(cbind(seeds.matrix[which(!is.na(seeds.matrix))], which(!is.na(seeds.matrix), arr.ind=TRUE))); aux <- aux[-1,]
  dim(aux) <- c(NROW(aux), NCOL(aux))
  seeds.grid <- data.frame(row=aux[1,], col=aux[2,], incept=seeds.rp@data$incept)
  if ("v0" %in% names(seeds.rp@data)) { seeds.grid$v0 <- seeds.rp@data$v0 }
  if ("off.path" %in% names(seeds.rp@data)) { seeds.grid$off.path <- seeds.rp@data$off.path }
  if ("horseback" %in% names(seeds.rp@data)) { seeds.grid$horseback <- seeds.rp@data$horseback }


  # Check if seeds are inside domain ----------------------------------------
  test <- raster::extract(dem, seeds.rp); length(which(test==0 | is.na(test)))>0
  if (length(which(test==0 | is.na(test)))>0) { stop('Seed(s) not inside valid domain. Please check and rerun.') }


  # Run Fast Hiking Method ------------------------------------------------
  fm <- gridFastHike(domain.grid, seeds.grid, spatial.res, fun)


  # Output ------------------------------------------------------------------
  TT <- raster::raster(fm$arrival.time, template=dem)

  fm$seeds <- seeds
  fm$domain <- dem
  fm$spatial.res <- spatial.res
  fm$arrival.time <- TT
  return(fm)
}
