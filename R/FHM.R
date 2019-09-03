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
#' @export
ModFastHiking <- function(domain, seeds, spatial.res, fun) {

  # Initializations ---------------------------------------------------------
  temp.res <- 1000    # temporal resolution
  Narrow_free <- 100000; Narrow_id <- 0; Narrow_list <- array(0,c(4,Narrow_free))   # Narrow Band arrays
  ne <- rbind(c(-1,0),c(1,0),c(0,-1),c(0,1))    # neighbour pixels
  clock <- 0    # external clock for seed inception


  # Seed cleanup ------------------------------------------------------------
  if ("v0" %in% names(seeds)) { v0 <- seeds$v0 } else { v0 <- rep(NULL,NROW(seeds)) }
  if ("off.path" %in% names(seeds)) { off.path <- seeds$off.path } else { off.path <- rep(F,NROW(seeds)) }
  if ("horseback" %in% names(seeds)) { horseback <- seeds$horseback } else { horseback <- rep(F,NROW(seeds)) }
  seeds <- t(seeds[,1:3])
  seeds[3,] <- seeds[3,]/temp.res
  incept <- matrix(c(1:(ncol(seeds)+1),c(seeds[3,],Inf)), nrow=2, byrow=TRUE) # incept times for seeds


  # Intialization of Fast Hiking grids ------------------------------------
  Map <- domain
  gridsize <- dim(Map)
  T <- array(0,c(gridsize[1],gridsize[2]))
  process <- array(0,c(gridsize[1],gridsize[2]))
  Tm2 <- array(0,c(1,4))
  Frozen <- Map==0


  # Fast Hiking Run -------------------------------------------------------
  F1 <- length(which(Frozen==1))
  pb <- utils::txtProgressBar(max=length(which(Frozen==0)), style=3)
  while(Narrow_id > -1) {

    if(clock >= min(incept[2,])) {
      # Incepts seeds to initiate processes at incept time ----------------------

      # choose seeds to incept
      seed.id <- incept[1, which.min(incept[2,] - clock)]
      seed.points <- floor( seeds[1:2,seed.id] );  dim(seed.points) <- c(2,length(seed.id))

      # incept seeds
      x <- seed.points[1,]; y <- seed.points[2,]
      for (kk in 1:length(seed.id)) {
        Frozen[x[kk],y[kk]] <- 1
        T[x[kk],y[kk]] <- seeds[3,seed.id]
        process[x[kk],y[kk]] <- seed.id

        # add neighbours of seeds to Narrow list
        for (k in 1:4) {
          i <- x[kk] + ne[k,1]; j <- y[kk] + ne[k,2]

          if ((i > 0) && (j > 0) && (i <= dim(T)[1]) && (j <= dim(T)[2]) && (Frozen[i,j]==0)) {
            process[i,j] <- process[x[kk], y[kk]]
            slope <- (Map[i,j] - Map[x[kk],y[kk]]) / spatial.res
            speed <- hiking.speed(slope, units='nounit', fun=fun, v0=v0[process[i,j]], off.path=off.path[process[i,j]], horseback=horseback[process[i,j]])*temp.res
            Tt <- seeds[3,seed.id] + spatial.res/(speed*1000)

            if (T[i,j] > 0) {
              Narrow_list[1,T[i,j]] <- min(Re(Tt),Narrow_list[1,T[i,j]])
              Narrow_list[4,T[i,j]] <- process[x[kk],y[kk]]
            } else {
              Narrow_id <- Narrow_id + 1
              if (Narrow_id > Narrow_free) { Narrow_free <- Narrow_free + 100000; Narrow_list[1,Narrow_free] <- 0 }
              Narrow_list[,Narrow_id] <- rbind( Re(Tt), i, j, process[x[kk],y[kk]] )
              T[i,j] <- Narrow_id
            }
          }
        }
      }
      d <- dim(incept); d[2] <- d[2] - 1
      incept <- incept[,-1]; dim(incept) <- d
    }

    # get the closest pixel to the wavefront and make it the current pixel
    t <- min(Narrow_list[1,1:Narrow_id]); index <- which.min(c(Narrow_list[1,1:Narrow_id])); clock <- t
    x <- Narrow_list[2,index]; y <- Narrow_list[3,index]
    Frozen[x,y] <- 1
    T[x,y] <- Narrow_list[1,index]
    process[x,y] <- Narrow_list[4,index]

    # replace min value with the last value in the array
    if(index < Narrow_id) {
      Narrow_list[,index] <- Narrow_list[,Narrow_id]
      x2 <- Narrow_list[2,index]; y2 <-  Narrow_list[3,index]
      T[x2,y2] <- index
    }
    Narrow_id <- Narrow_id - 1

    # loop through neighbours of current pixel
    for (k in 1:4) {
      i <- x + ne[k,1]; j <- y + ne[k,2]

      # check if neighbour is not yet frozen
      if((i > 0) && (j > 0) && (i <= dim(T)[1]) && (j <= dim(T)[2]) && (Frozen[i,j]==0) && (process[i,j]==0)) {
        Tpatch <- matrix(Inf,5,5)
        for (nx in -2:2) {
          for (ny in -2:2) {
            i.n <- i + nx; j.n <- j + ny
            if((i.n > 0) && (j.n > 0) && (i.n <= dim(T)[1]) && (j.n <= dim(T)[2]) && (Frozen[i.n,j.n]==1) && (process[i.n,j.n]==process[x,y])) {
              Tpatch[nx+3,ny+3] <- T[i.n,j.n]
            }
          }
        }

        # stores the derivative order of derivative to use
        Order <- array(0,c(1,4))

        # First Order derivatives in x-y and cross directions
        Tm <- array(0,c(1,4));
        Tm[1] <- min(Tpatch[2,3],Tpatch[4,3]); if(is.finite(Tm[1])) { Order[1] <- 1 }
        Tm[2] <- min(Tpatch[3,2],Tpatch[3,4]); if(is.finite(Tm[2])) { Order[2] <- 1 }
        Tm[3] <- min(Tpatch[2,2],Tpatch[4,4]); if(is.finite(Tm[3])) { Order[3] <- 1 }
        Tm[4] <- min(Tpatch[2,4],Tpatch[4,2]); if(is.finite(Tm[4])) { Order[4] <- 1 }

        # Second Order derivatives
        Tm2 <- array(0,c(1,4))
        ch1 <- (Tpatch[1,3] < Tpatch[2,3]) && is.finite(Tpatch[2,3])
        ch2 <- (Tpatch[5,3] < Tpatch[4,3]) && is.finite(Tpatch[4,3])
        if(ch1) { Tm2[1] <- (4*Tpatch[2,3] - Tpatch[1,3])/3; Order[1] <- 2 }
        if(ch2) { Tm2[1] <- (4*Tpatch[4,3] - Tpatch[5,3])/3; Order[1] <- 2 }
        if(ch1&&ch2) { Tm2[1] <- min((4*Tpatch[2,3]-Tpatch[1,3])/3, (4*Tpatch[4,3] - Tpatch[5,3])/3); Order[1] <- 2 }

        ch1 <- (Tpatch[3,1] < Tpatch[3,2]) && is.finite(Tpatch[3,2])
        ch2 <- (Tpatch[3,5] < Tpatch[3,4]) && is.finite(Tpatch[3,4])
        if(ch1) { Tm2[2] <- (4*Tpatch[3,2] - Tpatch[3,1])/3; Order[2] <- 2 }
        if(ch2) { Tm2[2] <- (4*Tpatch[3,4] - Tpatch[3,5])/3; Order[2] <- 2 }
        if(ch1&&ch2) { Tm2[2] <- min((4*Tpatch[3,2] - Tpatch[3,1])/3, (4*Tpatch[3,4] - Tpatch[3,5])/3); Order[2] <- 2}

        ch1 <- (Tpatch[1,1] < Tpatch[2,2]) && is.finite(Tpatch[2,2])
        ch2 <- (Tpatch[5,5] < Tpatch[4,4]) && is.finite(Tpatch[4,4])
        if(ch1) { Tm2[3] <- (4*Tpatch[2,2] - Tpatch[1,1])/3; Order[3] <- 2 }
        if(ch2) { Tm2[3] <- (4*Tpatch[4,4] - Tpatch[5,5])/3; Order[3] <- 2 }
        if(ch1&&ch2) { Tm2[3] <- min((4*Tpatch[2,2] - Tpatch[1,1])/3, (4*Tpatch[4,4] - Tpatch[5,5])/3); Order[3] <- 2}

        ch1 <- (Tpatch[1,5] < Tpatch[2,4]) && is.finite(Tpatch[2,4])
        ch2 <- (Tpatch[5,1] < Tpatch[4,2]) && is.finite(Tpatch[4,2])
        if(ch1) { Tm2[4] <- (4*Tpatch[2,4] - Tpatch[1,5])/3; Order[4] <- 2 }
        if(ch2) { Tm2[4] <- (4*Tpatch[4,2] - Tpatch[5,1])/3; Order[4] <- 2 }
        if(ch1&&ch2) { Tm2[4] <- min((4*Tpatch[2,4] - Tpatch[1,5])/3, (4*Tpatch[4,2] - Tpatch[5,1])/3); Order[4] <- 2}

        # calculates hiking speed to all neighbours
        slope.xy <- (Map[i,j] - Map[x,y]) / spatial.res
        speed.xy <- hiking.speed(slope.xy, units='nounit', fun=fun, v0=v0[process[x,y]], off.path=off.path[process[x,y]], horseback=horseback[process[x,y]])*1000*temp.res
        slope.cross <- (Map[i,j] - Map[x,y]) / (sqrt(2*spatial.res^2))
        speed.cross <- hiking.speed(slope.cross, units='nounit', fun=fun, v0=v0[process[x,y]], off.path=off.path[process[x,y]], horseback=horseback[process[x,y]])*1000*temp.res

        # calculates the distance using x-y and cross directions only
        Coeff <- c(0, 0, -1/(speed.xy^2))
        for (t in 1:2) { if(Order[t] > 0) { Coeff <- switch(Order[t], Coeff + c(1, -2*Tm[t], Tm[t]^2), Coeff + c(1, -2*Tm2[t], Tm2[t]^2)*9/4) }}
        Tt <- polyroot(rev(Coeff)); Tt <- max(Re(Tt))
        Coeff <- c(0, 0, -1/(speed.cross^2))
        # Coeff <- c(0, 0, -1/(speed.xy^2))
        for (t in 3:4) { if(Order[t] > 0) { Coeff <- switch(Order[t], Coeff + 0.5*c(1, -2*Tm[t], Tm[t]^2), Coeff + 0.5*c(1, -2*Tm2[t], Tm2[t]^2)*9/4) }}
        Tt2 <- polyroot(rev(Coeff))

        # Check for upwind condition
        if(length(Tt2) > 0) { Tt2 <- max(Re(Tt2)); Tt <- min(Tt,Tt2) }
        DirectNeigbInSol <- Tm[is.finite(Tm)]
        if(length(which(DirectNeigbInSol>=Tt)) > 0) { Tt <- max(DirectNeigbInSol) + 1/(speed.xy) }

        # Updates Narrow Band array
        if(T[i,j] > 0) {
          mT <- min(c(Tt, Narrow_list[1,T[i,j]])); ind <- which.min(c(Tt, Narrow_list[1,T[i,j]]))
          Narrow_list[1,T[i,j]] <- Re(mT)
          if (ind==1) { Narrow_list[4,T[i,j]] <- process[x,y] }
        } else {
          Narrow_id <- Narrow_id + 1
          if(Narrow_id > Narrow_free) { Narrow_free <- Narrow_free + 100000; Narrow_list[1,Narrow_free]<- 0 }
          Narrow_list[,Narrow_id] <- Re(c(Tt,i,j,process[x,y]))
          T[i,j] <- Narrow_id
        }
      }
    }
    utils::setTxtProgressBar(pb, length(which(Frozen==1))-F1)
  }


  # Output ------------------------------------------------------------------
  T <- T * temp.res
  cdist <- T
  for (i in 1:NROW(seeds[3,])) {
    aux <- which(process==i)
    cdist[aux] <- spatial.res*(cdist[aux] - seeds[3,i]*temp.res)
  }

  T[process==0] <- NA
  cdist[process==0] <- NA
  process[process==0] <- NA

  out <- list(domain = domain, seeds = seeds, spatial.res = spatial.res, arrival.time = T, process = process, cost.distance = cdist)
  class(out) <- 'fastmaRching'

  return(out)
}




#' Modified Fast Hiking Method on a gridded domain
#'
#' This function runs the Modified Fast Hiking Method on a gridded domain.
#' Output arrival time is in hours.
#' @param domain Grid (matrix) of chosen dimension with diffusivity values
#'  for every grid cell. Values above 1 will boost diffusivity, below 1 will
#'  inhibit it. Values of 0 should mark cells that block diffusion.
#' @param seeds A (4 x n) array containing the x-coordinate, y-coordinate,
#' incept time and  \emph{hiking.speed} parameters for each of the n seeds.
#' @param spatial.res (Optional) Spatial resolution of the grid in metre.
#' See example below. Defaults to 1 metre.
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



#' Modified Fast Hiking Method on a spatial domain
#'
#' This function runs the Modified Fast Hiking Method from \emph{sp} and
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
#' Defaults to that of the raster used for DEM
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
  pp <- raster::raster(fm$process, template=dem)
  cc <- raster::raster(fm$cost.distance, template=dem)

  fm$seeds <- seeds
  fm$domain <- dem
  fm$spatial.res <- spatial.res
  fm$arrival.time <- TT
  fm$process <- pp
  fm$cost.distance <- cc

  return(fm)
}
