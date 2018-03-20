#' Gridded Modified Fast Marching Method
#' @noRd
ModFastMarching <- function(domain, seeds, spatial.res=1) {

# Initializations ---------------------------------------------------------
  temp.res <- 1000    # temporal resolution
  Narrow_free <- 100000; Narrow_id <- 0; Narrow_list <- array(0,c(4,Narrow_free))   # Narrow Band arrays
  ne <- rbind(c(-1,0),c(1,0),c(0,-1),c(0,1))    # neighbour pixels
  clock <- 0    # external clock for seed inception


# Seed cleanup ------------------------------------------------------------
  dim(seeds) <- c(NROW(seeds), NCOL(seeds))
  seeds[3,] <- seeds[3,]/temp.res
  V <- seeds[4,]*(temp.res/spatial.res)
  incept <- matrix(c(1:(ncol(seeds)+1),c(seeds[3,],Inf)), nrow=2, byrow=TRUE) # incept times for seeds


# Intialization of Fast Marching grids ------------------------------------
  Map <- domain
  gridsize <- dim(Map)
  T <- array(0,c(gridsize[1],gridsize[2]))
  process <- array(0,c(gridsize[1],gridsize[2]))
  Tm2 <- array(0,c(1,4))
  Frozen <- Map==0


# Fast Marching Run -------------------------------------------------------
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
            Tt <- seeds[3,seed.id] + 1/(V[process[i,j]]+.Machine$double.eps)

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

        # calculates the distance using x-y and cross directions only
        Coeff <- c(0, 0, -1/((V[process[x,y]]*Map[x,y])^2));
        for (t in 1:2) { if(Order[t] > 0) { Coeff <- switch(Order[t], Coeff+c(1, -2*Tm[t], Tm[t]^2), Coeff+c(1, -2*Tm2[t], Tm2[t]^2)*9/4) }}
        Tt <- polyroot(rev(Coeff)); Tt <- max(Re(Tt))
        Coeff <- c(0, 0, -1/((V[process[x,y]]*Map[x,y])^2));
        for (t in 3:4) { if(Order[t] > 0) { Coeff <- switch(Order[t], Coeff+0.5*c(1, -2*Tm[t], Tm[t]^2), Coeff+0.5*c(1, -2*Tm2[t], Tm2[t]^2)*9/4) }}
        Tt2 <- polyroot(rev(Coeff))

        # Check for upwind condition
        if(length(Tt2) > 0) { Tt2 <- max(Re(Tt2)); Tt <- min(Tt,Tt2) }
        DirectNeigbInSol <- Tm[is.finite(Tm)]
        if(length(which(DirectNeigbInSol>=Tt)) > 0) { Tt <- min(DirectNeigbInSol) + 1/((V[process[x,y]]*Map[x,y])) }

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
  }


# Output ------------------------------------------------------------------
  T <- T * temp.res
  cdist <- T
  for (i in 1:NROW(seeds[3,])) {
    aux <- which(process==i)
    cdist[aux] <- spatial.res*V[i]*(cdist[aux] - seeds[3,i]*temp.res)
  }

  T[process==0] <- NA
  cdist[process==0] <- NA
  process[process==0] <- NA

  out <- list(domain = domain, seeds = seeds, spatial.res = spatial.res, arrival.time = T, process = process, cost.distance = cdist)
  class(out) <- 'fastmaRching'

  return(out)
}




#' Runs the grid version of the Modified Fast Marching Method
#'
#' This function runs the Modified Fast Marching Method of Silva and Steele
#' (2012,2014) on a gridded domain.
#' @param domain Grid (matrix) of chosen dimension with diffusivity values
#'  for every grid cell. Values above 1 will boost diffusivity, below 1 will
#'  inhibit it. Values of 0 should mark cells that block diffusion.
#' @param seeds A (4 x n) array containing the x-coordinate, y-coordinate,
#' incept time and rate-of-spread for each of the n seeds.
#' @param spatial.res (Optional) Spatial resolution of the grid, necessary only
#' to correct the rate-of-spread unit. See example below. Defaults to 1.
#' @references Sethian, J.A. (1996), A fast marching level set method for
#' monotonically advancing fronts, \emph{Proc. Natl. Acad. Sci.} 93 (4),
#' 1591-1595.
#' @references Silva, F. and Steele, J. (2012), Modeling Boundaries Between
#' Converging Fronts in Prehistory, \emph{Advances in Complex Systems},
#'  15(1-2), 1150005, <doi:10.1142/S0219525911003293>
#' @references Silva, F. and Steele, J. (2014), New methods for reconstructing
#' geographical effects on dispersal rates and routes from large-scale
#' radiocarbon databases, \emph{Journal of Archaeological Science} 52,
#' 609-620, <doi:10.1016/j.jas.2014.04.021>
#' @export
#' @examples
#' # Single process
#' grid <- matrix(1,10,10)
#' seed <- c(5,5,0,1)
#' fm <- gridFastMarch(grid, seed)
#' image(fm$arrival.time)
#'
#' # Two processes with same incept time
#' seeds <- cbind(c(7,7,0,1),c(2,2,0,1))
#' fm2 <- gridFastMarch(grid, seeds)
#' par(mfrow=c(1,3))
#' image(fm2$process, main='process')
#' image(fm2$arrival.time, main='arrival time')
#' image(fm2$cost.distance, main='distance')
#'
#' # Same as before but changing spatial.res parameter
#' fm3 <- gridFastMarch(grid, seeds, spatial.res=10)
#'
#' # Same as before but with a barrier in middle
#' grid[5,2:9] <- 0
#' fm4 <- gridFastMarch(grid, seeds, spatial.res=10)
#' par(mfrow=c(1,3))
#' image(fm4$process, main='process')
#' image(fm4$arrival.time, main='arrival time')
#' image(fm4$cost.distance, main='distance')
#'
#' # Same as before but with different incept times and speeds
#' seeds <- cbind(c(7,7,0,1),c(2,2,1,0.5))
#' fm5 <- gridFastMarch(grid, seeds, spatial.res=10)
#' par(mfrow=c(1,3))
#' image(fm5$process, main='process')
#' image(fm5$arrival.time, main='arrival time')
#' image(fm5$cost.distance, main='distance')
gridFastMarch <- function(domain, seeds, spatial.res=1) {
  compiler::enableJIT(3)
  gFM <- compiler::cmpfun(ModFastMarching)
  return(gFM(domain, seeds, spatial.res))
}



#' Runs the spatial version of the Modified Fast Marching Method
#'
#' This function runs the Modified Fast Marching Method of Silva and Steele
#' (2012,2014) from \emph{sp} and \emph{raster} objects and outputs results
#' in the same formats, makin it more convenient for (geo)spatial analyses
#' and simulation.
#' @param domain A \code{\link[raster]{raster}} object of chosen dimension
#' and resolution with diffusivity values for every cell. Values above 1 will
#'  boost diffusivity, below 1 will inhibit it. Values of 0 should mark cells
#'  that block diffusion.
#' @param seeds A \code{\link[sp]{SpatialPointsDataFrame}} object containing
#' the incept time and rate-of-spread for each of the n seeds in its data.frame,
#'  in columns named exactly \emph{incept} (for incept time) and \emph{speed} (
#'  for rate-of-spread).
#' This object will be automatically transformed to the projection of \emph{domain}.
#' @param spatial.res (Optional) Spatial resolution of the raster, necessary only
#' to correct the rate-of-spread unit. Defaults to that of the raster used for domain.
#' @references Sethian, J.A. (1996), A fast marching level set method for
#' monotonically advancing fronts, \emph{Proc. Natl. Acad. Sci.} 93 (4),
#' 1591-1595, doi:
#' @references Silva, F. and Steele, J. (2012), Modeling Boundaries Between
#' Converging Fronts in Prehistory, \emph{Advances in Complex Systems},
#'  15(1-2), 1150005, doi: 10.1142/S0219525911003293
#' @references Silva, F. and Steele, J. (2014), New methods for reconstructing
#' geographical effects on dispersal rates and routes from large-scale
#' radiocarbon databases, \emph{Journal of Archaeological Science} 52,
#' 609-620, doi: 10.1016/j.jas.2014.04.021
#' @export
#' @examples
#' library(raster); library(sp); library(rgdal)
#' domain <- raster(system.file("external/test.grd", package="raster")) # sample raster
#' domain <- domain > 0 # flattening elevation data
#' coords <- cbind(c(179000,181200), c(330000, 333000)) # coordinates for seeds
#' seed.df <- data.frame(incept=c(0,10), speed=c(.1,.1)) # incept time and speed for each seed
#' seeds <- SpatialPointsDataFrame(coords, seed.df, proj4string=crs(domain))
#'
#' fm <- spFastMarch(domain, seeds)
#' par(mfrow=c(1,3))
#' plot(fm$process, main='process')
#' plot(fm$arrival.time, main='arrival time')
#' plot(fm$cost.distance, main='distance')
spFastMarch <- function(domain, seeds, spatial.res) {
  # Convert Raster to Matrix ------------------------------------------------
  domain.grid <- raster::as.matrix(domain)
  domain.grid[is.na(domain.grid)] <- 0
  if (missing(spatial.res)) { spatial.res <- mean(raster::res(domain))/1000 } # in km


  # Convert SpatialPointsDataFrame to seeds array ---------------------------
  seeds.rp <- sp::spTransform(seeds, raster::crs(domain))
  seeds.matrix <- raster::as.matrix(raster::rasterize(seeds.rp, domain, background=NA)$ID)
  aux <- t(cbind(seeds.matrix[which(!is.na(seeds.matrix))], which(!is.na(seeds.matrix), arr.ind=TRUE)))
  aux <- rbind(aux, seeds.rp@data$incept, seeds.rp@data$speed); seeds.grid <- aux[-1,]


  # Run Fast Marching Method ------------------------------------------------
  fm <- gridFastMarch(domain.grid, seeds.grid, spatial.res)


  # Output ------------------------------------------------------------------
  TT <- raster::raster(fm$arrival.time, template=domain)
  pp <- raster::raster(fm$process, template=domain)
  cc <- raster::raster(fm$cost.distance, template=domain)

  fm$seeds <- seeds
  fm$domain <- domain
  fm$spatial.res <- spatial.res
  fm$arrival.time <- TT
  fm$process <- pp
  fm$cost.distance <- cc

  return(fm)
}
