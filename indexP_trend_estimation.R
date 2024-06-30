################################################################################
#
# We compute the trend and uncertainty estimates (95% CI and p-value) for 
# the Index P time series of each pixel for a given country or territory. 
# 
# Some of this code is adapted from the MannKendallTrends R package from 
# Collaud Coen et al. (https://doi.org/10.5194/amt-13-6945-2020)
#
################################################################################

################################################################################
# SET UP

# load the required packages
library(tidyverse)
library(raster)
library(trend)
library(terra)
library(parallel)
library(MannKendallTrends)

# define the method used to estimate trend: seasonal Mann-Kendall test with pre-whitening
trend_method <- "seasonal_mannkendall_PW"

# define the number of cores
ncores <- 4

# set the country/territory 
region <- "Togo"

# load the raster for the region of interest
FILE_PATH <- #INCLUDE FILE PATH
indexP_raster <- raster::brick()

indexP_raster <- raster::brick(file.path("indexP", "V2_indexP_denv_aegypti_maps", region, paste0(region, "_indexP_monthly_mean_rasters.tif")))


################################################################################
# HELPER FUNCTIONS

# compute_mk_stat takes as an argument a vector of index P values over time and 
# calculates various statistics based on the Mann-Kendall test and Sen's slope. 
compute_mk_stat <- function(x) {
  if (all(x==0)) return(list(P=NA, S=0, vari=0, Z=NA, slope=0, UCL=0, LCL=0))
  result <- list()
  # Mann-Kendall test
  mk_results <- mk.test(x)
  # Sen's slope
  sen_results <- sens.slope(x)
  
  # results
  result$P <- mk_results$p.value
  result$S <- mk_results$estimates[1]
  result$vari <- mk_results$estimates[2]
  result$Z <- mk_results$statistic
  result$slope <- sen_results$estimates
  result$UCL <- sen_results$conf.int[2]
  result$LCL <- sen_results$conf.int[1]
  return(result)
}

# prewrite.AR assesses for significant autocorrelation and then performs
# prewhitening. 
prewhite.AR <- function (data, alpha.ak = 95) {
  durb.levinson <- function(x, p = NULL) {
    fit <- function(acf, p) {
      ref <- numeric(p)
      g <- -acf[2]/acf[1]
      a <- g
      v <- Re((1 - g * Conj(g)) * acf[1])
      ref[1] <- g
      for (t in 2:p) {
        g <- -(acf[t + 1] + a %*% acf[seq(t, 2, by = -1)])/v
        a <- c((a + as.vector(g) * Conj(a[seq(t - 1, 1, -1)])), g)
        v <- v * (1 - Re(g * Conj(g)))
        ref[t] <- g
      }
      a <- c(1, a)
      return(list(a = a, v = v, ref = ref))
    }
    if ((!is.null(p) && (p != as.integer(p))) || (p < 2)) 
      stop("p must be integer >= 2.")
    if (is.numeric(x) && is.null(dim(x))) {
      lx <- length(x)
      if (is.null(p) || p >= lx) 
        p <- lx - 1
      r <- fit(x, p)
    }
    else {
      if (is.numeric(x) && !(is.null(dim(x)))) {
        lx <- dim(x)
        if (is.null(p) || p >= lx[1]) 
          p <- lx[1] - 1
        zr <- apply(x, 2, function(y) fit(y, p))
        zr <- matrix(unlist(zr), nrow = lx[2], byrow = TRUE)
        a <- zr[, 1:(p + 1), drop = FALSE]
        v <- zr[, p + 2]
        ref <- t(zr[, -(1:(p + 2)), drop = FALSE])
        r <- list(a = a, v = v, ref = ref)
      }
      else {
        stop("x must be a numeric vector or matrix.")
      }
    }
    return(r)
  }
  
  data[which(is.infinite(data) == TRUE)] <- NA
  data[which(is.nan(data) == TRUE)] <- NA
  p <- 5
  nblag <- 10
  if (nblag > length(data)/2) {
    nblag <- floor(length(data)/2)
  }
  if (p >= nblag) {
    p <- nblag - 1
  }
  x <- as.numeric(acf(data, lag.max=nblag, plot=FALSE)$acf)
  x1 <- x/length(na.omit(data))
  lev <- durb.levinson(x1, p)
  K <- lev$ref
  ak.coef <- -K
  uconf <- qnorm(1 - (1 - alpha.ak/100)/2, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)/sqrt(sum(complete.cases(data)))
  ak.lag <- x[2]
  if (abs(ak.coef[1]) - uconf < 0) {
    data.prewhite <- data
    ak.ss <- 0
    warning("no s.s. autocorrelation")
  } else {
    y <- rep(NA, length(data))
    y[-1] <- x[2] * data[-length(data)]
    data.prewhite <- data - y
    ak.ss <- alpha.ak
  }
  return(list(ak.lag = ak.lag, ak.ss = ak.ss, data.prewhite = data.prewhite))
}

# prewhite_ts calculates prewhitened time series using various algorithms as described in Collaud Coen et al.
# we use the trend-free prewhitening (TFPW.Y) and ignore the output for the other methods. 
# this code was sourced and adapted from the MannKendallTrends R package. 
prewhite_ts <- function (data.ts, column, resolution, alpha.ak = 95) {
  data <- data.ts[, column]
  data[is.infinite(abs(data))] <- NA
  epoch.time <- as.numeric(data.ts$Time)/3600/24
  t.time <- do.call(cbind.data.frame, lapply(c("%Y", "%m", "%d"), function(x) as.numeric(format(data.ts$Time, x))))
  names(t.time) <- c("year", "month", "day")
  nb.loop <- 0
  dataPW <- data.frame(time = data.ts[, 1])
  c <- data.frame(matrix(NA, nrow = 1, ncol = 0))
  out <- prewhite.AR(data = data, alpha.ak = alpha.ak)
  c$PW <- out$ak.lag
  dataARremoved <- as.numeric(out$data.prewhite)
  c$ss <- out$ak.ss
  rm(out)
  if (length(na.omit(dataARremoved)) > 0 & c$ss == alpha.ak & c$PW >= 0.05) {
    dataPW$PW <- dataARremoved
    dataPW$PW.cor <- dataARremoved/(1 - c$PW)
    ties <- Nb.tie(data = dataARremoved/(1 - c$PW), resolution = resolution)
    n <- S.test(data = dataARremoved/(1 - c$PW), t.time = t.time)$n
    vari <- Kendall.var(data = dataARremoved/(1 - c$PW), 
                        t = ties, n = n)
    b0.PW <- sen.slope(data = dataARremoved/(1 - c$PW), epoch.time = epoch.time, 
                       vari = vari, alpha.cl = 90)$slope
    t <- Nb.tie(data = data, resolution = resolution)
    n <- S.test(data = data, t.time = t.time)$n
    vari <- Kendall.var(data = data, t = ties, n = n)
    b0.or <- sen.slope(data = data, epoch.time = epoch.time, 
                       vari = vari)$slope
    dataDetrend.PW <- data - b0.PW * c(epoch.time - epoch.time[1])
    dataDetrend.or <- data - b0.or * c(epoch.time - epoch.time[1])
    out <- nanprewhite.AR(data = dataDetrend.or, alpha.ak = alpha.ak)
    c$VCTFPW <- out$ak.lag
    dataARremoved.or <- as.numeric(out$data.prewhite)
    c$ssVC <- out$ak.ss
    rm(out)
    out <- nanprewhite.AR(data = dataDetrend.PW, alpha.ak = alpha.ak)
    ak.PW <- out$ak.lag
    dataARremoved.PW <- as.numeric(out$data.prewhite)
    ss.PW <- out$ak.ss
    rm(out)
    if (sum(na.omit(dataARremoved.or) > 0)) {
      dataPW$TFPW.Y <- dataARremoved.or + b0.or * (epoch.time - epoch.time[1])
    } else {
      dataPW$TFPW.Y <- data
    }
    if (abs(ak.PW) >= 0.05 & ss.PW == 95) {
      c1 <- c$PW
      dataARremoved.PW <- c(data[1], (data[-1] - ak.PW * data[-length(data)])/(1 - ak.PW))
      t <- Nb.tie(data = dataARremoved.PW, resolution = resolution)
      n <- S.test(data = dataARremoved.PW, t.time = t.time)$n
      vari <- Kendall.var(dataARremoved.PW, t = ties, n = n)
      b1.PW <- sen.slope(data = dataARremoved.PW, epoch.time = epoch.time, vari = vari)$slope
      while (abs(ak.PW - c1) > 1e-04 & abs(b1.PW - b0.PW) > 1e-04) {
        if (ak.PW >= 0.05 & ss.PW == 95) {
          nb.loop <- nb.loop + 1
          dataDetrend.PW <- data - b1.PW * (epoch.time - epoch.time[1])
          c1 <- ak.PW
          b0.PW <- b1.PW
          out <- nanprewhite.AR(data = dataDetrend.PW, alpha.ak = alpha.ak)
          ak.PW <- out$ak.lag
          dataARremoved.2PW <- as.numeric(out$data.prewhite)
          ss.PW <- out$ak.ss
          rm(out)
          if (ak.PW > 0 & ss.PW == 95) {
            dataARremoved.2PW <- c(data[1], c(data[-1] - ak.PW * data[-length(data)])/(1 - ak.PW))
            t <- Nb.tie(data = dataARremoved.2PW, resolution = resolution)
            n <- S.test(data = dataARremoved.2PW, t.time = t.time)$n
            vari <- Kendall.var(data = dataARremoved.2PW, 
                                t = t, n = n)
            b1.PW <- sen.slope(data = dataARremoved.2PW, epoch.time = epoch.time, vari = vari)$slope
            dataARremoved.PW <- dataARremoved.2PW
            if (nb.loop > 10) 
              break
          }
          else {
            break
          }
        }
      }
    } else {
      b1.PW <- b0.PW
      ak.PW <- c$PW
    }
    if (sum(na.omit(dataARremoved.PW)) > 0) {
      dataPW$TFPW.WS <- dataARremoved.PW
      c$TFPW.WS <- ak.PW
    } else {
      dataPW$TFPW.WS <- data
    }
    var.data <- var(data, na.rm = TRUE)
    var.data.TFPW <- var(dataARremoved.or, na.rm = TRUE)
    dataARremoved.var <- dataARremoved.or * (var.data/var.data.TFPW)
    if (c$VCTFPW >= 0) {
      bVC <- b0.or/sqrt((1 + c$VCTFPW)/(1 - c$VCTFPW))
    } else if (c$VCTFPW < 0) {
      bVC <- b0.or
    }
    dataPW$VCTFPW <- dataARremoved.var + bVC * (epoch.time - epoch.time[1])
  } else {
    dataPW$PW <- data
    dataPW$PW.cor <- data
    dataPW$TFPW.Y <- data
    dataPW$TFPW.WS <- data
    dataPW$VCTFPW <- data
  }
  return(dataPW)
}

# MK_test_PW performs the seasonal Mann-Kendall test with monthly index P
# time series from 1979 to 2022. Sen's slope is used to estimate the linear 
# rate of change for each pixel. Prewhitening is performed to adjust for 
# temporal autocorrelation in the data. 
MK_test_PW <- function (data, seasons = NULL) {
  # organize data
  names(data)[1] <- "Time"

  dataPW <- prewhite_ts(data.ts = data, column = 2, resolution = 0.001, alpha.ak = 95)
  data.season <- split(x = dataPW, f = seasons)
  result <- data.frame(matrix(NA, nrow = length(data.season) + 1, ncol = 5))
  colnames(result) <- c("slope", "UCL", "LCL", "Z", "P")
  S.TFPW.Y <- 0
  vari.TFPW.Y <- 0
  for (m in 1:length(data.season)) {
    data.mois <- data.season[[m]]
    if (length(na.omit(data.mois$PW)) > 1) {
      out.TFPW.Y <- compute_mk_stat(x=na.omit(data.mois$TFPW.Y))
      result[m, c("slope", "UCL", "LCL", "Z", "P")] <- out.TFPW.Y[c("slope", "UCL", "LCL", "Z", "P")]
      S.TFPW.Y <- S.TFPW.Y+out.TFPW.Y$S
      vari.TFPW.Y <- vari.TFPW.Y + out.TFPW.Y$vari
    } else {
      result[m, c("slope", "UCL", "LCL", "Z", "P")] <- NA
    }
  }
  Ztot.TFPW.Y <- STD.normale.var(S.TFPW.Y, vari.TFPW.Y)
  Ptot.TFPW.Y <- 2 * (1 - pnorm(abs(Ztot.TFPW.Y), 0, 1))
  result[m + 1, c("Z", "P")] <- c(Ztot.TFPW.Y, Ptot.TFPW.Y)
  result$slope[m + 1] <- median(result$slope[1:m], na.rm = TRUE)
  result$UCL[m + 1] <- median(result$UCL[1:m], na.rm = TRUE)
  result$LCL[m + 1] <- median(result$LCL[1:m], na.rm = TRUE)
  return(result)
}

################################################################################
# ANALYSIS

# calc_slope is a wrapper function that calculates the seasonal MK test for 
# the prewhitened time series. 
calc_slope <- function(x) {
  # check if all NA
  if (all(is.na(x))) return(rep(NA, 5*13))
  
  # compute the trend using the Mann-Kendall statistic and Sen's slope
  x <- as.numeric(x)
  if (all(x==0)) return(c(rep(0, 13), rep(NA, 4*13)))
  dates <- seq(from=as.Date("1979/01/01"), as.Date("2022/12/01"), by="month")
  month <- rep(1:12, length=length(x))
  pixel_data <- data.frame(x=dates, y=x)
  result <- MK_test_PW(data=pixel_data, seasons=month)
  
  return(c(result$slope, result$UCL, result$LCL, result$Z, result$P))
}

# perform the analysis 
indexP_raster <- rast(indexP_raster)

cl <- parallel::makeCluster(ncores)
clusterEvalQ(cl, library(trend))
clusterEvalQ(cl, library(MannKendallTrends))
clusterExport(cl, list("MK_test_PW", "compute_mk_stat", "prewhite_ts", "prewhite.AR"))
trend_output <- app(indexP_raster, fun=calc_slope, cores=cl)
stopCluster(cl)
names(trend_output) <- c(paste0("slope", 1:13), paste0("UCL", 1:13), paste0("LCL", 1:13),  paste0("Z", 1:13), paste0("P", 1:13))
