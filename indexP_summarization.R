################################################################################
#
# We calculate the average Index P for different time periods. 
#
################################################################################

################################################################################
# SET UP

# load the required packages
library(tidyverse)
library(raster)
library(terra)
library(parallel)
library(sf)

# define the start and end of the required time period
args <- commandArgs(trailingOnly=TRUE)
start_year <- as.integer(args[1]) # start year 
if (!(start_year %in% 1979:2022)) stop("Start year must be between 1979-2022")
end_year <- as.integer(args[2]) # end year
if (!(end_year %in% 1979:2022)) stop("Start year must be between 1979-2022")

# define layers of rasters to subset
earliest_year <- 1979
subset_layers <- seq(from=((start_year-earliest_year)*12)+1, length.out=(end_year-start_year+1)*12)

# set the required region
region <- "Togo"

# load the raster for the region of interest
FILE_PATH <- #INCLUDE FILE PATH FOR THE MEDIAN INDEX P
indexP_raster_p50 <- raster::brick()
FILE_PATH <- #INCLUDE FILE PATH FOR THE Q3 INDEX P
indexP_raster_p75 <- raster::brick()
FILE_PATH <- #INCLUDE FILE PATH FOR THE Q1 INDEX P
indexP_raster_p25 <- raster::brick()

################################################################################
# COMPUTE THE MEAN INDEXP FOR A GIVEN INTERVAL

# define the distribution of sample means as a normal distribution using the CLT
indexP_raster_p50 <- rast(indexP_raster_p50)
indexP_raster_p50 <- app(indexP_raster_p50, fun="mean")
indexP_raster_p50 <- raster::raster(indexP_raster_p50)
indexP_raster_p25 <- rast(indexP_raster_p25)
indexP_raster_p25 <- app(indexP_raster_p25, fun="mean")
indexP_raster_p25 <- raster::raster(indexP_raster_p25)
indexP_raster_p75 <- rast(indexP_raster_p75)
indexP_raster_p75 <- app(indexP_raster_p75, fun="mean")
indexP_raster_p75 <- raster::raster(indexP_raster_p75)
indexP_raster_sd <- (indexP_raster_p75-indexP_raster_p25)/(2*0.6744898)
indexP_raster <- raster::brick(list(indexP_raster_p50, indexP_raster_sd))

# simulate the distribution
simulate_dist <- function(x, n_samples) {
  if (all(is.na(x))) return(c(NA, NA, NA))
  samples <- rnorm(n=n_samples, mean=x[1], sd=x[2])
  return(c(mean(samples), quantile(samples, probs=c(0.05, 0.95))))
}
total_samples <- 1000
indexP_raster <- rast(indexP_raster)
interval_raster <- app(indexP_raster, fun=simulate_dist, n_samples=total_samples)
interval_raster <- raster::brick(interval_raster)
names(interval_raster) <- c("mean", "lower.90", "upper.90")
