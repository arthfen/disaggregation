######### Supplementary Material - A scalable method for the estimation of spatial downscaling models
####### This file creates the simulated aggregated variable
######### 

rm(list=ls())
library(raster)
library(rgdal)

# Random seed for reproducibility
set.seed(1)
nc <- nr <- 1000
n <- nc * nr

# Loads the raster with explanatory variables
rx3 <- raster('0rx3.tif')
x3 <- rx3[]
x1 <- coordinates(rx3)[,1]
x2 <- coordinates(rx3)[,2]

# Simulates residuals, epsilon ~ N(0, I * sigma_e)
e <- rnorm(n, 0, sd=0.1)

# The (unobserved during fit) vector z is stored here as 'z.correct'
# Then, the next lines will only export the simulated datasets
z.correct <- -4.5 - 4*sin(x3)/(1 + x3^2) + 3*((x1 - 20)/100)^2*((x2 + 15)/100) - 2*((x2 + 25)/100)^2*((x1 - 5)/100)
bl <- readOGR('0rblock.gpkg')
bl <- rasterize(bl, rx3)
writeRaster(bl, '0rblock.tif', overwrite = TRUE)

Y <- data.frame(bl = bl[], y = (exp(z.correct) + e))
Y <- aggregate(y ~ bl, Y, FUN = sum)

bl[] <- Y[match(bl[], Y[,1]), 2]
writeRaster(bl, '0rY.tif', overwrite = TRUE)

colnames(Y) <- c('area', 'Y')
saveRDS(Y, '0bd.rds')

bl[] <- z.correct
writeRaster(bl, '0rz_correct.tif', overwrite = TRUE)

