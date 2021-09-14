######### Supplementary Material - A scalable method for the estimation of spatial downscaling models
####### This file is an analysis script to look at the results
######### 

rm(list=ls())
library(raster)
library(Rmpi)
library(mgcv)
library(snow)
library(parallel)

#Loads the observed variable Y
df.Y <- readRDS('0bd.rds')

#Stores the link function and derivative, and the explanatory variables, as used for the fit
g <- function(x) exp(x)
dg <- function(x) exp(x)
x3 <- raster('0rx3.tif')
x4 <- raster('0rx4.tif')
y.ok <- raster('0rz_correct.tif')

#Loads the reference smoothers and create penalty matrices
sm1 <- readRDS('out/1sm1.rds')
sm2 <- readRDS('out/1sm2.rds')
sm3 <- readRDS('out/1sm3.rds')
bd <- sapply(1:3, function(i) get(paste0('sm', i))[[1]]$df)
bd <- 1 + sum(bd)
S1 <- S2 <- S3 <- S4 <- matrix(0, bd, bd); idx <- 1
S1[idx + (1:sm1[[1]]$df), idx + (1:sm1[[1]]$df)] <- sm1[[1]]$S[[1]]
S2[idx + (1:sm1[[1]]$df), idx + (1:sm1[[1]]$df)] <- sm1[[1]]$S[[2]]; idx <- idx + sm1[[1]]$df
S3[idx + (1:sm2[[1]]$df), idx + (1:sm2[[1]]$df)] <- sm2[[1]]$S[[1]]; idx <- idx + sm2[[1]]$df
S4[idx + (1:sm3[[1]]$df), idx + (1:sm3[[1]]$df)] <- sm3[[1]]$S[[1]]

#Sample the explanatory variables to create the reference curves
set.seed(1)
m.x3 <- sampleRegular(x3, 1e5, asRaster = TRUE, useGDAL = TRUE)
m.x4 <- sampleRegular(x4, 1e5, asRaster = TRUE, useGDAL = TRUE)
m.y.ok <- sampleRegular(y.ok, 1e5, asRaster = TRUE, useGDAL = TRUE)

#Loads the solution to the iterative estimation
mtcs <- readRDS('out/mtcs_sol.rds')
Qty. <- mtcs[[1]]
R <- mtcs[[2]]
r <- mtcs[[3]]
N <- mtcs[[4]]
b <- mtcs[[5]]
w <- mtcs[[6]]
st <- mtcs[[7]]

#And generates the matrices necessary for prediction
rho <- exp(st)
Sl <- rho[2]*S1 + rho[3]*S2 + rho[4]*S3 + rho[5]*S4
V <- crossprod(R)/rho[1] + Sl
V.mod <- V + w*diag(ncol(V))
eigenV.mod <- eigen(V.mod)
sqrtinv <- sqrt(1/eigenV.mod$values)
V.mod.inv <- crossprod(sqrtinv * t(eigenV.mod$vectors))
hat.beta <- V.mod.inv %*% (crossprod(R, Qty.)/rho[1] + w*b)
hat.cov <- V.mod.inv

#Constructs an X matrix for predictions
tmp <- coordinates(m.x3)
tmp <- data.frame(x1 = tmp[,1], x2 = tmp[,2], x3 = m.x3[], x4 = m.x4[])
tmpgam.X1 <- PredictMat(sm1[[1]], tmp)
tmpgam.X2 <- PredictMat(sm2[[1]], tmp)
tmpgam.X3 <- PredictMat(sm3[[1]], tmp)
X <- cbind(1, tmpgam.X1, tmpgam.X2, tmpgam.X3)

#And generates credible intervals by simulation from the variance-covariance matrix
set.seed(1)
Cv <- chol(hat.cov, pivot = TRUE)
Cv <- Cv[,order(attr(Cv, 'pivot'))]
n.rep = 10000
nb <- dim(hat.cov)[1]
br <- t(Cv) %*% matrix(rnorm(n.rep*nb),nb,n.rep)
br <- as.numeric(hat.beta) + br

#Generates s(x3)
out <- X[,401:419] %*% br[401:419,]
out.mean <- apply(out, 1, mean)
out.sd <- apply(out, 1, sd)
om <- out.mean[order(tmp$x3)]
os <- out.sd[order(tmp$x3)]
y <- m.y.ok[order(tmp$x3)]
d <- m.x3[order(tmp$x3)]
rd <- sample(1:length(y), 5000)
y.s <- y[rd]
d.s <- d[rd]
y.tr <- -4*sin(d)/(1 + d^2)
## Here we zero-center the true smoother, and subtract the intercept from the data points
### 'om' is already centred by default
y.s <- y.s - (-4.5 + mean(y.tr))
y.tr <- y.tr - mean(y.tr)

### Export the true and estimated smoothers to separate pdfs
pdf('1rout_s(x3)-true.pdf')
plot(om ~ d, type = 'l', ylim = c(min(y.s), max(y.s)), xlim = c(min(d), max(d)), main = 's(x3) true', xlab = 'x3', ylab = 's(x3)', col = NA, cex.main = 1.5, cex.lab = 1.5)
lines(y.tr ~ d, xlim = c(min(d), max(d)), lty = 2, lwd = 2.5, col = 'red')
dev.off()

pdf(paste0('1rout_s(x3)-estimated.pdf'))
plot(om ~ d, type = 'l', ylim = c(min(y.s), max(y.s)), xlim = c(min(d), max(d)), main = 's(x3) estimated', xlab = 'x3', ylab = 's(x3)', cex.main = 1.5, cex.lab = 1.5)
points(d.s, y.s, type = 'p', col = rgb(0, 0, 0, 150, maxColorValue = 255), pch = '.', cex = 2)
polygon(c(d, rev(d)), c(om - os, rev(om + os)), col = rgb(80, 80, 80, 50, maxColorValue = 255), border = NA)
lines(om ~ d, col = 'red', lwd = 2.5)
dev.off()

#Generates s(x4)
out <- X[,420:438] %*% br[420:438,]
out.mean <- apply(out, 1, mean)
out.sd <- apply(out, 1, sd)
om <- out.mean[order(tmp$x4)]
os <- out.sd[order(tmp$x4)]
y <- m.y.ok[order(tmp$x4)]
d <- m.x4[order(tmp$x4)]
rd <- sample(1:length(y), 5000)
y.s <- y[rd]
d.s <- d[rd]
y.tr <- rep(0, length(d))
y.s <- y.s - (-4.5 + mean(y.tr))
y.tr <- y.tr - mean(y.tr)

## Export the true and estimated smoothers to separate pdfs
pdf('1rout_s(x4)-true.pdf')
plot(om ~ d, type = 'l', ylim = c(min(y.s), max(y.s)), xlim = c(min(d), max(d)), main = 's(x4) true', xlab = 'x4', ylab = 's(x4)', col = NA, cex.main = 1.5, cex.lab = 1.5)
abline(h = 0, col = 'red', lty = 2, lwd = 2.5)
dev.off()

pdf(paste0('1rout_s(x4)-estimated.pdf'))
plot(om ~ d, type = 'l', ylim = c(min(y.s), max(y.s)), xlim = c(min(d), max(d)), main = 's(x4) estimated', xlab = 'x4', ylab = 's(x4)', cex.main = 1.5, cex.lab = 1.5)
points(d.s, y.s, type = 'p', col = rgb(0, 0, 0, 150, maxColorValue = 255), pch = '.', cex = 2)
polygon(c(d, rev(d)), c(om - os, rev(om + os)), col = rgb(80, 80, 80, 50, maxColorValue = 255), border = NA)
lines(om ~ d, lwd = 2.5, col = 'red')
dev.off()

#Gere, we generate raster, containing: the response variable (y) predicted, and s(x1, x2) -> true, estimated, and s.d.
bss <- blockSize(x3, minblocks = 200)
print(paste0('bss: ', bss$n))
d <- bss$n
rout <- raster(x3)
out_mean <- writeStart(rout, '1rout_mean.tif', datatype = 'FLT4S', format = 'GTiff', overwrite = TRUE)
out_x12_mean <- writeStart(rout, '1rout_s(x1,x2)-mean.tif', datatype = 'FLT4S', format = 'GTiff', overwrite = TRUE)
out_x12_sd <- writeStart(rout, '1rout_s(x1,x2)-sd.tif', datatype = 'FLT4S', format = 'GTiff', overwrite = TRUE)
out_x12_true <- writeStart(rout, '1rout_s(x1,x2)-true.tif', datatype = 'FLT4S', format = 'GTiff', overwrite = TRUE)
for(i in 1:d){
	print(i)
	tmp <- data.frame(
		x3 = getValues(x3, row = bss$row[i], nrows = bss$nrow[i]),
		x4 = getValues(x4, row = bss$row[i], nrows = bss$nrow[i]))
	tmp <- data.frame(xyFromCell(x3, cellFromRow(x3, rownr = bss$row[i]:(bss$row[i] + bss$nrows[i] - 1))), tmp)
	colnames(tmp)[1:2] <- c('x1', 'x2')

	groups <- sort(unique(tmp$block))
	tmpgam.X1 <- PredictMat(sm1[[1]], data = tmp)
	tmpgam.X2 <- PredictMat(sm2[[1]], data = tmp)
	tmpgam.X3 <- PredictMat(sm3[[1]], data = tmp)
	tmpgam <- cbind(1, tmpgam.X1, tmpgam.X2, tmpgam.X3)
	rm(tmpgam.X1, tmpgam.X2, tmpgam.X3)

	out <- tmpgam %*% br
	out.mean <- sapply(1:nrow(out), function(j) mean(out[j,]))
	out_mean <- writeValues(out_mean, out.mean, bss$row[i])

	out <- tmpgam[,2:400] %*% br[2:400,]
	out.mean <- apply(out, 1, mean)
	out.sd <- apply(out, 1, sd)

	out_x12_mean <- writeValues(out_x12_mean, out.mean, bss$row[i])
	out_x12_sd <- writeValues(out_x12_sd, out.sd, bss$row[i])
	out_x12_true <- writeValues(out_x12_true, {
		with(tmp, 3*((x1 - 20)/100)^2*((x2 + 15)/100) - 2*((x2 + 25)/100)^2*((x1 - 5)/100))
	}, bss$row[i])
}
out_mean <- writeStop(out_mean)
out_x12_mean <- writeStop(out_x12_mean)
out_x12_sd <- writeStop(out_x12_sd)
out_x12_true <- writeStop(out_x12_true)

