######### Supplementary Material - A scalable method for the estimation of spatial downscaling models
####### This file contains the method for parameter estimation
######### 

rm(list=ls())
library(raster)
library(mgcv)
library(snow)
library(parallel)
library(Rmpi)

##########
## First section: some preliminary variables and functions for the method
#Read the aggregated observed variable Y
df.Y <- readRDS('0bd.rds')

#Defines the link function and its derivative
g <- function(x) exp(x)
dg <- function(x) exp(x)

#Function 'fit' is the likelihood for numerical optimization
fit <- function(rho, Qty., R, r, N, b0, type = 1){
	#rho[1] = variance, rho[2:4] = smoothing parameters, w = linearization weight

        rho <- exp(rho)

        if(any(!is.finite(rho))) return(NA)
        if(any(sapply(rho, all.equal, 0) == TRUE)) return(NA)

        Sl <- rho[2]*S1 + rho[3]*S2 + rho[4]*S3 + rho[5]*S4
	eigenSl <- eigen(Sl[-1,-1], symmetric = TRUE)
	ldetSl <- eigenSl$values
	if(any(ldetSl <= 0)) return(NA)
	ldetSl <- sum(log(ldetSl))

        V <- crossprod(R)/rho[1] + Sl
        eigenV <- eigen(V)
        ldetV <- eigenV$values
	if(any(ldetV <= 0)){
		return(NA)
	}
	ldetV <- sum(log(ldetV))

        V.mod <- V + w*diag(ncol(V))
        eigenV.mod <- eigen(V.mod)
	if(any(eigenV.mod$values <= 0)){
		return(NA)
	}
        sqrtinv <- sqrt(1/eigenV.mod$values)
        V.mod.inv <- crossprod(sqrtinv * t(eigenV.mod$vectors))

        hat.beta <- V.mod.inv %*% (crossprod(R, Qty.)/rho[1] + w*b0)

	logLik <- as.numeric(- crossprod(Qty. - R %*% hat.beta)/rho[1] - r/rho[1] - t(hat.beta) %*% Sl %*% hat.beta - w*crossprod(hat.beta - 
b0) - N * log(rho[1]) + ldetSl - ldetV)
        if(type == 2){
                attr(logLik, 'hat.beta') <- hat.beta
        }

	return(-logLik)
}

#Function 'grfit' is the derivative of the likelihood w.r.t. the vector [log(sigma^2), log(lambda)]
grfit <- function(rho, Qty., R, r, N, b0, type = 1){
	#rho[1] = variance, rho[2:4] = smoothing parameters, w = linearization weight

        rho <- exp(rho)
        if(any(!is.finite(rho))) return(NA)
        if(any(sapply(rho, all.equal, 0) == TRUE)) return(NA)

        Sl <- rho[2]*S1 + rho[3]*S2 + rho[4]*S3 + rho[5]*S4
	eigenSl <- eigen(Sl[-1,-1], symmetric = TRUE)
	ldetSl <- eigenSl$values
	if(any(ldetSl <= 0)) return(NA)
	ldetSl <- sum(log(ldetSl))
        sqrtinv <- sqrt(1/eigenSl$value)
        Sl.inv <- rbind(0, cbind(0, crossprod(sqrtinv * t(eigenSl$vectors))))

        V <- crossprod(R)/rho[1] + Sl
        eigenV <- eigen(V)
        ldetV <- eigenV$values
	if(any(ldetV <= 0)){
		return(NA)
	}
	ldetV <- sum(log(ldetV))
        sqrtinv <- sqrt(1/eigenV$value)
        V.inv <- crossprod(sqrtinv * t(eigenV$vectors))

        V.mod <- V + w*diag(ncol(V))
        eigenV.mod <- eigen(V.mod)
	if(any(eigenV.mod$values <= 0)){
		return(NA)
	}
        sqrtinv <- sqrt(1/eigenV.mod$value)
        V.mod.inv <- crossprod(sqrtinv * t(eigenV.mod$vectors))

        hat.beta <- V.mod.inv %*% (crossprod(R, Qty.)/rho[1] + w*b0)

	gr <- rep(NA, 5)
	#Below, tr(V.inv %*% R'R/rho[1]) = sum(V.inv * crossprod(R)/rho[1])
	gr[1] <- crossprod(Qty. - R %*% hat.beta)/rho[1] + r/rho[1] - N + sum(V.inv * crossprod(R)/rho[1])
	for(i in 2:5){
		gr[i] <- - t(hat.beta) %*% (rho[i] * get(paste0('S', i-1))) %*% hat.beta
		gr[i] <- gr[i] + sum(Sl.inv * rho[i] * get(paste0('S', i-1)))
		gr[i] <- gr[i] - sum(V.inv * rho[i] * get(paste0('S', i-1)))
	}

        return(-gr)
}

##########
## Second section: from here on, trying the algorithm
#Reads the spatial files
x3 <- raster('0rx3.tif')
x4 <- raster('0rx4.tif')
block <- raster('0rblock.tif')

#Builds a reference smoother based on a coarse version of the data
set.seed(1)
m.x3 <- sampleRegular(x3, 1e5, asRaster = TRUE, useGDAL = TRUE)
m.x4 <- sampleRegular(x4, 1e5, asRaster = TRUE, useGDAL = TRUE)
m.x1 <- coordinates(m.x3)[,1]
m.x2 <- coordinates(m.x3)[,2]
sm1 <- smoothCon(te(x1, x2, bs = 'cs', k = 20), data = data.frame(x1 = m.x1, x2 = m.x2), absorb.cons = TRUE, scale.penalty = TRUE); sm1[[1]]$X <- NULL
sm2 <- smoothCon(s(x3, bs = 'cs', k = 20), data = data.frame(x3 = m.x3[]), absorb.cons = TRUE, scale.penalty = TRUE); sm2[[1]]$X <- NULL
sm3 <- smoothCon(s(x4, bs = 'cs', k = 20), data = data.frame(x4 = m.x4[]), absorb.cons = TRUE, scale.penalty = TRUE); sm3[[1]]$X <- NULL
saveRDS(sm1, 'out/1sm1.rds')
saveRDS(sm2, 'out/1sm2.rds');
saveRDS(sm3, 'out/1sm3.rds')
#Important note: in an actual application, 'sampleRegular' may miss information on extreme values (e.g., max. or min.)
#This may generate poor smoothers with 'smoothCon', so must be addressed somehow.

#In this block, we generate the penalty matrices for smoothers
#In this simulated example, our matrix \Omega_{\vect \lambda}^{-1} = S1 * rho[2] + S2 * rho[3] + S3 * rho[4]
bd <- sapply(1:3, function(i) get(paste0('sm', i))[[1]]$df)
bd <- 1 + sum(bd)
S1 <- S2 <- S3 <- S4 <- matrix(0, bd, bd); idx <- 1
S1[idx + (1:sm1[[1]]$df), idx + (1:sm1[[1]]$df)] <- sm1[[1]]$S[[1]]
S2[idx + (1:sm1[[1]]$df), idx + (1:sm1[[1]]$df)] <- sm1[[1]]$S[[2]]; idx <- idx + sm1[[1]]$df
S3[idx + (1:sm2[[1]]$df), idx + (1:sm2[[1]]$df)] <- sm2[[1]]$S[[1]]; idx <- idx + sm2[[1]]$df
S4[idx + (1:sm3[[1]]$df), idx + (1:sm3[[1]]$df)] <- sm3[[1]]$S[[1]]

#Here we split the raster into chunks for efficient processing
bss <- blockSize(x3, minblocks = 200)
print(paste0('bss: ', bss$n))
d <- bss$n

#And we define a function to generate the aggregated matrices from each chunk of information
ag.m <- function(i, b){
	tmp <- data.frame(
		x3 = getValues(x3, row = bss$row[i], nrows = bss$nrow[i]),
		x4 = getValues(x4, row = bss$row[i], nrows = bss$nrow[i]),
		block = getValues(block, row = bss$row[i], nrows = bss$nrow[i]))
	tmp <- data.frame(xyFromCell(x3, cellFromRow(x3, rownr = bss$row[i]:(bss$row[i] + bss$nrows[i] - 1))), tmp)
	colnames(tmp)[1:2] <- c('x1', 'x2')

	#The ordering of groups in the aggregation function
	groups <- sort(unique(tmp$block))

	tmpgam.X1 <- PredictMat(sm1[[1]], data = tmp)
	tmpgam.X2 <- PredictMat(sm2[[1]], data = tmp)
	tmpgam.X3 <- PredictMat(sm3[[1]], data = tmp)
	tmpgam <- cbind(1, tmpgam.X1, tmpgam.X2, tmpgam.X3)
	rm(tmpgam.X1, tmpgam.X2, tmpgam.X3)
	nb <- length(b)

	#These are temporary vectors that will be summed to generate the aggregated matrices
	#They will be filled with data below. 'rowsum' is a faster alternative for the 'aggregate' function
	tmp_Agz. <- tmp_AGz. <- tmp_AGX <- tmp_MMt <- rep(list(NA), nb)

	for(j in 1:nb){
		z. <- tmpgam %*% b[[j]]

		gz. <- g(z.)
		v.out <- rowsum(gz., group = tmp$block)
		tmp_Agz.[[j]] <- v.out

		Gz. <- dg(z.) * z.
		v.out <- rowsum(Gz., group = tmp$block)
		tmp_AGz.[[j]] <- v.out

		v.out <- as.numeric(dg(z.)) * tmpgam
		v.out <- rowsum(v.out, group = tmp$block)
		tmp_AGX[[j]] <- v.out

		dg2 <- 0 * gz. + 1
		v.out <- rowsum(dg2, group = tmp$block)
		tmp_MMt[[j]] <- v.out
	}
	return(list(tmp_Agz., tmp_AGz., tmp_AGX, tmp_MMt, groups))
}

#Creates the cluster for parallel processing
nodes <- 4
cl <- makeCluster(nodes, type = 'MPI')
clusterEvalQ(cl, library(raster))
clusterEvalQ(cl, library(mgcv))
clusterExport(cl=cl, list('ag.m','x3','x4','block','g','dg','sm1','sm2','sm3','bd','bss','df.Y'))

#This is a function that received a guess for the betas (vector b)
	# and a guess of starting values for the numerical optimization (vector st)
out.fn <- function(b, st){
	b <- list(matrix(b, ncol = 1))
	nb <- length(b)

	#These are the N-dimension aggregated matrices that we fill iteratively by iterating over each chunk of the data
	Agz. <- rep(list(matrix(0, nrow = dim(df.Y)[1], ncol = 2)), nb)
	AGz. <- rep(list(matrix(0, nrow = dim(df.Y)[1], ncol = 2)), nb)
	for(j in 1:nb) AGz.[[j]][,1] <- Agz.[[j]][,1] <- df.Y[,1]
	X_ <- rep(list(matrix(0, nrow = dim(df.Y)[1], ncol = bd)), nb)
	diagMMt <- rep(list(matrix(0, nrow = dim(df.Y)[1], ncol = 2)), nb)
	Y_ <- rep(list(NA), nb)
	for(j in 1:(nb)) diagMMt[[j]][,1] <- df.Y[,1]

	#In this block of code, we will iterate over chunks of the dataset, doing:
		#1) Submitting the calculation of the N-dimensional matrices for chunk "i"
		#2) Summing the output of (1) into the matrices created above
	#After iterating over all of them, our matrices will then be used for the numerical optimization step
	for(i in 1:nodes){
		sendCall(cl[[i]], ag.m, list(i = i, b = b), tag = i)
	}
	for (i in 1:d) {
		k <- recvOneData(cl)
		if(!k$value$success) saveRDS(k, 'erro.rds')
    
		ni <- nodes + i
		if (ni <= d) {
			sendCall(cl[[k$node]], ag.m, list(i = ni, b = b), tag = ni)
		}

		tag <- k$value$tag
		if(tag %% 50 == 0) print(paste0('recebido: ', tag, ' - ag.m - ', Sys.time()))

		for(j in 1:nb){
			if(!is.null(k$value$value[[5]])){
				groups <- k$value$value[[5]]

				v.out <- k$value$value[[1]][[j]]
				v.out <- c(0, v.out)[match(df.Y[,1], c(NA, groups), nomatch = 1)]
				Agz.[[j]][,2] <- Agz.[[j]][,2] + v.out

				v.out <- k$value$value[[2]][[j]]
				v.out <- c(0, v.out)[match(df.Y[,1], c(NA, groups), nomatch = 1)]
				AGz.[[j]][,2] <- AGz.[[j]][,2] + v.out

				v.out <- k$value$value[[3]][[j]]
				v.out <- rbind(0, v.out)[match(df.Y[,1], c(NA, groups), nomatch = 1),]
				X_[[j]] <- X_[[j]] + as.matrix(v.out)

				v.out <- k$value$value[[4]][[j]]
				v.out <- c(0, v.out)[match(df.Y[,1], c(NA, groups), nomatch = 1)]
				diagMMt[[j]][,2] <- diagMMt[[j]][,2] + v.out

				rm(v.out)
			}
		}
		rm(k)
	}
	for(j in 1:nb) Y_[[j]] <- df.Y[,2] - Agz.[[j]][,2] + AGz.[[j]][,2]

	p <- nb
	diagMMt <- list(diagMMt[[p]])
	X_ <- list(X_[[p]])
	Y_ <- list(Y_[[p]])
	b <- b[[p]]

	#These are the auxiliary matrices that speed up likelihood calculation
	Qty. <- R <- r <- N <- rep(list(NA), 1)
	diagL <- 0 * diagMMt[[1]][,2]
	ql <- diagMMt[[1]][,2] > 0
	diagL[ql] <- sqrt(1/diagMMt[[1]][ql,2])
	y. <- diagL * Y_[[1]]
	X. <- diagL * X_[[1]]
	X.qr <- qr(X., LAPACK = TRUE)
	Qty.[[1]] <- qr.qty(X.qr, y.)[1:X.qr$rank,]
	R[[1]] <- qr.R(X.qr)[,sort.list(X.qr$pivot)]
	r[[1]] <- crossprod(y.) - crossprod(Qty.[[1]])
	N[[1]] <- sum(ql)

# A temporary print of the progress
#	er <- sapply(1:nb, function(j){
#		crossprod(diagL * (df.Y[,2] - Agz.[[j]][,2]))
#	})
#	print(er)

	return(list(Qty.[[1]], R[[1]], r[[1]], N[[1]], b))
}

##########
## Third section: now we finally define the method parameters and run
#'maxit' contains the maximum numbers of iterations to stop if the method hasn't converged
maxit <- 50
ml <- rep(NA, maxit)
mtcs <- rep(list(NA), maxit)
it <- 0

set.seed(1)
b0 <- rnorm(bd, 0, 0.1)
st <- rep(0, 5)
j <- 0
n.dif <- 1

#'w' contains the initial guess for the linearization weight
w <- 10
lin.k <- 2

#Now, iteratively we repeat:
stop <- FALSE
while(n.dif > 1e-4 & it <= maxit){
	it <- it + 1
	print(paste0('iteration #', it))

	#1) Generate the aggregated matrices in parallel
	ll.t1 <- out.fn(b=b0, st=st)

	#2) If we are in the first iteration,
	if(it == 1){
		#2.1) Then, we optimize the likelihood numerically
		l <- tryCatch({
			optim(st, fn=fit, gr=grfit, Qty.=ll.t1[[1]], R=ll.t1[[2]], r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]], method='BFGS')
		}, error = function(e){
			NA
		})

		ml[it] <- as.numeric(l$value)
		st <- l$par
		#2.2) Get a new guess for the betas
		b0 <- attr(fit(st, Qty.=ll.t1[[1]], R=ll.t1[[2]], r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]], type = 2), 'hat.beta')
		gg <- crossprod(grfit(st, Qty.=ll.t1[[1]], R=ll.t1[[2]], r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]]))
		mtcs[[it]] <- ll.t1
		mtcs[[it]][[6]] <- w
		mtcs[[it]][[7]] <- st
		mtcs[[it]][[8]] <- b0

		print(paste0('norm gr @ minimum: ', gg))
		print(paste0('ll: ', round(ml[it], 4)))
		print('n.dif > 0 (start value)')
		print(paste0('rhos: ', paste0(round(st, 4), collapse = ' ')))
		print(paste0('s2: ', round(exp(st[1]), 8)))
		print(paste0('new weight: ', round(w, 8)))

		#2.3) And update the weight
		w <- 1/lin.k * w
	} else if(it > 1) {
		#2.1) Otherwise, we first check if the likelihood has increased
		#Important observation: we are actually minimizing -1 times the likelihood, so we check for decrease
		ml[it] <- ml[it-1]

		#2.2) While the likelihood hasn't increased
		while(ml[it] >= ml[it-1]){
			#2.3) We increase the weight
			j <- j + 1
			w <- w * lin.k
			if(j > 1){
				print(paste0('Adjusting weight and repeating iteration: ', round(w, 8)))
			}

			if(j > 5){
				stop <- TRUE
				break
			}

			#2.4) Optimize again with the new weight
			l <- tryCatch({
				optim(st, fn=fit, gr=grfit, Qty.=ll.t1[[1]], R=ll.t1[[2]], r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]], 
method='BFGS')
			}, error = function(e){
				NA
			})

			#And return to 2.2)
			if(all(is.na(l))){
				ml[it] <- ml[it-1]
			} else {
				ml[it] <- as.numeric(fit(l$par, Qty.=ll.t1[[1]], R=ll.t1[[2]], r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]]))
			}
		}
		if(stop){
			ml[it] <- NA
			print('Finishing.')

			break
		}

		#3) When the likelihood increased
		n.dif <- ml[it-1] - ml[it]
		#3.1) We store new starting values for the numerical optimization and for the betas
		st <- l$par
		b0 <- attr(fit(st, Qty.=ll.t1[[1]], R=ll.t1[[2]], r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]], type = 2), 'hat.beta')
		gg <- crossprod(grfit(st, Qty.=ll.t1[[1]], R=ll.t1[[2]], r=ll.t1[[3]], N=ll.t1[[4]], b0=ll.t1[[5]]))

		mtcs[[it]] <- ll.t1
		mtcs[[it]][[6]] <- w
		mtcs[[it]][[7]] <- st
		mtcs[[it]][[8]] <- b0

		print(paste0('norm gr @ minimum: ', gg))
		print(paste0('ll: ', round(ml[it], 4)))
		print(paste0('n.dif: ', round(n.dif, 4)))
		print(paste0('rhos: ', paste0(round(st, 4), collapse = ' ')))
		print(paste0('s2: ', round(exp(st[1]), 8)))

		#3.2) Update the weight, and continue
		w <- w * 0.75
		print(paste0('new weight: ', round(w, 8)))

		w <- 1/lin.k * w
		j <- 0
	}
}
#At the end, we save the results in these files
saveRDS(mtcs, 'out/mtcs.rds')
saveRDS(mtcs[[it-1]], 'out/mtcs_sol.rds')

