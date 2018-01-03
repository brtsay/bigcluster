#!/usr/bin/Rscript

## Modify sandwich vcovCL to use bigglm and have reduced memory footprint

library(biglm)
library(bigmemory)
library(bigtabulate)
library(bootSVD)
library(data.table)
library(ff)
library(magrittr)
library(MASS)
library(Matrix)
library(parallel)

GetTheta <- function(model) {
    ## Gets the dispersion parameter for negbin models
    modelFamily <- model$family$family
    theta <- gsub("[\\(\\)]", "",
                  regmatches(modelFamily,
                             gregexpr("\\(.*?\\)", modelFamily))[[1]]) %>% as.numeric()
    return(theta)
}


CreateResid <- function(model, dataFrame, ycolName, xmat) {
    ## Create residuals for big models
    modelOffset <- 0
    if (!is.null(attr(model$terms, "offset"))) {
        offsetName <- attr(model$terms, "variables")[[attr(model$terms, "offset") + 1]] %>% as.character()
        if (grepl("log", offsetName[2])) {
            offsetVar <- gsub("[\\(\\)]", "",
                  regmatches(offsetName[2],
                             gregexpr("\\(.*?\\)", offsetName[2]))[[1]])
            modelOffset <- log(dataFrame[[offsetVar]])
        } else {
            modelOffset <- dataFrame[[offsetName[2]]]
        }
    }
    eta <- xmat %*% coef(model) + modelOffset
    mu <- exp(eta)
    Y <- dataFrame[[ycolName]]
    if (class(model)[1L] == "biglm") {
        xResid <- Y - eta
    } else if (grepl("Negative Binomial", model$family[[1]])) {
        theta <- GetTheta(model)
        xResid <- Y - mu * (Y + theta) / (mu + theta)        
    } else if (grepl("poisson", model$family[[1]])) {
        xResid <- Y - mu
    } else {
        stop("Model not supported")
    }
    return(as.vector(xResid))
}


EstfunBigGLM <- function(model, modelFormula, dataFrame) {
    ## Generate estimating functions
    xmat <- sparse.model.matrix(modelFormula, dataFrame) %>%
        naresid(model$na.action, .)
    if(any(alias <- is.na(coef(model)))) xmat <- xmat[, !alias, drop = FALSE]
    ## get residuals
    ycolName <- as.character(modelFormula[2])
    xResid <- CreateResid(model, dataFrame, ycolName, xmat)
    
    rval <- xResid * xmat
    attr(rval, "assign") <- NULL
    attr(rval, "contrasts") <- NULL
    return(rval)
}

UnscaledCov <- function(object, dispersion = NULL,  ...) {
    ## Return the unscaled covariance matrix
    ## equivalent to summary(model)$cov.unscaled
    ## taken from biglm
    if(!object$qr$checked)
        object$qr <- .Call("singcheckQR", object$qr)
    p <- length(object$qr$D)
    R <- diag(p)
    R[row(R)>col(R)] <- object$qr$rbar
    R <- t(R) * sqrt(object$qr$D)
    ok <- object$qr$D!=0
    R[ok,ok] <- chol2inv(R[ok,ok])
    R[!ok,] <- R[,!ok] <- NA
    dimnames(R) <- list(object$names, object$names)
    return(R)
}

ClusterIndexSum <- function(x, storage, ef, j) {
    interestIndices <- storage[index %in% x]
    rowIdx <- unlist(interestIndices$rowIdx)
    cutted <- findInterval(rowIdx, vec = ef@p,
                           rightmost.closed = FALSE, left.open = TRUE)
    group <- which(cutted == j)
    result <- sum(ef@x[rowIdx[group]])
    return(result)    
}

BigMeat <- function(bigmodel, clusters, type, dataFrame, sumMethod, multi0=FALSE, cadjust=TRUE) {
    ## Generates the "meat" of the clustered standard errors
    modelFormula <- formula(bigmodel$terms)
    ef <- EstfunBigGLM(bigmodel, modelFormula, dataFrame) # scores
    k <- ncol(ef)
    n <- nrow(ef)

    cluster <- data.frame(clusters)
    stopifnot(nrow(cluster) == n)
    p <- ncol(cluster)
    if (p > 1L) {
        ## generate cluster combinations
        clustCombs <- lapply(1L:p, function(j)
            combn(1L:p, j, simplify = FALSE))
        clustCombs <- unlist(clustCombs, recursive = FALSE)
        ## which clusters to subtract from others
        signs <- sapply(clustCombs, function(clustComb) (-1L)^(length(clustComb) + 1L))
        for (i in (p + 1L):length(clustCombs)) {
            ## these should be the combination clusters
            ## e.g. comb of cluster1 cluster2, etc
            stopifnot(length(clustCombs[[i]]) > 1)
            ## cluster <- cbind(cluster,
            ##                  Reduce(paste0, c(cluster[, clustCombs[[i]]])))
            combinedCluster <- apply(cluster[, clustCombs[[i]]],
                                     1, paste, collapse = "-")
            cluster <- cbind(cluster, combinedCluster)
        }
        if (multi0) cluster[[length(clustCombs)]] <- 1L:n        
    } else {
        clustCombs <- list(1)
        signs <- 1
    }


    ## number of clusters (and cluster interactions)
    ## e.g. if there are two cluster variables
    ## g will be of length 3
    ## 1st element: number cluster values for 1st cluster
    ## 2nd element: number cluster values for 2nd cluster
    ## 3rd element: number of unique cluster1-cluster2 combinations
    g <- sapply(1L:length(clustCombs), function(i) {
        if(is.factor(cluster[[i]])) {
            length(levels(cluster[[i]]))
        } else {
            length(unique(cluster[[i]]))
        }
    })
    gmin <- min(g[1L:p])
    if(FALSE) g[] <- gmin

    ## bias correction
    if (is.null(type)) {
        type <- ifelse(class(bigmodel)[1L] == "biglm", "HC1", "HC0")
    }
    
    rval <- matrix(0, nrow = k, ncol = k,
                   dimnames = list(colnames(ef), colnames(ef)))

    cluster <- as.big.matrix(cluster)

    ## add OPG for each cluster-aggregated estfun
    for (i in 1:length(clustCombs)) {
        ## estimating functions for aggregation by i-th clustering variable
        ## adjustments for bias correction
        adj <- if (multi0 & (i == length(clustCombs))) {
                   ifelse(type == "HC1", (n - k)/(n - 1L), 1)
               } else {
                   g[i] / (g[i] - 1L)
               }
        if (g[i] < n) {
            if (sumMethod == "slow") {
                storage <- data.table(index = ef@i,
                                      rowIdx = 1:length(ef@i))
                storage <- storage[, .(rowIdx = list(rowIdx)),
                                   by = index]
                clusterIndex <- bigsplit(cluster, i, splitcol = NA_real_)
                efi <- ff(vmode = "double", dim = c(length(clusterIndex), ncol = ncol(ef)))
                for (j in 1:ncol(ef)) {
                    efi[, j] <- mclapply(clusterIndex, function(x) {
                        rowIdx <- unlist(storage[x]$rowIdx, recursive = FALSE, use.names = FALSE)
                        cutted <- findInterval(rowIdx, vec = ef@p,
                                                  rightmost.closed = FALSE, left.open = TRUE)
                        group <- which(cutted == j)
                        result <- sum(ef@x[rowIdx[group]])
                        return(result)
                    })
                    if (j %% 10 == 0) {
                        message("Completed column ", j, "/", ncol(ef), " of cluster ", i, "/", length(clustCombs))
                    }
                }
            } else if (sumMethod == "fast") {
                efi <- apply(ef, 2L, tapply, cluster[,i], sum)
            }

        } else efi <- as.matrix(ef)
        crossproduct <- ffmatrixmult(efi, xt = TRUE, override.big.error = TRUE)
        crossproduct <- as.matrix(crossproduct[,])
        rval <- rval + signs[i] * adj * crossproduct/n
    }
    if (type == "HC1") rval <- (n - 1L)/(n - k) * rval
    if (!cadjust) rval <- (gmin - 1L)/gmin * rval
    return(rval)
}


BigCluster <- function(bigmodel, clusters, type = NULL, dataFrame, sumMethod = "slow") {
    ## Computes the cluster-robust covariance matrix for biglm objects

    ## Args:
    ##    bigmodel (biglm): Your model
    ##    clusters (data.frame): Cluster assignments for every observation
    ##    type (char): "HC0" or "HC1". Defaults to HC0 for bigglm and HC1 for biglm
    ##    dataFrame (data.frame): Your data
    ##    sumMethod (char): Can be "slow" or "fast". The fast option
    ##       is faster but uses more memory. The slow option is the
    ##       opposite.

    ## Returns:
    ##    Cluster-robust covariance matrix
    
    if (!match(class(bigmodel)[1L], c("biglm", "bigglm"))) stop("Model not supported")
    meat <- BigMeat(bigmodel, clusters, type, dataFrame, sumMethod)
    bread <- UnscaledCov(bigmodel) * bigmodel$n
    return(1/bigmodel$n * bread %*% meat %*% bread)
}

