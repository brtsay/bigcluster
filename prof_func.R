## Functions for profiling different cluster functions

library(sandwich)
source("bigcluster_sandwich.R")

GenerateX <- function(nRows, nRegressors, nClusterVars, nClusterCats) {
    ## Generates random data
    nClusterCols <- nClusterVars * (nClusterCats - 1)
    X <- matrix(rnorm(nRegressors * nRows),
                nrow = nRows, ncol = nRegressors)
    Beta <- runif(n = nRegressors + nClusterCols, min = -0.1, max = 0.1)

    clusterMatrix <- matrix(0, nrow = nRows, ncol = nClusterCols)
    clusterList <- list()
    startIdx <- 1
    for (i in 1:nClusterVars) {
        cluster <- as.factor(sample(1:nClusterCats,
                                    size = nRows,
                                    replace = TRUE))
        nClusterVal <- length(unique(cluster))
        clusterList[[i]] <- cluster
        endIdx <- startIdx + nClusterVal - 2
        clusterMatrix[1:nRows, startIdx:endIdx] <- model.matrix(~ cluster - 1)[, -1]
        startIdx <- endIdx + 1
    }
    
    X <- cbind(X, clusterMatrix)
    clusters <- do.call(cbind, clusterList)
    return(list(X, Beta, clusters))
}

GenerateData <- function(nRows, nRegressors, nClusterVars, nClusterCats) {
    y <- NA
    while (any(is.na(y))) {
        message('NA in y. Regenerating data')
        Xinfo <- GenerateX(nRows, nRegressors, nClusterVars, nClusterCats)
        errors <- rep(0, nRows)
        for (j in 1:ncol(Xinfo[[3]])) {
            scaled <- scale(Xinfo[[3]][, j])
            errors <- errors + rnorm(n = scaled,
                                     mean = scaled/100)
        }
        ## errors <- rnorm(nRows)
        ## print(mean(errors))
        mu <- exp(Xinfo[[1]] %*% Xinfo[[2]] + errors)
        theta <- runif(1, 0, 5)
        y <- rnegbin(nRows, mu, theta)
    }
    dataFrame <- data.frame(y, Xinfo[[1]])
    ## Xinfo[[3]] are the clusters
    return(list(dataFrame, Xinfo[[3]], theta, Xinfo[[2]]))
}

RunModels <- function(dataFrame, theta, startValues, modelFamily = "negbin") {
    modelFormula <- formula(paste("y ~", paste(names(dataFrame)[-1], collapse = "+")))
    models <- switch(modelFamily,
                     negbin = {
                         regmodel <- glm(modelFormula,
                                         data = dataFrame,
                                         start = startValues,
                                         family = negative.binomial(theta))
                         bigmodel <- bigglm(modelFormula,
                                            data = dataFrame,
                                            family = negative.binomial(theta),
                                            start = coef(regmodel),
                                            maxit = 200)
                         models <- list(regmodel, bigmodel)
                     },
                     poisson = {
                         regmodel <- glm(modelFormula,
                                         data = dataFrame,
                                         start = startValues,
                                         family = "quasipoisson")
                         bigmodel <- bigglm(modelFormula,
                                            data = dataFrame,
                                            family = quasipoisson(),
                                            start = coef(regmodel),
                                            maxit = 200)
                         models <- list(regmodel, bigmodel)
                     },
                     ols = {
                         regmodel <- lm(modelFormula,
                                        data = dataFrame)
                         bigmodel <- biglm(modelFormula,
                                           data = dataFrame)
                         models <- list(regmodel, bigmodel)
                     })
    return(models)
}


ProfileCluster <- function(type, ...) {
    startTime <- Sys.time()
    if (type == "reg") {
        profilename <- "profiles/reg.out"
        Rprofmem(filename = profilename, threshold = 1e5)                
        clustered <- vcovCL(...)
    } else if (type == "big") {
        profilename <- "profiles/big.out"
        Rprofmem(filename = profilename, threshold = 1e5)
        clustered <- BigCluster(...)
    }
    endTime <- Sys.time()
    timeUsed <- endTime - startTime
    Rprofmem(NULL)
    ## read Rprofmem output
    con <- file(profilename, open = "r")
    highestMem <- 0
    while (length(oneline <- noquote(readLines(con, n = 1))) > 0) {
        ## lines are always split by a colon
        splitted <- strsplit(oneline, ":")
        memUsed <- as.numeric(splitted[[1]][1])
        if (is.na(memUsed)) {}
        else if (memUsed > highestMem) highestMem <- memUsed
    }
    close(con)
    return(c(highestMem, timeUsed))
}
