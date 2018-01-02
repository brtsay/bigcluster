#!/usr/bin/Rscript

## Function to test equivalence of bigcluster and vcovCL

library(sandwich)
source("bigcluster_sandwich.R")
source("prof_func.R")

nRows <- 1000
nRegressors <- 5
nClusterVars <- 2
nSims <- 2
nClusterCats <- 9
## modelFamily options are poisson, negbin, or ols
modelFamily <- "ols"
errorsPerc <- list(rep(0, nSims))
errorsMag <- list(rep(0, nSims))
regvcovList <- list(rep(0, nSims))

for (nsim in 1:nSims) {
    dataClustersTheta <- GenerateData(nRows, nRegressors, nClusterVars, nClusterCats)
    models <- tryCatch(
    {
        models <- RunModels(dataClustersTheta[[1]], dataClustersTheta[[3]], startValues = c(0, dataClustersTheta[[4]]), modelFamily = modelFamily)
    },
    error = function(cond) {
        message(cond)
        message("Regenerating data and rerunning models")
        dataClustersTheta <- GenerateData(nRows, nRegressors, nClusterVars, nClusterCats)
        models <- RunModels(dataClustersTheta[[1]], dataClustersTheta[[3]], startValues = c(0, dataClustersTheta[[4]]), modelFamily = modelFamily)
    }
    )
    regvcov <- vcovCL(models[[1]], cluster = dataClustersTheta[[2]])
    bigvcov <- BigCluster(models[[2]],
                          clusters = dataClustersTheta[[2]],
                          dataFrame = dataClustersTheta[[1]],
                          sumMethod = "fast")
    randomRow <- sample(1:nrow(bigvcov), 1)
    randomCol <- sample(1:ncol(bigvcov), 1)
    print("Random entry in covariance matrix")
    print(paste(bigvcov[randomRow, randomCol], "(big)"))
    print(paste(regvcov[randomRow, randomCol], "(sandwich)"))
    errorsPerc[[nsim]] <- abs((bigvcov - regvcov)/regvcov)
    errorsMag[[nsim]] <- abs(bigvcov - regvcov)
    regvcovList[[nsim]] <- regvcov
    maxError <- which(errorsPerc[[nsim]] == max(errorsPerc[[nsim]]), arr.ind = TRUE)
    print("Entry where sandwich and big covariance matrix differ the most")
    print(paste(bigvcov[maxError[1], maxError[2]], "(big)"))
    print(paste(regvcov[maxError[1], maxError[2]], "(sandwich)"))
}

print("Mean error (percent)")
print(sapply(errorsPerc, mean))
print("Max error (magnitude)")
print(sapply(errorsMag, max))

## plot
errorsVec <- as.vector(unlist(errorsPerc))
regVec <- as.vector(unlist(regvcovList))
plot(log(abs(regVec)), errorsVec)

