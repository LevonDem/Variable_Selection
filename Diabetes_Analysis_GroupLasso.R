####################################################
## Description: Diabetes dataset variable selection
## using penalized regression (group lasso)
## Author: Levon Demirdjian
## Email: levondem@gmail.com
## Last updated: 05/30/2018
####################################################

library(grpregOverlap)

## Import design matrix and response
designMatrix <- read.table('~/Research/Ray-Bing/Diabetes_data/Metabolite_DATA_top1KCV.txt', header = TRUE)
response     <- read.table('~/Research/Ray-Bing/Diabetes_data/GSE25462_EXP.txt', sep = "\t", header = TRUE)

## Save design and response into matrices
X <- t(as.matrix(designMatrix[, 2:51]))
Y <- as.matrix(response$Sensitivity)

## Delete missing data
X <- X[which(!is.na(Y)), ]
Y <- as.matrix(Y[which(!is.na(Y))])
n <- nrow(X)

## Define variables
Xvariables  <- as.character(designMatrix[,1])
colnames(X) <- Xvariables

## Import group structure
group_data <- read.table('~/Research/Ray-Bing/Diabetes_data/Metabolite_KEGG_top1KCV.txt', sep = '\t', header = TRUE)
group_data <- group_data[as.character(group_data[,2]) != "", ]

## Eliminate variables that don't belong to a known pathway
invalid_vars <- !(Xvariables %in% as.character(group_data[,1]))
X            <- X[,-which(invalid_vars)]
Xvariables   <- colnames(X)
nvars        <- length(Xvariables)

## The jth element of groupStructure indicates which groups the jth covariate belongs to
uniqueGroups   <- c()
groupStructure <- list()
for(i in 1:nvars){
  groupStructure[[i]] <- unlist(strsplit(as.character(group_data[i,2]), '///'))
  uniqueGroups        <- c(uniqueGroups, groupStructure[[i]])
}
uniqueGroups <- unique(uniqueGroups) ## Keep this ordering for groups (i.e. group 1 is hsa03030)
nGroups1     <- length(uniqueGroups)


##############################
## Organize group structure ##
##############################

## "group" contains the group structure. 
## The jth element represents the variables in group j.

group        <- lapply(1:nGroups1, function(x) NA)
names(group) <- uniqueGroups
for(i in 1:nvars){
  curVar    <- Xvariables[i]
  curIndex  <- which(as.character(group_data[,1]) == curVar)
  curGroups <- (unlist(strsplit(as.character(group_data[curIndex,2]), '///')))
  for(j in 1:length(curGroups)){
    locIndex <- which(uniqueGroups == curGroups[j])
    group[[locIndex]] <- c(group[[locIndex]], curVar)
  }
}

## Delete empty groups
group   <- lapply(group, function(x) x[!is.na(x)])
group   <- group[sapply(group, length) != 0]
p       <- as.numeric(sapply(group, length))
nGroups <- length(p)

## Remove very large groups
group   <- group[p < 100]
p       <- p[p < 100]
nGroups <- length(p)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~ Try Penalized regression approach ~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

# Run penalized regression using cross validation to pick a model
cvfit1 <- cv.grpregOverlap(X, Y, group, penalty = 'gel')
cvfit2 <- cv.grpregOverlap(X, Y, group, penalty = 'grMCP')
cvfit3 <- cv.grpregOverlap(X, Y, group, penalty = 'grLasso')

## Look at results
summary(cvfit1)
summary(cvfit2)
summary(cvfit3)

## Find variables which are active
details1 <- cvfit1$fit
details2 <- cvfit2$fit
details3 <- cvfit3$fit

opt_lambda1  <- which.min(abs(details1$lambda - cvfit1$lambda.min))
opt_lambda2  <- which.min(abs(details2$lambda - cvfit2$lambda.min))
opt_lambda3  <- which.min(abs(details3$lambda - cvfit3$lambda.min))

latentGamma1 <- which(as.numeric(details1$beta.latent[-1,opt_lambda1])!=0)
latentGamma2 <- which(as.numeric(details2$beta.latent[-1,opt_lambda2])!=0)
latentGamma3 <- which(as.numeric(details3$beta.latent[-1,opt_lambda3])!=0)

## Save important groups and variables
gammahat_gel_final     <- unlist(group)[latentGamma1]
gammahat_grMCP_final   <- unlist(group)[latentGamma2]
gammahat_grLasso_final <- unlist(group)[latentGamma3]

