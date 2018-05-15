####################################################
## Description: Diabetes dataset variable selection
## Author: Levon Demirdjian
## Email: levondem@gmail.com
## Last updated: 05/14/2018
####################################################

library(grpregOverlap)
library(MASS)
library(pscl)
library(Rcpp)

## Import design matrix and response
designMatrix <- read.table('Metabolite_DATA_top1KCV.txt', header = TRUE)
response     <- read.table('GSE25462_EXP.txt', sep = "\t", header = TRUE)

## Save design and response into matrices
X <- t(as.matrix(designMatrix[, 2:51]))
Y <- as.matrix(response$Sensitivity)

## Delete missing data
X <- X[which(!is.na(Y)), ]
Y <- as.matrix(Y[which(!is.na(Y))])

## Define variables
Xvariables  <- as.character(designMatrix[,1])
colnames(X) <- Xvariables

## Import group structure
group_data <- read.table('Metabolite_KEGG_top1KCV.txt', sep = '\t', header = TRUE)
group_data <- group_data[as.character(group_data[,2]) != "", ]

## Eliminate variables that don't belong to a known pathway
valid_vars1 <- (Xvariables %in% as.character(group_data[,1]))
Xvariables  <- Xvariables[valid_vars1]
X           <- X[,valid_vars1]

## The jth element of groupStructure indicates which 
## groups the jth covariate belongs to
uniqueGroups   <- c()
groupStructure <- list()
for(i in 1:nrow(group_data)){
  groupStructure[[i]] <- (unlist(strsplit(as.character(group_data[i,2]), '///')))
  uniqueGroups <- c(uniqueGroups, groupStructure[[i]])
}
uniqueGroups <- unique(uniqueGroups) ## Keep this ordering for groups (i.e. group 1 is hsa03030)
nGroups1     <- length(unique(uniqueGroups))


################################
## Organize overlap structure ##
################################

## "variableIndicator" contains the group structure. 
## The jth element represents the variables in group j.

## Load in variables
variableIndicator        <- lapply(1:nGroups1, function(x) NA)
names(variableIndicator) <- uniqueGroups
for(i in 1:length(Xvariables)){
  curVar    <- Xvariables[i]
  curIndex  <- which(as.character(group_data[,1]) == curVar)
  curGroups <- (unlist(strsplit(as.character(group_data[curIndex,2]), '///')))
  for(j in 1:length(curGroups)){
    locIndex <- which(uniqueGroups == curGroups[j])
    variableIndicator[[locIndex]] <- c(variableIndicator[[locIndex]], curVar)
  }
}

## Delete empty groups
variableIndicator <- lapply(variableIndicator, function(x) x[!is.na(x)])
variableIndicator <- variableIndicator[sapply(variableIndicator, length) != 0]
p       <- as.numeric(sapply(variableIndicator, length))
nGroups <- length(p)


## Remove the very large group
variableIndicator <- variableIndicator[-2]
p       <- p[-2]
nGroups <- length(p)

X       <- X[,colnames(X) %in% unique(unlist(variableIndicator))]
Xvariables <- colnames(X)
  

#######################################
## Try Penalized regression approach ##
#######################################

# Run penalized regression
cvfit1 <- cv.grpregOverlap(X, Y, variableIndicator, penalty = 'gel')
cvfit2 <- cv.grpregOverlap(X, Y, variableIndicator, penalty = 'grMCP')

## Find variables which are active
details1 <- cvfit1$fit
details2 <- cvfit2$fit

opt_lambda1  <- which.min(abs(details1$lambda - cvfit1$lambda.min))
opt_lambda2  <- which.min(abs(details2$lambda - cvfit2$lambda.min))
latentGamma1 <- which(as.numeric(details1$beta.latent[-1,opt_lambda1])!=0)
latentGamma2 <- which(as.numeric(details2$beta.latent[-1,opt_lambda2])!=0)

## Save eta and gamma
gammahat_gel_final   <- unlist(variableIndicator)[latentGamma1]
gammahat_grMCP_final <- unlist(variableIndicator)[latentGamma2]


###########################
## Preliminary functions ##
###########################


my.grid <- function(x){
  as.matrix(expand.grid(rep(list(c(0,1)), x)))
}

## Used to compute the predicted and true values of beta
betaTotal <- function(betaList, variableIndicator, nVariables, nGroups, BurnIn){
  pred <- rep(0, nVariables)
  for(i in 1:nVariables){
    predSum <- 0
    for(j in 1:nGroups){
      if(Xvariables[i] %in% variableIndicator[[j]]){
        curBeta <- betaList[[j]][[which(variableIndicator[[j]] == Xvariables[i])]]
        l <- BurnIn * length(curBeta)
        predSum <- predSum + mean(tail(curBeta, l))
      }
    }
  pred[i] <- predSum
  }
  pred
}

## Read in Rcpp functions
sourceCpp('SampleBasedAlgorithm.cpp')

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~ Perform bi-level variable selection ~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

##############################
## Step 1: Define variables ##
##############################

n          <- nrow(X)
nVariables <- length(Xvariables)

## Organize X into groups
X.tilde <- list()
for(i in 1:nGroups){
  cur_vars <- match(variableIndicator[[i]], Xvariables)
  cur_vars <- cur_vars[!is.na(cur_vars)]
  X.tilde[[i]] <- cbind(X[, cur_vars])
  X.tilde[[i]] <- scale(X.tilde[[i]])
}

####################################
## Step 2: Define hyperparameters ##
####################################

gThresh <- 0.5
vThresh <- 0.5
BurnIn  <- 1/3   
S0      <- 1
S1      <- 25
S       <- 25

outerIter  <- 90    ## Iterations
CWiter     <- 1000  ## For the CW Gibbs sampler
numSim     <- 100   ## For the normalizing constant
burninSim  <- 1/3   ## For normalizing constant using f, not uniform
MRF        <- 0     ## 1 means yes, 0 means 0
MRFweight  <- 1


########################################
## Step 3: Perform bi-level selection ##
########################################


## Initialize parameter lists   
beta.tilde  <- list()
for(i in 1:nGroups){
  beta.tilde[[i]]  <- list()
  for(j in 1:length(variableIndicator[[i]])){
    beta.tilde[[i]][[j]] <- 0
  }
}

B <- list()
for(i in 1:nGroups){
  if(p[i] <= 10){
    B[[i]] <- my.grid(p[i])
  }else{
    B[[i]] <- diag(2)
  }
}

r <- rep(1, length(p))
s <- rep(1, length(p))

params <- list(CWiter = CWiter, outerIter = outerIter, 
               numSim = numSim, MRF = 0, BurnIn = 2/3, 
               S = S, S0 = S0, S1 = S1, 
               MRFweight = MRFweight, B = B)

## Begin analysis
timea <- Sys.time()
timea
res   <- outerLoopCS(X.tilde, Y, beta.tilde, variableIndicator, 
                    params, rbeta, rigamma, rnorm, r, s)
sigma.results <- res$sigmaTildeVec
timeb <- Sys.time()
timeb - timea


## Use only the last (1 - BurnIn)% of samples to handle the burn-in period
etahat.final        <- 1 * (res$etaFinal > 0.5) ## Median rule
gammahat.final_temp <- lapply(1:nGroups, function(x) 1 * (res$gammaFinal[[x]] > 0.5))

## MCSE for model variance
mcseSigma <- sd(tail(sigma.results, round(length(sigma.results) * (1 - BurnIn)))) / 
  sqrt(round(length(sigma.results) * (1 - BurnIn)))

## Find active variables
gammahat.final <- list()
for(i in 1:nGroups){
  active_var          <- variableIndicator[[i]][which(gammahat.final_temp[[i]] == 1)]
  gammahat.final[[i]] <- active_var
}
names(gammahat.final) <- names(variableIndicator)
gammahat.final        <- gammahat.final[etahat.final == 1]

## Compute the beta's
betahat.final     <- list()
for(i in 1:nGroups){
  betahat.final[[i]] <- 0
  cur_vars  <- which(gammahat.final_temp[[i]] == 1)
  if(length(cur_vars) == 0){
    next
  }
  for(j in 1:length(cur_vars)){
    cur_beta  <- beta.tilde[[i]][cur_vars[j]]
    betahat.final[[i]][j] <- mean(tail(cur_beta[[1]], BurnIn * length(cur_beta[[1]])))
  }
}
names(betahat.final) <- names(variableIndicator)
betahat.final        <- betahat.final[etahat.final == 1]

## Eliminate empty groups
finalGroups       <- names(variableIndicator[etahat.final == 1])
finalVariables    <- gammahat.final
finalCoefficients <- betahat.final

## Compute residual vector
predicted.beta <- betaTotal(res$betaNew, variableIndicator, nVariables, nGroups, BurnIn)

## Store results into a matrix
# RESULTS <- matrix(NA, nrow = 2*length(finalGroups), ncol = (max(sapply(finalVariables, function(x) length(x))) + 1))
# for(i in 1:(nrow(RESULTS)/2)){
#   RESULTS[((2*i)-1),1] <- finalGroups[i]
#   RESULTS[((2*i)-1),2:(length(finalVariables[[i]])+1)] <- finalVariables[[i]]
#   RESULTS[2*i,1] <- finalGroups[i]
#   RESULTS[2*i,2:(length(finalVariables[[i]])+1)] <- finalCoefficients[[i]]
# }
# 
# ## Write results to a .txt file
# write.table(RESULTS, file = "Diabetes_Results.txt", row.names=FALSE, 
#             quote = FALSE, na="", col.names=FALSE, sep="\t")

