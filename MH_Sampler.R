############################################
## Description: Simulation from categorical
## model, variable selection
## Author: Levon Demirdjian
## Email: levondem@gmail.com
## Last updated: 05/07/2018
############################################


## Function to generate simulation data according to statistical model
generate_data <- function(n, p, K, sigma0, sigma1, rho){

  ## Inputs
  ## n: Sample size
  ## p: Number of variables
  ## K: Number of classes
  ## sigma0: SD of beta component corresponding to gamma = 0
  ## sigma1: SD of beta component corresponding to gamma = 1
  ## rho: Gamma distribution's parameter
  
  X            <- matrix(rnorm(n * p), nrow = n, ncol = p)
  gamma_matrix <- matrix(0, nrow = p, ncol = 1)
  
  ## Make sure there is at least on active variable
  while(sum(gamma_matrix) == 0)
    gamma_matrix <- matrix(sample(0:1, p , TRUE, c(rho, 1 - rho)), nrow = p)
  
  tau_matrix  <- (sigma0)^(1 - gamma_matrix) * (sigma1)^(gamma_matrix) 
  beta_matrix <- sapply(1:(K-1), function(k) rnorm(p, 0, tau_matrix))
  exp_X_beta  <- cbind(exp(X %*% beta_matrix), 1)
  X_denom     <- rowSums(exp_X_beta)
  
  ## Generate Y
  prob_y <- matrix(0, nrow = n, ncol = K)
  Y      <- rep(0, n)
  for(i in 1:n){
    prob_y[i, ] <- sapply(1:K, function(k) exp_X_beta[,k][i]) / X_denom[i]
    Y[i] <- sample(1:K, 1, prob = prob_y[i,])
  }
  
  output <- list(gamma_true = gamma_matrix, 
                 tau_true = tau_matrix, 
                 beta_true = beta_matrix,
                 sigma0 = sigma0, 
                 sigma1 = sigma1,
                 rho = rho,
                 K = K,
                 X = X,
                 Y = Y)
  output
  
}

## Function to sample from P(gamma | beta, X, Y = k)
sample_gamma <- function(p, gamma_old, beta_old, sigma0, sigma1, rho, K){
  
  ## Inputs
  ## p: Number of variables
  ## gamma_old: px1 matrix of current gamma values
  ## beta_old: (K-1) x p matrix of current beta values
  ## sigma0: SD of beta component corresponding to gamma = 0
  ## sigma1: SD of beta component corresponding to gamma = 1
  ## rho: Gamma distribution's parameter
  ## K: Number of classes
  
  for(j in 1:p){
    
    ## Step 1: Pick one gamma to sample; fix the rest
    gamma_1 <- 1
    gamma_0 <- 0
    tau_1   <- (sigma0)^(1 - gamma_1) * (sigma1)^(gamma_1) 
    tau_0   <- (sigma0)^(1 - gamma_0) * (sigma1)^(gamma_0) 
    
    ## Step 2: Compute logP for gamma_1
    A1     <- -0.5 * (tau_1^(-2)) * sum(beta_old[j,]^2)
    A2     <- gamma_1 * log(1 - rho) + (1 - gamma_1) * log(rho)
    A3     <- (1 - K) * log(tau_1)
    log_P1 <- A1 + A2 + A3
    
    ## Step 3: Compute logP for gamma_0
    B1     <- -0.5 * (tau_0^(-2)) * sum(beta_old[j,]^2)
    B2     <- gamma_0 * log(1 - rho) + (1 - gamma_0) * log(rho)
    B3     <- (1 - K) * log(tau_0)
    log_P0 <- B1 + B2 + B3
    
    ## Step 4: Sample gamma
    p0 <- 1 / (1 + exp(log_P1 - log_P0))
    p1 <- 1 - p0
    
    new_gamma <- sample(1:0, 1, prob = c(p1, p0))
    
    ## Step 5: Update this gamma
    gamma_old[j] <- new_gamma
    
  }
  
  ## Output
  ## gamma_new: Vector of sampled gamma values
  gamma_new <- gamma_old
  gamma_new
  
}

## Function to sample from P(beta | gamma, X, Y = k)
sample_beta <- function(X, Y, p, beta_old, tau_old, K, eps){
  
  ## Inputs
  ## X: Data matrix
  ## Y: Vector of classes
  ## p: Number of variables
  ## beta_old: (K-1) x p matrix of current beta values
  ## tau_old: Matrix of tau values using current gamma
  ## K: Number of classes
  ## eps: Langevin step size
  
  n          <- nrow(X)
  accept_vec <- rep(0, p)
  
  for(j in 1:p){
    
    ## Step 1: Compute gradient of log probability
    cur_beta <- beta_old[j,]

    ## Compute dU
    beta_matrix_cur <- beta_old
    X_beta_cur      <- cbind(exp(X %*% beta_matrix_cur), 1)
    exp_X_beta      <- t(sapply(1:n, function(i) exp(X[i,] %*% beta_matrix_cur)))
    
    if(K == 2)
      exp_X_beta <- t(exp_X_beta)
    
    num <- t(sapply(1:n, function(i) X[i,j] * exp_X_beta[i,]))
    
    if(K == 2)
      num <- t(num)
    
    den <- rowSums(X_beta_cur)
    
    A1 <- cur_beta/(tau_old[j, 1]^2)
    A2 <- colSums(num/den)
    A3 <- -sapply(1:(K-1), function(k) sum(X[Y == k, j]))
    dU <- A1 + A2 + A3
    
    ## Step 2: Implement Langevin dynamics
    Z         <- rnorm(K-1, sd = 1)
    prop_beta <- cur_beta - 0.5 * (eps^2) * dU + eps * Z 
    
    ## Step 3a: MH acceptance probability, P(prop_beta)
    prop_beta_full     <- beta_old
    prop_beta_full[j,] <- prop_beta
    beta_matrix_prop   <- prop_beta_full
    
    a1 <- -sum(prop_beta^2) / (2 * tau_old[j, 1]^2) 
    a2 <- X[,j] %*% c(prop_beta, 0)[Y]
    a3 <- -sum(log(rowSums(cbind(exp(X %*% beta_matrix_prop), 1))))
    log_P_beta_prop <- as.numeric(a1 + a2 + a3)
    
    ## Step 3b: MH acceptance probability, P(cur_beta)
    b1 <- -sum(cur_beta^2) / (2 * tau_old[j, 1]^2) 
    b2 <- X[,j] %*% c(cur_beta, 0)[Y]
    #b3 <- -sum(log(rowSums(X_beta_cur)))
    b3 <- -sum(log(rowSums(cbind(exp(X %*% beta_old), 1))))
    log_P_beta_cur <- as.numeric(b1 + b2 + b3)
    
    ## Step 3: MH acceptance probability
    if(log_P_beta_prop > log_P_beta_cur){
      accept <- 1
    }else{
      p_rej  <- 1 / (1 + exp(log_P_beta_prop - log_P_beta_cur))
      p_acc  <- 1 - p_rej
      accept <- sample(1:0, 1, prob = c(p_acc, p_rej))
    }
    
    if(accept == 1){
      new_beta <- prop_beta
    }else{
      new_beta <- cur_beta
    }

    ######################
    ## COMMENT THIS OUT ##
    #new_beta <- prop_beta
    ######################
    
    ## Step 4: Update and save this beta
    beta_old[j, ] <- new_beta
    accept_vec[j] <- accept
    
    #for(k in 1:(K-1))
      #beta_history[[k]][j,iter] <- new_beta[k]
    
  }
  
  ## Output:
  ## accept: Binary vector indicating if proposal was accepted
  ## beta_new: Matrix of sample beta values
  beta_new <- beta_old
  output   <- list(accept = accept_vec, beta_new = beta_new)
  output
  
}

## Function to sample from P(beta, gamma | X, Y = k)
sample_posterior <- function(nIter, eps, X, Y, tau_true, sigma0, sigma1, rho, K){

  ## Inputs
  ## nIter: Iterations of sampling
  ## eps: Langevin step size
  ## X: Data matrix
  ## Y: Vector of classes
  ## tau_true: True tau matrix
  ## sigma0: SD of beta component corresponding to gamma = 0
  ## sigma1: SD of beta component corresponding to gamma = 1
  ## rho: Gamma distribution's parameter
  ## K: Number of classes
  
  p <- ncol(X)

  ## Randomly initialize our estimate of beta and gamma
  gamma_old     <- matrix(sample(0:1, p, TRUE, c(0.5, 0.5)), nrow = p)
  tau_old       <- (sigma0)^(1 - gamma_old) * (sigma1)^(gamma_old) 
  beta_old      <- sapply(1:(K-1), function(k) rnorm(p, 0, tau_old))
  beta_history  <- lapply(1:(K-1), function(k) matrix(NA, nrow = p, ncol = nIter))
  gamma_history <- matrix(NA, nrow = p, ncol = nIter)
  accept_matrix <- matrix(NA, nrow = p, ncol = nIter)

  ## Sample from P(beta, gamma | X, Y = k)
  for(iter in 1:nIter){
    
    ## Sample from P(gamma | beta, Y = k) directly
    new_gamma            <- sample_gamma(p, gamma_old, beta_old, sigma0, sigma1, rho, K)
    tau_old              <- (sigma0)^(1 - new_gamma) * (sigma1)^(new_gamma) 
    gamma_old            <- new_gamma
    gamma_history[,iter] <- new_gamma
    
    ## Sample from P(beta | gamma, Y = k) using Langevin dynamics
    beta_res <- sample_beta(X, Y, p, beta_old, tau_old, K, eps)
    accept   <- beta_res$accept
    new_beta <- beta_res$beta_new
    beta_old <- new_beta
    
    ## Save whether or not the Langevin proposal was accepted
    accept_matrix[,iter] <- accept

    for(j in 1:p){
      for(k in 1:(K-1)){
        beta_history[[k]][j,iter] <- new_beta[j, k]
      }
    }
    
  }
  
  ## Output:
  ## gamma_history: Matrix of sampled gammas
  ## beta_history: List of matrices of sample betas
  ## accept_matrix: Langevin acceptance counts
  output <- list(gamma_history = gamma_history, beta_history = beta_history, 
                 accept_matrix = accept_matrix)
  output

}

## Function to compute TCR, TPR, FPR, MSE
sim_accuracy <- function(gamma_true, gamma_pred, beta_true, beta_pred){
  
  ## Input
  ## gamma_true: True gamma
  ## gamma_pred: Predicted gamma
  ## beta_true: True beta
  ## beta_pred: Predicted beta
  
  ## Compute TPR, FPR, TCR
  TPR <- sum(gamma_true == 1 & gamma_pred == 1) / sum(gamma_true == 1)
  FPR <- sum(gamma_true == 0 & gamma_pred == 1) / sum(gamma_true == 0)
  TCR <- mean(gamma_true == gamma_pred)
  
  ## Compute MSE
  MSE <- mean((beta_true = beta_pred)^2)
  
  ## Output
  ## TCR: True classification rate
  ## TPR: True positive rate
  ## FPR: False positive rate
  ## MSE: Mean squared error
  output <- list(TCR = TCR, TPR = TPR, FPR = FPR, MSE = MSE)
  output
  
}

## Function to perform one simulation
simulation <- function(n, p, K, sigma0 = 0.01, sigma1 = 1, rho = 0.5, nIter = 1000, 
                       eps = 0.01, burn_percent = 0.33, stride_size = 5){
  
  
  ## Inputs
  ## n: Sample size
  ## p: Number of variables
  ## K: Number of classes
  ## sigma0: SD of beta component corresponding to gamma = 0
  ## sigma1: SD of beta component corresponding to gamma = 1
  ## rho: Gamma distribution's parameter
  ## nIter: Number of iterations of Gibbs sampling
  ## eps: Langevin step size
  ## burn_percent: Percent of samples to be burned out
  ## stride_size: Stride size for selecting samples
  
  ## Generate some data
  simulated_data <- generate_data(n, p, K, sigma0, sigma1, rho)
  
  ## Save each variable separately
  X          <- simulated_data$X
  Y          <- simulated_data$Y
  tau_true   <- simulated_data$tau_true
  gamma_true <- simulated_data$gamma_true
  beta_true  <- simulated_data$beta_true
  sigma0     <- simulated_data$sigma0
  sigma1     <- simulated_data$sigma1
  rho        <- simulated_data$rho
  K          <- simulated_data$K

  ## Sample from the posterior
  posterior_samples <- sample_posterior(nIter = nIter, eps = eps, X = X, Y = Y, 
                                        tau_true = tau_true, sigma0 = sigma0, sigma1 = sigma1, rho = rho, K = K)
  
  ## Save samples
  gamma_history <- posterior_samples$gamma_history
  beta_history  <- posterior_samples$beta_history
  accept_matrix <- posterior_samples$accept_matrix
  
  ##############################
  ## Apply burn in and stride ##
  ##############################
  
  ## Burn in percentage
  beta_burned  <- lapply(1:(K-1), function(k) beta_history[[k]][,floor(burn_percent*nIter):nIter])
  gamma_burned <- gamma_history[,floor(burn_percent*nIter):nIter]
  n_burned     <- ncol(beta_burned[[1]])
  
  ## Apply stride
  beta_strided  <- lapply(1:(K-1), function(k) beta_burned[[k]][,seq(1, n_burned, stride_size)])
  gamma_strided <- gamma_burned[,seq(1, n_burned, stride_size)]
  
  ## Save final values of parameter estimates
  beta_final  <- sapply(1:(K-1), function(k) rowSums(beta_strided[[k]])/ncol(beta_strided[[k]]))
  gamma_probs <- rowSums(gamma_strided)/ncol(gamma_strided)
  
  ## Apply median probability criterion
  gamma_pred <- 1 * (gamma_probs >= 0.5)

  ## Save simulation accuracy measures
  accuracy <- sim_accuracy(gamma_true, gamma_pred, beta_true, beta_final)
  TPR      <- accuracy$TPR
  FPR      <- accuracy$FPR
  TCR      <- accuracy$TCR
  MSE      <- accuracy$MSE
  
  ## Output
  ## TCR
  ## TPR
  ## FPR
  ## MSE
  output <- list(TCR = TCR, TPR = TPR, FPR = FPR, MSE = MSE)
  output

}


###################################
## Perform simulation many times ##
###################################

sample_sizes  <- c(10, 50, 100)
accuracy_list <- list()
numSim        <- 10
for(n_ind in sample_sizes){
  
  ## Print iteration number to screen
  cat('On n = ', n_ind, '...\n', sep = '')
  
  accuracy           <- matrix(0, nrow = numSim, ncol = 4)
  colnames(accuracy) <- c('TCR', 'TPR', 'FPR', 'MSE')
  for(i in 1:numSim){

    ## Perform simulation
    sim_res <- simulation(n = n_ind, p = 10, K = 3, eps = 10^(-2), nIter = 1000)
    
    ## Save results
    accuracy[i, ] <- unlist(sim_res)
  
  }

  ## Save results to list
  accuracy_list[[which(sample_sizes %in% n_ind)]] <- accuracy

}
names(accuracy_list) <- paste('n =', sample_sizes)


## Mean accuracies across all simulations for each sample size
mean_accuracy <- t(sapply(1:length(accuracy_list), function(k) colMeans(accuracy_list[[k]])))
rownames(mean_accuracy) <- paste('n=', sample_sizes, ':', sep = '')
