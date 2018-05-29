############################################
## Description: Simulation from categorical
## model, variable selection
## Author: Levon Demirdjian
## Email: levondem@gmail.com
## Last updated: 05/28/2018
############################################

## Function to generate simulation data according to statistical model
generate_data <- function(n, p, K, lambda = c(1e+2, 1), nu = c(1, 10), alpha = c(1, 1)){

  ## Inputs
  ## n: Sample size
  ## p: Number of variables
  ## K: Number of classes
  ## lambda: First components of IG distributions for sigma_0^2 and sigma_1^2
  ## nu: Second components of IG distributions for sigma_0^2 and sigma_1^2
  ## alpha: Parameter of Beta distribution for rho
  
  ## Sample values for hyperparameters
  sigma0 <- sqrt(1/rgamma(p, lambda[1], nu[1]))
  sigma1 <- sqrt(1/rgamma(p, lambda[2], nu[2]))
  rho    <- rbeta(p, alpha[1], alpha[2])
  
  ## Generate data
  X            <- matrix(rnorm(n * p), nrow = n, ncol = p)
  gamma_matrix <- matrix(0, nrow = p, ncol = 1)
  
  ## Make sure there is at least one active variable
  while(sum(gamma_matrix) == 0){
    for(j in 1:p)
      gamma_matrix[j, 1] <- sample(0:1, 1, 1 - rho[j]) 
  }

  tau_matrix  <- (1 - gamma_matrix) * sigma0 + gamma_matrix * sigma1
  beta_matrix <- as.matrix(sapply(1:(K-1), function(k) rnorm(p, 0, tau_matrix)))
  exp_X_beta  <- cbind(exp(X %*% beta_matrix), 1)
  X_denom     <- rowSums(exp_X_beta)
  
  ## Generate Y
  prob_y <- matrix(0, nrow = n, ncol = K)
  Y      <- rep(0, n)
  for(i in 1:n){
    prob_y[i, ] <- sapply(1:K, function(k) exp_X_beta[,k][i]) / X_denom[i]
    Y[i]        <- sample(1:K, 1, prob = prob_y[i,])
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

## Function to sample from P(gamma | beta, sigma_0, sigma_1, rho, X, Y = k)
sample_gamma <- function(p, gamma_old, beta_old, sigma0_old, sigma1_old, rho_old, K){
  
  ## Inputs
  ## p: Number of variables
  ## gamma_old: px1 matrix of current gamma values
  ## beta_old: p x (K - 1) matrix of current beta values
  ## sigma0: SD of beta components corresponding to gamma = 0
  ## sigma1: SD of beta components corresponding to gamma = 1
  ## rho_old: Current estimate of rho vector
  ## K: Number of classes
  
  for(j in 1:p){
    
    ## Step 1: Pick one gamma to sample; fix the rest
    gamma_1 <- 1
    gamma_0 <- 0
    tau_1   <- (1 - gamma_1) * sigma0_old[j] + gamma_1 * sigma1_old[j]
    tau_0   <- (1 - gamma_0) * sigma0_old[j] + gamma_0 * sigma1_old[j]
    
    ## Step 2: Compute logP for gamma_1
    A1     <- -0.5 * (tau_1^(-2)) * sum(beta_old[j,]^2)
    A2     <- gamma_1 * log(1 - rho_old[j]) + (1 - gamma_1) * log(rho_old[j])
    A3     <- (1 - K) * log(tau_1)
    log_P1 <- A1 + A2 + A3
    
    ## Step 3: Compute logP for gamma_0
    B1     <- -0.5 * (tau_0^(-2)) * sum(beta_old[j,]^2)
    B2     <- gamma_0 * log(1 - rho_old[j]) + (1 - gamma_0) * log(rho_old[j])
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

## Function to sample from P(beta | gamma, sigma_0, sigma_1, rho, X, Y = k)
sample_beta <- function(X, Y, p, beta_old, tau_old, K, eps){
  
  ## Inputs
  ## X: Data matrix
  ## Y: Vector of classes
  ## p: Number of variables
  ## beta_old: p X (K - 1) matrix of current beta values
  ## tau_old: Matrix of tau values using current gamma
  ## K: Number of classes
  ## eps: Langevin step size
  
  n  <- nrow(X)
  dU <- matrix(NA, nrow = p, ncol = K - 1)
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
    A3 <- -sapply(1:(K - 1), function(k) sum(X[Y == k, j]))
    dU[j, ] <- A1 + A2 + A3
    
  }
  
  ## Step 2: Implement Langevin dynamics
  Z         <- matrix(rnorm(p * (K - 1), sd = 1), nrow = p, ncol = K - 1)
  prop_beta <- beta_old - 0.5 * (eps^2) * dU + eps * Z
  
  ## Step 3a: Log MH transition probability, P(prop_beta | beta_old)
  old_to_prop <- sum(dnorm(x = beta_old - 0.5 * (eps^2) * dU, mean = 0, sd = eps, log = TRUE))
  old_to_prop <- 0 ## Comment out

  ## Step 3b: P(prop_beta)
  a1 <- -0.5 * sum(rowSums(prop_beta^2)/(tau_old^2))
  a2 <- sum(sapply(1:n, function(i) (X %*% cbind(prop_beta, 0))[i, Y[i]]))
  a3 <- -sum(log(rowSums(cbind(exp(X %*% prop_beta), 1))))
  log_P_beta_prop <- as.numeric(a1 + a2 + a3)

  ## Step 3c: Log MH transition probability, P(beta_old | prop_beta)
  dU2 <- matrix(NA, nrow = p, ncol = K - 1)
  for(j in 1:p){

    ## Step 1: Compute gradient of log probability
    cur_beta <- prop_beta[j,]

    ## Compute dU
    beta_matrix_cur <- prop_beta
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
    A3 <- -sapply(1:(K - 1), function(k) sum(X[Y == k, j]))
    dU2[j, ] <- A1 + A2 + A3

  }

  prop_to_old <- sum(dnorm(x = prop_beta - 0.5 * (eps^2) * dU2, mean = 0, sd = eps, log = TRUE))
  prop_to_old <- 0 ## Comment out

  ## Step 3d: P(beta_old)
  b1 <- -0.5 * sum(rowSums(beta_old^2)/(tau_old^2))
  b2 <- sum(sapply(1:n, function(i) (X %*% cbind(beta_old, 0))[i, Y[i]]))
  b3 <- -sum(log(rowSums(cbind(exp(X %*% beta_old), 1))))
  log_P_beta_cur <- as.numeric(b1 + b2 + b3)

  ## Step 3: MH acceptance probability
  log_prob <- log_P_beta_prop - log_P_beta_cur + prop_to_old - old_to_prop
  if(log_prob > 0){
    accept <- 1
  }else{
    p_rej  <- 1 / (1 + exp(log_prob))
    p_acc  <- 1 - p_rej
    accept <- sample(1:0, 1, prob = c(p_acc, p_rej))
  }

  if(accept == 1){
    new_beta <- prop_beta
  }else{
    new_beta <- beta_old
  }
  
  ## accept   <- 1
  ## new_beta <- prop_beta
  
  ## Step 4: Update and save this beta
  beta_old   <- new_beta
  accept_vec <- accept
  
  ## Output:
  ## accept: Binary vector indicating if proposal was accepted
  ## beta_new: Matrix of sample beta values
  output   <- list(accept = accept_vec, beta_new = beta_old)
  output
  
}

## Function to sample from P(sigma_0 | gamma, beta, sigma_1, rho, X, Y = k)
sample_sigma <- function(p, lambda, nu, K, gamma_old, beta_old){
  
  ## Inputs
  ## p: Number of variables
  ## lambda: First components of IG distributions for sigma_0^2 and sigma_1^2
  ## nu: Second components of IG distributions for sigma_0^2 and sigma_1^2
  ## K: Number of classes
  ## gamma_old: p x 1 matrix of current gamma values
  ## beta_old: p x (K - 1) matrix of current beta values

  ## First and second component of the posterior for sigma0
  a1 <- lambda[1] + (1 - gamma_old) * (K/2)
  a2 <- nu[1] + 0.5 * (1 - gamma_old) * rowSums(beta_old^2)
  
  ## First and second component of the posterior for sigma0
  b1 <- lambda[2] + gamma_old * (K/2)
  b2 <- nu[2] + 0.5 * gamma_old * rowSums(beta_old^2)
  
  ## Sample sigma0 and sigma1
  sigma0_new <- sqrt(1/rgamma(p, a1, a2))
  sigma1_new <- sqrt(1/rgamma(p, b1, b2))
  
  ## Output:
  ## sigma0_new: Vector of sampled sigma0 values
  ## sigma1_new: Vector of sampled sigma1 values
  output   <- list(sigma0_new = sigma0_new, sigma1_new = sigma1_new)
  output
  
}

## Function to sample from P(rho | gamma, beta, sigma_0, sigma_1, X, Y = k)
sample_rho <- function(alpha, gamma_old){
  
  ## Inputs
  ## alpha: Hyperparameters for the beta distribution
  ## gamma_old: px1 matrix of current gamma values

  p       <- length(gamma_old)
  rho_new <- rbeta(p, alpha[1] + 1 - gamma_old, alpha[2] + gamma_old)
  
  ## Output: 
  ## rho_new: Sampled value for rho
  rho_new 
  
}

## Function to sample from P(beta, gamma, sigma_0, sigma_1, rho | X, Y = k)
sample_posterior <- function(nIter, eps, X, Y, K){

  ## Inputs
  ## nIter: Iterations of sampling
  ## eps: Langevin step size
  ## X: Data matrix
  ## Y: Vector of classes
  ## sigma0: SD of beta component corresponding to gamma = 0
  ## sigma1: SD of beta component corresponding to gamma = 1
  ## rho: Gamma distribution's parameter
  ## K: Number of classes
  
  p <- ncol(X)

  ## Randomly initialize our estimate of beta, gamma, sigma0, sigma1, rho
  gamma_old     <- matrix(sample(0:1, p, TRUE, c(0.5, 0.5)), nrow = p)
  sigma0_old    <- sqrt(1/rgamma(p, 1e+2, 1))
  sigma1_old    <- sqrt(1/rgamma(p, 1, 10))
  rho_old       <- rbeta(p, 1, 1)
  tau_old       <- (1 - gamma_old) * sigma0_old + gamma_old * sigma1_old 
  beta_old      <- as.matrix(sapply(1:(K-1), function(k) rnorm(p, 0, tau_old)))
  
  beta_history   <- lapply(1:(K-1), function(k) matrix(NA, nrow = p, ncol = nIter))
  gamma_history  <- matrix(NA, nrow = p, ncol = nIter)
  sigma0_history <- matrix(NA, nrow = p, ncol = nIter)
  sigma1_history <- matrix(NA, nrow = p, ncol = nIter)
  rho_history    <- matrix(NA, nrow = p, ncol = nIter)
  accept_matrix  <- matrix(NA, nrow = 1, ncol = nIter)

  ## Sample from P(beta, gamma | X, Y = k)
  for(iter in 1:nIter){
    
    ## Sample from P(gamma | beta, Y = k) directly
    new_gamma            <- sample_gamma(p, gamma_old, beta_old, sigma0_old, sigma1_old, rho_old, K)
    tau_old              <- (1 - new_gamma) * sigma0_old + new_gamma * sigma1_old
    gamma_old            <- new_gamma
    gamma_history[,iter] <- new_gamma
    
    ## Sample from P(beta | gamma, Y = k) using Langevin dynamics
    beta_res <- sample_beta(X, Y, p, beta_old, tau_old, K, eps)
    accept   <- beta_res$accept
    new_beta <- beta_res$beta_new
    beta_old <- new_beta
    
    ## Save whether or not the Langevin proposal was accepted
    accept_matrix[iter] <- accept

    for(j in 1:p){
      for(k in 1:(K-1)){
        beta_history[[k]][j,iter] <- new_beta[j, k]
      }
    }
    
    ## Sample from P(sigma_0 | gamma, beta) and P(sigma_0 | gamma, beta)
    sampled_sigmas <- sample_sigma(p, lambda = c(1e+2, 1), nu = c(1,10), K, gamma_old, beta_old)
    sigma0_old     <- sampled_sigmas$sigma0_new
    sigma1_old     <- sampled_sigmas$sigma1_new
    
    sigma0_history[,iter] <- sigma0_old
    sigma1_history[,iter] <- sigma1_old
      
    ## Sample from P(rho | gamma)
    rho_old            <- sample_rho(alpha = c(1,1), gamma_old)
    rho_history[,iter] <- rho_old
      
  }
  
  ## Output:
  ## gamma_history: Matrix of sampled gammas
  ## beta_history: List of matrices of sample betas
  ## accept_matrix: Langevin acceptance counts
  output <- list(gamma_history = gamma_history, 
                 beta_history = beta_history, 
                 sigma0_history = sigma0_history, 
                 sigma1_history = sigma1_history,
                 rho_history = rho_history,
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
  MSE <- mean((beta_true - beta_pred)^2)
  
  ## Output
  ## TCR: True classification rate
  ## TPR: True positive rate
  ## FPR: False positive rate
  ## MSE: Mean squared error
  output <- list(TCR = TCR, TPR = TPR, FPR = FPR, MSE = MSE)
  output
  
}

## Function to perform one simulation
simulation <- function(n, p, K, nIter = 1000, eps = 0.01, burn_percent = 0.33, stride_size = 5){
  
  
  ## Inputs
  ## n: Sample size
  ## p: Number of variables
  ## K: Number of classes
  ## nIter: Number of iterations of Gibbs sampling
  ## eps: Langevin step size
  ## burn_percent: Percent of samples to be burned out
  ## stride_size: Stride size for selecting samples
  
  ## Generate some data
  simulated_data <- generate_data(n, p, K)
  
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
  posterior_samples <- sample_posterior(nIter = nIter, eps = eps, X = X, Y = Y, K = K)
  
  ## Save samples
  gamma_history  <- posterior_samples$gamma_history
  beta_history   <- posterior_samples$beta_history
  sigma0_history <- posterior_samples$sigma0_history
  sigma1_history <- posterior_samples$sigma1_history
  rho_history    <- posterior_samples$rho_history
  accept_matrix  <- posterior_samples$accept_matrix
  
  ##############################
  ## Apply burn in and stride ##
  ##############################
  
  ## Burn in percentage
  beta_burned   <- lapply(1:(K-1), function(k) beta_history[[k]][,floor(burn_percent*nIter):nIter])
  gamma_burned  <- gamma_history[,floor(burn_percent*nIter):nIter]
  sigma0_burned <- sigma0_history[,floor(burn_percent*nIter):nIter]
  sigma1_burned <- sigma1_history[,floor(burn_percent*nIter):nIter]
  rho_burned    <- rho_history[,floor(burn_percent*nIter):nIter]
  n_burned      <- ncol(beta_burned[[1]])
  
  ## Apply stride
  beta_strided   <- lapply(1:(K-1), function(k) beta_burned[[k]][,seq(1, n_burned, stride_size)])
  gamma_strided  <- gamma_burned[,seq(1, n_burned, stride_size)]
  sigma0_strided <- sigma0_burned[,seq(1, n_burned, stride_size)]
  sigma1_strided <- sigma1_burned[,seq(1, n_burned, stride_size)]
  rho_strided    <- rho_burned[,seq(1, n_burned, stride_size)]
  
  ## Save final values of parameter estimates
  beta_final   <- sapply(1:(K-1), function(k) rowSums(beta_strided[[k]])/ncol(beta_strided[[k]]))
  gamma_probs  <- rowSums(gamma_strided)/ncol(gamma_strided)
  sigma0_final <- rowMeans(sigma0_strided)
  sigma1_final <- rowMeans(sigma1_strided)
  rho_final    <- rowMeans(rho_strided)
  
  ## Apply median probability criterion
  gamma_pred <- 1 * (gamma_probs >= 0.5)

  ## Save simulation accuracy measures
  accuracy <- sim_accuracy(gamma_true, gamma_pred, beta_true, beta_final)
  TPR      <- accuracy$TPR
  FPR      <- accuracy$FPR
  TCR      <- accuracy$TCR
  MSE      <- accuracy$MSE
  
  ## Print output
  #accuracy
  #mean(accept_matrix)
  
  ## Output
  ## TCR
  ## TPR
  ## FPR
  ## MSE
  output <- list(TCR = TCR, TPR = TPR, FPR = FPR, MSE = MSE)
  output

}


########################
## Perform simulation ##
########################

sim_results <- simulation(n = 100, p = 3, K = 5)
