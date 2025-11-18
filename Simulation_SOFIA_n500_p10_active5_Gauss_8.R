



source("SOFIA_main.R")

#############################################################################

compute_snr <- function(Y, sigma_squared) {
  signal_variance <- var(Y)  # Directly use the variance of Y
  snr <- signal_variance / sigma_squared
  return(snr)
}


compute_noise_variance_from_snr <- function(Y, snr) {
  signal_variance <- var(Y)  # Compute variance of Y
  noise_variance <- signal_variance / snr  # Compute noise variance
  return(noise_variance)
}


###########################################################

###########################################################

#Function for simulation

simfx2_train <- function(n, p, tps, varx = rep(0, p), bx = 5, mx = 2*pi) {
  fx <- list()
  fobs <- list()
  means <- list()
  sds <- list()
  
  for (j in 1:p) {
    tmax <- max(tps[[j]]) # Maximum time point for the j-th coefficient function.
    
    fx[[j]] <- matrix(0, n, length(tps[[j]])) # Matrix for the true signal
    fobs[[j]] <- matrix(0, n, length(tps[[j]])) # Matrix for observed values (with noise)
    
    for (i in 1:n) {
      bij <- runif(5, 0, bx) # "bx" is the parameter of $a_{ijr}$ in (8).
      mij <- runif(5, 0, mx) # "mx" is the parameter $m_ijr}$ in (8).
      
      tfx <- function(tp) {
        (sum(bij * sin(tp * (5 - bij) * (2 * pi / tmax)) - mij) + 15) / 100
      }
      fx[[j]][i, ] <- sapply(tps[[j]], tfx) # Apply the function for all time points.
    }
    
    # Compute mean and standard deviation before standardization
    means[[j]] <- colMeans(fx[[j]])
    sds[[j]] <- apply(fx[[j]], 2, sd)
    
    # Standardize the matrix columns
    fx[[j]] <- scale(fx[[j]], center = means[[j]], scale = sds[[j]])
  }
  
  # Add noise to the functional observations.
  for (j in 1:p) {
    for (i in 1:n) {
      # Add noise based on the computed noise variance (varx[j])
      fobs[[j]][i, ] <- fx[[j]][i, ] + rnorm(length(tps[[j]]), 0, sqrt(varx[j])) 
    }
  }
  
  return(list("funx" = fx, "funcs" = fobs, "means" = means, "sds" = sds))
}





#######################################################

# Function for generating test matrices

simfx2_test <- function(n, p, tps, varx = rep(0, p), bx = 5, mx = 2*pi) {
  fx <- list()
  fobs <- list()
  
  for (j in 1:p) {
    tmax <- max(tps[[j]])
    
    fx[[j]] <- matrix(0, n, length(tps[[j]]))
    fobs[[j]] <- matrix(0, n, length(tps[[j]]))
    
    for (i in 1:n) {
      bij <- runif(5, 0, bx)
      mij <- runif(5, 0, mx)
      
      tfx <- function(tp) {
        (sum(bij * sin(tp * (5 - bij) * (2 * pi / tmax)) - mij) + 15) / 100
      }
      fx[[j]][i, ] <- sapply(tps[[j]], tfx)
    }
  }
  
  return(fx) # No standardization applied here!
}


#######################################################################

# Function to standardize test matrices using stored mean and std

standardize_test_matrices <- function(test_matrices, means, sds) {
  standardized_test <- list()
  for (j in 1:length(test_matrices)) {
    standardized_test[[j]] <- scale(test_matrices[[j]], center = means[[j]], scale = sds[[j]])
  }
  return(standardized_test)
}

###################################################################


# Prediction Function
predict_scalar_functional <- function(X_test, beta_est, y_bar, T_domain) {
  n <- nrow(X_test[[1]])
  p <- length(X_test)
  M_integ <- sapply(T_domain, function(td) length(td) / diff(range(td)))
  Y_pred <- rep(y_bar, n)
  
  for (i in 1:n) {
    for (j in 1:p) {
      Y_pred[i] <- Y_pred[i] + sum(X_test[[j]][i, ] * beta_est$data[, j]) / M_integ[j]
    }
  }
  
  return(Y_pred)
}


#####################################################



# 19- An Example with all kernels with 5 active predictors out of 10.

start_time <- Sys.time() # to compute the execution time

n_train <- 500
n_test<-100 
p <- 10
#length_out_values <- c(5, 10, 50 , 100, 150, 200)
length_out_values <- c(50)
#kernels <- c("gaussian", "Matern5/2", "Matern3/2", "exponential" )
kernels <- c("gaussian")
#params <- c( 1/2, 1, 2, 4 , 8, 16, 32)
params <- c(32)
snr_levels <- c(0.1, 0.5, 1, 10, 100)


results <- data.frame()
n_reps <- 10

for (length_out in length_out_values) {
  cat("Running analysis for length.out =", length_out, "\n")
  
  tp <- list()
  for (j in 1:p) { tp[[j]] <- seq(0, 1, length.out = length_out) }
  T_domain <- replicate(p, tp[[1]], simplify = FALSE)
  M_integ  <- length(T_domain[[1]]) /  (max(range(T_domain[[1]])) - min(range(T_domain[[1]])))  
  
  set.seed(123)
  bt <- dgamma(30 * tp[[1]], 3, 1/3)
  
  btrue <- list()
  for (j in 1:p) { btrue[[j]] <- 0 }
  
  btrue[[1]] <- (20 - 1 * 1) * bt
  btrue[[2]] <- (20 - 1 * 3) * bt
  btrue[[3]] <- (20 - 1 * 4) * bt
  btrue[[4]] <- 3.5 * dexp(0.02 * 300 * tp[[1]])
  btrue[[5]] <- 2.5 * dexp(0.01 * 300 * tp[[1]])
  
  btrue_matrix <- do.call(rbind, btrue)
  
  
  
  all_beta_values <- list()
  
  #signal_variance <- compute_signal_variance(btrue, T_domain)
  
  for (kernel in kernels) {
    for (param in params) {
      for (snr in snr_levels) {
        
        model_sizes <- numeric(n_reps)
        true_positives <- numeric(n_reps)
        false_positives <- numeric(n_reps)
        true_model_percentage <- numeric(n_reps)
        selection_counts <- matrix(0, nrow = p, ncol = n_reps)
        mse_y <- numeric(n_reps)
        beta_values <- vector("list", n_reps)
        RMSE_test_values <- numeric(n_reps)  
        
        
        
        key <- paste(length_out, kernel, param, snr, sep = "_")
        all_beta_values[[key]] <- vector("list", n_reps)
        
        for (rep in 1:n_reps) {
          set.seed(123 + rep )
          
          XX_train <- simfx2_train(n = n_train, p = p, tps = tp)
          X_train <- lapply(1:p, function(i) list(time_domain = T_domain[[i]], data = XX_train$funcs[[i]]))
          
          # Build standardized-space true betas on the time grid (p x T)
          btrue_std_matrix <- btrue_matrix
          for (j in 1:p) {
            btrue_std_matrix[j, ] <- btrue_matrix[j, ] * as.numeric(XX_train$sds[[j]])
          }
          
          
          Y_train <- numeric(n_train)
          for (i in 1:n_train) {
            for (j in 1:p) {
              Y_train[i] <- Y_train[i] + btrue_matrix[j, ] %*% XX_train$funcs[[j]][i, ] / M_integ
            }
          }
          sigma2 <- compute_noise_variance_from_snr(Y_train, snr)
          y_train <- Y_train + rnorm(n_train, mean = 0, sd = sqrt(sigma2))
          y_bar <- mean(y_train)
          y_train <- as.matrix(y_train) - mean(y_train)
          
          sim_5_active <- try(
            SOFIA(
              X = X_train,
              Y = y_train,
              type_kernel = kernel,
              param_kernel = param, 
              thres_eigen = 0.99,
              period_kernel = NULL,
              NoI = 10,
              thres_CD = 0.05,
              number_non_zeros = 10,
              ratio_lambda = 0.01,
              number_lambda = 100,
              proportion_training_set = 0.75,
              n_fold = 5,
              verbose = FALSE
            ),
            silent = TRUE
          )
          
          # Skip to next rep if there's an error
          if (inherits(sim_5_active, "try-error")) {
            all_beta_values[[key]][[rep]] <- list(error = TRUE)
            next  # skip this rep safely
          }
          
          
          selected_vars <- sim_5_active$predictors_adaptive
          
          
          model_sizes[rep] <- length(selected_vars)
          
          TP <- sum(selected_vars %in% 1:5)
          FP <- sum(selected_vars %in% 6:p)
          true_positives[rep] <- TP
          false_positives[rep] <- FP
          true_model_percentage[rep] <- ifelse(TP == 5 && FP == 0, 1, 0)
          selection_counts[selected_vars, rep] <- 1
          mse_y[rep] <- ifelse(is.null(sim_5_active$MSEY_adaptive), NA, sim_5_active$MSEY_adaptive)
          
          
          beta_values[[rep]] <- list(
            time_domain = T_domain[[1]],
            data = sim_5_active$beta_std_on_time_grid  # <-- standardized-space beta curves (T x p)
          )
          
          
          
          if (inherits(sim_5_active, "try-error")) {
            all_beta_values[[key]][[rep]] <- list(error = TRUE)
          } else {
            all_beta_values[[key]][[rep]] <- sim_5_active
          } 
          
          
          XX_test <- simfx2_test(n_test, p, tp)
          X_test <- standardize_test_matrices(XX_test, XX_train$means, XX_train$sds)
          
          y_pred <- predict_scalar_functional(X_test, beta_values[[rep]], y_bar, T_domain)
          y_true <- numeric(n_test)
          for (i in 1:n_test) {
            y_true[i] <- 0 
            for (j in 1:p) {
              y_true[i] <- y_true[i] + btrue_std_matrix[j, ] %*% X_test[[j]][i, ] / M_integ
            }
          }
          
          y_true <- y_true + rnorm(n_test, mean = 0, sd = sqrt(sigma2))
          #RMSE <-sqrt(mean((y_pred - y_true)^2))
          RMSE_test_values[rep] <- sqrt(mean((y_pred - y_true)^2))  # ✅ Store RMSE for each repetition
          
        }
        
        avg_model_size <- mean(model_sizes)
        avg_true_positive <- mean(true_positives)
        avg_false_positive <- mean(false_positives)
        percent_true_model <- mean(true_model_percentage) * 100
        selection_frequencies <- rowMeans(selection_counts) * 100
        avg_mse_y <- mean(mse_y)
        avg_RMSE_test <- mean(RMSE_test_values)
        
        result_row <- data.frame(
          Length_Out = length_out,
          Kernel = kernel,
          Param = param,
          SNR = snr,
          Sigma = sqrt(sigma2),
          Avg_Model_Size = avg_model_size,
          Avg_True_Positive = avg_true_positive,
          Avg_False_Positive = avg_false_positive,
          Percent_True_Model = percent_true_model,
          Avg_MSE_Y = avg_mse_y,
          #beta_values = beta_values,
          RMSE_test = avg_RMSE_test
        )
        
        for (i in 1:p) {
          result_row[paste0("Predictor_", i, "_Pct_Selected")] <- selection_frequencies[i]
        }
        
        results <- rbind(results, result_row)
      }
    }
  }
}

print(results)

write.csv(results, "SOFIA_n500_p10_active5_gauss_new.csv", row.names = FALSE)

saveRDS(all_beta_values, "SOFIA_n500_p10_active5_gauss_8.rds")

end_time <- Sys.time()

execution_time <- as.numeric(difftime(end_time, start_time, units = "mins"))

print(paste("Total execution time (minutes):", execution_time))



##################################################
################################################### 




