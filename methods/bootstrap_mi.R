## Code to compute the sandwich smooth mean estimator 
## and the pooled covariance estimator.

## Inputs: event_complete.RDS and nonevent_complete.RDS 
## Inputs: These are dataframes containing the relevent EDA
## Outputs: RDS files for event and nonevent dataframes:
## Outputs: (1) mean functions mu(t,s), 
## Outputs: (2) coefficients on the residual X(t,s) - mu(t,s)
## Outputs: (3) eigenvectors of the pooled covariance estimate 
## Outputs: (4) complete-case timestamps

## Note: one can check if reconstruction is good by taking
## hatX(t,s) = coefficients%*%eigenvectors + mu(t,s) 
## And seeing how close hatX(t,s) is to X(t,s).

args=(commandArgs(TRUE))

print(args)

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  current_bootstrap = 1 # CONVERT TO ACCEPT THIS 
}else{
  current_bootstrap = as.numeric(args[[1]])
}

print(paste("Current bootstrap:", current_bootstrap))

# num_bootstraps = 200
num_imputes = 2
library(refund); library(lubridate); library(dynr); library(Matrix)
# set.seed(1391307)
  
## WINDOWS
setwd("Z:/SI_data/")
## LINUX
# setwd('/mnt/turbo/SI_data/')

# set.seed(123193871)
set_of_types = c("eda", "acc")
for (type in set_of_types){
  event_complete = readRDS(paste(type, "_event_complete_HLP_2021-11-20.RDS", sep = ""))
  event_means = readRDS(file = paste(type, "_lin_event_means_2021-12-20.RDS", sep=""))
  event_times = readRDS(file = paste(type, "_lin_event_complete_case_timesincebaseline_2021-12-20.RDS", sep=""))
  event_Sigma = readRDS(file = paste(type, "_event_Sigma_2022-01-26.RDS", sep=""))
  sequence = seq(0,-30, length.out = ncol(event_complete) - 3); sensor_obs = 4:ncol(event_complete)
  bootstrap_rows = sample(1:nrow(event_complete), size = nrow(event_complete), replace = TRUE)
  bootstrap_event_complete = event_complete[bootstrap_rows,]
  impute_Sigma = event_Sigma
  eigen_Sigma = eigen(impute_Sigma)
  eigen_values = eigen_Sigma$values
  eigen_values[eigen_values < 0] = 0
  impute_Sigma = eigen_Sigma$vectors%*%diag(eigen_Sigma$values)%*%t(eigen_Sigma$vectors)
  ## 
  ## Mean: \mu_1 + Sigma_{12} (\Sigma_22^plus) (y - mu_2)
  ## COVARIANCE: \Sigma_11 - \Sigma{12} \sigma_{22}^plus \Sigma_{12}
  ## \Sigma_22 =T' T 
  ## \Sigma_22^\plus = T'(T' T)^{-2} T
  full_obs = apply(X = bootstrap_event_complete, MARGIN = 1, FUN = function(x){!any(is.na(x))})
  how_bad = apply(X = bootstrap_event_complete, MARGIN = 1, FUN = function(x){mean(is.na(x))})
  
  ## Remove the inflation factor
  inflation = 100
  impute_Sigma = event_Sigma/inflation^2
  
  for(current_mi in 1:num_imputes) {
    print(paste("Number of rows to impute:", length(which(!full_obs))))
    for(row in which(!full_obs)) {
      print(paste("Current rows:", which(row == which(!full_obs))))
      current_row = bootstrap_event_complete[row,]
      current_time = current_row[3]
      data = current_row[,4:length(current_row)]
      ## FIND MOST CLOSE MEAN FUNCTION
      distance = abs(current_time - event_times)
      current_mean = event_means[which(distance == min(distance))[1],]
      ## EXTRACT CONDITIONAL VARIANCE
      missing_spots = is.na(data)
      observed_data = data[!missing_spots]
      print(paste("Number of missing entries:", sum(!missing_spots)))
      Sigma_11 = as.matrix(impute_Sigma[missing_spots, missing_spots], nrow = sum(missing_spots))
      Sigma_12 = as.matrix(impute_Sigma[missing_spots, !missing_spots], nrow = sum(missing_spots))
      if(all(is.na(data))) {
        pseudo_inverse = 0
      } else {
        Sigma_22 = as.matrix(impute_Sigma[!missing_spots, !missing_spots])
        mu_1 = current_mean[missing_spots]
        mu_2 = current_mean[!missing_spots]
        eig_Sigma = eigen(Sigma_22)
        eig_vectors = eig_Sigma$vectors
        eig_values = eig_Sigma$values
        eig_values[eig_values < 0 ] = 0
        max_K = min(which(cumsum(eig_values)/sum(eig_values) > 0.9999))
        T = eig_vectors[,1:max_K]%*%diag(sqrt(eig_values[1:max_K]))
        pseudo_inverse = T%*%solve(t(T)%*%T) %*% solve(t(T)%*%T, t(T))
      }
      if(sum(missing_spots) > 1) {
        mu_conditional = mu_1 + Sigma_12%*%pseudo_inverse %*% (observed_data - mu_2)
        Sigma_conditional = Sigma_11 - Sigma_12%*%pseudo_inverse%*%t(Sigma_12)
        Sigma_conditional = (Sigma_conditional + t(Sigma_conditional))/2
        eig_Sigma = eigen(Sigma_conditional)
        eig_vectors = eig_Sigma$vectors
        eig_values = eig_Sigma$values
        eig_values[eig_values < 0 ] = 0
        if(sum(eig_values) == 0) {
          T = matrix(0, nrow = nrow(mu_conditional), ncol = ncol(eig_vectors))
        } else{
          max_K = min(which(cumsum(eig_values)/sum(eig_values) > 0.9999))
          if(max_K > 1) {
            T = eig_vectors[,1:max_K]%*%diag(sqrt(eig_values[1:max_K]))
          } else{
            T = eig_vectors[,1:max_K]%*%as.matrix(sqrt(eig_values[1:max_K]))
          }
        }
        imputed_obs = mu_conditional + T%*%rnorm(ncol(T))
      } else {
        mu_conditional = mu_1 + t(Sigma_12)%*%pseudo_inverse%*%(observed_data - mu_2)
        Sigma_conditional = Sigma_11 - t(Sigma_12)%*%pseudo_inverse%*%Sigma_12
        eig_Sigma = eigen(Sigma_conditional)
        eig_vectors = eig_Sigma$vectors
        eig_values = eig_Sigma$values
        eig_values[eig_values < 0 ] = 0
        if(sum(eig_values) == 0) {
          T = matrix(0, nrow = nrow(mu_conditional), ncol = ncol(eig_vectors))
        } else{
          max_K = min(which(cumsum(eig_values)/sum(eig_values) > 0.9999))
          if(max_K > 1) {
            T = eig_vectors[,1:max_K]%*%diag(sqrt(eig_values[1:max_K]))
          } else{
            T = eig_vectors[,1:max_K]%*%as.matrix(sqrt(eig_values[1:max_K]))
          }
        }
        imputed_obs = mu_conditional + T%*%rnorm(ncol(T))
      }
      data[missing_spots] = imputed_obs
      bootstrap_event_complete[row,4:length(current_row)] = data
    }
    
    full_obs = apply(X = bootstrap_event_complete, MARGIN = 1, FUN = function(x){!any(is.na(x))})
    Y = bootstrap_event_complete[full_obs,sensor_obs]
    Y = as.matrix(Y, nrow = nrow(Y), ncol = ncol(Y))
    x = as.numeric(bootstrap_event_complete[full_obs,3])
    z = sequence
    
    est <- fbps(Y,list(x=x,z=z))
    residual = inflation*(Y-est$Yhat)
    ## TAKE PAIRS OF ROWS and CALCULATE THE MSE
    print(paste("Made it to event, linear Sigma calc for", type, "under imputation #", current_mi))
    eig_Sigma_est = eigen(event_Sigma)
    eig_vectors = eig_Sigma_est$vectors
    eig_values = eig_Sigma_est$values
    eig_values[eig_values< 0] = 0
    K = 35
    print(paste("Variance explained: ", (cumsum(eig_values)/sum(eig_values))[K]))
    phi_vectors = eig_vectors[,1:K]
    coef = residual%*%phi_vectors
    ## To show the impact of K
    obs = 20
    optK = 10
    optphi_vectors = eig_vectors[,1:optK]
    optcoef = residual%*%optphi_vectors
    ## Timestamps
    timestamp = bootstrap_event_complete[full_obs,2]
    
    ## GENERATE SPLINES
    if (type == "eda") {
      sequence <- seq(-30,0, by = 1/60)
    } else if (type == "acc") {
      sequence <- seq(-30,0, by = 1/6)
    }
    K_b = 35
    num=K_b-3
    qtiles <- seq(0, 1, length = num + 2)[-c(1, num + 2)]
    knots <- quantile(sequence, qtiles)
    bb = cbind(1, sequence, sequence^2, sapply(knots, function(k) ((sequence - k > 0) * (sequence - k)) ^ 2))
    print("Generated Splines")
    
    ## GENERATE MODEL MATRIX
    print(paste("Building", type, "event model matrix"))
    event_eigen_vectors = phi_vectors
    event_coef = coef/inflation
    event_means = est$Yhat
    event_timesincebaseline = x
    event_timestamp = timestamp
    event_id = as.numeric(event_complete[full_obs,1])
    event_J_coef = event_coef%*%t(event_eigen_vectors)%*%bb
    event_J_means = event_means%*%bb
    event_model.matrix = cbind(event_id, event_timesincebaseline, hour(event_timestamp), event_J_coef + event_J_means)
    event_model.matrix = data.frame(event_model.matrix)
    names(event_model.matrix) = c("id", "timesincebaseline", "hour", paste("acc_X", 1:ncol(event_J_coef), sep = ""))
    
    saveRDS(event_model.matrix, file = paste("bootstrap_files/",type, "_event_modelmatrix_bootstrap_",current_bootstrap,"_mi_",current_mi,"_",today(),".RDS", sep = ""))
    
  }
  
}

## CLEAR WORKSPACE
rm(list=c("eig_Sigma", "eig_Sigma_est", "eig_vectors", "eigen_Sigma", "est", 
           "event_coef", "event_complete", "event_eigen_vectors", "event_J_coef",
           "event_means", "event_J_means", "event_model.matrix", "event_Sigma", "impute_Sigma",
           "pseudo_inverse", "residual", "Sigma_11", "Sigma_12", "Sigma_22", "Y"))
current_bootstrap = 1 # CONVERT TO ACCEPT THIS 
num_imputes = 2
library(refund); library(lubridate)

## WINDOWS
setwd("Z:/SI_data/")
## Linux
# setwd('/mnt/turbo/SI_data/')

set_of_types = c("eda", "acc")
for (type in set_of_types){
  nonevent_complete = readRDS(paste(type, "_nonevent_complete_HLP_2021-11-20.RDS", sep=""))
  nonevent_means = readRDS(file = paste(type, "_lin_nonevent_means_2021-12-20.RDS", sep=""))
  nonevent_times = readRDS(file = paste(type, "_lin_nonevent_complete_case_timesincebaseline_2021-12-20.RDS", sep=""))
  nonevent_Sigma = readRDS(file = paste(type, "_nonevent_Sigma_2022-01-26.RDS", sep=""))
  sequence = seq(0,-30,length.out = ncol(nonevent_complete) - 3); sensor_obs = 4:ncol(nonevent_complete)
  # plot(sequence, nonevent_complete[4,sensor_obs])
  
  bootstrap_rows = sample(1:nrow(nonevent_complete), size = nrow(nonevent_complete), replace = TRUE)
  bootstrap_nonevent_complete = nonevent_complete[bootstrap_rows,]
  impute_Sigma = nonevent_Sigma
  eigen_Sigma = eigen(impute_Sigma)
  eigen_values = eigen_Sigma$values
  eigen_values[eigen_values < 0] = 0
  impute_Sigma = eigen_Sigma$vectors%*%diag(eigen_Sigma$values)%*%t(eigen_Sigma$vectors)
  ## 
  ## Mean: \mu_1 + Sigma_{12} (\Sigma_22^plus) (y - mu_2)
  ## COVARIANCE: \Sigma_11 - \Sigma{12} \sigma_{22}^plus \Sigma_{12}
  ## \Sigma_22 =T' T 
  ## \Sigma_22^\plus = T'(T' T)^{-2} T
  full_obs = apply(X = bootstrap_nonevent_complete, MARGIN = 1, FUN = function(x){!any(is.na(x))})
  how_bad = apply(X = bootstrap_nonevent_complete, MARGIN = 1, FUN = function(x){mean(is.na(x))})
  
  ## Remove the inflation factor
  inflation = 100
  impute_Sigma = nonevent_Sigma/inflation^2
  
  for(current_mi in 1:num_imputes) {
    for(row in which(!full_obs)) {
      current_row = bootstrap_nonevent_complete[row,]
      current_time = current_row[3]
      data = current_row[,4:length(current_row)]
      ## FIND MOST CLOSE MEAN FUNCTION
      distance = abs(current_time - nonevent_times)
      current_mean = nonevent_means[which(distance == min(distance))[1],]
      ## EXTRACT CONDITIONAL VARIANCE
      missing_spots = is.na(data)
      observed_data = data[!missing_spots]
      Sigma_11 = as.matrix(impute_Sigma[missing_spots, missing_spots], nrow = sum(missing_spots))
      Sigma_12 = as.matrix(impute_Sigma[missing_spots, !missing_spots], nrow = sum(missing_spots))
      Sigma_22 = as.matrix(impute_Sigma[!missing_spots, !missing_spots])
      mu_1 = current_mean[missing_spots]
      mu_2 = current_mean[!missing_spots]
      eig_Sigma = eigen(Sigma_22)
      eig_vectors = eig_Sigma$vectors
      eig_values = eig_Sigma$values
      eig_values[eig_values < 0 ] = 0
      max_K = min(which(cumsum(eig_values)/sum(eig_values) > 0.9999))
      T = eig_vectors[,1:max_K]%*%diag(sqrt(eig_values[1:max_K]))
      pseudo_inverse = T%*%solve(t(T)%*%T) %*% solve(t(T)%*%T, t(T))
      if(sum(missing_spots) > 1) {
        mu_conditional = mu_1 + Sigma_12%*%pseudo_inverse %*% (observed_data - mu_2)
        Sigma_conditional = Sigma_11 - Sigma_12%*%pseudo_inverse%*%t(Sigma_12)
        Sigma_conditional = (Sigma_conditional + t(Sigma_conditional))/2
        eig_Sigma = eigen(Sigma_conditional)
        eig_vectors = eig_Sigma$vectors
        eig_values = eig_Sigma$values
        eig_values[eig_values < 0 ] = 0
        if(sum(eig_values) == 0) {
          T = matrix(0, nrow = nrow(mu_conditional), ncol = ncol(eig_vectors))
        } else{
          max_K = min(which(cumsum(eig_values)/sum(eig_values) > 0.9999))
          if(max_K > 1) {
            T = eig_vectors[,1:max_K]%*%diag(sqrt(eig_values[1:max_K]))
          } else{
            T = eig_vectors[,1:max_K]%*%as.matrix(sqrt(eig_values[1:max_K]))
          }
        }
        imputed_obs = mu_conditional + T%*%rnorm(ncol(T))
      } else {
        mu_conditional = mu_1 + t(Sigma_12)%*%pseudo_inverse%*%(observed_data - mu_2)
        Sigma_conditional = Sigma_11 - t(Sigma_12)%*%pseudo_inverse%*%Sigma_12
        eig_Sigma = eigen(Sigma_conditional)
        eig_vectors = eig_Sigma$vectors
        eig_values = eig_Sigma$values
        eig_values[eig_values < 0 ] = 0
        if(sum(eig_values) == 0) {
          T = matrix(0, nrow = nrow(mu_conditional), ncol = ncol(eig_vectors))
        } else{
          max_K = min(which(cumsum(eig_values)/sum(eig_values) > 0.9999))
          if(max_K > 1) {
            T = eig_vectors[,1:max_K]%*%diag(sqrt(eig_values[1:max_K]))
          } else{
            T = eig_vectors[,1:max_K]%*%as.matrix(sqrt(eig_values[1:max_K]))
          }
        }
        imputed_obs = mu_conditional + T%*%rnorm(ncol(T))
      }
      data[missing_spots] = imputed_obs
      bootstrap_nonevent_complete[row,4:length(current_row)] = data
    }
    
    full_obs = apply(X = bootstrap_nonevent_complete, MARGIN = 1, FUN = function(x){!any(is.na(x))})
    ## print(paste(length(which(!full_obs)), "out of", length(full_obs), "missing at least 1 entry"))
    Y = bootstrap_nonevent_complete[,sensor_obs]
    Y = as.matrix(Y, nrow = nrow(Y), ncol = ncol(Y))
    x = as.numeric(bootstrap_nonevent_complete[,3])
    z = sequence
    id = bootstrap_nonevent_complete[,1]
    ## Timestamps
    timestamp = as_datetime(bootstrap_nonevent_complete[,2])
    nonevent_id = as.numeric(bootstrap_nonevent_complete[full_obs,1])
    rm("nonevent_complete", "nonevent_means") # Remove the full data to free up memory
    
    est <- fbps(Y,list(x=x,z=z))
    ## COMPUTE ESTIMATED MARGINAL COVARIANCE 
    inflation = 100
    residual = inflation*(Y-est$Yhat)
    nonevent_means = est$Yhat
    rm("Y"); rm("est")
    ## TAKE PAIRS OF ROWS and CALCULATE THE MSE
    print(paste("Made it to nonevent, linear Sigma calc for", type))
    Sigma = matrix(nrow = ncol(residual), ncol = ncol(residual))
    for (i in 1:ncol(residual)) {
      for (j in 1:ncol(residual)) {
        Sigma[i,j] = sum(residual[,i]*residual[,j])/nrow(residual)
      }
    }
    Sigma_est <- fbps(Sigma)
    final_Sigma = (Sigma_est$Yhat + t(Sigma_est$Yhat))/2
    eig_Sigma_est = eigen(final_Sigma)
    eig_vectors = eig_Sigma_est$vectors
    eig_values = eig_Sigma_est$values
    eig_values[eig_values<0] = 0
    
    K = 35
    print(paste("Variance explained: ", (cumsum(eig_values)/sum(eig_values))[K]))
    phi_vectors = eig_vectors[,1:K]
    coef = residual%*%phi_vectors
    
    ## GENERATE SPLINES
    if (type == "eda") {
      sequence <- seq(-30,0, by = 1/60)
    } else if (type == "acc") {
      sequence <- seq(-30,0, by = 1/6)
    }
    K_b = 35
    num=K_b-3
    qtiles <- seq(0, 1, length = num + 2)[-c(1, num + 2)]
    knots <- quantile(sequence, qtiles)
    bb = cbind(1, sequence, sequence^2, sapply(knots, function(k) ((sequence - k > 0) * (sequence - k)) ^ 2))
    print("Generated Splines")
    
    ## GENERATE MODEL MATRIX
    print(paste("Building", type, "nonevent model matrix"))
    nonevent_eigen_vectors = phi_vectors
    nonevent_coef = coef/inflation
    nonevent_timesincebaseline = x
    nonevent_timestamp = timestamp
    nonevent_J_coef = nonevent_coef%*%t(nonevent_eigen_vectors)%*%bb
    nonevent_J_means = nonevent_means%*%bb
    nonevent_model.matrix = cbind(nonevent_id, nonevent_timesincebaseline, hour(nonevent_timestamp), nonevent_J_coef + nonevent_J_means)
    nonevent_model.matrix = data.frame(nonevent_model.matrix)
    names(nonevent_model.matrix) = c("id", "timesincebaseline", "hour", paste("acc_X", 1:ncol(nonevent_J_coef), sep = ""))
    
    saveRDS(nonevent_model.matrix, file = paste("bootstrap_files/",type, "_nonevent_modelmatrix_bootstrap_",current_bootstrap,"_mi_",current_mi,"_",today(),".RDS", sep = ""))
    
  }
}
