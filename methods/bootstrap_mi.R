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

num_bootstraps = 200
num_imputes = 2
set.seed(1391307)
  
setwd("Z:/SI_data/")
set_of_types = c("eda", "acc")
for (type in set_of_types){
  event_complete = readRDS(paste(type, "_event_complete_HLP_2021-11-20.RDS", sep = ""))
  event_means = readRDS(file = paste(type, "_lin_event_means_2021-12-20.RDS", sep=""))
  event_times = readRDS(file = paste(type, "_lin_event_complete_case_timesincebaseline_2021-12-20.RDS", sep=""))
  event_Sigma = readRDS(file = paste(type, "_event_Sigma_2022-01-26.RDS", sep=""))
  sequence = seq(0,-30, length.out = ncol(event_complete) - 3); sensor_obs = 4:ncol(event_complete)
  bootstrap_rows = sample(1:nrow(event_complete), size = nrow(event_complete), replace = TRUE)
  bootstrap_event_complete = event_complete[bootstrap_rows,]
  
  
  
  library(refund); library(lubridate)
  full_obs = apply(X = bootstrap_event_complete, MARGIN = 1, FUN = function(x){!any(is.na(x))})
  ## Remove the inflation factor
  inflation = 100
  impute_Sigma = event_Sigma/inflation^2
  for(row in which(!full_obs)) {
    current_row = bootstrap_event_complete[row,]
    current_time = current_row[3]
    data = current_row[,4:length(current_row)]
    ## FIND MOST CLOSE MEAN FUNCTION
    distance = abs(current_time - event_times)
    current_mean = event_means[which(distance == min(distance))[1],]
    ## EXTRACT CONDITIONAL VARIANCE
    missing_spots = is.na(data)
    observed_data = data[!missing_spots]
    Sigma_11 = as.matrix(impute_Sigma[missing_spots, missing_spots], nrow = sum(missing_spots))
    Sigma_12 = as.matrix(impute_Sigma[missing_spots, !missing_spots], nrow = sum(missing_spots))
    Sigma_22 = impute_Sigma[!missing_spots, !missing_spots]
    mu_1 = current_mean[missing_spots]
    mu_2 = current_mean[!missing_spots]
    pseudo_inverse = ginv(Sigma_22)
    mu_conditional = mu_1 + t(Sigma_12)%*%pseudo_inverse%*%(observed_data - mu_2)
    Sigma_conditional = Sigma_11 - t(Sigma_12)%*%pseudo_inverse%*%Sigma_12
  }
  
  
  Y = event_complete[full_obs,sensor_obs]
  Y = as.matrix(Y, nrow = nrow(Y), ncol = ncol(Y))
  x = as.numeric(event_complete[full_obs,3])
  z = sequence
  
  est <- fbps(Y,list(x=x,z=z))
  residual = inflation*(Y-est$Yhat)
  ## TAKE PAIRS OF ROWS and CALCULATE THE MSE
  print(paste("Made it to event, linear Sigma calc for", type))
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
  timestamp = event_complete[full_obs,2]
  ## SAVE THE MEAN, COEFFICIENT MATRIX, AND FIRST K EIGEN VECTORS
  saveRDS(object = est$Yhat, file = paste(type, "_lin_event_means_",today(),".RDS", sep=""))
  saveRDS(object = coef/inflation, file = paste(type, "_lin_event_coef_matrix_",today(),".RDS", sep=""))
  saveRDS(object = phi_vectors, file = paste(type, "_lin_event_eigen_vectors_",today(),".RDS", sep=""))
  saveRDS(object = x, file = paste(type, "_lin_event_complete_case_timesincebaseline_",today(),".RDS", sep=""))
  saveRDS(object = timestamp, file = paste(type, "_lin_event_complete_case_timestamp_",today(),".RDS", sep=""))
  saveRDS(object = as.numeric(event_complete[full_obs,1]), 
          file = paste(type, "_lin_event_complete_case_ids_",today(),".RDS", sep=""))
  
}

## CLEAR WORKSPACE
rm(list=ls()) 
## WINDOWS
setwd("Z:/SI_data/")
## Linux
set_of_types = c("eda", "acc")
for (type in set_of_types){
  nonevent_complete = readRDS(paste(type, "_nonevent_complete_HLP_2021-11-20.RDS", sep=""))
  sequence = seq(0,-30,length.out = ncol(nonevent_complete) - 3); sensor_obs = 4:ncol(nonevent_complete)
  # plot(sequence, nonevent_complete[4,sensor_obs])
  
  library(refund); library(lubridate)
  full_obs = apply(X = nonevent_complete, MARGIN = 1, FUN = function(x){!any(is.na(x))})
  ## print(paste(length(which(!full_obs)), "out of", length(full_obs), "missing at least 1 entry"))
  Y = nonevent_complete[full_obs,sensor_obs]
  Y = as.matrix(Y, nrow = nrow(Y), ncol = ncol(Y))
  mu = rowMeans(Y)
  quad_Y = sweep(Y,1,mu)^2
  x = as.numeric(nonevent_complete[full_obs,3])
  z = sequence
  id = nonevent_complete[full_obs,1]
  ## Timestamps
  timestamp = nonevent_complete[full_obs,2]
  rm("nonevent_complete") # Remove the full data to free up memory
  
  est <- fbps(Y,list(x=x,z=z))
  quad_est <- fbps(quad_Y,list(x=x,z=z))
  
  ## COMPUTE ESTIMATED MARGINAL COVARIANCE 
  inflation = 100
  residual = inflation*(Y-est$Yhat)
  quad_residual = inflation*(quad_Y-quad_est$Yhat)
  saveRDS(object = est$Yhat, file = paste(type, "_lin_nonevent_means_",today(),".RDS", sep=""))
  saveRDS(object = quad_est$Yhat, file = paste(type, "_quad_nonevent_means_",today(),".RDS", sep=""))
  rm("Y"); rm("est"); rm("quad_Y"); rm("quad_est")
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
  saveRDS(object = final_Sigma, file = paste(type, "_nonevent_Sigma_",today(),".RDS", sep=""))
  eig_Sigma_est = eigen(final_Sigma$Yhat)
  eig_vectors = eig_Sigma_est$vectors
  eig_values = eig_Sigma_est$values
  
  eig_values[eig_values<0] = 0
  
  K = 35
  print(paste("Variance explained: ", (cumsum(eig_values)/sum(eig_values))[K]))
  phi_vectors = eig_vectors[,1:K]
  coef = residual%*%phi_vectors
  
  
  ## SAVE THE MEAN, COEFFICIENT MATRIX, AND FIRST K EIGEN VECTORS
  saveRDS(object = coef/inflation, file = paste(type, "_lin_nonevent_coef_matrix_",today(),".RDS", sep=""))
  saveRDS(object = phi_vectors, file = paste(type, "_lin_nonevent_eigen_vectors_",today(),".RDS", sep=""))
  saveRDS(object = x, file = paste(type, "_lin_nonevent_complete_case_timesincebaseline_",today(),".RDS", sep=""))
  saveRDS(object = timestamp, file = paste(type, "_lin_nonevent_complete_case_timestamp_",today(),".RDS", sep=""))
  saveRDS(object = id, 
          file = paste(type, "_nonevent_complete_case_ids_",today(),".RDS", sep=""))
  
  ## TAKE PAIRS OF ROWS and CALCULATE THE MSE
  print(paste("Made it to nonevent, quadratic Sigma calc for", type))
  Sigma = matrix(nrow = ncol(quad_residual), ncol = ncol(quad_residual))
  for (i in 1:ncol(quad_residual)) {
    for (j in 1:ncol(quad_residual)) {
      Sigma[i,j] = sum(quad_residual[,i]*quad_residual[,j])/nrow(quad_residual)
    }
  }
  
  Sigma_est <- fbps(Sigma)
  Sigma_est$Yhat[1:10,1:10]
  eig_Sigma_est = eigen(Sigma_est$Yhat)
  
  eig_vectors = eig_Sigma_est$vectors
  eig_values = eig_Sigma_est$values
  
  eig_values[eig_values<0] = 0
  
  K = 35
  print(paste("Variance explained: ", (cumsum(eig_values)/sum(eig_values))[K]))
  phi_vectors = eig_vectors[,1:K]
  coef = residual%*%phi_vectors

  ## SAVE THE MEAN, COEFFICIENT MATRIX, AND FIRST K EIGEN VECTORS
  saveRDS(object = coef/inflation, file = paste(type, "_quad_nonevent_coef_matrix_",today(),".RDS", sep=""))
  saveRDS(object = phi_vectors, file = paste(type, "_quad_nonevent_eigen_vectors_",today(),".RDS", sep=""))
}
