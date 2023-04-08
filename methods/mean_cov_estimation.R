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

## WINDOWS
setwd("Z:/SI_data/")
## LINUX
# setwd("/mnt/turbo/SI_data/")

set_of_types = c("eda", "acc")
Delta = 15
for (type in set_of_types){
  event_complete = readRDS(paste(type, "_event_complete_Delta_",Delta,"_HLP_2023-04-06.RDS", sep = ""))
  sequence = seq(0,-Delta, length.out = ncol(event_complete) - 3); sensor_obs = 4:ncol(event_complete)
  # plot(sequence, event_complete[3,sensor_obs], type = "l")

  library(refund); library(lubridate)
  full_obs = apply(X = event_complete, MARGIN = 1, FUN = function(x){!any(is.na(x))})
  Y = event_complete[full_obs,sensor_obs]
  Y = as.matrix(Y, nrow = nrow(Y), ncol = ncol(Y))
  x = as.numeric(event_complete[full_obs,3])
  z = sequence
  est <- fbps(Y,list(x=x,z=z))
  
  ## COMPUTE ESTIMATED MARGINAL COVARIANCE
  inflation = 100
  residual = inflation*(Y-est$Yhat)
  ## TAKE PAIRS OF ROWS and CALCULATE THE MSE
  print(paste("Made it to event, linear Sigma calc for", type))
  Sigma = matrix(nrow = ncol(residual), ncol = ncol(residual))
  for (i in 1:ncol(residual)) {
    for (j in 1:ncol(residual)) {
      Sigma[i,j] = sum(residual[,i]*residual[,j])/nrow(residual)
    }
  }
  Sigma_est <- fbps(Sigma)
  final_Sigma = (Sigma_est$Yhat + t(Sigma_est$Yhat))/2
  saveRDS(object = final_Sigma, file = paste(type, "_event_Sigma_Delta_",Delta,"_",today(),".RDS", sep=""))
  eig_Sigma_est = eigen(final_Sigma)

  eig_vectors = eig_Sigma_est$vectors
  eig_values = eig_Sigma_est$values

  eig_values[eig_values<0] = 0
  K = 35
  print(paste("Variance explained: ", (cumsum(eig_values)/sum(eig_values))[K]))
  phi_vectors = eig_vectors[,1:K]
  coef = residual%*%phi_vectors
  ## To show the impact of K
  obs = 20
  optK = 10
  optphi_vectors = eig_vectors[,1:optK]
  optcoef = residual%*%optphi_vectors
  if(type == "eda") {legend_name = "EDA"} else{legend_name = "ACC"}
  par(mar = c(4,3,2,1)+ 0.1)
  png(filename = paste("~/GitHub/fda-recurrentevents/figures/",
  # png(filename = paste("~/Documents/github/fda-recurrentevents/figures/",
                       type, "_fda_Delta_",Delta,"_fitplot.png", sep = ""),
      width = 600, height = 480, units = "px", pointsize = 12)
  plot(sequence, phi_vectors%*%coef[obs,]/inflation + est$Yhat[obs,], type = "l", axes = FALSE, xlab = "Time until button press", 
       ylab = "", col = "red", lwd = 2, ylim = c(min(Y[obs,]), max(Y[obs,]))) # , ylim = c(0.02,0.08))
  axis(side = 1); axis(side = 2)
  lines(sequence, Y[obs,], col = "black", lwd = 1, lty = 2)
  lines(sequence, optphi_vectors%*%optcoef[obs,]/inflation + est$Yhat[obs,], col = "red", lty = 2)
  legend(-10,max(Y[obs,])*0.8, c(paste("Observed", legend_name), "K=10", "K=35"), col = c("black", "red", "red"), lty = c(1,1,2), cex = 0.75, bty = "n")
  dev.off()
  
  ## Timestamps
  timestamp = event_complete[full_obs,2]
  ## SAVE THE MEAN, COEFFICIENT MATRIX, AND FIRST K EIGEN VECTORS
  saveRDS(object = est$Yhat, file = paste(type, "_lin_event_means_Delta_",Delta,"_",today(),".RDS", sep=""))
  saveRDS(object = coef/inflation, file = paste(type, "_lin_event_coef_matrix_Delta_",Delta,"_",today(),".RDS", sep=""))
  saveRDS(object = phi_vectors, file = paste(type, "_lin_event_eigen_vectors_Delta_",Delta,"_",today(),".RDS", sep=""))
  saveRDS(object = x, file = paste(type, "_lin_event_complete_case_timesincebaseline_Delta_",Delta,"_",today(),".RDS", sep=""))
  saveRDS(object = timestamp, file = paste(type, "_lin_event_complete_case_timestamp_Delta_",Delta,"_",today(),".RDS", sep=""))
  saveRDS(object = as.numeric(event_complete[full_obs,1]), 
          file = paste(type, "_lin_event_complete_case_ids_Delta_",Delta,"_",today(),".RDS", sep=""))
  
}

## CLEAR WORKSPACE
rm(list=ls()) 
## WINDOWS
# setwd("Z:/SI_data/")
## Linux
setwd("/mnt/turbo/SI_data/")
set_of_types = c("eda", "acc")
Delta = 15
for (type in set_of_types){
  nonevent_complete = readRDS(paste(type, "_nonevent_complete_Delta_",Delta,"_HLP_2023-04-06.RDS", sep = ""))
  sequence = seq(0,-Delta,length.out = ncol(nonevent_complete) - 3); sensor_obs = 4:ncol(nonevent_complete)
  # plot(sequence, nonevent_complete[4,sensor_obs])
  
  library(refund); library(lubridate)
  full_obs = apply(X = nonevent_complete, MARGIN = 1, FUN = function(x){!any(is.na(x))})
  ## print(paste(length(which(!full_obs)), "out of", length(full_obs), "missing at least 1 entry"))
  Y = nonevent_complete[full_obs,sensor_obs]
  Y = as.matrix(Y, nrow = nrow(Y), ncol = ncol(Y))
  mu = rowMeans(Y)
  x = as.numeric(nonevent_complete[full_obs,3])
  z = sequence
  id = nonevent_complete[full_obs,1]
  ## Timestamps
  timestamp = nonevent_complete[full_obs,2]
  rm("nonevent_complete") # Remove the full data to free up memory
  
  est <- fbps(Y,list(x=x,z=z))
  
  ## COMPUTE ESTIMATED MARGINAL COVARIANCE 
  inflation = 100
  residual = inflation*(Y-est$Yhat)
  saveRDS(object = est$Yhat, file = paste(type, "_lin_nonevent_means_Delta_",Delta,"_",today(),".RDS", sep=""))
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
  saveRDS(object = final_Sigma, file = paste(type, "_nonevent_Sigma_Delta_",Delta,"_",today(),".RDS", sep=""))
  eig_Sigma_est = eigen(final_Sigma)
  eig_vectors = eig_Sigma_est$vectors
  eig_values = eig_Sigma_est$values
  
  eig_values[eig_values<0] = 0
  
  K = 35
  print(paste("Variance explained: ", (cumsum(eig_values)/sum(eig_values))[K]))
  phi_vectors = eig_vectors[,1:K]
  coef = residual%*%phi_vectors
  
  
  ## SAVE THE MEAN, COEFFICIENT MATRIX, AND FIRST K EIGEN VECTORS
  saveRDS(object = coef/inflation, file = paste(type, "_lin_nonevent_coef_matrix_Delta_",Delta,"_",today(),".RDS", sep=""))
  saveRDS(object = phi_vectors, file = paste(type, "_lin_nonevent_eigen_vectors_Delta_",Delta,"_",today(),".RDS", sep=""))
  saveRDS(object = x, file = paste(type, "_lin_nonevent_complete_case_timesincebaseline_Delta_",Delta,"_",today(),".RDS", sep=""))
  saveRDS(object = timestamp, file = paste(type, "_lin_nonevent_complete_case_timestamp_Delta_",Delta,"_",today(),".RDS", sep=""))
  saveRDS(object = id, 
          file = paste(type, "_nonevent_complete_case_ids_Delta_",Delta,"_",today(),".RDS", sep=""))
}
