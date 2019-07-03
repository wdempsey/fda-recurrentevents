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

setwd("/Volumes/murphy_lab/users/wdempsey/data-for-fda/data/")
event_complete = readRDS("./event_complete_2019-07-01.RDS")
sequence = seq(-30,0,by =1/60); sensor_obs = 4:1804
plot(sequence, event_complete[3,sensor_obs], type = "l")

library(refund)
full_obs = apply(X = event_complete, MARGIN = 1, FUN = function(x){!any(is.na(x))})
greater_than_ten_percent_missing = apply(X = event_complete, MARGIN = 1, FUN = function(x){mean(!is.na(x[sensor_obs])) > 0.1})
greater_than_thirty_percent_missing = apply(X = event_complete, MARGIN = 1, FUN = function(x){mean(!is.na(x[sensor_obs])) > 0.3})
greater_than_fifty_percent_missing = apply(X = event_complete, MARGIN = 1, FUN = function(x){mean(!is.na(x[sensor_obs])) > 0.5})
Y = event_complete[full_obs,sensor_obs]
x = event_complete[full_obs,2]
z = sequence
est <- fbps(Y,list(x=x,z=z))
obs = 100
plot(sequence,est$Yhat[obs,]- mean(est$Yhat[obs,]))
lines(sequence, Y[obs,]-mean(Y[obs,]))

## COMPUTE ESTIMATED MARGINAL COVARIANCE 

inflation = 100
residual = inflation*(Y-est$Yhat)
## TAKE PAIRS OF ROWS and CALCULATE THE MSE
Sigma = matrix(nrow = ncol(residual), ncol = ncol(residual))
for (i in 1:ncol(residual)) {
  for (j in 1:ncol(residual)) {
    Sigma[i,j] = sum(residual[,i]*residual[,j])/nrow(residual)
  }
}

Sigma_est <- fbps(Sigma)
Sigma_est$Yhat[1:10,1:10]
eig_Sigma_est = eigen(Sigma_est$Yhat)

eig_vectors = eig_Sigma_est$vectors
eig_values = eig_Sigma_est$values

eig_values[eig_values<0] = 0

(cumsum(eig_values)/sum(eig_values))[1:35]
K = 35
phi_vectors = eig_vectors[,1:K]
coef = residual%*%phi_vectors
plot(sequence, phi_vectors%*%coef[obs,]/inflation + est$Yhat[obs,], type = "l")
lines(sequence, Y[obs,], col = "red", lwd = 2)

## SAVE THE MEAN, COEFFICIENT MATRIX, AND FIRST K EIGEN VECTORS
saveRDS(object = est$Yhat, file = paste("event_means_",today(),".RDS", sep=""))
saveRDS(object = coef/inflation, file = paste("event_coef_matrix_",today(),".RDS", sep=""))
saveRDS(object = phi_vectors, file = paste("event_eigen_vectors_",today(),".RDS", sep=""))
saveRDS(object = x, file = paste("event_complete_case_times_",today(),".RDS", sep=""))

## CLEAR WORKSPACE
rm(list=ls()) 
setwd("/Volumes/murphy_lab/users/wdempsey/data-for-fda/data/")
nonevent_complete = readRDS("./nonevent_complete_2019-07-01.RDS")
sequence = seq(-30,0,by =1/60); sensor_obs = 4:1804
plot(sequence, nonevent_complete[3,sensor_obs], type = "l")

library(refund)
full_obs = apply(X = nonevent_complete, MARGIN = 1, FUN = function(x){!any(is.na(x))})
Y = nonevent_complete[full_obs,sensor_obs]
x = nonevent_complete[full_obs,2]
z = sequence
rm("nonevent_complete") # Remove the full data to free up memory

est <- fbps(Y,list(x=x,z=z))
obs = 30
plot(sequence,est$Yhat[obs,]- mean(est$Yhat[obs,]))
lines(sequence, Y[obs,]-mean(Y[obs,]))

## COMPUTE ESTIMATED MARGINAL COVARIANCE 
inflation = 100
residual = inflation*(Y-est$Yhat)
saveRDS(object = est$Yhat, file = paste("nonevent_means_",today(),".RDS", sep=""))
rm("Y"); rm("est")
## TAKE PAIRS OF ROWS and CALCULATE THE MSE
Sigma = matrix(nrow = ncol(residual), ncol = ncol(residual))
for (i in 1:ncol(residual)) {
  for (j in 1:ncol(residual)) {
    Sigma[i,j] = sum(residual[,i]*residual[,j])/nrow(residual)
  }
}

Sigma_est <- fbps(Sigma)
Sigma_est$Yhat[1:10,1:10]
eig_Sigma_est = eigen(Sigma_est$Yhat)

eig_vectors = eig_Sigma_est$vectors
eig_values = eig_Sigma_est$values

eig_values[eig_values<0] = 0

(cumsum(eig_values)/sum(eig_values))[1:35]
K = 35
phi_vectors = eig_vectors[,1:K]
coef = residual%*%phi_vectors
plot(sequence, phi_vectors%*%coef[obs,]/inflation + est$Yhat[obs,], type = "l")
lines(sequence, Y[obs,], col = "red", lwd = 2)

## SAVE THE MEAN, COEFFICIENT MATRIX, AND FIRST K EIGEN VECTORS
saveRDS(object = coef/inflation, file = paste("nonevent_coef_matrix_",today(),".RDS", sep=""))
saveRDS(object = phi_vectors, file = paste("nonevent_eigen_vectors_",today(),".RDS", sep=""))
saveRDS(object = x, file = paste("nonevent_complete_case_times_",today(),".RDS", sep=""))
