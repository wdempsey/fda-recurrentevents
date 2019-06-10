setwd("/Volumes/murphy_lab/users/wdempsey/fda-recurrent-events/methods/")
event_complete = readRDS("./event_complete.RDS")
sequence = seq(-30,0,by =1/60)
plot(sequence, event_complete[3,3:1803], type = "l")

library(refund)
full_obs = apply(X = event_complete, MARGIN = 1, FUN = function(x){!any(is.na(x))})
Y = event_complete[full_obs,3:1803]
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
saveRDS(object = est$Yhat, file = "event_means.RDS")
saveRDS(object = coef/inflation, file = "event_coef_matrix.RDS")
saveRDS(object = phi_vectors, file = "event_eigen_vectors.RDS")
saveRDS(object = x, file = "event_complete_case_times.RDS")

## CLEAR WORKSPACE
rm(list=ls()) 
setwd("/Volumes/murphy_lab/users/wdempsey/fda-recurrent-events/methods/")
nonevent_complete = readRDS("./nonevent_complete.RDS")
sequence = seq(-30,0,by =1/60)
plot(sequence, nonevent_complete[3,3:1803], type = "l")

library(refund)
full_obs = apply(X = nonevent_complete, MARGIN = 1, FUN = function(x){!any(is.na(x))})
Y = nonevent_complete[full_obs,3:1803]
x = nonevent_complete[full_obs,2]
z = sequence
rm("nonevent_complete") # Remove the full data to free up memory

est <- fbps(Y,list(x=x,z=z))
obs = 100
plot(sequence,est$Yhat[obs,]- mean(est$Yhat[obs,]))
lines(sequence, Y[obs,]-mean(Y[obs,]))

## COMPUTE ESTIMATED MARGINAL COVARIANCE 
inflation = 100
residual = inflation*(Y-est$Yhat)
saveRDS(object = est$Yhat, file = "nonevent_means.RDS")
rm(c("Y", "est"))
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
saveRDS(object = coef/inflation, file = "nonevent_coef_matrix.RDS")
saveRDS(object = phi_vectors, file = "nonevent_eigen_vectors.RDS")
saveRDS(object = x, file = "nonevent_complete_case_times.RDS")
