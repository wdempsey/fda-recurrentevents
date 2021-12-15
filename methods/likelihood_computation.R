## Code to compute the approximate likelihood function.
## Recall the approximation is design unbiased
## That is, E [ loglik (\theta; subssampled data)] = loglik(\theta; full data)
## where the expecation is over the random subsampling

## Inputs: RDS files for event and nonevent dataframes:
## Inputs: (1) mean functions mu(t,s), 
## Inputs: (2) coefficients on the residual X(t,s) - mu(t,s)
## Inputs: (3) eigenvectors of the pooled covariance estimate 
## Inputs: (4) complete-case timestamps
## Output: RDS file (mle_output.RDS) that contains 
## Output: (1) the approximate MLE (mle_output$estimate),
## Output: (2) the standard errors (mle_output$stderr), and 
## Output: (2) subsampling component of variance (mle_output$var_subsampl), and 
## Output: (4) recurrent event component of variance (mle_output$var_event). 
## Relative efficiency for any functional \Var (c^\prime \hat theta) can be computed via 
## numerator =  c^prime %*% mle_output$var_event %*% c 
## denominator =  c^prime %*% (mle_output$var_subsampl + mle_output$var_event ) %*% c 
## rel eff = numerator / denominator

## LIBRARIES
library(splines)
library(lubridate)
## WINDOWS
setwd("Z:/SI_data/")
## LINUX 
setwd('/mnt/turbo/SI_data/')

## PULL IN EVENT RDS FILES
set_of_types = c("acc", "eda")
## GENERATE SPLINES
if (type == "eda") {
  sequence <- seq(-30,0, by = 1/60)
} else if (type == "acc") {
  sequence <- seq(-30,0, by = 1/6)
}
# knots <- seq(-40, 10, 1)  # 10 => 10-4 = 6 Basis splines
# knots <- c(-32,seq(-30,0,1), 2)  # 10 => 10-4 = 6 Basis splines
# x <- seq(-30, 0, by = 1/60)
# bb <- splineDesign(knots, x = x, outer.ok = FALSE, ord = 2)
x <- seq(-30,-5, by = 5)
bb <- matrix(nrow = length(sequence), ncol = length(x) + 1)
bb[,1] = 1
for (i in 1:length(x)) {
  bb[,i+1] = unlist(lapply(X = sequence, function(y) {max(y - x[i],0)}))
}
print("Generated Splines")

if(!file.exists(paste(type, "_event_modelmatrix_2021-11-23.RDS", sep = ""))) {
  event_eigen_vectors = readRDS(paste(type, "_lin_event_eigen_vectors_2021-11-22.RDS", sep = ""))
  event_coef = readRDS(paste(type, "_lin_event_coef_matrix_2021-11-22.RDS", sep = ""))
  event_means = readRDS(paste(type, "_lin_event_means_2021-11-22.RDS", sep = ""))
  event_times = readRDS(paste(type, "_lin_event_complete_case_times_2021-11-22.RDS", sep = ""))
  event_id = readRDS(paste(type, "_lin_event_complete_case_ids_2021-11-22.RDS", sep = ""))
  ## BUILD outer product of b-spline and event eigen
  event_J_coef = event_coef%*%t(event_eigen_vectors)%*%bb
  event_J_means = event_means%*%bb
  event_model.matrix = cbind(hour(as_datetime(event_times)), event_J_coef + event_J_means)
  
  rm("event_means"); rm("event_coef"); rm("event_eigen_vectors"); rm("event_times")
  rm("event_J_coef"); rm("event_J_means")
  
  saveRDS(event_model.matrix, file = paste(type, "_event_modelmatrix_",today(), ".RDS", sep = ""))
} else {
  event_model.matrix = readRDS(paste(type,"_event_modelmatrix_2021-12-14.RDS",sep =""))
}
## PULL IN NONEVENT RDS FILES
if(!file.exists("eda_nonevent_modelmatrix_2021-11-23.RDS")) {
  nonevent_eigen_vectors = readRDS(paste(type, "_lin_nonevent_eigen_vectors_2021-11-22.RDS", sep =""))
  nonevent_coef = readRDS(paste(type, "_lin_nonevent_coef_matrix_2021-11-22.RDS", sep = ""))
  nonevent_J_coef = nonevent_coef%*%t(nonevent_eigen_vectors)%*%bb
  rm("nonevent_eigen_vectors"); rm("nonevent_coef")
  
  nonevent_ids = readRDS(paste(type, "_nonevent_complete_case_ids_2021-11-22.RDS", sep =""))
  nonevent_means = readRDS(paste(type, "_lin_nonevent_means_2021-11-22.RDS", sep = ""))
  nonevent_J_means = nonevent_means%*%bb
  rm("nonevent_means")
  ## BUILD outer product of b-spline and event eigen
  ## BUILD outer product of b-spline and event eigen
  nonevent_times = readRDS(paste(type, "_lin_nonevent_complete_case_times_2021-11-22.RDS", sep =""))
  nonevent_model.matrix = cbind(hour(as_datetime(nonevent_times)), nonevent_J_coef + nonevent_J_means)
  rm("nonevent_J_coef"); rm("nonevent_J_means")
  saveRDS(nonevent_model.matrix, file = paste(type, "_nonevent_modelmatrix_",today(), ".RDS", sep = ""))
} else {
  nonevent_model.matrix = readRDS("acc_nonevent_modelmatrix_2021-12-14.RDS")
}

## PULL IN PI_IDS
log_sampling_rate = log(0.5)
print("RDS Files Reading correctly")

## Bring in glmnet
if(!require("glmnet")){install.packages('glmnet')}
library('glmnet')

Y = c(rep(1,nrow(event_model.matrix)), rep(0, nrow(nonevent_model.matrix)))

model.matrix = rbind(event_model.matrix[,-1], nonevent_model.matrix[,-1])
hour.info = c(event_model.matrix[,1], nonevent_model.matrix[,1])
daytime_obs = (hour.info > 9) & (hour.info < 20)

temp = glmnet(Y~as.factor(hour.info) + model.matrix - 1, family = "binomial", offset = rep(log_sampling_rate,length(Y)))
summary(temp)

temp_daytime = glmnet(Y[daytime_obs]~model.matrix[daytime_obs,] - 1, family = "binomial", offset = rep(log_sampling_rate,length(Y[daytime_obs])))
summary(temp_daytime)

saveRDS(temp, file = "linear_edaonly_alldata_fit.RDS")
saveRDS(temp_daytime, file = "linear_edaonly_daytimedata_fit.RDS")

saveRDS(temp_daytime_plusacc, file = "linear_edaacc_daytimedata_fit.RDS")

beta_obs = 25:(25+ncol(bb)-1)
beta = temp$coefficients[beta_obs]
beta[is.na(beta)] = 0

png("~/Downloads/linearfit.png", width = 720, 
    height = 480, units = "px", pointsize = 12)
par(mar = c(4,4,1,1) + 0.1)
obs = sequence > -30 & sequence < 0
plot(sequence[obs], (bb%*%beta)[obs], type= "l", 
     axes = FALSE, xlab = "Time until event", 
     ylab = expression(paste(beta, "(s)")))
axis(side = 1); axis(side = 2)
Sigma = vcov(temp)[beta_obs, beta_obs]
Sigma[is.na(Sigma)] = 0
stderr = sqrt(diag(bb%*%Sigma%*%t(bb)))

lines(sequence[obs], (bb%*%beta + 1.96 * stderr)[obs], col = "red")
lines(sequence[obs], (bb%*%beta - 1.96 * stderr)[obs], col = "red")
abline(h = 0, col = "blue")
dev.off()

sig_obs_pos = bb%*%beta + 1.96 * stderr > 0 & bb%*%beta - 1.96 * stderr > 0
sig_obs_neg = bb%*%beta + 1.96 * stderr < 0 & bb%*%beta - 1.96 * stderr < 0
sig_obs = sig_obs_pos | sig_obs_neg
which(sig_obs)
lines(sequence[1:80], (bb%*%beta)[1:80], col = "red", lwd = 4)
lines(sequence[1412:1530], (bb%*%beta)[1412:1530], col = "red", lwd = 4)
lines(sequence[1762:1768], (bb%*%beta)[1762:1768], col = "red", lwd = 4)



# saveRDS(object = output, file = "wedidit2.RDS")

