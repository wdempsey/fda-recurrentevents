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

## GENERATE SPLINES
library(splines)
sequence <- seq(-30,0, by = 1/60)
# knots <- seq(-40, 10, 1)  # 10 => 10-4 = 6 Basis splines
# knots <- c(-32,seq(-30,0,1), 2)  # 10 => 10-4 = 6 Basis splines
# x <- seq(-30, 0, by = 1/60)
# bb <- splineDesign(knots, x = x, outer.ok = FALSE, ord = 2)
x <- seq(-30, -5, by = 5)
bb <- matrix(nrow = length(sequence), ncol = length(x) + 1)
bb[,1] = 1
for (i in 1:length(x)) {
  bb[,i+1] = unlist(lapply(X = sequence, function(y) {max(y - x[i],0)}))
}
# knots <- seq(0, 30, 1)  # 10 => 10-4 = 6 Basis splines
# x <- seq(0, 30, by = 1/60)
# bb <- splineDesign(knots, x = x, outer.ok = TRUE)
print("Generated Splines")

library(lubridate)
# setwd("/Volumes/murphy_lab/users/wdempsey/data-for-fda/data/")
setwd("/n/murphy_lab/users/wdempsey/data-for-fda/data/")
## PULL IN EVENT RDS FILES
event_eigen_vectors = readRDS("eda_lin_event_eigen_vectors_2019-07-27.RDS")
event_coef = readRDS("eda_lin_event_coef_matrix_2019-07-27.RDS")
event_means = readRDS("eda_lin_event_means_2019-07-27.RDS")
event_times = readRDS("eda_lin_event_complete_case_times_2019-07-27.RDS") 
# event_id = readRDS("event_complete_case_ids_2019-07-03.RDS")  ## NEED TO FIX
## BUILD outer product of b-spline and event eigen
event_J_coef = event_coef%*%t(event_eigen_vectors)%*%bb
event_J_means = event_means%*%bb
event_model.matrix = cbind(hour(as_datetime(event_times)), event_J_coef + event_J_means)

rm("event_means"); rm("event_coef"); rm("event_eigen_vectors"); rm("event_times")
rm("event_J_coef"); rm("event_J_means")

## PULL IN NONEVENT RDS FILES
nonevent_eigen_vectors = readRDS("eda_lin_nonevent_eigen_vectors_2019-07-27.RDS")
nonevent_coef = readRDS("eda_lin_nonevent_coef_matrix_2019-07-27.RDS")
nonevent_J_coef = nonevent_coef%*%t(nonevent_eigen_vectors)%*%bb
rm("nonevent_eigen_vectors"); rm("nonevent_coef")

nonevent_means = readRDS("eda_lin_nonevent_means_2019-07-27.RDS")
# nonevent_ids = readRDS("nonevent_complete_case_ids_2019-07-03.RDS") ## NEED TO FIX
nonevent_J_means = nonevent_means%*%bb
rm("nonevent_means")
## BUILD outer product of b-spline and event eigen
## BUILD outer product of b-spline and event eigen
nonevent_times = readRDS("eda_lin_nonevent_complete_case_times_2019-07-27.RDS")
nonevent_model.matrix = cbind(hour(as_datetime(nonevent_times)), nonevent_J_coef + nonevent_J_means)
rm("nonevent_J_coef"); rm("nonevent_J_means")

## PULL IN PI_IDS
log_sampling_rate = log(0.25)
print("RDS Files Readin correctly")

Y = c(rep(1,nrow(event_model.matrix)), rep(0, nrow(nonevent_model.matrix)))

model.matrix = rbind(event_model.matrix[,-1], nonevent_model.matrix[,-1])
hour.info = c(event_model.matrix[,1], nonevent_model.matrix[,1])

daytime_obs = (hour.info > 9) & (hour.info < 20)

temp = glm(Y~as.factor(hour.info) + model.matrix - 1, family = "binomial", offset = rep(log_sampling_rate,length(Y)))
summary(temp)

temp_daytime = glm(Y[daytime_obs]~as.factor(hour.info[daytime_obs]) + model.matrix[daytime_obs,] - 1, family = "binomial", offset = rep(log_sampling_rate,length(Y[daytime_obs])))
summary(temp_daytime)

saveRDS(temp, file = "linear_edaonly_alldata_fit.RDS")
saveRDS(temp_daytime, file = "linear_edaonly_daytimedata_fit.RDS")

beta_obs = 25:(25+ncol(bb)-1)
beta = temp$coefficients[beta_obs]
beta[is.na(beta)] = 0

png("~/Downloads/linearfit.png", width = 720, 
    height = 480, units = "px", pointsize = 12)
par(mar = c(4,4,1,1) + 0.1)
obs = sequence > -30 & sequence < 0
plot(sequence[obs], (bb%*%beta)[obs], type= "l", 
     axes = FALSE, xlab = "Time until event", 
     ylab = expression(paste(beta, "(s)")), ylim = c(-0.02,0.02))
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
