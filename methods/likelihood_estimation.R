## WINDOWS
setwd("Z:/SI_data/")
## LINUX 
setwd('/mnt/turbo/SI_data/')

## Read in acc and eda model.matrix
## ACC
acc_event_model.matrix = readRDS("acc_event_modelmatrix_2021-12-16.RDS")
acc_nonevent_model.matrix = readRDS("acc_nonevent_modelmatrix_2021-12-16.RDS")
## EDA
eda_event_model.matrix = readRDS("eda_event_modelmatrix_2021-12-16.RDS")
eda_nonevent_model.matrix = readRDS("eda_nonevent_modelmatrix_2021-12-16.RDS")

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

n.tmp = length(Y)
p.tmp = ncol(model.matrix)
subsample_offset = rep(log_sampling_rate,nrow(dataset))
p.fac = rep(1, ncol(model.matrix))
p.fac[1:4] = 0 #no penalty on the first 4 variables
lambda_max <- 1/n.tmp
epsilon <- .0001
K <- 30
lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon), 
                            length.out = K)), digits = 10)

# set.seed("97139817")
# Start the clock!
ptm <- proc.time()
ridge.fit.cv <- cv.glmnet(model.matrix[daytime_obs,], Y[daytime_obs], alpha = 0, intercept = TRUE, 
                          penalty.factor = p.fac, standardize = FALSE,
                          lambda = lambdapath, nfolds = 20,
                          family = "binomial")
# Stop the clock
runtime = proc.time() - ptm

ridge.fit.lambda <- ridge.fit.cv$lambda.min
# plot(ridge.fit.cv)

# Extract coefficient values for lambda.1se (without intercept)
ridge.coef <- (coef(ridge.fit.cv, s = ridge.fit.lambda))[-1]
intercept <- (coef(ridge.fit.cv, s = ridge.fit.lambda))[1]

betaHat.net <- bb %*% ridge.coef

length(betaHat.net)

plot(sequence[1:1500],betaHat.net[1:1500])

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

