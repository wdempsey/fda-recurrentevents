## WINDOWS
setwd("Z:/SI_data/")
## LINUX 
setwd('/mnt/turbo/SI_data/')

## LIBRARIES
library(ggplot2)
library('dplyr')

## Read in acc and eda model.matrix
## ACC
acc_event_model.matrix = readRDS("acc_event_modelmatrix_2021-12-21.RDS")
acc_nonevent_model.matrix = readRDS("acc_nonevent_modelmatrix_2021-12-21.RDS")
## EDA
eda_event_model.matrix = readRDS("eda_event_modelmatrix_2021-12-21.RDS")
eda_nonevent_model.matrix = readRDS("eda_nonevent_modelmatrix_2021-12-21.RDS")

## GENERATE SPLINES
K_b = 35
num=K_b-3
qtiles <- seq(0, 1, length = num + 2)[-c(1, num + 2)]
## EDA
eda_sequence <- seq(-30,0, by = 1/60)
eda_knots <- quantile(eda_sequence, qtiles)
eda_bb = cbind(1, eda_sequence, eda_sequence^2, sapply(eda_knots, function(k) ((eda_sequence - k > 0) * (eda_sequence - k)) ^ 2))
print("Generated EDA Splines")
## ACC 
acc_sequence <- seq(-30,0, by = 1/6)
acc_knots <- quantile(acc_sequence, qtiles)
acc_bb = cbind(1, acc_sequence, acc_sequence^2, sapply(acc_knots, function(k) ((acc_sequence - k > 0) * (acc_sequence - k)) ^ 2))
print("Generated ACC Splines")
## PULL IN PI_IDS
log_sampling_rate = log(0.5)
print("RDS Files Reading correctly")

## Bring in glmnet
if(!require("glmnet")){install.packages('glmnet')}
library('glmnet')

## Model per type
acc_Y = c(rep(1,nrow(acc_event_model.matrix)), rep(0, nrow(acc_nonevent_model.matrix)))
acc_model.matrix = rbind(acc_event_model.matrix[,-c(1:3)], acc_nonevent_model.matrix[,-c(1:3)])
## HOUR WAS IN UTC BUT NEED TO BE TRANSLATED TO ETS
## THIS IS DONE BY -5 HOURS FUNCTION
hour.info = (c(acc_event_model.matrix[,3], acc_nonevent_model.matrix[,3]) - 5)%%24
daytime_obs = (hour.info > 9) & (hour.info < 20)
n.tmp = length(acc_Y)
p.tmp = ncol(acc_model.matrix)
subsample_offset = rep(log_sampling_rate,nrow(acc_model.matrix))
p.fac = rep(1, ncol(acc_model.matrix))
p.fac[1:3] = 0 #no penalty on the first 3 variables ## HOW IS INTERCEPT HANDLED?
# set.seed("97139817")
# Start the clock!
acc_model.matrix = as.matrix(acc_model.matrix)
ptm <- proc.time()
ridge.fit.cv <- cv.glmnet(acc_model.matrix[daytime_obs,], acc_Y[daytime_obs], alpha = 0, 
                          intercept = TRUE, penalty.factor = p.fac, standardize = F,
                          offset = rep(log_sampling_rate,length(acc_Y[daytime_obs])),
                          family = "binomial")
# Stop the clock
runtime = proc.time() - ptm

ridge.fit.lambda <- ridge.fit.cv$lambda.min
# Extract coefficient values for lambda.1se (without intercept)
ridge.coef <- (coef(ridge.fit.cv, s = ridge.fit.lambda))[-1]
intercepts <- (coef(ridge.fit.cv, s = ridge.fit.lambda))[1]

X = cbind(1,acc_model.matrix[daytime_obs,])
totals = X%*%coef(ridge.fit.cv, s = ridge.fit.lambda) + log_sampling_rate
probs = 1/(1+exp(-totals))
W = diag(as.vector(probs * (1-probs)))
fisher_info = t(X)%*%W%*%X
ridge_penalty = diag(ridge.fit.lambda*c(1,p.fac))
Sigma = solve(fisher_info+ridge_penalty)
updated_Sigma = Sigma[-1,-1]
stderr = sqrt(diag(acc_bb%*%updated_Sigma%*%t(acc_bb)))
betaHat.net <- acc_bb %*% ridge.coef

df_acc_summary <- data.frame(sequence = acc_sequence+30, estimate = acc_betaHat.net,
                             lowerCI = acc_betaHat.net - 1.96 * acc_stderr,
                             upperCI = acc_betaHat.net + 1.96 * acc_stderr)

acc_xmax = max(df_acc_summary$sequence[which(df_acc_summary$lowerCI > 0)])
acc_xmin = min(df_acc_summary$sequence[which(df_acc_summary$lowerCI > 0)])

png("C:/Users/Balthazar/Documents/GitHub/fda-recurrentevents/figures/acc_coef.png",
    width = 480, height = 480, units = "px", pointsize = 16)
ggplot(df_acc_summary, aes(x=sequence, y=estimate)) +
  geom_line(size=1, alpha=0.8) +
  geom_ribbon(aes(ymin=lowerCI, ymax=upperCI) ,fill="blue", alpha=0.2) +
  annotate("rect",xmin=acc_xmin,xmax=acc_xmax,ymin=-Inf,ymax=Inf, alpha=0.1, fill="black") +
  xlab("Time until Button Press") + ylab(expression(paste(beta, "(s)")))
dev.off()

### EDA
eda_Y = c(rep(1,nrow(eda_event_model.matrix)), rep(0, nrow(eda_nonevent_model.matrix)))
eda_model.matrix = rbind(eda_event_model.matrix[,-c(1:3)], eda_nonevent_model.matrix[,-c(1:3)])
hour.info = c(eda_event_model.matrix[,3], eda_nonevent_model.matrix[,3])
daytime_obs = (hour.info > 9) & (hour.info < 20)

n.tmp = length(eda_Y)
p.tmp = ncol(eda_model.matrix)
subsample_offset = rep(log_sampling_rate,nrow(eda_model.matrix))
p.fac = rep(1, ncol(eda_model.matrix))
p.fac[1:3] = 0 #no penalty on the first 3 variables

# set.seed("97139817")
# Start the clock!
eda_model.matrix = as.matrix(eda_model.matrix)
ptm <- proc.time()
ridge.fit.cv <- cv.glmnet(eda_model.matrix[daytime_obs,], eda_Y[daytime_obs], 
                          alpha = 0, intercept = TRUE, 
                          penalty.factor = p.fac, standardize = FALSE,
                          # lambda = lambdapath, nfolds = 20,
                          family = "binomial", 
                          offset = rep(log_sampling_rate,length(eda_Y[daytime_obs])))

# Stop the clock
runtime = proc.time() - ptm

ridge.fit.lambda <- ridge.fit.cv$lambda.min
# Extract coefficient values for lambda.1se (without intercept)
ridge.coef <- (coef(ridge.fit.cv, s = ridge.fit.lambda))[-1]
intercept <- (coef(ridge.fit.cv, s = ridge.fit.lambda))[1]

## STDERR
X = cbind(1,eda_model.matrix[daytime_obs,])
totals = X%*%coef(ridge.fit.cv, s = ridge.fit.lambda) + log_sampling_rate
probs = 1/(1+exp(-totals))
W = diag(as.vector(probs * (1-probs)))
fisher_info = t(X)%*%W%*%X
ridge_penalty = diag(ridge.fit.lambda*c(1,p.fac))
Sigma = solve(fisher_info+ridge_penalty)
updated_Sigma = Sigma[-1,-1]
stderr = sqrt(diag(eda_bb%*%updated_Sigma%*%t(eda_bb)))
betaHat.net <- eda_bb %*% ridge.coef

df_eda_summary <- data.frame(sequence = eda_sequence+30, estimate = eda_betaHat.net,
                             lowerCI = eda_betaHat.net - 1.96 * eda_stderr,
                             upperCI = eda_betaHat.net + 1.96 * eda_stderr)

# eda_xmax = max(df_eda_summary$sequence[which(df_eda_summary$lowerCI > 0)])
# eda_xmin = min(df_eda_summary$sequence[which(df_eda_summary$lowerCI > 0)])

png("C:/Users/Balthazar/Documents/GitHub/fda-recurrentevents/figures/eda_coef.png",
    width = 480, height = 480, units = "px", pointsize = 16)
ggplot(df_eda_summary, aes(x=sequence, y=estimate)) +
  geom_line(size=1, alpha=0.8) +
  geom_ribbon(aes(ymin=lowerCI, ymax=upperCI) ,fill="blue", alpha=0.2) +
  # annotate("rect",xmin=acc_xmin,xmax=acc_xmax,ymin=-Inf,ymax=Inf, alpha=0.1, fill="black") +
  xlab("Time until Button Press") + ylab(expression(paste(beta, "(s)")))
dev.off()


### MERGE DATASETS
names(eda_event_model.matrix)[4:length(eda_event_model.matrix)] = paste("eda_X", 1:35, sep ="")
names(eda_nonevent_model.matrix)[4:length(eda_nonevent_model.matrix)] = paste("eda_X", 1:35, sep ="")

all_event_model.matrix = left_join(acc_event_model.matrix, eda_event_model.matrix)
all_nonevent_model.matrix = left_join(acc_nonevent_model.matrix, eda_nonevent_model.matrix)

isna_event = apply(X = all_event_model.matrix, MARGIN = 1, FUN = function(x){any(is.na(x))})
isna_nonevent = apply(X = all_nonevent_model.matrix, MARGIN = 1, FUN = function(x){any(is.na(x))})
all_event_model.matrix = all_event_model.matrix[!isna_event,]
all_nonevent_model.matrix = all_nonevent_model.matrix[!isna_nonevent,]

all_Y = c(rep(1,nrow(all_event_model.matrix)), rep(0, nrow(all_nonevent_model.matrix)))
all_model.matrix = rbind(all_event_model.matrix[,-c(1:3)], all_nonevent_model.matrix[,-c(1:3)])
hour.info = c(all_event_model.matrix[,3], all_nonevent_model.matrix[,3])
daytime_obs = (hour.info > 9) & (hour.info < 20)

n.tmp = length(all_Y)
p.tmp = ncol(all_model.matrix)
subsample_offset = rep(log_sampling_rate,nrow(all_model.matrix))
p.fac = rep(1, ncol(all_model.matrix))
p.fac[c(1:3)] = 0 #no penalty on the first 4 variables

# set.seed("97139817")
# Start the clock!
all_model.matrix = as.matrix(all_model.matrix)
ptm <- proc.time()
ridge.fit.cv <- cv.glmnet(all_model.matrix[daytime_obs,], all_Y[daytime_obs], 
                          alpha = 0, intercept = TRUE, 
                          penalty.factor = p.fac, standardize = FALSE,
                          # lambda = lambdapath, nfolds = 20,
                          offset = rep(log_sampling_rate,length(all_Y[daytime_obs])),
                          family = "binomial")
# Stop the clock
runtime = proc.time() - ptm

ridge.fit.lambda <- ridge.fit.cv$lambda.min
# plot(ridge.fit.cv)

# Extract coefficient values for lambda.1se (without intercept)
ridge.coef <- (coef(ridge.fit.cv, s = ridge.fit.lambda))[-1]
intercept <- (coef(ridge.fit.cv, s = ridge.fit.lambda))[1]

acc_betaHat.net <- acc_bb %*% ridge.coef[1:35]
eda_betaHat.net <- eda_bb %*% ridge.coef[36:70]

plot(acc_sequence, acc_betaHat.net)
plot(eda_sequence, eda_betaHat.net)

## STDERR
X = cbind(1,all_model.matrix[daytime_obs,])
totals = X%*%coef(ridge.fit.cv, s = ridge.fit.lambda) + log_sampling_rate
probs = 1/(1+exp(-totals))
W = diag(as.vector(probs * (1-probs)))
fisher_info = t(X)%*%W%*%X
ridge_penalty = diag(ridge.fit.lambda*c(1,p.fac))
Sigma = solve(fisher_info+ridge_penalty)
updated_Sigma = Sigma[-1,-1]
acc_updated_Sigma = updated_Sigma[1:35,1:35]
eda_updated_Sigma = updated_Sigma[36:70,36:70]
acc_stderr = sqrt(diag(acc_bb%*%acc_updated_Sigma%*%t(acc_bb)))
eda_stderr = sqrt(diag(eda_bb%*%eda_updated_Sigma%*%t(eda_bb)))

## JOINT FIGURES
df_acc_summary <- data.frame(sequence = acc_sequence+30, estimate = acc_betaHat.net,
                 lowerCI = acc_betaHat.net - 1.96 * acc_stderr,
                 upperCI = acc_betaHat.net + 1.96 * acc_stderr)

acc_xmax = max(df_acc_summary$sequence[which(df_acc_summary$lowerCI > 0)])
acc_xmin = min(df_acc_summary$sequence[which(df_acc_summary$lowerCI > 0)])

png("C:/Users/Balthazar/Documents/GitHub/fda-recurrentevents/figures/acc_coef_joint.png",
    width = 480, height = 480, units = "px", pointsize = 16)
ggplot(df_acc_summary, aes(x=sequence, y=estimate)) +
  geom_line(size=1, alpha=0.8) +
  geom_ribbon(aes(ymin=lowerCI, ymax=upperCI) ,fill="blue", alpha=0.2) +
  annotate("rect",xmin=acc_xmin,xmax=acc_xmax,ymin=-Inf,ymax=Inf, alpha=0.1, fill="black") +
  xlab("Time until Button Press") + ylab(expression(paste(beta, "(s)")))
dev.off()
  # geom_line(data=df_tidy, aes(x=Time, y=Ratio, group=Cell), color="grey") +

df_eda_summary <- data.frame(sequence = eda_sequence+30, estimate = eda_betaHat.net,
                         lowerCI = eda_betaHat.net - 1.96 * eda_stderr,
                         upperCI = eda_betaHat.net + 1.96 * eda_stderr)

# eda_xmax = max(df_eda_summary$sequence[which(df_eda_summary$lowerCI > 0)])
# eda_xmin = min(df_eda_summary$sequence[which(df_eda_summary$lowerCI > 0)])

png("C:/Users/Balthazar/Documents/GitHub/fda-recurrentevents/figures/eda_coef_joint.png",
    width = 480, height = 480, units = "px", pointsize = 16)
ggplot(df_eda_summary, aes(x=sequence, y=estimate)) +
  geom_line(size=1, alpha=0.8) +
  geom_ribbon(aes(ymin=lowerCI, ymax=upperCI) ,fill="blue", alpha=0.2) +
  # annotate("rect",xmin=acc_xmin,xmax=acc_xmax,ymin=-Inf,ymax=Inf, alpha=0.1, fill="black") +
  xlab("Time until Button Press") + ylab(expression(paste(beta, "(s)")))
dev.off()

