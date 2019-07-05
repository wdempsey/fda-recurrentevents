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

## Source the underlying functions
source("./likelihood_functions.R")

library(doParallel)
cl = Sys.getenv("SLURM_NTASKS_PER_NODE")
cl

registerDoParallel(cores = (cl))

# Shows the number of Parallel Workers to be used
getDoParWorkers()
print("Parallel setup complete")
## GENERATE SPLINES
library(splines)
sequence <- seq(-30,0, by = 1/60)
# knots <- seq(-36, 6, 2)  # 10 => 10-4 = 6 Basis splines
# x <- seq(-30, 0, by = 1/60)
# bb <- splineDesign(knots, x = x, outer.ok = FALSE)
knots <- seq(0, 30, 1)  # 10 => 10-4 = 6 Basis splines
x <- seq(0, 30, by = 1/60)
bb <- splineDesign(knots, x = x, outer.ok = TRUE)
print("Generated Splines")

setwd("/Volumes/murphy_lab/users/wdempsey/data-for-fda/data/")
# setwd("/n/murphy_lab/users/wdempsey/data-for-fda/data/")
## PULL IN EVENT RDS FILES
event_eigen_vectors = readRDS("event_eigen_vectors_2019-07-03.RDS")
event_coef = readRDS("event_coef_matrix_2019-07-03.RDS")
event_means = readRDS("event_means_2019-07-03.RDS")
event_times = readRDS("event_complete_case_times_2019-07-03.RDS") 
# event_id = readRDS("event_complete_case_ids_2019-07-03.RDS")  ## NEED TO FIX
## BUILD outer product of b-spline and event eigen
event_J_coef = event_coef%*%t(event_eigen_vectors)%*%bb
event_J_means = event_means%*%bb
event_model.matrix = cbind(hour(as_datetime(event_times)), event_J_coef + event_J_means)

rm("event_means"); rm("event_coef"); rm("event_eigen_vectors"); rm("event_times")
rm("event_J_coef"); rm("event_J_means")

## PULL IN NONEVENT RDS FILES
nonevent_eigen_vectors = readRDS("nonevent_eigen_vectors_2019-07-03.RDS")
nonevent_coef = readRDS("nonevent_coef_matrix_2019-07-03.RDS")
nonevent_J_coef = nonevent_coef%*%t(nonevent_eigen_vectors)%*%bb
rm("nonevent_eigen_vectors"); rm("nonevent_coef")

nonevent_means = readRDS("nonevent_means_2019-07-03.RDS")
# nonevent_ids = readRDS("nonevent_complete_case_ids_2019-07-03.RDS") ## NEED TO FIX
nonevent_J_means = nonevent_means%*%bb
rm("nonevent_means")
## BUILD outer product of b-spline and event eigen
## BUILD outer product of b-spline and event eigen
nonevent_times = readRDS("nonevent_complete_case_times_2019-07-03.RDS")
nonevent_model.matrix = cbind(hour(as_datetime(nonevent_times)), nonevent_J_coef + nonevent_J_means)
rm("nonevent_J_coef"); rm("nonevent_J_means")

## PULL IN PI_IDS
log_sampling_rate = log(0.5)
print("RDS Files Readin correctly")

Y = c(rep(1,nrow(event_model.matrix)), rep(0, nrow(nonevent_model.matrix)))

model.matrix = rbind(event_model.matrix, nonevent_model.matrix)

temp = glm(Y~model.matrix-1, family = binomial, offset = rep(log_sampling_rate,length(Y)))

beta = temp$coefficients[-1]

par(mar = c(4,4,2,1) + 0.1)
plot(sequence, bb%*%beta, type= "l")

stderr = sqrt(diag(bb%*%vcov(temp)[2:28, 2:28]%*%t(bb)))

lines(sequence, bb%*%beta + 1.96 * stderr, col = "red")
lines(sequence, bb%*%beta - 1.96 * stderr, col = "red")

# ## Bring in rootSolve library for 
# ## using Newton Rhapson method on 
# ## collected data
# library("rootSolve")
# 
# llik <- llik_computation(event_coef, event_means, event_J,
#                          nonevent_coef, nonevent_means, nonevent_J,
#                          nonevent_ids, pi_ids) 
# 
# llik_deriv <- llik_computation_derivative(event_coef, event_means, event_J,
#                                           nonevent_coef, nonevent_means, nonevent_J,
#                                           nonevent_ids, pi_ids) 
# temp = readRDS("wedidit.RDS")
# init_beta = temp$par ## init_beta = c(1, rnorm(ncol(nonevent_J), sd = 0.1))
# 
# print("Optimization begins")
# 
# system.time(print(llik(init_beta)))
# 
# output <- optim(init_beta, fn = llik, gr = llik_deriv, method = "BFGS", control=list(trace=TRUE, maxit = 1000))
# #method = "BFGS", lower = rep(-10,length(init_beta)), upper = rep(10, length(init_beta)))
# 
# print("Optimization ends")
# 
# print(output$par)

saveRDS(object = output, file = "wedidit2.RDS")
