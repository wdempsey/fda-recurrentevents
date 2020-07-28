## Source the functions and necessary packages
source("./basic_simulation_functions.R")
if(!require("glmnet")){install.packages('glmnet')}
library('glmnet')

### Matern covariance parameters
nu = 1/2 # controls smoothness
sigma_sq = 1 # controls overall variance
rho = 0.3 # controls correlation
constant = sigma_sq / (gamma(nu) * 2^(nu-1))

### User-day information
### Each user day is t \in [0,1]
### We use M dense equally spaced timepoints
M = 1000 # Num of time points
times = seq(0,1, length.out = M) # Equally spaced timepoints in [0,1]
gap = diff(times)[1] # What is assumed gap times
abs_diff = abs(outer(times, times, "-")) # Absolute distance 
x = sqrt(2 * nu) * abs_diff/ rho 
Cov_X = constant * x^nu * besselK(x, nu) # Matern Covariance
diag(Cov_X) = 1
C = chol(Cov_X) # Cholesky decomposition

## Generate corresponding event times
## Simulation assumes we want events near
## (A) Want Up then Down pattern trends -> beta_1t
## (B) Overall high mean in the window -> beta_2t
# beta_1t = 30*c(rep(-1, 11), rep(1, 22), rep(-1, 11))
# beta_1t = 50*sin(1:44/44*2*pi-pi/2)
# beta_1t = 20*exp(0.1*1:44)/mean(exp(0.1*1:44))
beta_1t = 100*sin(1:44/44*2*pi-pi/2)
# beta_2t = 50*rep(1, 44)
# plot(times[1:44],beta_1t+beta_2t)

## Collect coefficients for first 35
## Take covariance in the 44 timeslots
Sigma = Cov_X[1:length(beta_1t),1:length(beta_1t)]
eig_Sigma = eigen(Sigma)
K_x = 35  # Pick first 35 eigen-vectors
cumsum(eig_Sigma$values)[K_x]/sum(eig_Sigma$values) # Explains most of the variation

# set.seed("131312")
id = runif(1, min = 0, max = 100000)
base_num_events = 5
base_rate = logit(5/length(times))
max_sampling_rate = 12*4/1000 # Max sampling rate is 1 obs every 1/2 hour
dataset = generate_complete_data(N = 500, Cov_X, C, times, eig_Sigma, beta_1t, max_sampling_rate, base_rate)
rm(abs_diff, C, Cov_X, x)

dataset =data.frame(dataset)
names(dataset)[1] = "userday"
names(dataset)[2] = "Y"

agg_results = aggregate(Y~userday,dataset, sum)
agg_summary = c(mean(agg_results[,2]),var(agg_results[,2])) # Report the mean and variance of number of events per day 
write.table(agg_summary, file = "aggregate_summary.csv", append = TRUE)

thinning_rates = c(1/1, 1/2, 1/4, 1/8)

for(rates in thinning_rates) {
  print(paste("On rate", rates*max_sampling_rate))
  subdataset = subsample_dataset(dataset, thinning_rate = rates) 
  intermediate_step = construct_J(times, eig_Sigma, subdataset)
  output = runglmnet(max_sampling_rate*rates, subdataset, intermediate_step$w, intermediate_step$Basis, epsilon = 0.0001)
  write.table(t(output$betahat), file = "betahat.csv", row.names = F, col.names = F, append = T, sep = ",")
  write.table(output$runtime, file = "runtime.csv", row.names = F, col.names = F, append = T, sep = ",")
  write.table(rates, file = "samplingrate.csv", row.names = F, col.names = F, append = T, sep = ",")
  write.table(id, file = "ids.csv", row.names = F, col.names = F, append = T, sep = ",")
}

