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
beta_1t = 30*exp(0.2*1:44)/mean(exp(0.2*1:44))
setting = "exponential"
# beta_1t = 100*sin(1:44/44*2*pi-pi/2)
# setting = "sine"

# test if there is at least one argument: if not, return an error
print("Made it to window length")
args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  print("Using default value!")
  window_length = 44
} else if (length(args)>=1) {
  # default output file
  window_length = as.numeric(args)
}
cat(paste0("Window length is ", window_length, "\n"))

## Collect coefficients for first 35
## Take covariance in the 44 timeslots
Sigma = Cov_X[1:length(beta_1t),1:length(beta_1t)]
eig_Sigma = eigen(Sigma)
K_x = 35  # Pick first 35 eigen-vectors
cumsum(eig_Sigma$values)[K_x]/sum(eig_Sigma$values) # Explains most of the variation

arrayid=Sys.getenv("SLURM_ARRAY_TASK_ID")
# arrayid = 10
print(paste("Current ARRAY TASK ID", arrayid))
#your array of seeds in x
allseeds = readRDS("basic_simulation_seeds.RDS")
# select the ith index from the seed array
seed=allseeds[as.numeric(arrayid)]

set.seed(seed)
id = runif(1, min = 0, max = 100000)
base_num_events = 5
base_rate = logit(5/length(times))
max_sampling_rate = 12*12/1000 # Max sampling rate is 1 obs every 5 minutes
dataset = generate_complete_data(N = 500, Cov_X, C, times, eig_Sigma, beta_1t, max_sampling_rate, base_rate)
rm(abs_diff, C, Cov_X, x)

dataset =data.frame(dataset)
names(dataset)[1] = "userday"
names(dataset)[2] = "Y"

agg_results = aggregate(Y~userday,dataset, sum)
agg_summary = matrix(c(mean(agg_results[,2]),var(agg_results[,2])), nrow = 1, ncol = 2) # Report the mean and variance of number of events per day
# write.table(agg_summary, file = "./output_csv/aggregate_summary.csv",  row.names = F, col.names = F, append = T, sep = ",")

thinning_rates = c(1/1, 1/3, 1/3* c(1/2, 1/4, 1/8))
col_names = as.vector(c("ids", "runtime", "rates", "mean_count", "var_count", paste0(rep("beta",44), as.character(1:44))))

for(rates in thinning_rates) {
  print(paste("On rate", rates*max_sampling_rate))
  subdataset = subsample_dataset(dataset, thinning_rate = rates)
  intermediate_step = construct_J(times, eig_Sigma, subdataset)
  output = runglmnet(max_sampling_rate*rates, subdataset, intermediate_step$w, intermediate_step$Basis, epsilon = 0.0001)
  results = matrix(c(id, output$runtime, rates, agg_summary, output$beta), nrow = 1)
  colnames(results) = col_names
  if(!file.exists(paste("./output_csv/", setting, "_results.csv", sep = ""))) {
    write.table(results, file = paste("./output_csv/", setting, "_results.csv", sep = ""), row.names = F, col.names = T, append = T, sep = ",")
  } else {
    write.table(results, file = paste("./output_csv/", setting, "_results.csv", sep = ""), row.names = F, col.names = F, append = T, sep = ",")
  }
}

