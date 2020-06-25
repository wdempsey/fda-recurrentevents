## Source the functions
source("./basic_simulation_functions.R")

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

# set.seed("134163")
# set.seed("131312")
base_num_events = 5
base_rate = logit(5/length(times))
sampling_rate = 10/1000
dataset = generate_complete_data(N = 4000, Cov_X, C, times, eig_Sigma, beta_1t, sampling_rate, base_rate)
rm(abs_diff, C, Cov_X, x)

dataset =data.frame(dataset)
names(dataset)[1] = "userday"
names(dataset)[2] = "Y"

agg_results = aggregate(Y~userday,dataset, sum)
table(agg_results[,2])

## Construct the J matrix
K_b = 35
local_times = times[1:44]
num=K_b-3
qtiles <- seq(0, local_times[length(local_times)], length = num + 2)[-c(1, num + 2)]
knots <- quantile(local_times, qtiles)
## Basis = bs(t, kb)
Basis = cbind(1, local_times, local_times^2, sapply(knots, function(k) ((local_times - k > 0) * (local_times - k)) ^ 2))
Psi = t(eig_Sigma$vectors[,1:K_x])
Phi = Basis
J = Psi%*%Phi
model.matrix = dataset[,3:ncol(dataset)]
new.model.matrix = as.matrix(model.matrix)%*%J*gap
w = cbind(new.model.matrix)
X=as.matrix(w[,1:3])
Z=as.matrix(w[,4:dim(w)[2]])

n.tmp = length(dataset$Y)
p.tmp = ncol(w)

library("glmnet")
subsample_offset = rep(log(10/1000),nrow(dataset))
p.fac = rep(1, ncol(w))
p.fac[1:4] = 0 #no penalty on the first 4 variables
lambda_max <- 1/n.tmp
epsilon <- .0001
K <- 100
lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon), 
                            length.out = K)), digits = 10)

# set.seed("97139817")
# Start the clock!
ptm <- proc.time()

ridge.fit.cv <- cv.glmnet(w, dataset$Y, alpha = 0, intercept = TRUE, 
                          penalty.factor = p.fac, standardize = FALSE,
                          lambda = lambdapath, nfolds = 20,
                          family = "binomial")
# Stop the clock
proc.time() - ptm

ridge.fit.lambda <- ridge.fit.cv$lambda.min
# plot(ridge.fit.cv)

# Extract coefficient values for lambda.1se (without intercept)
ridge.coef <- (coef(ridge.fit.cv, s = ridge.fit.lambda))[-1]

betaHat.net <- Basis %*% ridge.coef
par(mar = c(4,2.5,1,1)+0.1)
plot(max(local_times)-local_times, betaHat.net, col = "blue", type = "l", 
     ylim = range(beta_1t,betaHat.net),
     axes = FALSE, xlab = "", ylab = "")
axis(side = 1); axis(side = 2, labels = FALSE)
# lines(local_times, betaHat.net.star, col = "red")
lines(local_times, beta_1t)
mtext("Time until event", side = 1, line = 2)
mtext(expression(paste(beta, "(s)")),side = 2, line = 1)

mise_calc(beta_1t, betaHat.net)
