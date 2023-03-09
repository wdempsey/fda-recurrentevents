source("./basic_simulation_functions.R")
setting = "sine"
num_simulations = 500
seq_windowlengths = seq(36,52,by=4)

df36 = read.csv(x = df, file = paste("./output_csv/",setting,"_simdelta_allresults_wl_36.csv", sep = ""))
df40 = read.csv(x = df, file = paste("./output_csv/",setting,"_simdelta_allresults_wl_40.csv", sep = ""))
df44 = read.csv(x = df, file = paste("./output_csv/",setting,"_simdelta_allresults_wl_44.csv", sep = ""))
df48 = read.csv(x = df, file = paste("./output_csv/",setting,"_simdelta_allresults_wl_48.csv", sep = ""))
df52 = read.csv(x = df, file = paste("./output_csv/",setting,"_simdelta_allresults_wl_52.csv", sep = ""))

if(setting == "sine") {
  beta_truth = 100*sin(1:44/44*2*pi-pi/2)
} else if (setting == "exponential") {
  beta_truth = 30*exp(0.2*1:44)/mean(exp(0.2*1:44))
}

df = df[order(df$ids, df$rates),]

M = 1000 # Num of time points
times = seq(0,1, length.out = M) # Equally spaced timepoints in [0,1]
gap = diff(times)[1] # What is assumed gap times

true_variance = mise_calc(beta_truth, beta_truth*0)

unique_rates = unique(df$rates)
unique_rates = sort(unique_rates, decreasing = T)
betameans = matrix(nrow = length(unique_rates), ncol = length(6:ncol(df)))
for (i in 1:length(unique_rates)) {
  betameans[i,] = colMeans(df[df$rates == unique_rates[i],6:ncol(df)], na.rm = TRUE)
}

## Simple visualization to show that means are same across subsamples
## This is a good thing!  Means there is no bias due to subsampling 
## The bias inherent in even most data case is due to small n and penalization
par(mar = c(3,3,2,1)+0.1)
plot(1:length(beta_truth), beta_truth, type = "l", axes = F, xlab = "Time", ylab = "")
axis(side = 1, labels = F)
axis(side = 2, labels = F)
for(i in 1:length(unique_rates)) {
  lines(as.numeric(betameans[i,]), type = "l", col = "darkgrey")
}

sqbias = vector(length = length(unique_rates))
for(i in 1:length(unique_rates)) {
  sqbias[i] = mise_calc(beta_1t = beta_truth, betaHat.net = betameans[i,])
}

avg_runtime = rep(0,length(unique_rates))
for(i in 1:length(unique_rates)) {
  avg_runtime[i] = mean(df[df$rates == unique_rates[i],]$runtime)
}

variance = rep(0,length(unique_rates))
mise = rep(0,length(unique_rates))
for(i in 1:length(unique_rates)) {
  subsetdf = df[df$rates == unique_rates[i],6:ncol(df)]
  for (row in 1:nrow(subsetdf)) {
    variance[i] = variance[i] + as.numeric(mise_calc(beta_1t = betameans[i,], betaHat.net = subsetdf[row,])/nrow(subsetdf))
    mise[i] = mise[i] + as.numeric(mise_calc(beta_1t = beta_truth, betaHat.net = subsetdf[row,])/nrow(subsetdf))
  }
}

subsample_variance = rep(0,length(unique_rates))
subsample_variance_count = rep(0,length(unique_rates))
unique_ids = unique(df$ids)
for (id in unique_ids) {
  subset_df = df[df$ids == id,]
  if(any(subset_df$rates==1.00)) {
    beta_baseline = subset_df[subset_df$rates==1.00, 6:ncol(subset_df)]
    for (row in 1:nrow(subset_df)){
      which_rate = which(subset_df[row,]$rates == unique_rates)
      temp_error = as.numeric(mise_calc(beta_1t = beta_baseline, betaHat.net = subset_df[row,6:ncol(subset_df)])/nrow(subsetdf))
      subsample_variance_count[which_rate] = subsample_variance_count[which_rate] + 1
      subsample_variance[which_rate] = temp_error/subsample_variance_count[which_rate] + (subsample_variance_count[which_rate]-1)/subsample_variance_count[which_rate] * subsample_variance[which_rate]
    }
  }
}

table = rbind(unique_rates, subsample_variance/true_variance, variance/true_variance, sqbias/true_variance, mise/true_variance, avg_runtime)

round(table, 5)
mise/mise[1]