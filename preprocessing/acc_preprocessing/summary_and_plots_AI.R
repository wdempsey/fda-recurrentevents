# Summary and plots on AI and pause time

## ---- Please place this file in a folder named "methods" 
##      under the directory where the scaled AI data is stored. ----

## Recent update: May 27, 2021, by Xinrui Wu

library(tidyverse)
library(lubridate)

setwd('..')

## ------------------ Summary for AI and Pause time------------------ ##
n_obs = 91  # the largest number of observations (files)
epch = 10   # the length of the time window (in seconds) used to compute AI

# summary of AI for each person
summary_AI = data.frame(matrix(0, nrow = n_obs, ncol = 14))
colnames(summary_AI) = c('ID', 'AI_mean', 'AI_sd', 'AI_max', 'AI_min', 
                         'AI.25', 'AI.50', 'AI.75', 'AI.95', 'AI.99', 
                         'Days', 'Continuous Periods', 
                         'Pause_max(hour)', 'Pause_min(hour)')
# pause and period (after the pause) length in seconds
all_pause = data.frame(ID=integer(), pause_num=integer(), pause=double(), period=double()) 


none_ind = 0
for (i in 1:n_obs){
  id = 1000 + i 
  input_name = paste0("./", id, "_AI.rds")
  
  ## check the progress of the loop##
  if (i%%10 == 0){
    print(i)
  }
  
  if (file.exists(input_name)){
    acc_ai = readRDS(input_name)
    summary_AI[i, 'ID'] = id
    summary_AI[i, 'AI_mean'] = mean(acc_ai$AI)
    summary_AI[i, 'AI_sd'] = sd(acc_ai$AI)
    summary_AI[i, 'AI_max'] = max(acc_ai$AI)
    summary_AI[i, 'AI_min'] = min(acc_ai$AI)
    summary_AI[i, 'AI_max'] = max(acc_ai$AI)
    summary_AI[i, c('AI.25', 'AI.50', 'AI.75', 'AI.95', 'AI.99')] = 
      quantile(acc_ai$AI, probs = c(0.25, 0.5, 0.75, 0.95, 0.99))
    summary_AI[i, 'Days'] = length(unique(as.Date(acc_ai$timestamp)))
    summary_AI[i, 'Continuous Periods'] = length(unique(acc_ai$Period))
    
    # pause information
    time_diff = round(diff(as.numeric(acc_ai$timestamp)))
    cum_time = cumsum(time_diff)[(c(which(time_diff > epch), length(time_diff)+1)-1)] - 
      c(0, cumsum(time_diff)[which(time_diff > epch)])
    n_pause = length(which(time_diff > epch))
    all_pause = rbind(all_pause, data.frame(ID=id, pause_num=0, pause=0, period=cum_time[1]))
    if (n_pause > 0){
      sub_pause = time_diff[which(time_diff > epch)]
      all_pause = rbind(all_pause, data.frame(ID=rep(id, n_pause), pause_num=1:n_pause, 
                                              pause=sub_pause, period=cum_time[2:length(cum_time)]))
      summary_AI[i, 'Pause_max(hour)'] = round(max(sub_pause)/3600, 2)
      summary_AI[i, 'Pause_min(hour)'] = round(min(sub_pause)/3600, 2)
    }
    else{
      summary_AI[i, 'Pause_max(hour)'] = 0
      summary_AI[i, 'Pause_min(hour)'] = 0
    }
  }
  
  else{
    none_ind = c(none_ind, i)
  }
}

# remove rows in "summary_AI" where there are no observation data
none_ind = none_ind[-1]
if(length(none_ind) > 0) {summary_AI = summary_AI[-none_ind,]}

saveRDS(summary_AI, file = './methods/summary_data/summary_AI.rds')
saveRDS(all_pause, file = './methods/summary_data/pauses_AI.rds')

# summary of pause time for each observation
summary_pause = all_pause %>%
  group_by(ID) %>%
  summarise(num_less_than_1min = sum(pause_num > 0 & pause < 60), 
            num_1min_to_30min = sum(pause_num > 0 & pause >= 60 & pause < 30*60), 
            num_30min_to_1day = sum(pause_num > 0 & pause >= 30*60 & pause < 60*60*24),
            num_more_than_1day = sum(pause_num > 0 & pause >= 60*60*24), 
            num_total = sum(pause_num > 0))
saveRDS(summary_pause, file = './methods/summary_data/summary_pause_AI.rds')


## ---------- Frequency of missing AI data in the 30 minuetes before button presses ---------- ##
all_pause = readRDS('./methods/summary_data/pauses_AI.rds')

time_with_miss_30min = transmute(all_pause, ID, pause_num, period) %>%
  mutate(time_with_miss = period * (period <= 30*60 & pause_num > 0) + 
           rep(30*60, dim(all_pause)[1]) * (period > 30*60 & pause_num > 0))

# note: "frequency of missing" = 
#       "real number of data entries" / "number of data entries there should be in 30 minuetes"
# frequency of missing for each observation

prob_miss_30min_obs = time_with_miss_30min %>%
  group_by(ID) %>%
  summarize(prob = sum(time_with_miss) / sum(period))

# frequency of missing for all observations
prob_miss_30min_all = sum(time_with_miss_30min$time_with_miss) / sum(time_with_miss_30min$period)
print(prob_miss_30min_all)


## ---------- Summary of devices ---------- ##
num_devices = data.frame(ID=integer(), num_device=integer()) # number of devices for each observation

n_obs = 91
for (i in 1:n_obs){
  id = 1000 + i 
  input_name = paste0("./", id, "_AI.rds")
  
  if (file.exists(input_name)){
    acc_ai = readRDS(input_name)
    temp = data.frame(ID = id, num_device = length(unique(acc_ai$Device)))
    num_devices = rbind(num_devices, temp)
  }
}
saveRDS(num_devices, file = './methods/summary_data/num_devices_AI.rds')

table(num_devices$num_device)




## ------------------ Plots ------------------ ##
## summary of AI ##
#### Histogram of AI mean ####
hist(summary_AI$AI_mean, breaks = 25, axes = F, xlab = '', ylab = '', 
     main = 'Histogram of mean(AI) across Observations')
axis(side = 1)
axis(side = 2)


#### AI quantiles for each obs. ####
par(xpd=T, mar=c(5,5,5,7))
ymax = max(summary_AI$AI_max)
ymin = min(summary_AI$AI_min)
plot(x = (summary_AI$ID - 1000), y = summary_AI$AI_mean, type = 'l', axes = F, ylim = c(ymin-50, ymax+50), 
     xlab = 'Observations', ylab = 'AI', main = 'Summary of AI acorss Observations')
axis(side = 1)
axis(side = 2)
points(x = (summary_AI$ID - 1000), y = summary_AI$AI_max, type = 'l', lty = 1, col = 'grey')
points(x = (summary_AI$ID - 1000), y = summary_AI$AI_min, type = 'l', lty = 1, col = 'grey')
points(x = (summary_AI$ID - 1000), y = summary_AI$AI.25, type = 'l', lty = 3, col = 'blue')
points(x = (summary_AI$ID - 1000), y = summary_AI$AI.50, type = 'l', lty = 3, col = 'red')
points(x = (summary_AI$ID - 1000), y = summary_AI$AI.75, type = 'l', lty = 2, col = 'blue')
points(x = (summary_AI$ID - 1000), y = summary_AI$AI.95, type = 'l', lty = 4, col = 'red')
points(x = (summary_AI$ID - 1000), y = summary_AI$AI.99, type = 'l', lty = 4, col = 'blue')
text(x = (summary_AI$ID - 1000), y = rep(ymax+20, dim(summary_AI)[1]), 
     labels = summary_AI$Days, col = "dimgrey", cex = .5)
text(x = (summary_AI$ID - 1000), y = rep(ymax+45, dim(summary_AI)[1]), 
     labels = summary_AI$`Continuous Periods`, col = "pink", cex = .5)
legend(95,900, legend=c("AI_max", "AI_min", "AI_mean", 'Days', 'Periods'), 
       col=c('red', 'blue', 'black', 'dimgrey', 'pink'), box.col = "white", 
       horiz=F, lty=c(2,2,1,1,1), cex=0.8)
legend(91,900, legend=c("mean", "maximum", "minimum", 
                        '25% quantile', '50% quantile', '75% quantile', 
                        '95% quantile', '99% quantile'), 
       col=c('black', 'grey', 'grey','blue','red', 'blue','red','blue'), box.col = "white", 
       horiz=F, lty=c(1,1,1,3,3,2,4,4), cex=0.8)
legend(91,100, legend=c("mean", 
                        '25% quantile', '50% quantile', '75% quantile'), 
       col=c('black', 'red', 'blue','dimgrey'), box.col = "white", 
       horiz=F, lty=c(1,3,3,3), cex=0.8)
dev.off()



#### AI histogram for each obs. ####
par(mfrow = c(15, 6), mar=c(1.5,2,1.5,0.5))
n_obs = 91
for (i in 1:n_obs){
  id = 1000 + i 
  input_name = paste0("./", id, "_AI.rds")
  
  if (file.exists(input_name)){
    acc_ai = readRDS(input_name)
    hist(acc_ai$AI, breaks=50, freq=F, xlab='AI', main = id)
    abline(v = mean(acc_ai$AI), col='blue')
    abline(v = quantile(acc_ai$AI, 0.95), col='red')
    abline(v = quantile(acc_ai$AI, 0.99), col='grey')
  }
}
dev.off()




## summary of pause time ##
#### boxplot ####
hist(all_pause$pause)
boxplot((pause/3600)~ID, data = all_pause, ylab='Pause Time / hour')

#### barplot ####
df = t(as.matrix(summary_pause[,2:5]))
barplot(df, col=c('red', 'blue', 'yellow', 'grey'),
        legend = c('< 1min', '[1min, 30mins)', '[30mins, 1day)', '>= 1day'), 
        xlab = 'ID', ylab = 'Count', main = "Number of Pauses")
par(mfrow = c(2, 2), mgp = c(1,1,0))
barplot(df[1,], col='red', xlab = 'ID', main = 'Pauses < 1 min')
barplot(df[2,], col='blue', xlab = 'ID', main = 'Pauses in [1min, 30mins)')
barplot(df[3,], col='yellow', xlab = 'ID', main = 'Pauses in [30mins, 1day)')
barplot(df[4,], col='grey', xlab = 'ID', main = 'Pauses > 1 day')
dev.off()


## probability of having missing data in the past 30 mins for each observed time ##
barplot(prob_miss_30min_obs$prob, type = 'l', xlab='ID', ylab='Probability', 
        main = 'Prob. of having missing data in the past 30 mins')
hist(prob_miss_30min_obs$prob, xlab='Probability', breaks = seq(from=0, to=0.15, by = 0.01), 
     main = 'Histogram of the Prob. of having missing data in the past 30 mins')
