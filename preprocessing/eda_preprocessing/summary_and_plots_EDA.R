## Summary and plots on EDA and its pause time
library(tidyverse)
library(lubridate)

setwd('/mnt/turbo/SI_data/R21_Study - EDA - scaled/')

## ------------------ Summary for Pause time in EDA------------------ ##
n_obs = 91  # the largest number of observations (files)
td = 0.25 # time difference between 2 consecutive data

# summary of AI for each person
summary_EDA = data.frame(matrix(0, nrow = n_obs, ncol = 14))
colnames(summary_EDA) = c('ID', 'EDA_mean', 'EDA_sd', 'EDA_max', 'EDA_min', 
                         'EDA.25', 'EDA.50', 'EDA.75', 'EDA.95', 'EDA.99', 
                         'Days', 'Continuous Periods', 
                         'Pause_max(hour)', 'Pause_min(hour)')

#  pause and period (after the pause) length in seconds
all_pause = data.frame(ID=integer(), pause_num=integer(), pause_len=double(), 
                       period_len=double(), num_na=integer()) 

for (i in 1:n_obs){
  id = 1000 + i 
  print(paste("At patient ", id))
  input_name = paste0("./", id, "_EDA_scaled.rds")
  
  ## check the progress of the loop and save temporary result ##
  if (i%%10 == 0){
    print(i)
    saveRDS(all_pause, file = '../summary_data/pauses_EDA.rds')
  }
  
  if (file.exists(input_name)){
    eda = readRDS(input_name)
    summary_EDA[i, 'ID'] = id
    summary_EDA[i, 'EDA_mean'] = mean(eda$EDA_scaled)
    summary_EDA[i, 'EDA_sd'] = sd(eda$EDA_scaled)
    summary_EDA[i, 'EDA_max'] = max(eda$EDA_scaled)
    summary_EDA[i, 'EDA_min'] = min(eda$EDA_scaled)
    summary_EDA[i, 'EDA_max'] = max(eda$EDA_scaled)
    summary_EDA[i, c('EDA.25', 'EDA.50', 'EDA.75', 'EDA.95', 'EDA.99')] = 
      quantile(eda$EDA_scaled, probs = c(0.25, 0.5, 0.75, 0.95, 0.99))
    summary_EDA[i, 'Days'] = length(unique(as.Date(eda$timestamp)))
    # pause information
    timestamps = eda$timestamp[order(eda$timestamp)]
    time_diff = round(diff(as.numeric(timestamps)),3)
    cum_time = cumsum(time_diff)[(c(which(time_diff > td), length(time_diff)+1)-1)] - 
      c(0, cumsum(time_diff)[which(time_diff > td)])
    n_pause = length(which(time_diff > td))
    all_pause = rbind(all_pause, data.frame(ID=id, pause_num=0, pause=0, period=cum_time[1]))
    if (n_pause > 0){
      sub_pause = time_diff[which(time_diff > td)]
      all_pause = rbind(all_pause, data.frame(ID=rep(id, n_pause), pause_num=1:n_pause, 
                                              pause=sub_pause, period=cum_time[2:length(cum_time)]))
      summary_EDA[i, 'Pause_max(hour)'] = round(max(sub_pause)/3600, 4)
      summary_EDA[i, 'Pause_min(hour)'] = round(min(sub_pause)/3600, 4)
      summary_EDA[i, 'Continuous Periods'] = n_pause+1
    }
    else{
      summary_EDA[i, 'Pause_max(hour)'] = 0
      summary_EDA[i, 'Pause_min(hour)'] = 0
    }
  } else{
    none_ind = c(none_ind, i)
  }
}
saveRDS(all_pause, file = '../summary_data/pauses_EDA.rds')

# summary of pause time for each observation
summary_pause_EDA = all_pause %>%
  group_by(ID) %>%
  summarise(num_less_than_1min = sum(pause_num > 0 & pause_len < 60), 
            num_1min_to_30min = sum(pause_num > 0 & pause_len >= 60 & pause_len < 30*60), 
            num_30min_to_1day = sum(pause_num > 0 & pause_len >= 30*60 & pause_len < 60*60*24),
            num_more_than_1day = sum(pause_num > 0 & pause_len >= 60*60*24), 
            num_total = sum(pause_num > 0))
saveRDS(summary_pause_EDA, file = './methods/summary_data/summary_pause_EDA.rds')


## ---------- Frequency of missing EDA data in the 30 minuetes before button presses ---------- ##

time_with_miss_30min = all_pause %>%
  mutate(time_with_miss = period_len * (period_len <= 30*60 & pause_num > 0) + 
           rep(30*60, dim(all_pause)[1]) * (period_len > 30*60 & pause_num > 0), 
         ave_pause=0)

for (i in 1:dim(time_with_miss_30min)[1]){
  time_with_miss = time_with_miss_30min$time_with_miss[i]
  pause = time_with_miss_30min$pause_len[i]
  if(time_with_miss > 0){
    if(pause>=1800){
      time_with_miss_30min$ave_pause[i] = (1800+(1800-time_with_miss))/2
    }
    else{
      if((1800-pause)>=time_with_miss){
        time_with_miss_30min$ave_pause[i] = pause
      }
      else{
        time_with_miss_30min$ave_pause[i] = ((1800-pause)*pause + 
                                               (pause+(1800-time_with_miss))*(time_with_miss+pause-1800)/2)/time_with_miss
      }
    }
  }
}

ave_pause = time_with_miss_30min %>%
  group_by(ID) %>%
  summarise(miss_ave = sum(ave_pause*time_with_miss)/max(sum(time_with_miss),1), 
            all_ave = sum(ave_pause*time_with_miss)/sum(period_len)) %>%
  mutate(miss_ave = round(miss_ave/60,2), all_ave = round(all_ave/60,2))

prob_miss_30min_obs = time_with_miss_30min %>%
  group_by(ID) %>%
  summarize(prob = sum(time_with_miss) / sum(period_len))
prob_miss_30min_all = sum(time_with_miss_30min$time_with_miss) / sum(time_with_miss_30min$period_len)
prob_miss_30min_all



## ------------------ Pause in AI and EDA------------------ ##
pauses_ai = readRDS('../R21_Study - ACC - AI/methods/summary_data/pauses_AI.rds')
pauses_eda = readRDS('./methods/summary_data/pauses_EDA.rds')

par(mfrow=c(17, 7), mar = c(1,1,1,1))
for (i in 1:n_obs){
  id = 1000 + i
  if (length(which(pauses_ai$ID==id)) > 0){
    sub_ai = pauses_ai[which(pauses_ai$ID==id),] %>%
      transmute(pause, period)
    ai = cumsum(matrix(t(as.matrix(sub_ai)), ncol = 1))
    sub_eda = pauses_eda[which(pauses_eda$ID==id),] %>%
      transmute(pause=pause_len, period=period_len)
    eda = cumsum(matrix(t(as.matrix(sub_eda)), ncol = 1))
    
    plot(x=1:ai[2], y=rep(1, ai[2]), type='l', col='red', lwd=1, 
         xlim = c(0, max(c(max(ai), max(eda)))), ylim=c(0.8, 2.2), axes = F, 
         xlab='time with data (secs)', ylab='', main=id)
    #axis(side = 1, labels = F)
    if (dim(sub_ai)[1] > 1){
      for (j in 2:dim(sub_ai)[1]){
        points(x=(ai[2*j-1]+1):ai[2*j], y=rep(1, (ai[2*j]-ai[2*j-1])), type='l', col='red', lwd=1)
      }
    }
    points(x=1:eda[2], y=rep(2, eda[2]), type='l', col='blue', lwd=1)
    if (dim(sub_eda)[1] > 1){
      for (j in 2:dim(sub_eda)[1]){
        points(x=(eda[2*j-1]+1):eda[2*j], y=rep(2, (eda[2*j]-eda[2*j-1])), type='l', col='blue', lwd=1)
      }
    }
  }
}



## ------------------ Plots ------------------ ##

## pause time ##
hist(all_pause$pause_len)
boxplot((pause_len/3600)~ID, data = all_pause, ylab='Pause Time (EDA) / hour')


df = t(as.matrix(summary_pause_EDA[,2:5]))
barplot(df, col=c('red', 'blue', 'yellow', 'grey'),
        legend = c('< 1min', '[1min, 30mins)', '[30mins, 1day)', '>= 1day'), 
        xlab = 'ID', ylab = 'Count', main = "Number of Pauses (EDA)")
par(mfrow = c(2, 2), mgp = c(1,1,0))
barplot(df[1,], col='red', xlab = 'ID', main = 'Pauses < 1 min (EDA)')
barplot(df[2,], col='blue', xlab = 'ID', main = 'Pauses in [1min, 30mins) (EDA)')
barplot(df[3,], col='yellow', xlab = 'ID', main = 'Pauses in [30mins, 1day) (EDA)')
barplot(df[4,], col='grey', xlab = 'ID', main = 'Pauses > 1 day (EDA)')
dev.off()


## probability of having missing data in the past 30 mins for each observed time ##
barplot(prob_miss_30min_obs$prob, type = 'l', xlab='ID', ylab='Probability', 
        main = 'Prob. of having missing data in the past 30 mins (EDA)')
hist(prob_miss_30min_obs$prob, xlab='Probability', breaks = seq(from=0, to=0.15, by = 0.01), 
     main = 'Histogram of the Prob. of having missing data in the past 30 mins (EDA)')

## average missing time ##
par(mfrow=c(1,2))
hist(ave_pause$miss_ave, breaks = 50, xlab='Average missing time/(min) (across individuals)', 
     main='Avg. from time points with missing-EDA')
hist(ave_pause$all_ave, breaks = 50, xlab='Average missing time/(min) (across individuals)', 
     main='Avg. from all observed time points-EDA')
dev.off()







