## Summary and plots on EDA and its pause time
library(tidyverse)
library(lubridate)
library(ggplot2)

## Windows
setwd('Z:/SI_data/R21_Study - EDA - scaled/')
## Linux
setwd('/mnt/turbo/SI_data/R21_Study - EDA - scaled/')

## ------------------ Summary for Pause time in EDA------------------ ##
n_obs = 91  # the largest number of observations (files)
td = 0.50 # time difference between 2 consecutive data

# summary of AI for each person
summary_EDA = data.frame(matrix(0, nrow = n_obs, ncol = 14))
colnames(summary_EDA) = c('ID', 'EDA_mean', 'EDA_sd', 'EDA_max', 'EDA_min', 
                         'EDA.25', 'EDA.50', 'EDA.75', 'EDA.95', 'EDA.99', 
                         'Days', 'Continuous Periods', 
                         'Pause_max(hour)', 'Pause_min(hour)')

#  pause and period (after the pause) length in seconds
all_pause = data.frame(ID=integer(), pause_num=integer(), pause_len=double(), 
                       period_len=double(), num_na=integer()) 
none_ind = 0
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
      quantile(eda$EDA_scaled, probs = c(0.25, 0.5, 0.75, 0.95, 0.99), na.rm = T)
    summary_EDA[i, 'Days'] = length(unique(as.Date(eda$timestamp)))
    # pause information
    timestamps = eda$timestamp[order(eda$timestamp)]
    time_diff = round(diff(as.numeric(timestamps)),3)
    cumulative_timediff = cumsum(time_diff)
    gap_obs = which(time_diff > td)
    cum_time = cumulative_timediff[(c(gap_obs, length(time_diff)+1)-1)] - 
      c(0,cumulative_timediff[gap_obs])
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
  summarise(num_less_than_1min = sum(pause_num > 0 & pause < 60), 
            num_1min_to_30min = sum(pause_num > 0 & pause >= 60 & pause < 30*60), 
            num_30min_to_1day = sum(pause_num > 0 & pause >= 30*60 & pause < 60*60*24),
            num_more_than_1day = sum(pause_num > 0 & pause >= 60*60*24), 
            num_total = sum(pause_num > 0))
saveRDS(summary_pause_EDA, file = '../summary_data/summary_pause_EDA.rds')


## ---------- Frequency of missing EDA data in the 30 minutes before button presses ---------- ##
time_with_miss_30min = transmute(all_pause, ID, pause_num, period) %>%
  mutate(time_with_miss = period * (period <= 30*60 & pause_num > 0) + 
           rep(30*60, dim(all_pause)[1]) * (period > 30*60 & pause_num > 0))

time_with_miss_30min$frac_missing = (time_with_miss_30min$time_with_miss / time_with_miss_30min$period)
time_with_miss_30min = subset(time_with_miss_30min, !is.na(frac_missing))

# note: "frequency of missing" = 
#       "real number of data entries" / "number of data entries there should be in 30 minuetes"
# frequency of missing for each observation
prob_miss_30min_obs = time_with_miss_30min %>%
  group_by(ID) %>%
  summarize(prob = sum(time_with_miss) / sum(period))


# frequency of missing for all observations
prob_miss_30min_all = sum(time_with_miss_30min$time_with_miss) / sum(time_with_miss_30min$period)
print(prob_miss_30min_all)

#### EDA quantiles for each obs. ####
summary_EDA = summary_EDA[summary_EDA$ID > 0,]
png("C:/Users/Balthazar/Documents/GitHub/fda-recurrentevents/figures/eda_summary.png",
    width = 720, height = 480, units = "px", pointsize = 18)
par(xpd=T, mar=c(4,4,1,7))
ymax = max(summary_EDA$EDA.99)
ymin = min(summary_EDA$EDA.25)
plot(x = (summary_EDA$ID - 1000), y = summary_EDA$EDA_mean, type = 'l', axes = F, ylim = c(ymin*0.8, ymax*1.2), 
     xlab = '', ylab = '')
axis(side = 1)
axis(side = 2)
points(x = (summary_EDA$ID - 1000), y = summary_EDA$EDA.25, type = 'l', lty = 3, col = 'blue')
points(x = (summary_EDA$ID - 1000), y = summary_EDA$EDA.50, type = 'l', lty = 3, col = 'red')
points(x = (summary_EDA$ID - 1000), y = summary_EDA$EDA.75, type = 'l', lty = 2, col = 'blue')
points(x = (summary_EDA$ID - 1000), y = summary_EDA$EDA.95, type = 'l', lty = 4, col = 'red')
points(x = (summary_EDA$ID - 1000), y = summary_EDA$EDA.99, type = 'l', lty = 4, col = 'blue')
text(x = (summary_EDA$ID - 1000), y = rep(ymax*1.1, dim(summary_EDA)[1]), 
     labels = summary_EDA$Days, col = "black", cex = .5)
text(x = (summary_EDA$ID - 1000), y = rep(ymax*1.15, dim(summary_EDA)[1]), 
     labels = summary_EDA$`Continuous Periods`, col = "grey", cex = .5)
# legend(95,900, legend=c("AI_max", "AI_min", "AI_mean", 'Days', 'Periods'), 
# col=c('red', 'blue', 'black', 'dimgrey', 'pink'), box.col = "white", 
# horiz=F, lty=c(2,2,1,1,1), cex=0.8)
legend(91,4, legend=c("mean", "maximum", "minimum", 
                        '25% quantile', '50% quantile', '75% quantile', 
                        '95% quantile', '99% quantile'), 
       col=c('black', 'grey', 'grey','blue','red', 'blue','red','blue'), box.col = "white", 
       horiz=F, lty=c(1,1,1,3,3,2,4,4), cex=0.8)
mtext("EDA", side = 2, line = 2)
mtext("User ID", side = 1, line = 2)
dev.off()

## probability of having missing data in the past 30 mins for each observed time ##
png("C:/Users/Balthazar/Documents/GitHub/fda-recurrentevents/figures/missing_data_eda.png",
    width = 480, height = 480, units = "px", pointsize = 12)
par(mar = c(4, 4, 1, 1), mfrow = c(1,1))
qplot(prob_miss_30min_obs$prob, 
      geom="histogram", xlab = "Patient Average Fraction of Missing EDA data in 30-minutes",
      bins = 20) 
dev.off()
