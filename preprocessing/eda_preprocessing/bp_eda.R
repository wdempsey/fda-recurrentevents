## Code for the mean EDA plots backward from the button press

## ---- Please place this file in a folder named "methods" 
##      under the directory where the scaled EDA data is stored. ----

## Recent update: May 27, 2021, by Xinrui Wu

library(tidyverse)
library(lubridate)


setwd('..')


# read in the button press data and transform the data type of the timestamp
bp = readRDS('../R21_Study - tags/button_presses.RDS')
bp = transmute(bp, ID, time=as_datetime(ts))


#### Match button presses and the previous EDA (scaled) data ####
bp_eda = data.frame(ID=NULL, BP_num=NULL, BP_time=NULL, Device=NULL, 
                    EDA_time=NULL, EDA=NULL, EDA_num=NULL)

n_obs = 91 # the largest number of observations (files)
for (i in 1:n_obs){
  id = 1000 + i
  input_name = paste0('./', id, '_EDA_scaled.rds')
  
  ## check the progress of the loop and save temporary result ##
  if (i%%10 == 0){
    print(i)
    saveRDS(bp_eda, file = './methods/summary_data/bp_eda.rds')
  }
  
  ## match bp and EDA for this observation
  ## observations 1010 and 1012 should be deleted since they do not have ACC data
  if (file.exists(input_name) & i!=10 & i!=12){ 
    sub_bp = filter(bp, ID==id)
    eda = readRDS(input_name)
    for (j in 1:dim(sub_bp)[1]){
      bp_time = sub_bp$time[j]
      # EDA_time in [BP_time-30min, BP_time)
      sub_eda = filter(eda, timestamp<bp_time, timestamp>=(bp_time-30*60))
      if (dim(sub_eda)[1] > 0){
        sub_bp_eda = sub_eda %>%
          transmute(ID=rep(id,n()), BP_num=rep(j,n()), BP_time=rep(bp_time, n()), 
                    Device, EDA_time=timestamp, EDA=EDA_scaled) %>%
          group_by(Device) %>%
          mutate(EDA_num = 1:n()) %>%
          ungroup() %>%
          arrange(Device, EDA_num)
        bp_eda = rbind(bp_eda, sub_bp_eda)
      }
    }
  }
}
saveRDS(bp_eda, file = './methods/summary_data/bp_eda.rds')

#### Check the matched data ####
bp_eda_read = readRDS('./methods/summary_data/bp_eda.rds')
bp_eda_summary = bp_eda_read %>%
  group_by(ID, BP_num, Device) %>%
  summarise(EDA_entries = n())
max(bp_eda_summary$EDA_entries) # should be no larger than 7200 = 30*60*4


#### Draw the mean plot and the boxplot ####
bp_eda = bp_eda_read %>%
  group_by(ID, BP_num, Device) %>%
  mutate(time_to_bp = max(EDA_num)+1-EDA_num) %>%
  ungroup() %>%
  arrange(time_to_bp)
eda_mean_plot =  bp_eda %>%
  group_by(time_to_bp) %>%
  summarize(EDA_mean = mean(EDA), EDA_sd = sd(EDA), num_bp = n()) %>%
  ungroup()


# mean plot
plot(x=eda_mean_plot$time_to_bp, y=eda_mean_plot$EDA_mean, axes=F, cex = 0.25, 
     main = "Mean of EDA across observations", ylab = 'mean(EDA)', 
     xlab = "time before a button press / 0.25 secs")
axis(side = 1)
axis(side = 2)


# boxplot for time points in the first minute before a button press
sub_bp_eda = filter(bp_eda, ID!=1009, time_to_bp<=240)  # 1009 has trouble raw EDA, at about 240
boxplot(EDA~time_to_bp, data=sub_bp_eda, cex=0.08, axes = F, 
        main = "Boxplots of EDA across observations", ylab = 'EDA',
        xlab = "time before a button press / 0.25 secs")
axis(side = 1)
axis(side = 2)
axis(side = 2, at=10, col = 'blue')
abline(h = 10, col = 'blue')
