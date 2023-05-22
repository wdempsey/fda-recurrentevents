## Code for the mean AI plots 30 minuets backward from the button press
library(tidyverse)
library(lubridate)

setwd('Z:/SI_data/R21_Study - ACC - AI/')

# read in the button press data and transform the data type of the time stamp
bp = readRDS('../R21_Study - tags/button_presses.RDS')
bp = transmute(bp, ID, time=as_datetime(ts))


#### Match button presses and the previous AI data ####
# data frame to store the matched bp and AI
if(!file.exists("../summary_data/bp_ai.rds")) {
  bp_ai = data.frame(ID=NULL, BP_num=NULL, BP_time=NULL, 
                     Device=NULL, AI_num=NULL, AI_time=NULL, AI=NULL)
  
  n_obs = 91 # the largest number of observations (files)
  for (i in 1:n_obs){
    print(paste("Current participant", i))
    id = 1000 + i
    input_name = paste0('./', id, '_AI.rds')
    
    ## check the progress of the loop and save temporary result ##
    if (i%%10 == 0){
      print(i)
      saveRDS(bp_ai, file = './methods/summary_data/bp_ai.rds')
    }
    
    ## match bp and AI for this observation
    if (file.exists(input_name)){
      sub_bp = filter(bp, ID==id)
      ai = readRDS(input_name)
      device = unique(ai$Device)
      
      # process the data from different devices seperately due to the timestamp overlap 
      for (d in 1:length(device)){
        sub_ai =  filter(ai, Device==device[d])
        for (j in 1:dim(sub_bp)[1]){
          bp_time = sub_bp$time[j]
          # AI_time in [BP_time-30min, BP_time)
          ai_per_device_bp = filter(sub_ai, timestamp<bp_time, timestamp>=(bp_time-30*60))
          if (dim(ai_per_device_bp)[1] > 0){
            sub_bp_ai = transmute(ai_per_device_bp, ID=rep(id,n()), 
                                  BP_num=rep(j,n()), BP_time=rep(bp_time, n()), 
                                  Device=rep(device[d],n()), AI_num=1:n(), AI_time=timestamp, AI)
            bp_ai = rbind(bp_ai, sub_bp_ai)
          }
        }
      }
    }
  }
  saveRDS(bp_ai, file = '../summary_data/bp_ai.rds')
} else {
  bp_ai = readRDS(file = '../summary_data/bp_ai.rds')
}

#### Check the matched data ####
bp_ai_summary = bp_ai %>%
  group_by(ID, BP_num, Device) %>%
  summarise(AI_entries = n())
max(bp_ai_summary$AI_entries) # should be no larger than 180 = 30*60/epch

#### Draw the mean plot and the boxplot ####
bp_ai = bp_ai %>%
  group_by(ID, BP_num, Device) %>%
  mutate(time_to_bp = max(AI_num)+1-AI_num) %>%
  ungroup() %>%
  arrange(time_to_bp)
ai_mean_plot =  bp_ai %>%
  group_by(time_to_bp) %>%
  summarize(AI_mean = mean(AI), num_bp = n()) %>%
  ungroup()

ai_mean_plot$time_to_bp = ai_mean_plot$time_to_bp * 10/60 # switch to minutes

png("../figures/smoothed_ai.png",
    width = 480, height = 480, units = "px", pointsize = 12)

par(mfrow = c(1,1), mar = c(5,4,1,1)+0.1)
# mean plot
plot(x=ai_mean_plot$time_to_bp, y=ai_mean_plot$AI_mean, axes=F, cex = 0.5, 
     # main = "Mean of AI across observations", 
     ylab = 'Average Activity Index (AI)', 
     xlab = "Time before Button Press (in minutes)",
     ylim = c(55,100))
axis(side = 1)
axis(side = 2)
smoothed_mean = loess(AI_mean ~ time_to_bp, data = ai_mean_plot, span = 1/2)
lines(ai_mean_plot$time_to_bp, smoothed_mean$fitted, col= "red", lwd = 2)

# boxplot
# boxplot(AI~time_to_bp, data=bp_ai, cex=0.08, axes = F, 
#         main = "Boxplots of AI across observations", ylab = 'AI',
#         xlab = "time before a button press / 10 secs")
# axis(side = 1)
# axis(side = 2)

dev.off()
