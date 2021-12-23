rm(list=ls())
setwd("Z:/SI_data/")
library("lubridate")
col_names = c("id", "start_ts", "end_ts")
buttonpresses = readRDS("./R21_Study - tags/button_presses.RDS")
ID_list = unique(buttonpresses$ID)
buttonpresses$timestamp = as_datetime(buttonpresses$ts)
buttonpresses$timesincebaseline = NA
buttonpresses$endtime = NA
if(!file.exists("./R21_Study - tags/button_presses_extrainfo.RDS")) {
  for(id in ID_list) {
    print(paste("On participant",id))
    if (file.exists(paste("./R21_Study - EDA - scaled3/",id,"_EDA_scaled.rds", sep = ""))) {
      edaX <- readRDS(paste("./R21_Study - EDA - scaled3/",id,"_EDA_scaled.rds", sep = ""))
      eda_start_ts = min(edaX$timestamp)
      eda_end_ts = max(edaX$timestamp)
      rm(edaX)
    } else {
      eda_start_ts = NA
      eda_end_ts = NA
    }
    if (file.exists(paste("./R21_Study - ACC - AI/",id,"_AI.rds", sep = ""))) {
      accX <- readRDS(paste("./R21_Study - ACC - AI/",id,"_AI.rds", sep = ""))
      acc_start_ts = min(accX$timestamp)
      acc_end_ts = max(accX$timestamp)
      rm(accX)
    } else {
      acc_start_ts = NA
      acc_end_ts = NA
    }
    start_ts = as_datetime(min(acc_start_ts, eda_start_ts, na.rm= T))
    end_ts = as_datetime(max(acc_end_ts, eda_end_ts, na.rm= T))
    id_bps = subset(buttonpresses, ID == id)
    buttonpresses$timesincebaseline[buttonpresses$ID == id] = difftime(id_bps$timestamp, start_ts, units = "hours")
    buttonpresses$endtime[buttonpresses$ID == id] = difftime(end_ts, start_ts, units = "hours")
  }
  saveRDS(buttonpresses, "./R21_Study - tags/button_presses_extrainfo.RDS")
} else {
  buttonpresses = readRDS("./R21_Study - tags/button_presses_extrainfo.RDS")
}

## PLOT BPS Over TIME
library(ggplot2)
point_df = data.frame(id = buttonpresses$ID-1000,
                      timesincebaseline = buttonpresses$timesincebaseline,
                      endtime = buttonpresses$endtime)

point_df$id = as.factor(point_df$id)
levels(point_df$id) = 1:length(unique(point_df$id))
unique_ids = unique(point_df$id)
png("C:/Users/Balthazar/Documents/GitHub/fda-recurrentevents/figures/recurrentevents_plot.png",
    width = 640, height = 480, units = "px", pointsize = 16)
ggplot(point_df, aes(x=timesincebaseline, y = id)) +
  geom_hline(yintercept=unique_ids, color = "white") +
  geom_point(colour = "red", size = 1.5) +
  xlim(0,720) + 
  ylab("Participants") + xlab("Time Since Study Entry (in hours)") +
  geom_point(aes(x=endtime, y = id), colour = "black", size = 1.5, shape = 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_blank(), axis.ticks.y =element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))
dev.off()

## BUILD DAILY BP DATASET
## UP UNTIL LAST MEASUREMENT TIME
buttonpresses$timestamp = with_tz(buttonpresses$timestamp, tzone = "America/New_York")
count_df = data.frame(id = NA, days = NA, count = NA)
for(id in unique(buttonpresses$ID)) {
  subset_df = subset(buttonpresses, ID == id)
  end_time = max(subset_df$endtime)
  how_many_days = ceiling(end_time/24)
  start_of_first_day = floor_date(min(subset_df$timestamp), "day")
  which_day = ceiling(difftime(subset_df$timestamp, start_of_first_day, units = "days"))
  temp_df = data.frame(id = id, days = 1:how_many_days, count = NA)
  for(day in 1:how_many_days) {
    temp_df$count[temp_df$days == day] = sum(which_day == day)
  }
  count_df = rbind(count_df, temp_df)
}

count_df$id = as.factor(count_df$id)
## RESTRICT TO FIRST 30 DAYS
count_df = subset(count_df, days <= 30)
count_df$days = as.factor(count_df$days)
two.way = aov(count ~ id + days, data = count_df) 
summary(two.way)

## MOMENTS FOR PAPER
mean(count_df$count)
sd(count_df$count)
