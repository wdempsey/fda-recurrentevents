rm(list=ls())
setwd("Z:/SI_data/")
library("lubridate")
col_names = c("id", "start_ts", "end_ts")
buttonpresses = readRDS("./R21_Study - tags/button_presses.RDS")
ID_list = unique(buttonpresses$ID)
buttonpresses$timestamp = as_datetime(buttonpresses$ts)
buttonpresses$timesincebaseline = NA
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
    start_ts = min(acc_start_ts, eda_start_ts, na.rm= T)
    id_bps = subset(buttonpresses, ID == id)
    buttonpresses$timesincebaseline[buttonpresses$ID == id] = id_bps$timestamp - start_ts
  }
  saveRDS(buttonpresses, "./R21_Study - tags/button_presses_extrainfo.RDS")
} else {
  buttonpresses = readRDS("./R21_Study - tags/button_presses_extrainfo.RDS")
}

## PLOT BPS Over TIME
library(ggplot2)
point_df = data.frame(id = buttonpresses$ID-1000,
                      timesincebaseline = buttonpresses$timesincebaseline/60)

point_df$id = as.factor(point_df$id)
levels(point_df$id) = 1:length(unique(point_df$id))
unique_ids = unique(point_df$id)
png("C:/Users/Balthazar/Documents/GitHub/fda-recurrentevents/figures/recurrentevents_plot.png",
    width = 480, height = 480, units = "px", pointsize = 16)
ggplot(subset(point_df, timesincebaseline < 250), aes(x=timesincebaseline, y = id)) +
  geom_hline(yintercept=unique_ids, color = "white") +
  geom_point(colour = "red", size = 1) +
  ylab("Participants") + xlab("Time Since Study Entry (in hours)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_blank(), axis.ticks.y =element_blank())
dev.off()
