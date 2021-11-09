# Code to scale EDA for files #

## ---- Please place this file in a folder named "methods" 
##      under the directory where the scaled EDA data is stored. ----

## Recent update: May 27, 2021, by Xinrui Wu

#### For each user:
####     EDA_scaled = (EDA - mean_by_date) / std_by_date

library(tidyverse)
library(lubridate)

setwd('Z:/SI_data/R21_Study - EDA - scaled/')

## ------------------ Scale EDA for all observations ------------------ ##

n_obs = 91  # the largest number of observations (files)
eda_dir = "../R21_Study - EDA" # the directory where the raw EDA data is
for (i in 1:n_obs){
  id = 1000 + i 
  input_name = paste0(eda_dir, "/", id, "_EDA.rds")
  output_name = paste0("./", id, "_EDA_scaled.rds")
  
  if (file.exists(input_name)){
    eda_raw = readRDS(input_name)
    eda_raw = subset(eda_raw, EDA_HighLowPass < 20)
    eda_scaled = mutate(eda_raw, time=as_datetime(ts/1000)) %>%  # change data type of the timestamp
      mutate(date=as_date(time)) %>%  # add date column
      transmute(date, time, Device=E4_serial, EDA=EDA_HighLowPass) %>%
      drop_na() %>%
      arrange(Device, time) %>%
      group_by(Device) %>%  # data on different devices may have time overlap
      mutate(time_diff = c(-1, diff(time))) %>%
      ungroup() %>%
      filter(time_diff != 0) %>%  # drop duplicate data
      group_by(date) %>%
      mutate(EDA_scaled = scale(EDA)) %>%
      ungroup() %>%
      arrange(Device, time) %>%
      transmute(timestamp=time, Device, EDA_scaled)
    saveRDS(eda_scaled, file = output_name)
  }
}


## drop trouble data (EDA_HighLowPass = 0) in observation 1077
eda_scaled = readRDS('./1077_EDA_scaled.rds')
eda_scaled = drop_na(eda_scaled)
saveRDS(eda_scaled, file = './1077_EDA_scaled.rds')
