# Code to compute AI for files #
## ---- Please place this file in a folder named "methods" 
##      under the directory where the AI data is stored. ----

## Recent update: May 27, 2021, by Xinrui Wu

library(tidyverse)
library(lubridate)
source('./acc_measure_single.R')


## ------------------ Compute AI for all observations ------------------ ##
setwd("..")

epch = 10 # integer, the length of the time window (in seconds) to compute AI
n_obs = 91  # the largest number of observations (files)
acc_dir = "../R21_Study - ACC" # the directory where the accelerometer data is
for (i in 1:n_obs){
  id = 1000 + i 
  input_name = paste0(acc_dir, "/", id, "_acc.rds")
  output_name = paste0("./", id, "_AI.rds")
  
  if (file.exists(input_name)){
    acc = readRDS(input_name)
    acc$time = as_datetime(acc$ts/1000) # transform timestamp from numeric to required data type
    rate = round(1000 / (acc$ts[2] - acc$ts[1]))
    acc_ai = acc_AI(acc, rate, epch)
    saveRDS(acc_ai, file = output_name)
  }
}






