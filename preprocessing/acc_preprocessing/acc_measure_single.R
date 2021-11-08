# Function to compute AI #

## -------------------- Note -------------------- ##
## ---- Please place this file in a folder named "methods" 
##      under the directory where the AI data is stored. ----

## Recent update: May 27, 2021, by Xinrui Wu

## Function name:
##    acc_AI: compute AI for a single obs.

## Install the required package to compute Activity Index from Accelerometer data:
install.packages("devtools")
devtools::install_github("javybai/ActivityIndex")

## The function:
acc_AI = function(data, rate, epch){
  # Input:
  #   data: accelerometry data for one person with variables 
  #         time: record time in timestamp data type ("POSIXct", "POSIXt")
  #         acc_x, acc_y, acc_z: position variables
  #         Device: device name
  #   rate: the rate of the accelerometry data
  #   epch: integer, the length of the time window (in seconds) to compute AI
  # Output:
  #   measure: dataframe for AI, including Priod (which continuous period the record is in), Device,
  #                                        timestamp (position of the record in that period) and AI
  
  require(lubridate) # arrange timestamp
  require("ActivityIndex")
  require(tidyverse)
  
  acc_all = data.frame(transmute(data, time, Device=E4_serial, acc_x, acc_y, acc_z))
  devices = unique(acc_all$Device)
  measure = data.frame(Period=NULL, Device=NULL, timestamp=NULL, AI=NULL)
  
  # process the data from different devices seperately due to the timestamp overlap 
  for (d in 1:length(devices)){ 
    acc = filter(acc_all, Device == devices[d]) %>%
      arrange(time)
    
    # view loss of record larger than 1/32 (about 0.032) second as discontinuity.
    cutpoint = c(0, which(diff(acc$time) > 0.032), dim(acc)[1])
    
    # process the data in each consecutive period, 
    # because the function "computeActivityIndex" does not recognize time
    for (i in 1:(length(cutpoint) - 1)){
      # drop a few data on the tail of the period to make sure 
      # the number of data entries is a multiple of rate*epch,
      # as rate*epch entries in the acc data will be used to compute an AI entry
      nrow = (cutpoint[i+1] - cutpoint[i]) %/% (rate*epch) * (rate*epch)
      sub = data.frame(cbind(1:(cutpoint[i+1] - cutpoint[i]), 
                             acc$acc_x[(cutpoint[i]+1):cutpoint[i+1]], 
                             acc$acc_y[(cutpoint[i]+1):cutpoint[i+1]], 
                             acc$acc_z[(cutpoint[i]+1):cutpoint[i+1]])) [1:nrow, ]
      colnames(sub) = c("Index", "X", "Y", "Z")
      
      sub_ai = computeActivityIndex(sub, sigma0 = 1e-07, epoch = epch, hertz = rate)
      
      # keep the first timestamp of rate*epch entries in the acc data as the time of 
      # the corresponding AI entry
      time_ind = seq(from = cutpoint[i]+1, to = cutpoint[i]+nrow, by = rate*epch)
      sub_time = acc[time_ind, 'time'] 
      
      temp = data.frame(Period = rep(i, dim(sub_ai)[1]), 
                        Device = rep(devices[d], dim(sub_ai)[1]), 
                        timestamp = sub_time, AI = sub_ai$AI / 1e+07)
      measure = rbind(measure, temp)
    }
  }
  return(measure)
}
