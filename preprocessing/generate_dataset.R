## Code to generate the dataframes that will be used
## to construct functional covariates 
## Requires subsampling non-event times 
## and the pooled covariance estimator.

## Inputs: RDS files for participant level eda and button_press
## Outputs: event_complete.RDS and nonevent_complete.RDS
source('./gendata_functions.R')
sampling_rate = 2 ## Every 30 minutes

# test if there is at least one argument: if not, return an error
print("Made it to window length")
args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  print("Using default value!")
  Delta = 30
} else if (length(args)>=1) {
  # default output file
  Delta = as.numeric(args)
}
cat(paste0("Window length is ", Delta, "\n"))


require(lubridate)

## Windows
setwd("Z:/SI_data/")
## Linux
setwd("/mnt/turbo/SI_data/")
## GREAT LAKES
setwd("/nsf/turbo/sph-wdem/SI_data")

buttonpress = readRDS("./R21_Study - tags/button_presses.RDS")

library(doParallel)
set.seed(87514)
sequence_eda = seq(0,30,1/60);
sequence_acc = seq(0,30,10/60);
eda_event_complete = eda_nonevent_complete = matrix(nrow = 0, ncol = 3+length(sequence))
acc_event_complete = acc_nonevent_complete = matrix(nrow = 0, ncol = 3+length(sequence))

cl = Sys.getenv("SLURM_NTASKS_PER_NODE")
cl

registerDoParallel(cores = (cl))

# Shows the number of Parallel Workers to be used
getDoParWorkers()

for (id in 1001:1091) {
  print(paste("Current participant is", id))
  
  id_bp = subset(buttonpress, ID == id)
  
  eda_file_name = paste("./R21_Study - EDA - scaled/",id,"_EDA_scaled.rds", sep = "")
  acc_file_name = paste("./R21_Study - ACC - AI/",id,"_AI.rds", sep = "")
  if(!(file.exists(eda_file_name) & file.exists(acc_file_name))) {
    print(paste("No EDA and/or ACC file with id", id))
  } else{
    eda = readRDS(eda_file_name)
    print("Made it to sampling nonevent times")
    sampled_times = generate_noneventtimes(eda$timestamp, sampling_rate, max.iters = 20000)
    
    print("And now onto EDA calculations")
    eda_output_event = foreach(iter=1:nrow(id_bp), .combine=rbind) %dopar% approximate_eda_apply(iter, sequence_eda, id_bp, Delta, eda)
    eda_output_nonevent = foreach(iter=1:length(sampled_times), .combine=rbind) %dopar% nonevent_approximate_eda_apply(iter, sequence_eda, sampled_times, Delta, eda)
    
    print("And now onto ACC calculations")
    acc = readRDS(acc_file_name)
    acc_output_event = foreach(iter=1:nrow(id_bp), .combine=rbind) %dopar% approximate_acc_apply(iter, sequence_acc, id_bp, Delta, acc)
    acc_output_nonevent = foreach(iter=1:length(sampled_times), .combine=rbind) %dopar% nonevent_approximate_acc_apply(iter, sequence_acc, sampled_times, Delta, acc)
    
    ## Define time as since minimum time in EDA or ID_BP
    id_bp$timestamp = as_datetime(id_bp$ts)
    base_time = min(eda$timestamp, id_bp$timestamp, acc$timestamp)
    event_times_mins_since_base = interval(base_time, id_bp$timestamp) %/% seconds(1) / 60 ## Minutes since Base Time
    nonevent_times_mins_since_base = interval(base_time, as_datetime(sampled_times)) %/% seconds(1) / 60 ## Minutes since Base Time
    
    print("Building Event Data Frame")
    if (nrow(id_bp) > 1) {
      eda_event_temp = data.frame(id = id, timestamp = id_bp$timestamp, 
                                  timesincebaseline = event_times_mins_since_base, 
                                  eda_output_event)
      acc_event_temp = data.frame(id = id, timestamp = id_bp$timestamp, 
                                  timesincebaseline = event_times_mins_since_base, 
                                  acc_output_event)
    } else {
      eda_event_temp = data.frame(id = id, timestamp = id_bp$timestamp, 
                                  timesincebaseline = event_times_mins_since_base, 
                                  t(unlist(eda_output_event)))     
      acc_event_temp = data.frame(id = id, timestamp = id_bp$timestamp, 
                                  timesincebaseline = event_times_mins_since_base, 
                                  t(unlist(acc_output_event)))     
    }
    
    print("Building Non-Event Data Frame")
    if (length(sampled_times) > 1) {
      eda_nonevent_temp = data.frame(id = id, timestamp = sampled_times, 
                                     timesincebaseline = nonevent_times_mins_since_base, 
                                     eda_output_nonevent)
      acc_nonevent_temp = data.frame(id = id, timestamp = sampled_times, 
                                     timesincebaseline = nonevent_times_mins_since_base, 
                                     acc_output_nonevent)
    } else {
      eda_nonevent_temp = data.frame(id = id, timestamp = sampled_times, 
                                  timesincebaseline = nonevent_times_mins_since_base, 
                                  t(unlist(eda_output_nonevent)))    
      acc_nonevent_temp = data.frame(id = id, timestamp = sampled_times, 
                                     timesincebaseline = nonevent_times_mins_since_base, 
                                     t(unlist(acc_output_nonevent)))    
    }
    print("Appending EDA data")
    eda_event_complete = rbind(eda_event_complete, eda_event_temp)
    eda_nonevent_complete = rbind(eda_nonevent_complete, eda_nonevent_temp)    
    print(paste("Current number of rows for sampled times:", nrow(eda_nonevent_complete)))
    
    print("Appending ACC data")
    acc_event_complete = rbind(acc_event_complete, acc_event_temp)
    acc_nonevent_complete = rbind(acc_nonevent_complete, acc_nonevent_temp)    
    print(paste("Current number of rows for sampled times:", nrow(acc_nonevent_complete)))
  }
}  
saveRDS(object = eda_event_complete, file = paste("eda_event_complete_Delta_",Delta, "_HLP_", today(), ".RDS", sep = ""))
saveRDS(object = eda_nonevent_complete, file = paste("eda_nonevent_complete_Delta_",Delta, "_HLP_", today(), ".RDS", sep = ""))

saveRDS(object = acc_event_complete, file = paste("acc_event_complete_Delta_",Delta, "_HLP_", today(), ".RDS", sep = ""))
saveRDS(object = acc_nonevent_complete, file = paste("acc_nonevent_complete_Delta_",Delta, "_HLP_", today(), ".RDS", sep = ""))