## Code to generate the dataframes that will be used
## to construct functional covariates 
## Requires subsampling non-event times 
## and the pooled covariance estimator.

## Inputs: RDS files for participant level eda and button_press
## Outputs: event_complete.RDS and nonevent_complete.RDS
source('./gendata_functions.R')
source('../inputs.R')
require(lubridate)

setwd("/n/murphy_lab/users/wdempsey/data-for-fda/data/")
# setwd("/Volumes/murphy_lab/users/wdempsey/data-for-fda/data/")
# setwd("/Users/walterdempsey/Harvard University/R21_Study - raw_data/")
buttonpress = readRDS("./R21_Study - tags/button_presses.RDS")

# bp_ts = as_datetime(buttonpress$ts)
# temp = aggregate(rep(1, length(bp_ts))~month(bp_ts) + day(bp_ts) + as.factor(buttonpress$ID), FUN = sum)
# summary(temp$`rep(1, length(bp_ts))`)

library(doParallel)
set.seed(87514)
sequence = seq(-30,0,1/60);
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
  
  eda_file_name = paste("./R21_Study - EDA/",id,"_EDA.rds", sep = "")
  acc_file_name = paste("./R21_Study - ACC/",id,"_acc.rds", sep = "")
  if(!(file.exists(eda_file_name) & file.exists(acc_file_name))) {
    print(paste("No EDA file with id", id))
  } else{
    eda = readRDS(eda_file_name)
    acc = readRDS(acc_file_name); max_acc = 127 * 0.5773503 * 3
    acc$g = sqrt((acc$acc_x/max_acc)^2 + (acc$acc_y/max_acc)^2 + (acc$acc_z/max_acc)^2)
    print("Made it to sampling nonevent times")
    datetime_ts = as_datetime(eda$ts/1000)
    sampled_times = generate_noneventtimes(datetime_ts, sampling_rate, max.iters = 20000)
    eda_output_event = foreach(iter=1:nrow(id_bp), .combine=rbind) %dopar% approximate_eda_apply(iter, sequence, id_bp, eda)
    acc_output_event = foreach(iter=1:nrow(id_bp), .combine=rbind) %dopar% approximate_acc_apply(iter, sequence, id_bp, acc)
    eda_output_nonevent = foreach(iter=1:length(sampled_times), .combine=rbind) %dopar% nonevent_approximate_eda_apply(iter, sequence, sampled_times, eda)
    acc_output_nonevent = foreach(iter=1:length(sampled_times), .combine=rbind) %dopar% nonevent_approximate_acc_apply(iter, sequence, sampled_times, acc)
    
    ## Define time as since minimum time in EDA or ID_BP
    base_time = min(eda$ts, id_bp$ts_ms, acc$ts)
    event_times_mins_since_base = (id_bp$ts_ms-base_time)/1000/60
    nonevent_times_mins_since_base = (sampled_times*1000-base_time)/1000/60
    
    if (nrow(id_bp) > 1) {
      eda_event_temp = cbind(id, id_bp$ts, event_times_mins_since_base, eda_output_event)
      acc_event_temp = cbind(id, id_bp$ts, event_times_mins_since_base, acc_output_event)
    } else {
      eda_event_temp = c(id, id_bp$ts, event_times_mins_since_base, eda_output_event)
      acc_event_temp = c(id, id_bp$ts, event_times_mins_since_base, acc_output_event)
    }
    
    if (length(sampled_times) > 1) {
      eda_nonevent_temp = cbind(id, sampled_times, nonevent_times_mins_since_base, eda_output_nonevent)
      acc_nonevent_temp = cbind(id, sampled_times, nonevent_times_mins_since_base, acc_output_nonevent)
    } else {
      eda_nonevent_temp = c(id, sampled_times, nonevent_times_mins_since_base, eda_output_nonevent)
      acc_nonevent_temp = c(id, sampled_times, nonevent_times_mins_since_base, acc_output_nonevent)
    }
    eda_event_complete = rbind(eda_event_complete, eda_event_temp)
    eda_nonevent_complete = rbind(eda_nonevent_complete, eda_nonevent_temp)    
    
    acc_event_complete = rbind(acc_event_complete, acc_event_temp)
    acc_nonevent_complete = rbind(acc_nonevent_complete, acc_nonevent_temp)    
  }
}  
saveRDS(object = eda_event_complete, file = paste("eda_event_complete_HLP_", today(), ".RDS", sep = ""))
saveRDS(object = eda_nonevent_complete, file = paste("eda_nonevent_complete_HLP_", today(), ".RDS", sep = ""))

saveRDS(object = acc_event_complete, file = paste("acc_event_complete_HLP_", today(), ".RDS", sep = ""))
saveRDS(object = acc_nonevent_complete, file = paste("acc_nonevent_complete_HLP_", today(), ".RDS", sep = ""))