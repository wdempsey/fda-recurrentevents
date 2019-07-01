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
buttonpress = readRDS("./R21_Study - tags/button_presses.RDS")

library(doParallel)
set.seed(87514)
sequence = seq(-30,0,1/60);
event_complete = nonevent_complete = matrix(nrow = 0, ncol = 3+length(sequence))

cl = Sys.getenv("SLURM_NTASKS_PER_NODE")
cl

registerDoParallel(cores = (cl))

# Shows the number of Parallel Workers to be used
getDoParWorkers()

for (id in 1001:1091) {
  print(paste("Current participant is", id))
  
  id_bp = subset(buttonpress, ID == id)
  
  file_name = paste("./R21_Study - EDA/",id,"_EDA.rds", sep = "")
  if(!file.exists(file_name)) {
    print(paste("No EDA file with id", id))
  } else{
    eda = readRDS(file_name)
    print("Made it to sampling nonevent times")
    datetime_ts = as_datetime(eda$ts/1000)
    sampled_times = generate_noneventtimes(datetime_ts, sampling_rate, max.iters = 20000)
    output_event = foreach(iter=1:nrow(id_bp), .combine=rbind) %dopar% approximate_eda_apply(iter, sequence, id_bp, eda)
    output_nonevent = foreach(iter=1:length(sampled_times), .combine=rbind) %dopar% nonevent_approximate_eda_apply(iter, sequence, sampled_times, eda)
    
    ## Define time as since minimum time in EDA or ID_BP
    base_time = min(eda$ts, id_bp$ts_ms)
    event_times_mins_since_base = (id_bp$ts_ms-base_time)/1000/60
    nonevent_times_mins_since_base = (sampled_times*1000-base_time)/1000/60
    
    if (nrow(id_bp) > 1) {
      event_temp = cbind(id, id_bp$ts, event_times_mins_since_base, output_event)
    } else {
      event_temp = c(id, id_bp$ts, event_times_mins_since_base, output_event)
    }
    
    if (length(sampled_times) > 1) {
      nonevent_temp = cbind(id, sampled_times, nonevent_times_mins_since_base, output_nonevent)
    } else {
      nonevent_temp = c(id, sampled_times, nonevent_times_mins_since_base, output_nonevent)
    }
    event_complete = rbind(event_complete, event_temp)
    nonevent_complete = rbind(nonevent_complete, nonevent_temp)    
  }
}  
saveRDS(object = event_complete, file = paste("event_complete_", today(), ".RDS", sep = ""))
saveRDS(object = nonevent_complete, file = paste("nonevent_complete_", today(), ".RDS", sep = ""))