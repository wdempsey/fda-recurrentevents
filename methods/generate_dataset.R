source('./estimation_functions.R')
source('./inputs.R')
require(lubridate)

setwd("/n/murphy_lab/users/wdempsey/data-for-fda/data/")
buttonpress = readRDS("./R21_Study - tags/button_presses.RDS")

library(doParallel)
set.seed(87514)
sequence = seq(-30,0,1/60);
event_complete = nonevent_complete = matrix(nrow = 0, ncol = 2+length(sequence))

cl = Sys.getenv("SLURM_NTASKS_PER_NODE")
cl

registerDoParallel(cores = (cl))

# Shows the number of Parallel Workers to be used
getDoParWorkers()

for (id in 1001:1091) {
  print(paste("Current participant is", id))
  
  id_bp = subset(buttonpress, ID == id)
  
  file_name = paste("./R21_Study - EDA/",id,"_EDA.rds", sep = "")
  if(file.exists(file_name)) {
    eda = readRDS(file_name)
  } else{
    print(paste("No EDA file with id", id))
    break
  }
  print("Made it to sampling nonevent times")
  datetime_ts = as_datetime(eda$ts/1000)
  sampled_times = generate_noneventtimes(datetime_ts, sampling_rate, max.iters = 20000)
  output_event = foreach(iter=1:nrow(id_bp), .combine=rbind) %dopar% approximate_eda_apply(iter, sequence, id_bp, eda)
  output_nonevent = foreach(iter=1:length(sampled_times), .combine=rbind) %dopar% nonevent_approximate_eda_apply(iter, sequence, sampled_times, eda)
  
  ## Define time as since minimum time in EDA or ID_BP
  base_time = min(eda$ts, id_bp$ts_ms)
  event_times_mins_since_base = (id_bp$ts_ms-base_time)/1000/60
  nonevent_times_mins_since_base = (sampled_times*1000-base_time)/1000/60
  
  event_temp = cbind(id, event_times_mins_since_base, output_event)
  nonevent_temp = cbind(id, nonevent_times_mins_since_base, output_nonevent)
  
  event_complete = rbind(event_complete, event_temp)
  nonevent_complete = rbind(nonevent_complete, nonevent_temp)
  
}

saveRDS(object = event_complete, file = "event_complete.RDS")
saveRDS(object = nonevent_complete, file = "nonevent_complete.RDS")