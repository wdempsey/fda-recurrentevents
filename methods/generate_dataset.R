source('./estimation_functions.R')

setwd("/Users/walterdempsey/Harvard University/R21_Study - EDA")
buttonpress = readRDS("../R21_Study - tags/button_presses.RDS")

library(doParallel)
set.seed(87514)
sequence = seq(-30,0,1/60); num.iters = 1000; 
event_complete = nonevent_complete = matrix(nrow = 0, ncol = 2+length(sequence))

cl = Sys.getenv("SLURM_NTASKS_PER_NODE")
cl

registerDoParallel(cores = (cl))

# Shows the number of Parallel Workers to be used
getDoParWorkers()

for (id in 1001:1091) {
  print(paste("Current participant is", id))
  
  id_bp = subset(buttonpress, ID == id)
  
  file_name = paste(id,"_EDA.rds", sep = "")
  if(file.exists(file_name)) {
    eda = readRDS(paste(id,"_EDA.rds", sep = ""))
  } else{
    print(paste("No EDA file with id", id))
    break
  }
  
  sampled_times = sample(eda$ts, size = num.iters, replace = FALSE)
  output_event = foreach(iter=1:nrow(id_bp), .combine=rbind) %dopar% approximate_eda_apply(iter, sequence, id_bp, eda)
  output_nonevent = foreach(iter=1:num.iters, .combine=rbind) %dopar% nonevent_approximate_eda_apply(iter, sequence, sampled_times, eda)
  
  ## Define time as since minimum time in EDA or ID_BP
  base_time = min(eda$ts, id_bp$ts_ms)
  event_times_mins_since_base = (id_bp$ts_ms-base_time)/1000/60
  nonevent_times_mins_since_base = (sampled_times-base_time)/1000/60
  
  event_temp = cbind(id, event_times_mins_since_base, output_event)
  nonevent_temp = cbind(id, nonevent_times_mins_since_base, output_nonevent)
  
  event_complete = rbind(event_complete, event_temp)
  nonevent_complete = rbind(nonevent_complete, nonevent_temp)
  
}

saveRDS(object = event_complete, file = "event_complete.RDS")
saveRDS(object = nonevent_complete, file = "nonevent_complete.RDS")