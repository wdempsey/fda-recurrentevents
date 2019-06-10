#First read in the arguments listed at the command line
args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  ## Supply default participant id
  id = 1040
}else{
  id = as.numeric(args[[1]])
}

print(paste("Current participant is", id))
source('./functions.R')

setwd("/Users/walterdempsey/Harvard University/R21_Study - EDA")
buttonpress = readRDS("../R21_Study - tags/button_presses.RDS")

file_name = paste(id,"_EDA.rds", sep = "")
if(file.exists(file_name)) {
  eda = readRDS(paste(id,"_EDA.rds", sep = ""))
} else{
  print(paste("No EDA file with id", id))
  break
}

id_bp = subset(buttonpress, ID == id)
library(doParallel)
library(doRNG)

cl = Sys.getenv("SLURM_NTASKS_PER_NODE")
cl

registerDoParallel(cores = (cl))

# Shows the number of Parallel Workers to be used
getDoParWorkers()

sequence = seq(-30,10, by = 0.25) # sequence along which we compute

event_times_output = foreach(bp_obs=1:nrow(id_bp), .combine = 'cbind') %dopar% 
  event_kernel_apply(bp_obs, eda, id_bp, current_bp_ts, sequence)

event_numerator_output = matrix(event_times_output[1,], ncol = length(sequence), byrow = TRUE)
event_denominator_output = matrix(event_times_output[2,], ncol = length(sequence), byrow = TRUE)

num.iters = 5000 # Number of non-event samples per participant
nonevent_times_output = foreach(bp_obs=1:num.iters, .combine = 'cbind', .options.RNG=197311) %dopar% 
  nonevent_kernel_apply(eda, id_bp, sequence)

nonevent_numerator_output = matrix(nonevent_times_output[1,], ncol = length(sequence), byrow = TRUE)
nonevent_denominator_output = matrix(nonevent_times_output[2,], ncol = length(sequence), byrow = TRUE)

full_data = list("event_numerator_output"=event_numerator_output, 
                 "event_denominator_output"=event_denominator_output,
                 "nonevent_numerator_output"=nonevent_numerator_output,
                 "nonevent_numerator_output"=nonevent_numerator_output)

setwd("/Volumes/murphy_lab/users/wdempsey/data-for-fda/data/")
saveRDS(object = full_data, file = paste("kerneloutput_id_",id,".RDS", sep=""))

