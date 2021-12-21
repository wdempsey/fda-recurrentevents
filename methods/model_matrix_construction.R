## Code to compute the approximate likelihood function.
## Recall the approximation is design unbiased
## That is, E [ loglik (\theta; subssampled data)] = loglik(\theta; full data)
## where the expecation is over the random subsampling

## Inputs: RDS files for event and nonevent dataframes:
## Inputs: (1) mean functions mu(t,s), 
## Inputs: (2) coefficients on the residual X(t,s) - mu(t,s)
## Inputs: (3) eigenvectors of the pooled covariance estimate 
## Inputs: (4) complete-case timestamps
## Output: RDS file (mle_output.RDS) that contains 
## Output: (1) the approximate MLE (mle_output$estimate),
## Output: (2) the standard errors (mle_output$stderr), and 
## Output: (2) subsampling component of variance (mle_output$var_subsampl), and 
## Output: (4) recurrent event component of variance (mle_output$var_event). 
## Relative efficiency for any functional \Var (c^\prime \hat theta) can be computed via 
## numerator =  c^prime %*% mle_output$var_event %*% c 
## denominator =  c^prime %*% (mle_output$var_subsampl + mle_output$var_event ) %*% c 
## rel eff = numerator / denominator

## LIBRARIES
library(splines)
library(lubridate)
## WINDOWS
setwd("Z:/SI_data/")
## LINUX 
setwd('/mnt/turbo/SI_data/')

## PULL IN EVENT RDS FILES
set_of_types = c("acc", "eda")
for (type in set_of_types) {
  ## GENERATE SPLINES
  if (type == "eda") {
    sequence <- seq(-30,0, by = 1/60)
  } else if (type == "acc") {
    sequence <- seq(-30,0, by = 1/6)
  }
  K_b = 35
  num=K_b-3
  qtiles <- seq(0, 1, length = num + 2)[-c(1, num + 2)]
  knots <- quantile(sequence, qtiles)
  bb = cbind(1, sequence, sequence^2, sapply(knots, function(k) ((sequence - k > 0) * (sequence - k)) ^ 2))
  print("Generated Splines")
  
  if(!file.exists(paste(type, "_event_modelmatrix_2021-12-21.RDS", sep = ""))) {
    print(paste("Building", type, "event model matrix"))
    event_eigen_vectors = readRDS(paste(type, "_lin_event_eigen_vectors_2021-12-20.RDS", sep = ""))
    event_coef = readRDS(paste(type, "_lin_event_coef_matrix_2021-12-20.RDS", sep = ""))
    event_means = readRDS(paste(type, "_lin_event_means_2021-12-20.RDS", sep = ""))
    event_timesincebaseline = readRDS(paste(type, "_lin_event_complete_case_timesincebaseline_2021-12-20.RDS", sep = ""))
    event_timestamp = readRDS(paste(type, "_lin_event_complete_case_timestamp_2021-12-20.RDS", sep = ""))
    event_id = readRDS(paste(type, "_lin_event_complete_case_ids_2021-12-20.RDS", sep = ""))
    ## BUILD outer product of b-spline and event eigen
    event_J_coef = event_coef%*%t(event_eigen_vectors)%*%bb
    event_J_means = event_means%*%bb
    event_model.matrix = cbind(event_id, event_timesincebaseline, hour(event_timestamp), event_J_coef + event_J_means)
    event_model.matrix = data.frame(event_model.matrix)
    names(event_model.matrix) = c("id", "timesincebaseline", "hour", paste("acc_X", 1:ncol(event_J_coef), sep = ""))
    
    rm("event_means"); rm("event_coef"); rm("event_eigen_vectors"); rm("event_timestamp")
    rm("event_timesincebaseline"); rm("event_J_coef"); rm("event_J_means"); rm("event_id")
    
    saveRDS(event_model.matrix, file = paste(type, "_event_modelmatrix_",today(), ".RDS", sep = ""))
  } else {
    print(paste("Already built", type, "event"))
    event_model.matrix = readRDS(paste(type,"_event_modelmatrix_2021-12-21.RDS",sep =""))
  }
  ## PULL IN NONEVENT RDS FILES
  if(!file.exists(paste(type, "_nonevent_modelmatrix_2021-12-21.RDS", sep = ""))) {
    print(paste("Bulding", type, "non-event model matrix"))
    nonevent_eigen_vectors = readRDS(paste(type, "_lin_nonevent_eigen_vectors_2021-12-20.RDS", sep =""))
    nonevent_coef = readRDS(paste(type, "_lin_nonevent_coef_matrix_2021-12-20.RDS", sep = ""))
    nonevent_J_coef = nonevent_coef%*%t(nonevent_eigen_vectors)%*%bb
    rm("nonevent_eigen_vectors"); rm("nonevent_coef")
      
    nonevent_ids = readRDS(paste(type, "_nonevent_complete_case_ids_2021-12-20.RDS", sep =""))
    nonevent_means = readRDS(paste(type, "_lin_nonevent_means_2021-12-20.RDS", sep = ""))
    nonevent_J_means = nonevent_means%*%bb
    rm("nonevent_means")
    ## BUILD outer product of b-spline and event eigen
    ## BUILD outer product of b-spline and event eigen
    nonevent_timesincebaseline = readRDS(paste(type, "_lin_nonevent_complete_case_timesincebaseline_2021-12-20.RDS", sep = ""))
    nonevent_timestamp = readRDS(paste(type, "_lin_nonevent_complete_case_timestamp_2021-12-20.RDS", sep = ""))
    nonevent_model.matrix = cbind(nonevent_ids, nonevent_timesincebaseline, hour(as_datetime(nonevent_timestamp)), nonevent_J_coef + nonevent_J_means)
    nonevent_model.matrix = data.frame(nonevent_model.matrix)
    names(nonevent_model.matrix) = c("id", "timesincebaseline", "hour", paste("acc_X", 1:ncol(nonevent_J_coef), sep = ""))
    
    rm("nonevent_J_coef"); rm("nonevent_J_means")
    saveRDS(nonevent_model.matrix, file = paste(type, "_nonevent_modelmatrix_",today(), ".RDS", sep = ""))
  } else {
    print(paste("Already built", type, "nonevent"))
    nonevent_model.matrix = readRDS(paste(type, "_nonevent_modelmatrix_2021-12-21.RDS",  sep =""))
  }
}

## BUILDING THE MODEL MATRICES IS COMPLETE 