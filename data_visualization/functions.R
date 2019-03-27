kernel_outcome <- function(eda, current_bp_ts, window_length = 3, bandwidth = 20/60) {
  # Sets up the dataset and current time-stamp to compare to
  # Choice of window_length is mostly to speed up computation
  # Choice of bandwidth
  ## NOTE: If all NAs near 
  inside_fn <- function(time_until_event) {
    # Computes numerator and denominator of the kernel estimates
    diff = eda$ts - current_bp_ts
    minute_diff = diff/1000/60 
    good_obs = (time_until_event - window_length < minute_diff) & (minute_diff <= time_until_event + window_length) & !is.na(eda$EDA_filtered)
    weights = dnorm(abs(minute_diff[good_obs]-time_until_event)^2, sd = bandwidth)
    return(list("numerator_weight"=sum((weights*eda$EDA_FeatureScaled_Filtered[good_obs]), na.rm = TRUE), 
                "denominator_weight"=sum(weights)))
  }
  return(inside_fn)
}

event_kernel_apply <- function(bp_obs, eda, id_bp, current_bp_ts, sequence) {
  current_bp_ts = id_bp$ts_ms[bp_obs]
  temp_output <- matrix(unlist(lapply(sequence, kernel_outcome(eda,current_bp_ts))), nrow = 2)
  return(temp_output)
}

nonevent_kernel_apply <- function(eda, id_bp, sequence) {
  ## Throw away if within 30 minutes of button press
  current_ts = sample(eda$ts, size = 1)
  if(min(abs(current_ts-id_bp$ts_ms)/1000/60) < 30){bad_choice = TRUE} else{bad_choice=FALSE}
  while(bad_choice) {
    current_ts = sample(eda$ts, size = 1)
    if(min(abs(current_ts-id_bp$ts_ms)/1000/60) < 30){bad_choice = TRUE} else{bad_choice=FALSE}
  }
  temp_output <- matrix(unlist(lapply(sequence, kernel_outcome(eda,current_ts))), nrow = 2)
  return(temp_output)
}
