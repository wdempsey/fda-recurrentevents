approximate_eda <- function(x, y, tol = 1/60) {
  ## Input: x = minute_diff from current_time; y = corresponding EDA observations
  ## Output: Function that takes in time-until-event and outputs closest observed EDA
  inside_fn <-function(time_until_event) {
    ## Input: Time until event = s
    ## Output: Finds the observations within 1 second of t-s; 
    ## Output: outputs mean of these if multiple min(gap) 
    gap = abs(x-time_until_event)
    which_obs = which(gap == min(gap))
    if( min(gap) < tol ) {
      return (mean(y[which_obs]))  
    } else {
      return (NA)
    }
  }
  return(inside_fn)
}

approximate_eda_apply <- function(iter, sequence, id_bp, Delta, eda) {
  ## Input: Current iter, sequence of times until event, set of event times, EDA data
  ## Output: Outputs corresponding EDA at time-until-event
  current_ts = as_datetime(id_bp$ts[iter])
  minute_diff = interval(eda$timestamp,current_ts) %/% seconds(1) / 60
  keep_obs = which((Delta > minute_diff) & (minute_diff >= 0))
  output = unlist(lapply(sequence, approximate_eda(minute_diff[keep_obs],eda$EDA_scaled[keep_obs], tol = 1/60)))
  return(output)
}

approximate_acc_apply <- function(iter, sequence, id_bp, Delta, acc) {
  ## Input: Current iter, sequence of times until event, set of event times, EDA data
  ## Output: Outputs corresponding EDA at time-until-event
  current_ts = as_datetime(id_bp$ts[iter])
  minute_diff = interval(acc$timestamp,current_ts) %/% seconds(1) / 60
  keep_obs = which((Delta > minute_diff) & (minute_diff >= 0))
  output = unlist(lapply(sequence, approximate_eda(minute_diff[keep_obs],acc$AI[keep_obs], tol = 10/60)))
  return(output)
}

nonevent_approximate_eda_apply <- function(iter, sequence, sampled_times, Delta, eda) {
  ## Input: Current iter, sequence of times until event, set of sampled times, EDA data
  ## Output: Outputs corresponding EDA at time-until-event
  current_ts = as_datetime(sampled_times[iter])
  minute_diff = interval(eda$timestamp,current_ts) %/% seconds(1) / 60
  keep_obs = which((Delta > minute_diff) & (minute_diff >= 0))
  output = unlist(lapply(sequence, approximate_eda(minute_diff[keep_obs],eda$EDA_scaled[keep_obs])))
  return(output)
}

nonevent_approximate_acc_apply <- function(iter, sequence, sampled_times, Delta, acc) {
  ## Input: Current iter, sequence of times until event, set of sampled times, EDA data
  ## Output: Outputs corresponding EDA at time-until-event
  current_ts = as_datetime(sampled_times[iter])
  minute_diff = interval(acc$timestamp,current_ts) %/% seconds(1) / 60
  keep_obs = which((Delta > minute_diff) & (minute_diff >= 0))
  output = unlist(lapply(sequence, approximate_eda(minute_diff[keep_obs],acc$AI[keep_obs], tol = 10/60)))
  return(output)
}

generate_noneventtimes <- function(datetime_ts, sampling_rate, max.iters) {
  ## Code to generate sequence of sampled times
  ## When sample is not near a point, we move forward to next time
  ## This is because of gaps in coverage; luckily exponentials are memoryless
  ## Input: Sequence of datetimes, sampling rate, and max number of samples
  ## Output: Sequence of random sampled datetimes
  sampled_times = c()
  current_time = datetime_ts[1]
  for(iter in 1:max.iters) {
    # print(iter)
    random.next.time <- current_time + rexp(1,rate = sampling_rate)*60*60
    if(random.next.time > max(datetime_ts)) {break}
    diff = datetime_ts-random.next.time
    if (seconds(min(abs(diff))) < 1) {
      ## If the sampled time is < 1 minute from
      ## sensor times, then keep as a sampled time
      sampled_obs = which(abs(diff) == min(abs(diff)))
      sampled_times = c(sampled_times, datetime_ts[sampled_obs[1]])
      current_time = datetime_ts[sampled_obs[1]]
    } else{
      ## Else, set current.time to min time with diff > 0
      print("ELSE")
      new_obs = min(which(diff > 0))
      current_time = datetime_ts[new_obs[1]]
    }
    # print(current_time)
  }
  return(sampled_times)
}
