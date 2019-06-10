approximate_eda <- function(x, y) {
  inside_fn <-function(time_until_event) {
    gap = abs(x-time_until_event)
    which_obs = which(gap == min(gap))
    if( gap[which_obs] < 1/60 ) {
      return (y[which_obs])  
    } else {
      return (NA)
    }
  }
  return(inside_fn)
}

approximate_eda_apply <- function(iter, sequence, id_bp, eda) {
  current_ts = id_bp$ts_ms[iter]
  diff = eda$ts - current_ts
  minute_diff = diff/1000/60 
  keep_obs = (-30 < minute_diff) & (minute_diff <= 0)
  output = unlist(lapply(sequence, approximate_eda(minute_diff[keep_obs],eda$EDA_FeatureScaled_Filtered[keep_obs])))
  return(output)
}

nonevent_approximate_eda_apply <- function(iter, sequence, sampled_times, eda) {
  current_ts = sampled_times[iter]
  diff = eda$ts - current_ts
  minute_diff = diff/1000/60 
  keep_obs = (-30 < minute_diff) & (minute_diff <= 0)
  output = unlist(lapply(sequence, approximate_eda(minute_diff[keep_obs],eda$EDA_FeatureScaled_Filtered[keep_obs])))
  return(output)
}

generate_noneventtimes <- function(datetime_ts, sampling_rate, max.iters) {
  sampled_times = c()
  current_time = datetime_ts[1]
  for(iter in 1:max.iters) {
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
  }
  return(sampled_times)
}
