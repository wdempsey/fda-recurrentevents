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
