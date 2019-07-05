eventtime_llik <- function(obs, event_model.matrix, beta) {
  return(event_model.matrix[obs,]%*%beta)
}

eventtime_llik_derivative <- function(obs, event_model.matrix, beta) {
  return(event_model.matrix[obs,])
}

noneventtime_llik <- function(obs, nonevent_model.matrix, nonevent_ids, pi_ids, beta) {
  current_pi = pi_ids[nonevent_ids[obs]-1000]
  return(exp(nonevent_model.matrix[obs,]%*%beta)/current_pi)
}

noneventtime_llik_derivative <- function(obs, nonevent_model.matrix, nonevent_ids, pi_ids, beta) {
  current_pi = pi_ids[nonevent_ids[obs]-1000]
  return(exp(nonevent_model.matrix[obs,]%*%beta)/current_pi*nonevent_model.matrix[obs,])
}

# score_computation <- function(event_coef, event_means, event_J,
#                               nonevent_coef, nonevent_means, nonevent_J,
#                               nonevent_ids, pi_ids) {
#   inside_fn <- function(beta) {
#     term1 = foreach(obs=1:nrow(event_means), .combine = "+") %dopar% eventtime_llik(obs, event_coef, event_means, event_J, beta)    
#     term2 = foreach(obs=1:nrow(nonevent_means), .combine = "+") %dopar% noneventtime_llik(obs, nonevent_coef, nonevent_means, nonevent_J, nonevent_ids, pi_ids, beta)
#     return(-term1+term2)
#   }
#   return(inside_fn)
# }

llik_computation <- function(event_coef, event_means, event_J,
                             nonevent_coef, nonevent_means, nonevent_J,
                             nonevent_ids, pi_ids) {
  inside_fn <- function(beta) {
    term1 = sum(event_model.matrix %*% beta)    
    term2 = foreach(obs=1:nrow(nonevent_means), .combine = "+") %dopar% noneventtime_llik(obs, nonevent_model.matrix, nonevent_ids, pi_ids, beta)
    return(-term1+term2)
  }
  return(inside_fn)
}

llik_computation_derivative <- function(event_coef, event_means, event_J,
                                        nonevent_coef, nonevent_means, nonevent_J,
                                        nonevent_ids, pi_ids) {
  inside_fn <- function(beta) {
    term1 = foreach(obs=1:nrow(event_means), .combine = "+") %dopar% eventtime_llik_derivative(obs, event_coef, event_means, event_J, beta)    
    term2 = foreach(obs=1:nrow(nonevent_means), .combine = "+") %dopar% noneventtime_llik_derivative(obs, nonevent_coef, nonevent_means, nonevent_J, nonevent_ids, pi_ids, beta)
    return(-term1+term2)
  }
  return(inside_fn)
}
