generate_userday <- function(Cov_X, C, times, eig_Sigma, beta_1t) {
  ## Generate X(t) for t \in [0,1]
  X = t(C)%*%rnorm(M)  # Use Chol Decomp and Random Normals
  # plot(times, X)
  
  ## Compute the integrals at each time t
  result = rep(0,0)
  event_nodes = (length(beta_1t)+1):length(times)
  event_times = times[event_nodes]
  for (t in event_nodes) {
    temp_X = X[(t+1-length(beta_1t)):t]
    temp_Xb = beta_1t*temp_X
    term2 = (temp_Xb[1] + temp_Xb[length(temp_Xb)])/2 * gap + sum(temp_Xb[2:(length(temp_Xb)-1)]*gap)
    result = c(result, term2)
  }
  
  ## Calculate the logit based probabilities
  probs = expit(logit(10/length(times))+result)
  # plot(event_times, probs)
  Y = sapply(probs,function(p) rbinom(n=1,size = 1,p))
  # abline(v = event_times[Y==1])
  events = event_times[Y==1]
  # plot(event_times, probs, type = 'l')
  
  ## Given time stamp I want to generate the summary
  nonevents = event_times[rbinom(n = length(event_times), size = 1, prob = 10/1000)==1]
  nonevents = nonevents[!is.element(nonevents, events)]
  
  data = rep(0,0)
  # print("Events at times")
  # print(events)
  if(length(events) > 0) {
    for (t in events) {
      t_loc = which(is.element(times,t))
      X_Delta = X[(t_loc+1-length(beta_1t)):(t_loc)] # Past window of X
      C_x = rep(0, K_x)  
      for (i in 1:length(C_x)) {
        C_x[i] = eig_Sigma$vectors[,i]%*%X_Delta
      }
      row_data = c(1,C_x)
      data = rbind(data, row_data)
    }
  }
  
  # print("Nonevents at times")
  # print(nonevents)
  if(length(nonevents) > 0) {
    for (t in nonevents) {
      t_loc = which(is.element(times,t))
      X_Delta = X[(t_loc+1-length(beta_1t)):t_loc] # Past window of X
      C_x = rep(0, K_x)  
      for (i in 1:length(C_x)) {
        C_x[i] = eig_Sigma$vectors[,i]%*%X_Delta
      }
      row_data = c(0,C_x)
      data = rbind(data, row_data)
    }
  }
  return(data)
}

generate_complete_data <- function(N = 500, Cov_X, C, times, eig_Sigma, beta_1t) {
  full_data = rep(0,0)
  for (i in 1:N) {
    print(paste("On userday", i))
    userday_data = cbind(i, generate_userday(Cov_X, C, times, eig_Sigma, beta_1t))
    full_data = rbind(full_data, userday_data)
  }
  return(full_data)
}


expit <- function(x) {exp(x)/(1+exp(x))}

logit <- function(p) {log(p/(1-p))}

mise_calc <- function(beta_1t, betaHat.net.list) {
  temp_Xb = (beta_1t - betaHat.net.list[i] - mean(betaHat.net.list))^2
  term2 = (temp_Xb[1] + temp_Xb[length(temp_Xb)])/2 * gap + sum(temp_Xb[2:(length(temp_Xb)-1)]*gap) 
}