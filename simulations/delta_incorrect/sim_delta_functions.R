generate_userday <- function(Cov_X, C, times, eig_Sigma, beta_1t, sampling_rate, base_rate, window_length = NA) {
  if (is.na(window_length)) { window_length = length(beta_1t) }
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
  probs = expit(base_rate+result)
  # plot(event_times, probs)
  Y = sapply(probs,function(p) rbinom(n=1,size = 1,p))
  # abline(v = event_times[Y==1])
  events = event_times[Y==1]
  # plot(event_times, probs, type = 'l')
  
  ## Given time stamp I want to generate the summary
  nonevents = event_times[rbinom(n = length(event_times), size = 1, prob = sampling_rate)==1]
  nonevents = nonevents[!is.element(nonevents, events)]
  
  data = rep(0,0)
  ## Construct event times
  if(length(events) > 0) {
    for (t in events) {
      t_loc = which(is.element(times,t))
      min_loc = max(t_loc+1-window_length, 1)
      X_Delta = X[(min_loc):(t_loc)] # Past window of X
      if(length(X_Delta) < window_length) {
        X_Delta = c(rep(0,window_length - length(X_Delta)), X_Delta)
      }
      C_x = rep(0, K_x)  
      for (i in 1:length(C_x)) {
        C_x[i] = eig_Sigma$vectors[,i]%*%X_Delta
      }
      row_data = c(1,C_x)
      data = rbind(data, row_data)
    }
  }
  
  ## Construct non-event times
  if(length(nonevents) > 0) {
    for (t in nonevents) {
      t_loc = which(is.element(times,t))
      min_loc = max(t_loc+1-window_length, 1)
      X_Delta = X[(min_loc):(t_loc)] # Past window of X
      if(length(X_Delta) < window_length) {
        X_Delta = c(rep(0,window_length - length(X_Delta)), X_Delta)
      } 
      C_x = rep(0, K_x)  
      for (i in 1:length(C_x)) {
        C_x[i] = eig_Sigma$vectors[,i]%*%X_Delta
      }
      row_data = c(0,C_x)
      data = rbind(data, row_data)
    }
  }
  return(list("data" = data, "event_times" = events, "nonevent_times" = nonevents))
}

generate_complete_data <- function(N = 500, Cov_X, C, times, eig_Sigma, beta_1t, sampling_rate, base_rate, window_length = NA) {
  full_data = rep(0,0)
  full_events = rep(0,0)
  full_nonevents = rep(0,0)
  for (i in 1:N) {
    # print(paste("On userday", i))
    userday_data = generate_userday(Cov_X, C, times, eig_Sigma, beta_1t, sampling_rate, base_rate, window_length)
    full_data = rbind(full_data, cbind(i,userday_data$data))
    if(length(userday_data$event_times) > 0) {
      full_events = rbind(full_events, cbind(i, userday_data$event_times))
    }
    if(length(userday_data$event_times) > 0) {
      full_nonevents = rbind(full_nonevents, cbind(i, userday_data$nonevent_times))
    }
  }
  # return(list("data" = full_data, "events" = full_events, "nonevents" = full_nonevents))
  return("data" = full_data)
}

expit <- function(x) {exp(x)/(1+exp(x))}

logit <- function(p) {log(p/(1-p))}

mise_calc <- function(beta_1t, betaHat.net) {
  diff_delta = length(beta_1t) - length(betaHat.net)
  new_beta_1t = c(rep(0,max(-diff_delta,0)), beta_1t)
  new_betaHat.net = c(rep(0,max(diff_delta,0)), as.numeric(betaHat.net))
  temp_Xb = (new_beta_1t - new_betaHat.net)^2
  term2 = (temp_Xb[1] + temp_Xb[length(temp_Xb)])/2 * gap + sum(temp_Xb[2:(length(temp_Xb)-1)]*gap) 
  return(term2)
}

partialmise_calc <- function(beta_1t, betaHat.net) {
  diff_delta = length(beta_1t) - length(betaHat.net)
  new_beta_1t = beta_1t[max(diff_delta+1,1):length(beta_1t)]
  new_betaHat.net = as.numeric(betaHat.net[max(-diff_delta+1,1):length(betaHat.net)])
  temp_Xb = (new_beta_1t - new_betaHat.net)^2
  term2 = (temp_Xb[1] + temp_Xb[length(temp_Xb)])/2 * gap + sum(temp_Xb[2:(length(temp_Xb)-1)]*gap) 
  return(term2 * length(beta_1t)/min(length(beta_1t), length(betaHat.net)))
}

subsample_dataset <- function(dataset, thinning_rate) {
  keep_obs = rbinom(nrow(dataset), size = 1, prob = thinning_rate*(1-dataset$Y) + dataset$Y)
  return(dataset[keep_obs==1,])
}

construct_J <- function(times, eig_Sigma, dataset, window_length) {
  ## Construct the J matrix
  K_b = 35
  local_times = times[1:window_length]
  num=K_b-2
  qtiles <- seq(0, 1, length = num + 2)[-c(1, num + 2)]
  knots <- quantile(local_times, qtiles)
  ## Basis = bs(t, kb)
  Basis = cbind(1, local_times, sapply(knots, function(k) ((local_times - k > 0) * (local_times - k))))
  Psi = t(eig_Sigma$vectors[,1:K_x])
  Phi = Basis
  J = Psi%*%Phi
  model.matrix = dataset[,3:ncol(dataset)]
  new.model.matrix = as.matrix(model.matrix)%*%J*gap
  w = cbind(new.model.matrix)
  return(list("Basis" = Basis, "w" = w))
}

runglmnet <- function(sampling_rate, dataset, w, Basis, epsilon = 0.0001) {
  n.tmp = length(dataset$Y) 
  p.tmp = ncol(w)
  subsample_offset = rep(log(sampling_rate),nrow(dataset))
  p.fac = rep(1, ncol(w))
  p.fac[1:3] = 0 #no penalty on the first 4 variables
  lambda_max <- 1/n.tmp^2 # THESE ARE FOR SIN SETTING 
  epsilon <- .000001 # THESE ARE FOR SIN SETTING 
  # lambda_max <- 1/n.tmp # THESE ARE FOR EXP SETTING 
  # epsilon <- .0001 # THESE ARE FOR EXP SETTING 
  K <- 50
  lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon), 
                              length.out = K)), digits = 15)
  
  # set.seed("97139817")
  # Start the clock!
  ptm <- proc.time()
  ridge.fit.cv <- cv.glmnet(w, dataset$Y, alpha = 0, intercept = TRUE, 
                            penalty.factor = p.fac, standardize = FALSE,
                            lambda = lambdapath, nfolds = 20,
                            family = "binomial")
  # Stop the clock
  runtime = proc.time() - ptm
  
  ridge.fit.lambda <- ridge.fit.cv$lambda.min
  # plot(ridge.fit.cv)
  
  # Extract coefficient values for lambda.min (without intercept)
  ridge.coef <- (coef(ridge.fit.cv, s = ridge.fit.lambda))[-1]
  
  betaHat.net <- Basis %*% ridge.coef
  
  # Compute the deviance for lambda.min
  bindeviance = deviance(ridge.fit.cv$glmnet.fit)[lambdapath == ridge.fit.lambda]
  
  return(list("runtime" = runtime[3], "betahat" = betaHat.net, "dev" = bindeviance))
}

makeplot <- function(local_times, betaHat.net, beta_1t) {
  par(mar = c(4,2.5,1,1)+0.1)
  plot(max(local_times)-local_times, betaHat.net, col = "blue", type = "l",
       ylim = range(beta_1t,betaHat.net),
       axes = FALSE, xlab = "", ylab = "")
  axis(side = 1); axis(side = 2, labels = FALSE)
  # lines(local_times, betaHat.net.star, col = "red")
  lines(local_times, beta_1t)
  mtext("Time until event", side = 1, line = 2)
  mtext(expression(paste(beta, "(s)")),side = 2, line = 1)
  mise_calc(beta_1t, betaHat.net)
}
