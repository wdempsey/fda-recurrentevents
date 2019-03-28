sequence = seq(-30,10, by = 0.25) # sequence along which we compute
event_total_numerator = event_total_denominator = rep(NA, length(sequence))
nonevent_total_numerator = nonevent_total_denominator = rep(NA, length(sequence))
id_counter = NA
setwd("/Volumes/murphy_lab/users/wdempsey/fda-recurrent-events/data")

for(id in 1001:1091) {
  print(paste("On participant ",id))
  temp = readRDS(file = paste("kerneloutput_id_",id,".RDS",sep=""))
  
  id_counter = c(id_counter, id)
  event_total_numerator = rbind(event_total_numerator,colSums(temp$event_numerator_output, na.rm = TRUE))
  event_total_denominator = rbind(event_total_denominator,colSums(temp$event_denominator_output, na.rm = TRUE))
  nonevent_total_numerator = rbind(nonevent_total_numerator,colSums(temp$event_numerator_output, na.rm = TRUE))
  nonevent_total_denominator = rbind(nonevent_total_denominator,colSums(temp$event_denominator_output, na.rm = TRUE))
  
  event_mean = colSums(temp$event_numerator_output, na.rm = TRUE)/colSums(temp$event_denominator_output, na.rm = TRUE)
  nonevent_mean = colSums(temp$nonevent_numerator_output, na.rm = TRUE)/colSums(temp$nonevent_denominator_output, na.rm = TRUE)
  
  print(paste("Number of events is", nrow(temp$event_numerator_output)))
  
  filename = paste("~/Documents/github/fda-recurrentevents/data_visualization/figs/participant_",id,"_meancurves.png", sep="")
  png(filename, width = 480, height = 480, units = "px", pointsize = 12)
  par(mfrow = c(2,1), mar = c(3,3,1,1)+0.1)
  plot(sequence, event_mean-mean(event_mean), 
       type = "l", xlim = c(min(sequence), max(sequence)), 
       axes = FALSE, ylab = "", xlab = "")
  mtext("EDA (mean-centered)", side = 2, line = 2, cex = 0.75)
  #mtext("Time until event", side = 1, line = 2, cex = 0.75)
  mtext(paste("Participant", id), side = 3, line = -1)
  axis(side =1, cex.axis = 0.75); axis(side=2,  cex.axis = 0.75)
  lines(sequence[plot_obs], nonevent_mean-mean(nonevent_mean), col = "red")
  
  max.obs = max(c(event_mean, nonevent_mean), na.rm = TRUE); min.obs = min(c(event_mean, nonevent_mean), na.rm = TRUE)
  plot(sequence, event_mean, 
       type = "l", xlim = c(min(sequence), max(sequence)),
       ylim = c(min.obs*0.95, max.obs*1.05),
       axes = FALSE, ylab = "", xlab = "")
  mtext("EDA", side = 2, line = 2, cex = 0.75)
  mtext("Time until event", side = 1, line = 2, cex = 0.75)
  axis(side =1, cex.axis = 0.75); axis(side=2,  cex.axis = 0.75)
  lines(sequence[plot_obs], nonevent_mean, col = "red")
  dev.off()
}