## SETWD
setwd("Z:/SI_data/bootstrap_files_Delta/estimates/")

number_bootstraps = 200
number_mi = 2
seq_of_Delta = c(5,15,30)

current_Delta = 5

acc_size = length(seq(-current_Delta,0, by = 1/6))
eda_size = length(seq(-current_Delta,0, by = 1/60))
acc_theta_bootstrap = matrix(ncol = acc_size, nrow = number_bootstraps)
eda_theta_bootstrap = matrix(ncol = eda_size, nrow = number_bootstraps) 
acc_msw = matrix(0,nrow = acc_size, ncol = acc_size)
eda_msw = matrix(0,nrow = eda_size, ncol = eda_size)



for (current_bootstrap in 1:number_bootstraps) {
    acc_current_file_mi1 = paste("acc_bootstrap_Delta_",current_Delta, "_bs_", current_bootstrap, "_mi_1_2023-05-08.RDS", sep = "")
    eda_current_file_mi1 = paste("eda_bootstrap_Delta_",current_Delta, "_bs_", current_bootstrap, "_mi_1_2023-05-08.RDS", sep = "")
    acc_current_file_mi2 = paste("acc_bootstrap_Delta_",current_Delta, "_bs_", current_bootstrap, "_mi_2_2023-05-08.RDS", sep = "")
    eda_current_file_mi2 = paste("eda_bootstrap_Delta_",current_Delta, "_bs_", current_bootstrap, "_mi_2_2023-05-08.RDS", sep = "")
    if(all(file.exists(acc_current_file_mi1, acc_current_file_mi2))) {
      print(paste("Bootstrap", current_bootstrap, "has ACC data"))
      acc_temp_mi1 = readRDS(acc_current_file_mi1)
      acc_temp_mi2 = readRDS(acc_current_file_mi2)
      acc_temp_mi = (acc_temp_mi1$estimate + acc_temp_mi2$estimate)/2
      error1 = acc_temp_mi1$estimate - acc_temp_mi
      error2 = acc_temp_mi2$estimate - acc_temp_mi
      acc_msw = acc_msw + outer(error1, error1) + outer(error2, error2)
      acc_theta_bootstrap[current_bootstrap,] = acc_temp_mi
    }
    if(all(file.exists(eda_current_file_mi1, eda_current_file_mi2))) {
      print(paste("Bootstrap", current_bootstrap, "has EDA data"))
      eda_temp_mi1 = readRDS(eda_current_file_mi1)
      eda_temp_mi2 = readRDS(eda_current_file_mi2)
      eda_temp_mi = (eda_temp_mi1$estimate + eda_temp_mi2$estimate)/2
      error1 = eda_temp_mi1$estimate - eda_temp_mi
      error2 = eda_temp_mi2$estimate - eda_temp_mi
      eda_msw = eda_msw + outer(error1, error1) + outer(error2, error2)
      eda_theta_bootstrap[current_bootstrap,] = eda_temp_mi
    }
}

acc_theta_bs_final = colMeans(acc_theta_bootstrap, na.rm = TRUE)
eda_theta_bs_final = colMeans(eda_theta_bootstrap, na.rm = TRUE)
acc_msb = matrix(0, nrow = acc_size, ncol = acc_size)
eda_msb = matrix(0, nrow = eda_size, ncol = eda_size)
acc_B = 0
eda_B = 0

for(current_bootstrap in 1:number_bootstraps) {
  print(paste("At BootStrap:", current_bootstrap))
  current_acc = acc_theta_bootstrap[current_bootstrap,]
  current_eda = eda_theta_bootstrap[current_bootstrap,]
  if(!any(is.na(current_acc))) {
    acc_B = acc_B + 1 
    acc_error = current_acc - acc_theta_bs_final
    acc_msb = acc_msb + outer(acc_error, acc_error)
  }
  if(!any(is.na(current_eda))) {
    eda_B = eda_B + 1 
    eda_error = current_eda - eda_theta_bs_final
    eda_msb = eda_msb + outer(eda_error, eda_error)
  }
}

acc_msw = acc_msw/acc_B*(2-1)
acc_msb = acc_msb/(acc_B-1)
accSigma_BM = (acc_msb * (1+1/acc_B) - acc_msw)/2
acc_stderr = sqrt(diag(accSigma_BM))

eda_msw = eda_msw/eda_B*(2-1)
eda_msb = eda_msb/(eda_B-1)
edaSigma_BM = (eda_msb * (1+1/eda_B) - eda_msw)/2
eda_stderr = sqrt(diag(edaSigma_BM))

## CALCULATE Satterthwaite
numerator = ((diag(acc_msb)*(1+1/acc_B) - diag(acc_msw))/2)^2
denominator = diag(acc_msb)^2 * (1+1/acc_B)^2/4 / (acc_B -1) + diag(acc_msw)^2/(acc_B * 4 * (2-1))
satt_df = numerator/denominator
acct_stat = qt(0.975, satt_df)

numerator = ((diag(eda_msb)*(1+1/eda_B) - diag(eda_msw))/2)^2
denominator = diag(eda_msb)^2 * (1+1/eda_B)^2/4 / (eda_B -1) + diag(eda_msw)^2/(eda_B * 4 * (2-1))
satt_df = numerator/denominator
edat_stat = qt(0.975, satt_df)


ymax = max(acc_theta_bs_final + acc_stderr * acct_stat)
ymin = min(acc_theta_bs_final - acc_stderr * acct_stat)
plot(acc_temp_mi1$sequence, acc_theta_bs_final, ylim = c(ymin, ymax), type = "l", col = "red")
lines(acc_temp_mi1$sequence, acc_theta_bs_final - acc_stderr * acct_stat, col = "red", lty = 2)
lines(acc_temp_mi1$sequence, acc_theta_bs_final + acc_stderr * acct_stat, col = "red", lty = 2)
abline(h=0, lty = 2)

ymax = max(eda_theta_bs_final + eda_stderr * edat_stat)
ymin = min(eda_theta_bs_final - eda_stderr * edat_stat)
plot(eda_temp_mi1$sequence, eda_theta_bs_final, ylim = c(ymin, ymax), type = "l", col = "red")
lines(eda_temp_mi1$sequence, eda_theta_bs_final - eda_stderr * edat_stat, col = "red", lty = 2)
lines(eda_temp_mi1$sequence, eda_theta_bs_final + eda_stderr * edat_stat, col = "red", lty = 2)
abline(h=0, lty = 2)

## PLOT OF BOOTSTRAP
library(ggplot2)
df_acc_summary <- data.frame(sequence = acc_temp_mi1$sequence,
                             estimate = acc_theta_bs_final,
                             lowerCI = acc_theta_bs_final - acct_stat * acc_stderr,
                             upperCI = acc_theta_bs_final + acct_stat * acc_stderr)

acc_xmax = max(df_acc_summary$sequence[which(df_acc_summary$lowerCI > 0)])
acc_xmin = min(df_acc_summary$sequence[which(df_acc_summary$lowerCI > 0)])

png(paste("C:/Users/Balthazar/Documents/GitHub/fda-recurrentevents/figures/acc_coef_joint_bootstrap_Delta_",current_Delta,".png", sep=""),
    width = 480, height = 480, units = "px", pointsize = 16)
ggplot(df_acc_summary, aes(x=sequence, y=estimate)) +
  geom_line(size=1, alpha=0.8) +
  geom_ribbon(aes(ymin=lowerCI, ymax=upperCI) ,fill="blue", alpha=0.2) +
  annotate("rect",xmin=acc_xmin,xmax=acc_xmax,ymin=-Inf,ymax=Inf, alpha=0.1, fill="black") +
  xlab("Time until Button Press") + ylab(expression(paste(beta, "(s)"))) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))
dev.off()

df_eda_summary <- data.frame(sequence = eda_temp_mi1$sequence,
                             estimate = eda_theta_bs_final,
                             lowerCI = eda_theta_bs_final - edat_stat * eda_stderr,
                             upperCI = eda_theta_bs_final + edat_stat * eda_stderr)

png(paste("C:/Users/Balthazar/Documents/GitHub/fda-recurrentevents/figures/eda_coef_joint_bootstrap_Delta_",current_Delta,".png", sep=""),
    width = 480, height = 480, units = "px", pointsize = 16)
ggplot(df_eda_summary, aes(x=sequence, y=estimate)) +
  geom_line(size=1, alpha=0.8) +
  geom_ribbon(aes(ymin=lowerCI, ymax=upperCI) ,fill="blue", alpha=0.2) +
  xlab("Time until Button Press") + ylab(expression(paste(beta, "(s)"))) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))
dev.off()
