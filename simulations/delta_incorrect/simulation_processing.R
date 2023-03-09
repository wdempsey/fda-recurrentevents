settings = c("sine", "exponential")
num_simulations = 500
seq_windowlengths = seq(36,52,by=4)

## BUILD THE DATASET
for (setting in settings) {
  for (windowlength in seq_windowlengths) {
    print(paste("On setting =", setting, "and window length =", windowlength))
    df = rep(0,0)
    for (iter in 2:num_simulations) {
      if(file.exists(paste("./output_csv/",setting,"_simdelta_results_wl_",windowlength,"_arrayid_",iter,".csv", sep = ""))) {
        temp_df = data.frame(read.table(paste("./output_csv/",setting,"_simdelta_results_wl_",windowlength,"_arrayid_",iter,".csv", sep = ""), sep = ",", header = T))
        if( nrow(temp_df)!= 4) {print(paste("Iteration", iter, "has wrong number of rows"))}
      } else {
        print(paste("Iteration", iter, "did not generate file"))
      }
      df = rbind(df, temp_df)
    }  
    write.csv(x = df, file = paste("./output_csv/",setting,"_simdelta_allresults_wl_",windowlength,".csv", sep = ""))
  }
}
