setting = "exponential"
test36 = read.csv(paste("./output_csv/", setting, "_simdelta_results_36.csv", sep = ""))
test40 = read.csv(paste("./output_csv/", setting, "_simdelta_results_40.csv", sep = ""))
test44 = read.csv(paste("./output_csv/", setting, "_simdelta_results_44.csv", sep = ""))
test48 = read.csv(paste("./output_csv/", setting, "_simdelta_results_48.csv", sep = ""))
test52 = read.csv(paste("./output_csv/", setting, "_simdelta_results_52.csv", sep = ""))

overlapping_ids = Reduce(intersect, list(test36$ids, test40$ids, test44$ids, test48$ids, test52$ids))
overlapping_ids = overlapping_ids[order(overlapping_ids)]
whichmin = function(x) {which(x == min(x))}
result = rep(0,0)
for(id in overlapping_ids) {
  
  temp = cbind(subset(test36, ids == id)[,3],
               subset(test40, ids == id)[,3],
               subset(test44, ids == id)[,3],
               subset(test48, ids == id)[,3],
               subset(test52, ids == id)[,3])
  result_id = cbind(id, apply(temp, 1, whichmin), subset(test36, ids == id)$rates)
  result = rbind(result, result_id)

}

result = data.frame(result)
names(result) = c("id", "choice", "rate")
unique(result$rate)
table(subset(result, rate == 1.0)[,2])
table(subset(result, rate == 0.5)[,2])
table(subset(result, rate == 0.25)[,2])
table(subset(result, rate == 0.125)[,2])