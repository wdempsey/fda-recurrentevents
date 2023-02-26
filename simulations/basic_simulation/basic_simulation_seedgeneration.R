set.seed(4183)
generate_seeds = sample(1:100000, size = 500, replace = FALSE) 
saveRDS(generate_seeds, file = "basic_simulation_seeds.RDS")
