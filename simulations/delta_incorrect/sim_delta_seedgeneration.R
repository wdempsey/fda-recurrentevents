set.seed(31245)
generate_seeds = sample(1:100000, size = 500, replace = FALSE) 
saveRDS(generate_seeds, file = "sim_delta_seeds.RDS")
