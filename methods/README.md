# Methods #

This directory contains all the files you need to implement the
subsampling methodology. These methods are: 

## Code description ##
Steps to run code are as follows:
1. [Mean and covariance estimation file](/mean_cov_estimation.R) can
   be run to generate `mean_estimates` and `pooled_estimates` from the pre-processed datasets.
2. [Likelihood file](likelihood_computation.R) can be run to fit the penalized logistic regression.
3. [The penalized mixed-effects likelihood file]
4. [Imputation](imputation.R) and [Bootstrap](bootstrap.R) takes input files

The order of operations is then
- Run `mean_cov_estimation.R` 
- Run `likelihood_computation.R` to get fixed effect model. 
- Run this for each choice of Delta and choose that which maximizes AIC. This completes the complete-case analysis.
- For optimal Delta, then run [bootstrap](bootstrap.R) $200$ times.
- For each boostrapped sample, impute missing rows using  `imputation.R`.
- Then re-run main analysis to get fixed effects estimates per dataset.

