# Recurrent event analysis in the presence of functional covariates via random subsampling  #
Code for project on combining functional data analysis (FDA) with
random subsampling methods for fast fitting of recurrent event models
in the presence of functional covariates (i.e., real-time sensor-based
measures).

## Project Description ##
This project includes the code needed to reproduce the results in the
paper: "Recurrent event analysis in the presence of functional
covariates via random subsampling". If using this code please cite the
paper using the following bibtex:

```tex
@InProceedings{Dempsey:2019,
author = {Dempsey, Walter},
title = {Recurrent event analysis in the presence of functional
covariates via random subsampling},
year = {2019}}
```

The goal of this project is to build subsampling-unbiased
likelihood-based estimators, that tradeoff statistical efficiency for
computational efficiency.
Standard error estimates can be decomposed to understand the relative
statistical efficiency, and whether further subsampling may be
required.
Methods allow for fast model fitting to help scientists interested in
exploratory data analysis for recurrent event models in the presence
of functional covariates.

## Code Description ##

Below we list important details in order use this repository:

0. **Dependencies**: All code is written in R. The code currently does
   not have dependencies beyond core R packages. The code was run on
   the research computing system managed by SLURM.
1. **Data**: Data is from 91 suicidal-inpatients at Franciscan
   Children's Hospital. 
* Data is stored on Harvard's Sharepoint system; access to data is
  managed by [Evan Kleiman](https://kleimanlab.org)
* Data used in this paper consists of three data frames per participant: 
  * **Electrodermal activity (EDA)**: Timestamped measurements of
  EDA. Each row consists of Participant ID, E4 wearable ID, timestamp,
  raw EDA measurement, and 5 different processed EDA measurements.
  * **Accelerometer (ACC)**: Timestamped measurements of x-,y-, and
  z-coordinates for ACC as well as generated features. Each row
  consists of Participant ID, E4 wearable ID, timestamp, raw x,y,z
  accelerometer measurements, and processed ACC mesurements 
  * **Button presses (BP)**: Timestamps at which the participant
  pressed the button to indicate a moment of distress.  Each row
  consists of Participant ID, timestamp (both in seconds and
  milliseconds).
2. **[Data visualization](/visualization)**: Directory with code to
   generate plots per participant
3. **[Methods](/methods)**: Directory with code implementation of the
   random subsampling methodology
* **Inputs file**: Fill in the [`inputs.R`](/methods/inputs.R)
  file. File contains spots for:
  * `eda_dir`: Directory for EDA participant data (1 file per
  participant)
  * `acc_dir`: Directory for ACC participant data (1 file per
  participant)
  * `button_press_dir`: Directory for button press data (1 file
  containing all button press data)
  * `participant_ids`: List of all participant ids
  * `feature_generator`: User-specified function that takes in
  timestamp, participant id, and the button press data and
  outputs a feature vector. One simple example is number of button
  presses in prior 30-minutes.
  * `partial_pooling_indicator`: User-specified indicator of whether
  partial pooling should be used. If 1, fit hierarchical model; if 0,
  fit population model.
  * `sampling-rate`: User-specified sampling rate (per hour) for the
  non-event times. Default is set to 1 sample every 30-minutes.
  * `Delta`: user-specified window-length (per hour). Default is set
    to 30-minute window-length.
* **Methods run file**: The following code snippet can be run once
`inputs.R` has been completed: ```Rscript fda-recurrent.R```
* **Intermediate output files**:
  * `sampled_data.RDS`: data frame where each row is participant id,
  sampled_timestamp, minute-by-minute smoothed EDA for prior
  * `sandwich-mean-output`: Sandwich smooth estimate of mu(s,t) at
  sampled_timestamps
  * `pooled-covariance`: Pooled sample covariance estimate
  * `eigenfunctions` and `eigenvalues`: Vector of eigenfunctions and
    eigenvalues
* **Outputs file**: All output is saved in `outputs.RDS`. Contains:
  * `sample-ts`: List of all sampled time points
  * `theta-hat`: Penalized maximum likelihood estimate
  * `std-err-theta`: Standard errors
  * `var-sample-comp`: Component of the variance related to
  subsampling
  * `var-process-comp`: Component of the variance related to
  recurrent event process
4. **[Evaluate](/evaluation)**
* Code to evaluate method on simulation studies [simulation subdirectory](/evaluation/simulationstudies)
* Code to evaluate method on Franciscan dataset [case study subdirectory](/evaluation/casestudy)
5. **Write-up**: Corresponding statistical writeup

## Contact me ##

Feel free to contact [Walter Dempsey](mailto:wdem@umich.edu) with
questions, comments, and advice.
