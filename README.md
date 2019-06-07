# Recurrent event analysis in the presence of functional covariates via random subsampling  #
Code for combining functional data analysis (FDA) with random
subsampling methods for fast fitting of recurrent event models in the
presence of functional covariates (i.e., real-time sensor-based
measures)

## Project Description ##
This project includes the code needed to reproduce the results in the
paper: "Recurrent event analysis in the presence of functional
covariates via random subsampling". If using this code please cite the
paper using the following bibtex:

```tex
@InProceedings{Dempsey:2019,
author = {Dempsey, Walter},
title = {Recurrent event analysis in the presence of functional covariates via random subsampling},
year = {2019}}
```

The goal of this project is to do. 

## Code Description ##

Below we list important details in order use this repository:

0. **Dependencies**: All code is written in R. The code currently does
   not have dependencies beyond core R packages. The code was run on
   the research computing system managed by SLURM.
1. **Data**: Data is from 91 suicidal-inpatients at Franciscan
Children's Hospital. 
* Data is stored on Harvard's Sharepoint system; access to data is
managed by [Evan Kleiman](https://kleimanlab.org)
* Data consists of two data frames per participant:
  * **EDA**: Timestamped measurements of EDA. Each row consists of Participant ID,
    E4 wearable ID, timestamp, raw EDA measurement, and 5 different
    processed EDA measurements.
  * **Button presses**: Timestamps at which the participant pressed
    the button to indicate a moment of distress.  Each row consists of
    Participant ID, timestamp (both in seconds and milliseconds). 
2. **[Data visualization](/visualization)**: Directory with code to
generate plots per participant
3. **[Methods](/methods)**: Directory with code implementation of the
random subsampling methodology
* **Inputs file**: Fill in the `inputs.json` file
  * `eda-dir`: Directory for EDA participant data (1 EDA file per
  participant)
  * `button-press-dir`: Directory for button press data (1 file
  containing all button press data)
  * `feature-generator-dir`: User-completed `feature-generator.R` file that
  takes in timestamp, participant id, and the button press data and
  outputs a feature vector. One simple example is number of button
  presses in prior 30-minutes.
  * `partial-pooling-indicator`: User-specified indicator of whether
  partial pooling should be used.
  * `sampling-rate`: User-specified sampling rate (per hour) for the non-event
  times. Default is set to 1 sample every 30-minutes *(need to check we
  cannot improve it)*.
* **Methods run file**: The following code snippet can be run once
`inputs.json` has been completed:

```Rscript fda-recurrent inputs.json```

* **Outputs file**: All output is saved in `outputs.RDS`. Contains:
  * `sample-ts`: List of all sampled time points
  * `theta-hat`: Penalized maximum likelihood estimate
  * Pooled Sample Covariance estimate
  * Eigenfunctions and eigenvalues
  * Parameter estimates along with standard errors
4. **Evaluate**
* Code to evaluate method on simulation studies [simulation subdirectory](/evaluation/simulationstudies)
* Code to evaluate method on Franciscan dataset [case study subdirectory](/evaluation/casestudy)
5. **Write-up**
* Corresponding statistical writeup can be found in
2. Run methods. Point people towards the folder with methods. [the methods directory](/methods)

## Contact me ##

Feel free to contact [Walter Dempsey](mailto:wdem@umich.edu) with
questions, comments, and advice.
