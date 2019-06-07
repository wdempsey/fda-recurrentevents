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
2. **Data visualization**: Generates
*alfja
3. **Methods**:
*
4. **Evaluate**
* Code to evaluate method on simulation studies [simulation subdirectory](/evaluation/simulationstudies)
* Code to evaluate method on Franciscan dataset [case study subdirectory](/evaluation/casestudy)
5. **Write-up**
* Corresponding statistical writeup can be found in
2. Run methods. Point people towards the folder with methods. [the methods directory](/methods)

# Contact us #

Feel free to contact [Walter Dempsey](mailto:wdem@umich.edu) with
questions, comments, and advice.
