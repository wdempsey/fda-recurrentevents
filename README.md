# Recurrent event analysis in the presence of functional covariates via random subsampling  #
Code for project on combining functional data analysis (FDA) with random subsampling methods for fast fitting of recurrent event models in the presence of functional covariates (i.e., real-time sensor-based measures).  

## Project Description ##
This project includes the code needed to reproduce the results in the paper: "Recurrent event analysis in the presence of functional covariates via random subsampling". If using this code please cite the paper using the following bibtex:

```tex
@misc{dempsey2022recurrent,
      title={Recurrent event analysis in the presence of real-time high frequency data via random subsampling}, 
      author={Walter Dempsey},
      year={2022},
      eprint={2204.06632},
      archivePrefix={arXiv},
      primaryClass={stat.ME}
}
```

The goal of this project is to build subsampling-unbiased likelihood-based estimators that trade off statistical and  computational efficiency.
Standard error estimates can be decomposed to understand the relative statistical efficiency and whether further subsampling may be required.
Methods allow for fast model fitting to help scientists interested in exploratory data analysis for recurrent event models in the presence of functional covariates.

## Code Description ##

Below we list important details in order use this repository:

0. **Dependencies**: All code is written in R. The code currently does not have dependencies beyond core R packages. The code was run on the research computing system managed by SLURM.
1. **Data**: Data is from 91 suicidal-inpatients at Franciscan Children's Hospital.
* Data is stored on Harvard's Sharepoint system; access to data is managed by [Evan Kleiman](https://kleimanlab.org)
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
2. **[Simulations](/simulations**: Directory with code to generate and analyze simulated datasets
3. **[Preprocessing](/preprocessing)**: Directory with code to
   pre-process case study data (AI and EDA).
4. **[Methods](/methods)**: Directory with code implementation of the
   random subsampling methodology
5. **Write-up**: Latex file containing the [statistical writeup](/write-up/fda-recurrent.tex)

## Contact me ##

Feel free to contact [Walter Dempsey](mailto:wdem@umich.edu) with questions, comments, and advice.
