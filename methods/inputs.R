# Directory for EDA participant data
eda_dir <- "Z:/SI_data/R21_Study - EDA - scaled"

# Directory for ACC participant data
acc_dir <- "Z:/SI_data/R21_Study - ACC - AI"

# Directory for BP data
button_press_dir <- "Z:/SI_data/R21_Study - tags"

# List of all participant ids
participant_ids <- seq(1001,1091)

# Feature Generator Function
feature_generator <- function(timestamp, id, bp_data) {
  # Generates Features from BP data
  # Inputs: Timestamp, participant id, button press data
  # Output: Feature vector
  # Currently an indicator of BP in past hour
  
  # Take subset of bp data with participant ID == id
  bp_data_participant_subset <- subset(bp_data, ID==id)
  # Compute time since timestamp given
  temp <- (timestamp - bp_data_participant_subset$ts)/(60*60)
  # Return indicator of <1hr since last BP
  return(as.numeric(any(temp[temp > 0.0] < 1.0)))
}

## User-specified indicator of whether partial
## pooling should be used. If 1, fit hierarchical model; 
## if 0, fit population model.
partial_pooling_indicator <- 0


## User-specified sampling rate (per hour)
## for the non-event times
sampling_rate = 2 # 2 times per hour (i.e., every 30-minutes)
