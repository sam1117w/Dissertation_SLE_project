cat("\014")
rm(list = ls())

source("Model_input.R")
source("regression models.R")

# Initialise simulation variables
cycles <- 100
trace.all.treatment1 <- vector("list", n)
trace.all.treatment2 <- vector("list", n)

# Pre-generate a list of random numbers for each patient and cycle
set.seed(123)
random_numbers <- list()
for (i in 1:n) {
  random_numbers[[i]] <- matrix(runif(cycles * 16), nrow = cycles, ncol = 16)  # Pre-generate 15 random numbers per cycle
}

# Initialize storage for results with additional columns for cost breakdown
results <- tibble(
  treatment = rep(NA, n * 2),  # Adjusted to accommodate both treatments
  cumulative_costs_discounted = rep(NA, n * 2),
  cumulative_QALY_discounted = rep(NA, n * 2),
  cumulative_costs_undiscounted = rep(NA, n * 2),
  cumulative_QALY_undiscounted = rep(NA, n * 2),
  life_years = rep(NA, n * 2),  # New column to record life years
  SLICC = rep(NA, n * 2),  # New column to record SLICC score
  AMS = rep(NA, n * 2),  # New column to record SS score
  # New columns for cost breakdown
  cumulative_drug_costs_discounted = rep(NA, n * 2),
  cumulative_organ_damage_costs_discounted = rep(NA, n * 2),
  cumulative_disease_activity_costs_discounted = rep(NA, n * 2),
  cumulative_drug_costs_undiscounted = rep(NA, n * 2),
  cumulative_organ_damage_costs_undiscounted = rep(NA, n * 2),
  cumulative_disease_activity_costs_undiscounted = rep(NA, n * 2)
)

# Function to incorporate life table with SMR
get_life_table_qx <- function(age, gender, life_table) {
  if (gender == "Male") {
    qx_value <- life_table$Male[life_table$Age == age]
  } else if (gender == "Female") {
    qx_value <- life_table$Female[life_table$Age == age]
  } else {
    stop("Invalid gender specified. Gender should be 'Male' or 'Female'.")
  }
  
  # Handle cases where age exceeds the life table
  if (length(qx_value) == 0 || is.na(qx_value)) {
    return(1)  # Certain death if age exceeds the table or value is missing
  }
  
  # Retrieve the SMR for the given age
  smr_value <- life_table$SMR[life_table$Age == age]
  
  # Adjust the qx value by the SMR
  adjusted_qx <- qx_value * smr_value
  
  return(adjusted_qx)
}

# Initialize the progress bar
pb <- txtProgressBar(min = 0, max = n * 2, style = 3)

# Run the simulation for both treatments
for (treatment in 1:2) {
  
  # Simulation loop
  for (i in 1:n) 
  { 
    trace.i <- tibble(
      cycle = 0:cycles, 
      SLEDAI_score = NA, 
      AMS = NA, 
      steroid_current = NA, 
      steroid_cumulative = NA,
      SLICC_score = NA, 
      organ_musculoskeletal_score = NA, 
      organ_ocular_score = NA, 
      organ_neuro_score = NA,
      organ_peripheral_vascular_score = NA, 
      organ_gastro_score = NA, 
      organ_gonadal_score = NA, 
      organ_diabetes_score = NA,
      organ_malignancy_score = NA, 
      organ_cardio_score = NA, 
      organ_skin_score = NA, 
      organ_renal_score = NA, 
      organ_pulmonary_score = NA,
      cumulative_costs_discounted = 0,  # Initialize to zeros
      cumulative_costs_undiscounted = 0,  # Initialize to zeros
      cumulative_QALY_discounted = 0,  # Initialize to zeros
      cumulative_QALY_undiscounted = 0,  # Initialize to zeros
      clean_utility = NA,
      current_utility = NA,
      organ_multiplier = NA,  # column to store the organ multiplier to make sure the model is running as intended
      discontinued = NA,  # To track discontinuation status
      responder = NA,  # To track responder status
      death = NA,
      # New columns to store costs
      drug_costs_discounted = 0,
      organ_damage_costs_discounted = 0,
      disease_activity_costs_discounted = 0,
      drug_costs_undiscounted = 0,
      organ_damage_costs_undiscounted = 0,
      disease_activity_costs_undiscounted = 0
    )
    
    # Access the random numbers for this patient
    patient_randoms <- random_numbers[[i]]
    
    # Initialize response status
    response <- NA
    discontinued <- FALSE
    
    # setting up starting SLICC, SLEDAI scores, and steroid use
    trace.i$SLICC_score[1] <- baseline$SLICC_baseline[i]
    trace.i$SLEDAI_score[1] <- baseline$SLEDAI_baseline[i]
    trace.i$steroid_cumulative[1] <- 0
    trace.i$steroid_current[1] <- baseline$steroid_baseline[i] * 30
    trace.i$AMS[1] <- baseline$SLEDAI_baseline[i]
    current_age <- baseline$age[i]
    current_disease_duration <- baseline$disease_duration[i]
    trace.i$death[1] <- 0
    
    # Setting up organ damage score at baseline
    trace.i$organ_musculoskeletal_score[1] <- baseline$musculoskeletal_score[i]
    trace.i$organ_ocular_score[1] <- baseline$ocular_score[i]
    trace.i$organ_neuro_score[1] <- baseline$neuropsychiatric_score[i]
    trace.i$organ_peripheral_vascular_score[1] <- baseline$peripheral_vascular_score[i]
    trace.i$organ_gastro_score[1] <- baseline$gastrointestinal_score[i]
    trace.i$organ_gonadal_score[1] <- baseline$premature_gonadal_failure_score[i]
    trace.i$organ_diabetes_score[1] <- baseline$diabetes_score[i]
    trace.i$organ_malignancy_score[1] <- baseline$malignancy_score[i]
    trace.i$organ_cardio_score[1] <- baseline$cardiovascular_score[i]
    trace.i$organ_skin_score[1] <- baseline$skin_score[i]
    trace.i$organ_renal_score[1] <- baseline$renal_score[i]
    trace.i$organ_pulmonary_score[1] <- baseline$pulmonary_score[i]
    
    for (c in 2:(cycles + 1)) { 
      
      # Get the random numbers for this cycle
      random_cycle <- patient_randoms[c, ]
      
      # Update l and u
      l <- c - 2
      u <- c - 1
      
      # Calculate mortality probability and update death status
      death_prob_regression <- mortality_probability(
        baseline$age_at_diagnosis[i],
        trace.i$AMS[c-1], 
        baseline$haematological_inv[i],
        baseline$renal_inv[i], 
        baseline$HCQ[i], 
        trace.i$organ_cardio_score[c-1], 
        trace.i$organ_renal_score[c-1], 
        trace.i$organ_musculoskeletal_score[c-1], 
        trace.i$organ_peripheral_vascular_score[c-1], 
        trace.i$organ_gastro_score[c-1], 
        trace.i$organ_malignancy_score[c-1], 
        baseline$ethnicity[i], 
        baseline$cholesterol[i], 
        baseline$infection[i], 
        l, u
      )
      
      # Retrieve qx from the life table based on the updated age and gender
      patient_gender <- ifelse(baseline$gender[i] == 0, "Female", "Male")  
      qx_life_table <- get_life_table_qx(trunc(current_age), patient_gender, life_table_wide)
      
      # Take the minimum of the two probabilities
      death_prob <- max(death_prob_regression, qx_life_table)
      
      if (random_cycle[1] < death_prob) {
        trace.i$death[(c):(cycles + 1)] <- 1
        break
      } else {
        trace.i$death[c] <- 0
      }
      
      
      # Determine responder status and check for discontinuation
      if (c == 2 && !discontinued) {  # First cycle, if not discontinued
        if (treatment == 1) {  # Belimumab
          if (random_cycle[2] < responder_belimumab) {  # If patient responds to Belimumab
            responder <- TRUE
          } else {  # If patient does not respond
            responder <- FALSE
          }
        } else if (treatment == 2) {  # TAK279
          if (random_cycle[2] < responder_TAK279) {  # If patient responds to TAK279
            responder <- TRUE
          } else {  # If patient does not respond
            responder <- FALSE
          }
        }
      }
      
      # Record responder status
      trace.i$responder[c] <- responder
      
      # Check for discontinuation based on responder status
      if (!discontinued && responder == TRUE) {  # Only allow discontinuation if the patient is a responder
        if (treatment == 1) {  # Belimumab
          if (c == 2 && random_cycle[3] < discontinue_belimumab_y1) {  # Year 1 discontinuation
            discontinued <- TRUE
          } else if (c > 2 && random_cycle[4] < discontinue_belimumab_y2) {  # Year 2 onwards discontinuation
            discontinued <- TRUE
          }
        } else if (treatment == 2) {  # TAK279
          if (c == 2 && random_cycle[3] < discontinue_TAK279_y1) {  # Year 1 discontinuation
            discontinued <- TRUE
          } else if (c > 2 && random_cycle[4] < discontinue_TAK279_y2) {  # Year 2 onwards discontinuation
            discontinued <- TRUE
          }
        }
      }
      
      # Record discontinuation status
      trace.i$discontinued[c] <- discontinued
      
      # Update age and disease duration
      current_age <- current_age + 1
      current_disease_duration <- current_disease_duration + 1
      
      # Predict SLEDAI score
      if (!discontinued) {
        if (responder == TRUE) {  # If patient is a responder
          if (c == 2) {  # Special condition for the second cycle
            if (treatment == 1) {  # Apply treatment-specific effect for responders in cycle 2
              trace.i$SLEDAI_score[c] <- max(0, trace.i$SLEDAI_score[c - 1] * (1 + Y1_SS_score$Tx_1 + Y1_SS_score$Tx_1_responder))
            } else if (treatment == 2) {
              trace.i$SLEDAI_score[c] <- max(0, trace.i$SLEDAI_score[c - 1] * (1 + Y1_SS_score$Tx_2 + Y1_SS_score$Tx_2_responder))
            }
          } else {  # For other cycles
            if (treatment == 1) {  # Apply treatment-specific absolute effect
              trace.i$SLEDAI_score[c] <- max(0, predict_sledai(
                trace.i$SLEDAI_score[c - 1],
                baseline$gender[i],
                current_age,
                baseline$renal_inv[i],
                baseline$ethnicity[i],
                baseline$dna_binding[i],
                baseline$complement[i],
                baseline$haematological_inv[i],
                baseline$anaemia[i]
              ) * (1 + Y1_SS_score$Tx_1+Y1_SS_score$Tx_1_responder - Y1_SS_score$SOC))  # Apply absolute treatment effect for treatment 1
            } else if (treatment == 2) {
              trace.i$SLEDAI_score[c] <- max(0, predict_sledai(
                trace.i$SLEDAI_score[c - 1],
                baseline$gender[i],
                current_age,
                baseline$renal_inv[i],
                baseline$ethnicity[i],
                baseline$dna_binding[i],
                baseline$complement[i],
                baseline$haematological_inv[i],
                baseline$anaemia[i]
              ) * (1 + Y1_SS_score$Tx_2 + Y1_SS_score$Tx_2_responder - Y1_SS_score$SOC))  # Apply absolute treatment effect for treatment 2
            }
          }
        } else {  # If patient is a non-responder or SOC
          if (c == 2) {  # Special condition for the second cycle for non-responders
            trace.i$SLEDAI_score[c] <- max(0, trace.i$SLEDAI_score[c - 1] * (1 + Y1_SS_score$SOC))  # Apply SOC effect for non-responders in the second cycle
          } else {
            # Apply the original logic for non-responders in other cycles
            trace.i$SLEDAI_score[c] <- max(0, predict_sledai(
              trace.i$SLEDAI_score[c - 1],
              baseline$gender[i],
              current_age,
              baseline$renal_inv[i],
              baseline$ethnicity[i],
              baseline$dna_binding[i],
              baseline$complement[i],
              baseline$haematological_inv[i],
              baseline$anaemia[i]
            ))
          }
        }
      } else {  # Natural history (discontinued)
        trace.i$SLEDAI_score[c] <- max(0, predict_sledai(
          trace.i$SLEDAI_score[c - 1],
          baseline$gender[i],
          current_age,
          baseline$renal_inv[i],
          baseline$ethnicity[i],
          baseline$dna_binding[i],
          baseline$complement[i],
          baseline$haematological_inv[i],
          baseline$anaemia[i]
        ))
      }
      
      # Calculate AMS
      trace.i$AMS[c] <- calculate_ams(trace.i$SLEDAI_score[1:c], trace.i$cycle[1:c])
      
      # Calculate clean utility
      trace.i$clean_utility[c] <- clean_u(current_age, baseline$ethnicity[i], trace.i$SLEDAI_score[c])
      
      # Predict current steroid use
      trace.i$steroid_current[c] <- predict_steroid(
        trace.i$SLEDAI_score[c - 1], 
        baseline$gender[i], 
        current_age, 
        baseline$renal_inv[i], 
        baseline$ethnicity[i], 
        baseline$dna_binding[i], 
        baseline$complement[i], 
        baseline$haematological_inv[i], 
        baseline$anaemia[i]
      )
      
      # Update cumulative steroid use
      trace.i$steroid_cumulative[c] <- sum(c(trace.i$steroid_current[2:(c)]))/((c-1))
      
      # Organ damage calculations with calibration factors
      if (!discontinued && responder == TRUE && c <= calibration_duration) {  # Only apply within first 5 years, for responders who have not discontinued
        if (treatment == 1) {
          calibration_factor <- calibration_belimumab
        } else if (treatment == 2) {
          calibration_factor <- calibration_TAK279
        }
      } else {
        calibration_factor <- 1  # No calibration after 5 years, if discontinued, or if not a responder
      }
      
      # Musculoskeletal damage
      if (trace.i$organ_musculoskeletal_score[c - 1] == 0) {  
        musculoskeletal_prob <- musculoskeletal_probability(
          current_age, 
          trace.i$steroid_cumulative[c-1], 
          baseline$cytotoxic_Tx[i],
          trace.i$SLICC_score[c - 1], 
          l, u
        ) * calibration_factor
        if (!is.na(musculoskeletal_prob) && random_cycle[5] < musculoskeletal_prob) {
          trace.i$organ_musculoskeletal_score[c] <- SLICC_score$Musculoskeletal
        } else {
          trace.i$organ_musculoskeletal_score[c] <- trace.i$organ_musculoskeletal_score[c - 1]
        }
      } else {
        trace.i$organ_musculoskeletal_score[c] <- trace.i$organ_musculoskeletal_score[c - 1]
      }
      
      # Ocular damage
      if (trace.i$organ_ocular_score[c - 1] == 0) {
        ocular_prob <- ocular_probability(
          current_age, 
          baseline$neuropsychiatric_inv[i], 
          trace.i$steroid_cumulative[c-1], 
          baseline$hypertension[i], 
          l, u
        ) * calibration_factor
        if (!is.na(ocular_prob) && random_cycle[6]  < ocular_prob) {
          trace.i$organ_ocular_score[c] <- SLICC_score$Ocular
        } else {
          trace.i$organ_ocular_score[c] <- trace.i$organ_ocular_score[c - 1]
        }
      } else {
        trace.i$organ_ocular_score[c] <- trace.i$organ_ocular_score[c - 1]
      }
      
      # Neuropsychiatric damage
      if (trace.i$organ_neuro_score[c - 1] == 0) {
        neuro_prob <- neuro_probability(
          current_age, 
          baseline$neuropsychiatric_inv[i], 
          trace.i$steroid_cumulative[c-1], 
          baseline$cholesterol[i], 
          baseline$hypertension[i], 
          l, u
        ) * calibration_factor
        if (!is.na(neuro_prob) && random_cycle[7]  < neuro_prob) {
          trace.i$organ_neuro_score[c] <- SLICC_score$Neuropsychiatric
        } else {
          trace.i$organ_neuro_score[c] <- trace.i$organ_neuro_score[c - 1]
        }
      } else {
        trace.i$organ_neuro_score[c] <- trace.i$organ_neuro_score[c - 1]
      }
      
      # Peripheral Vascular damage
      if (trace.i$organ_peripheral_vascular_score[c - 1] == 0) {
        peripheral_vascular_prob <- probability_peri_vascular(
          baseline$smoker[i], 
          baseline$cholesterol[i], 
          baseline$anticoagulant_positive[i], 
          l, u
        ) * calibration_factor
        if (!is.na(peripheral_vascular_prob) && random_cycle[8]  < peripheral_vascular_prob) {
          trace.i$organ_peripheral_vascular_score[c] <- SLICC_score$Peripheral_vascular
        } else {
          trace.i$organ_peripheral_vascular_score[c] <- trace.i$organ_peripheral_vascular_score[c - 1]
        }
      } else {
        trace.i$organ_peripheral_vascular_score[c] <- trace.i$organ_peripheral_vascular_score[c - 1]
      }
      
      # Gastrointestinal damage
      if (trace.i$organ_gastro_score[c - 1] == 0) {
        gastro_prob <- probability_gastro(
          baseline$vasculitis[i], 
          trace.i$steroid_cumulative[c-1], 
          l, u
        ) * calibration_factor
        if (!is.na(gastro_prob) && random_cycle[9]  < gastro_prob) {
          trace.i$organ_gastro_score[c] <- SLICC_score$Gastrointestinal
        } else {
          trace.i$organ_gastro_score[c] <- trace.i$organ_gastro_score[c - 1]
        }
      } else {
        trace.i$organ_gastro_score[c] <- trace.i$organ_gastro_score[c - 1]
      }
      
      # Gonadal damage
      if (trace.i$organ_gonadal_score[c - 1] == 0) {
        gonadal_prob <- probability_gonadal(
          trace.i$steroid_cumulative[c-1], 
          baseline$cytotoxic_Tx[i], 
          baseline$HCQ[i], 
          baseline$cholesterol[i], 
          l, u
        ) * calibration_factor
        if (!is.na(gonadal_prob) && random_cycle[10]  < gonadal_prob) {
          trace.i$organ_gonadal_score[c] <- SLICC_score$gonadal_failure
        } else {
          trace.i$organ_gonadal_score[c] <- trace.i$organ_gonadal_score[c - 1]
        }
      } else {
        trace.i$organ_gonadal_score[c] <- trace.i$organ_gonadal_score[c - 1]
      }
      
      # Diabetes damage
      if (trace.i$organ_diabetes_score[c - 1] == 0) {
        diabetes_prob <- probability_diabetes(
          current_age, 
          baseline$renal_inv[i], 
          trace.i$steroid_cumulative[c-1], 
          baseline$ethnicity[i], 
          l, u
        ) * calibration_factor
        if (!is.na(diabetes_prob) && random_cycle[11]  < diabetes_prob) {
          trace.i$organ_diabetes_score[c] <- SLICC_score$Diabetes
        } else {
          trace.i$organ_diabetes_score[c] <- trace.i$organ_diabetes_score[c - 1]
        }
      } else {
        trace.i$organ_diabetes_score[c] <- trace.i$organ_diabetes_score[c - 1]
      }
      
      # Malignancy damage
      if (trace.i$organ_malignancy_score[c - 1] == 0) {
        malignancy_prob <- probability_malignancy(
          current_disease_duration, 
          baseline$age_at_diagnosis[i],
          trace.i$SLICC_score[c - 1], 
          baseline$cholesterol[i], 
          l, u
        ) * calibration_factor
        if (!is.na(malignancy_prob) && random_cycle[12]  < malignancy_prob) {
          trace.i$organ_malignancy_score[c] <- SLICC_score$Malignancy
        } else {
          trace.i$organ_malignancy_score[c] <- trace.i$organ_malignancy_score[c - 1]
        }
      } else {
        trace.i$organ_malignancy_score[c] <- trace.i$organ_malignancy_score[c - 1]
      }
      
      # Cardiovascular damage
      if (trace.i$organ_cardio_score[c - 1] == 0) {
        cardio_prob <- probability_cardio(
          current_disease_duration, 
          baseline$age_at_diagnosis[i], 
          trace.i$AMS[c-1],
          baseline$serositis_inv[i], 
          trace.i$steroid_cumulative[c-1], 
          baseline$cholesterol[i], 
          baseline$hypertension[i], 
          l, u
        ) * calibration_factor
        if (!is.na(cardio_prob) && random_cycle[13]  < cardio_prob) {
          trace.i$organ_cardio_score[c] <- SLICC_score$Cardiovascular
        } else {
          trace.i$organ_cardio_score[c] <- trace.i$organ_cardio_score[c - 1]
        }
      } else {
        trace.i$organ_cardio_score[c] <- trace.i$organ_cardio_score[c - 1]
      }
      
      # Skin damage
      if (trace.i$organ_skin_score[c - 1] == 0) {
        skin_prob <- probability_skin( 
          baseline$skin_inv[i], 
          baseline$HCQ[i], 
          baseline$ethnicity[i],
          baseline$smoker[i], 
          l, u
        ) * calibration_factor
        if (!is.na(skin_prob) && random_cycle[14]  < skin_prob) {
          trace.i$organ_skin_score[c] <- SLICC_score$Skin
        } else {
          trace.i$organ_skin_score[c] <- trace.i$organ_skin_score[c - 1]
        }
      } else {
        trace.i$organ_skin_score[c] <- trace.i$organ_skin_score[c - 1]
      }
      
      # Renal damage
      if (trace.i$organ_renal_score[c - 1] == 0) {
        renal_prob <- probability_renal( 
          trace.i$AMS[c-1],
          baseline$renal_inv[i], 
          baseline$cholesterol[i], 
          l, u
        ) * calibration_factor
        if (!is.na(renal_prob) && random_cycle[15]  < renal_prob) {
          trace.i$organ_renal_score[c] <- SLICC_score$Renal
        } else {
          trace.i$organ_renal_score[c] <- trace.i$organ_renal_score[c - 1]
        }
      } else {
        trace.i$organ_renal_score[c] <- trace.i$organ_renal_score[c - 1]
      }
      
      # Pulmonary damage
      if (trace.i$organ_pulmonary_score[c - 1] == 0) {
        pulmonary_prob <- probability_pulmonary( 
          current_age,
          baseline$renal_inv[i], 
          baseline$serositis_inv[i],
          baseline$dna_binding[i],
          trace.i$SLICC_score[c - 1],
          baseline$aCL[i],
          l, u
        ) * calibration_factor
        if (!is.na(pulmonary_prob) && random_cycle[16]  < pulmonary_prob) {
          trace.i$organ_pulmonary_score[c] <- SLICC_score$Pulmonary
        } else {
          trace.i$organ_pulmonary_score[c] <- trace.i$organ_pulmonary_score[c - 1]
        }
      } else {
        trace.i$organ_pulmonary_score[c] <- trace.i$organ_pulmonary_score[c - 1]
      }
      
      # Calculate the utility multipliers based on organ damage
      organ_multipliers <- numeric(12)  # Store multipliers for each organ
      organs <- c("musculoskeletal", "ocular", "neuro", "peripheral_vascular", "gastro", "gonadal", "diabetes", "malignancy", "cardio", "skin", "renal", "pulmonary")
      
      for (j in 1:12) {
        organ <- organs[j]
        if (trace.i[[paste0("organ_", organ, "_score")]][c] > 0) {
          if (c == 2) {  # First year of damage
            organ_multipliers[j] <- organ_damage_multiplier_y1[[organ]]
          } else {  # Subsequent years
            organ_multipliers[j] <- organ_damage_multiplier_y2[[organ]]
          }
        } else {
          organ_multipliers[j] <- 1  # No damage, no utility loss
        }
      }
      
      # Apply the largest utility multiplier loss to clean utility and store it
      trace.i$organ_multiplier[c] <- min(organ_multipliers)
      trace.i$current_utility[c] <- trace.i$clean_utility[c] * trace.i$organ_multiplier[c]
      
      # Calculate the QALY for this cycle (assuming each cycle is 1 year)
      qaly_for_cycle <- trace.i$current_utility[c]
      
      # Accumulate the QALYs for both discounted and undiscounted
      qaly_for_cycle_discounted <- trace.i$current_utility[c] * (1 / ((1 + qaly_d) ^ (c - 1)))
      trace.i$cumulative_QALY_discounted[c] <- trace.i$cumulative_QALY_discounted[c - 1] + qaly_for_cycle_discounted
      trace.i$cumulative_QALY_undiscounted[c] <- trace.i$cumulative_QALY_undiscounted[c - 1] + qaly_for_cycle
      
      # Update SLICC score
      trace.i$SLICC_score[c] <- 
        trace.i$organ_musculoskeletal_score[c] +
        trace.i$organ_ocular_score[c] +
        trace.i$organ_neuro_score[c] +
        trace.i$organ_peripheral_vascular_score[c] +
        trace.i$organ_gastro_score[c] +
        trace.i$organ_gonadal_score[c] +
        trace.i$organ_diabetes_score[c] +
        trace.i$organ_malignancy_score[c] +
        trace.i$organ_cardio_score[c] +
        trace.i$organ_skin_score[c] +
        trace.i$organ_renal_score[c] +
        trace.i$organ_pulmonary_score[c]
      
    }
    
    # Store the trace for the current patient based on treatment
    if (treatment == 1) {
      trace.all.treatment1[[i]] <- trace.i  # Store in the list for Treatment 1
    } else if (treatment == 2) {
      trace.all.treatment2[[i]] <- trace.i  # Store in the list for Treatment 2
    }
    # Update the progress bar after each patient simulation
    setTxtProgressBar(pb, (treatment - 1) * n + i)
  }
}

# Initialize cumulative costs for each patient and update them based on organ damage occurrence
for (treatment in 1:2) {  # Iterate over both treatments
  if (treatment == 1) {
    trace.all <- trace.all.treatment1
  } else if (treatment == 2) {
    trace.all <- trace.all.treatment2
  }
  
  for (i in 1:n) {
    # Loop through each cycle to accumulate costs
    for (c in 2:(cycles + 1)) {
      # Discount factor for the current cycle
      discount_factor <- 1 / ((1 + cost_d) ^ (c - 1))
      
      # Check if the patient is alive before proceeding
      if (trace.all[[i]]$death[c] == 1) break  # Stop if the patient has died
      
      # Initialize per-cycle costs to zero
      per_cycle_drug_cost_discounted <- 0
      per_cycle_organ_damage_cost_discounted <- 0
      per_cycle_disease_activity_cost_discounted <- 0
      per_cycle_drug_cost_undiscounted <- 0
      per_cycle_organ_damage_cost_undiscounted <- 0
      per_cycle_disease_activity_cost_undiscounted <- 0
      
      # 1. Add disease activity-related costs based on SLEDAI_score
      sle_score <- min(ceiling(trace.all[[i]]$SLEDAI_score[c]), 13)
      if (!is.na(sle_score)) {
        # Get the cost based on SLEDAI score
        if (sle_score == 13) {
          disease_activity_cost <- disease_activity_costs$SS_13_20
        } else {
          disease_activity_cost <- disease_activity_costs[[paste0("SS_", sle_score)]]
        }
        
        # Store per-cycle disease activity costs
        per_cycle_disease_activity_cost_discounted <- disease_activity_cost * discount_factor
        per_cycle_disease_activity_cost_undiscounted <- disease_activity_cost
      }
      
      # 2. Apply drug costs if the patient is a responder and has not discontinued
      if (!is.na(trace.all[[i]]$discontinued[c]) && !trace.all[[i]]$discontinued[c] &&
          !is.na(trace.all[[i]]$responder[c]) && trace.all[[i]]$responder[c]) {
        if (treatment == 1) {
          drug_cost <- belimumab_cost
        } else if (treatment == 2) {
          drug_cost <- TAK279_cost
        }
        # Store per-cycle drug costs
        per_cycle_drug_cost_discounted <- drug_cost * discount_factor
        per_cycle_drug_cost_undiscounted <- drug_cost
      }
      
      # 3. Add organ damage-related costs for each organ system (both discounted and undiscounted)
      for (organ in c("Musculoskeletal", "Ocular", "Neuro", "Peripheral_vascular", "Gastro", "Gonadal", "Diabetes", "Malignancy", "Cardio", "Skin", "Renal", "Pulmonary")) {
        # Check if the organ damage occurred in this cycle (i.e., damage occurred now but not in previous cycle)
        damage_occurred <- trace.all[[i]][[paste0("organ_", tolower(organ), "_score")]][c] > 0
        previous_damage <- trace.all[[i]][[paste0("organ_", tolower(organ), "_score")]][c-1] == 0
        
        if (damage_occurred && previous_damage) {
          # First year costs (when damage first occurs)
          organ_damage_cost <- organ_damage_costs_y1[[organ]]
        } else if (damage_occurred) {
          # Ongoing costs for subsequent years
          organ_damage_cost <- organ_damage_costs_y2[[organ]]
        } else {
          organ_damage_cost <- 0
        }
        
        # Add organ damage cost to per-cycle organ damage costs
        per_cycle_organ_damage_cost_discounted <- per_cycle_organ_damage_cost_discounted + organ_damage_cost * discount_factor
        per_cycle_organ_damage_cost_undiscounted <- per_cycle_organ_damage_cost_undiscounted + organ_damage_cost
      }
      
      # Now, store per-cycle costs in trace.all[[i]]
      trace.all[[i]]$drug_costs_discounted[c] <- per_cycle_drug_cost_discounted
      trace.all[[i]]$organ_damage_costs_discounted[c] <- per_cycle_organ_damage_cost_discounted
      trace.all[[i]]$disease_activity_costs_discounted[c] <- per_cycle_disease_activity_cost_discounted
      trace.all[[i]]$drug_costs_undiscounted[c] <- per_cycle_drug_cost_undiscounted
      trace.all[[i]]$organ_damage_costs_undiscounted[c] <- per_cycle_organ_damage_cost_undiscounted
      trace.all[[i]]$disease_activity_costs_undiscounted[c] <- per_cycle_disease_activity_cost_undiscounted
      
      # Update cumulative costs
      trace.all[[i]]$cumulative_costs_discounted[c] <- trace.all[[i]]$cumulative_costs_discounted[c - 1] + per_cycle_drug_cost_discounted + per_cycle_organ_damage_cost_discounted + per_cycle_disease_activity_cost_discounted
      trace.all[[i]]$cumulative_costs_undiscounted[c] <- trace.all[[i]]$cumulative_costs_undiscounted[c - 1] + per_cycle_drug_cost_undiscounted + per_cycle_organ_damage_cost_undiscounted + per_cycle_disease_activity_cost_undiscounted
      
      # QALY calculation (both discounted and undiscounted)
      qaly_for_cycle_discounted <- trace.all[[i]]$current_utility[c] * (1 / ((1 + qaly_d) ^ (c - 1)))
      qaly_for_cycle_undiscounted <- trace.all[[i]]$current_utility[c]
      
      # Accumulate the QALY (both discounted and undiscounted)
      trace.all[[i]]$cumulative_QALY_discounted[c] <- trace.all[[i]]$cumulative_QALY_discounted[c - 1] + qaly_for_cycle_discounted
      trace.all[[i]]$cumulative_QALY_undiscounted[c] <- trace.all[[i]]$cumulative_QALY_undiscounted[c - 1] + qaly_for_cycle_undiscounted
    }
    
    # Store the trace with updated cumulative costs back to the appropriate treatment list
    if (treatment == 1) {
      trace.all.treatment1[[i]] <- trace.all[[i]]
    } else if (treatment == 2) {
      trace.all.treatment2[[i]] <- trace.all[[i]]
    }
    
  }
  
  # Close the progress bar after simulation completes
  close(pb)
  
  # Store the results of the simulation for each patient
  for (i in 1:n) {
    # Identify the cycle where the patient dies (if any)
    death_cycle <- which(trace.all[[i]]$death == 1)[1]
    
    if (!is.na(death_cycle)) {
      # If the patient dies, stop accumulating costs and QALYs beyond the death cycle
      final_cumulative_costs_discounted <- trace.all[[i]]$cumulative_costs_discounted[death_cycle-1]
      final_cumulative_QALY_discounted <- trace.all[[i]]$cumulative_QALY_discounted[death_cycle-1]
      final_cumulative_costs_undiscounted <- trace.all[[i]]$cumulative_costs_undiscounted[death_cycle-1]
      final_cumulative_QALY_undiscounted <- trace.all[[i]]$cumulative_QALY_undiscounted[death_cycle-1]
      life_years <- death_cycle - 1  # Record the life years (death cycle - 1 as cycle starts at 0)
      final_SLICC <-  trace.all[[i]]$SLICC_score[death_cycle-1]
      final_AMS <-  trace.all[[i]]$SLEDAI_score[death_cycle-1]
      # Extract cumulative cost components up to death_cycle -1
      final_cumulative_drug_costs_discounted <- sum(trace.all[[i]]$drug_costs_discounted[2:(death_cycle-1)])
      final_cumulative_organ_damage_costs_discounted <- sum(trace.all[[i]]$organ_damage_costs_discounted[2:(death_cycle-1)])
      final_cumulative_disease_activity_costs_discounted <- sum(trace.all[[i]]$disease_activity_costs_discounted[2:(death_cycle-1)])
      final_cumulative_drug_costs_undiscounted <- sum(trace.all[[i]]$drug_costs_undiscounted[2:(death_cycle-1)])
      final_cumulative_organ_damage_costs_undiscounted <- sum(trace.all[[i]]$organ_damage_costs_undiscounted[2:(death_cycle-1)])
      final_cumulative_disease_activity_costs_undiscounted <- sum(trace.all[[i]]$disease_activity_costs_undiscounted[2:(death_cycle-1)])
    } else {
      # If the patient does not die, use the final cycle's cumulative costs and QALYs
      final_cumulative_costs_discounted <- trace.all[[i]]$cumulative_costs_discounted[cycles + 1]
      final_cumulative_QALY_discounted <- trace.all[[i]]$cumulative_QALY_discounted[cycles + 1]
      final_cumulative_costs_undiscounted <- trace.all[[i]]$cumulative_costs_undiscounted[cycles + 1]
      final_cumulative_QALY_undiscounted <- trace.all[[i]]$cumulative_QALY_undiscounted[cycles + 1]
      life_years <- cycles  # If alive, set life years to the maximum cycle
      final_SLICC <-  trace.all[[i]]$SLICC_score[cycles + 1]
      final_AMS <-  trace.all[[i]]$SLEDAI_score[cycles + 1]
      # Extract cumulative cost components up to the last cycle
      final_cumulative_drug_costs_discounted <- sum(trace.all[[i]]$drug_costs_discounted[2:(cycles+1)])
      final_cumulative_organ_damage_costs_discounted <- sum(trace.all[[i]]$organ_damage_costs_discounted[2:(cycles+1)])
      final_cumulative_disease_activity_costs_discounted <- sum(trace.all[[i]]$disease_activity_costs_discounted[2:(cycles+1)])
      final_cumulative_drug_costs_undiscounted <- sum(trace.all[[i]]$drug_costs_undiscounted[2:(cycles+1)])
      final_cumulative_organ_damage_costs_undiscounted <- sum(trace.all[[i]]$organ_damage_costs_undiscounted[2:(cycles+1)])
      final_cumulative_disease_activity_costs_undiscounted <- sum(trace.all[[i]]$disease_activity_costs_undiscounted[2:(cycles+1)])
    }
    
    # Store the treatment type (1 or 2) for the current patient in the results tibble
    results$treatment[(treatment - 1) * n + i] <- treatment
    
    # Store the cumulative costs and QALYs (both discounted and undiscounted)
    results$cumulative_costs_discounted[(treatment - 1) * n + i] <- final_cumulative_costs_discounted
    results$cumulative_QALY_discounted[(treatment - 1) * n + i] <- final_cumulative_QALY_discounted
    results$cumulative_costs_undiscounted[(treatment - 1) * n + i] <- final_cumulative_costs_undiscounted
    results$cumulative_QALY_undiscounted[(treatment - 1) * n + i] <- final_cumulative_QALY_undiscounted
    results$life_years[(treatment - 1) * n + i] <- life_years  # Store the life years
    results$SLICC[(treatment - 1) * n + i] <- final_SLICC # Store the SLICC scores
    results$AMS[(treatment - 1) * n + i] <- final_AMS
    
    # Store the cumulative cost components
    results$cumulative_drug_costs_discounted[(treatment - 1) * n + i] <- final_cumulative_drug_costs_discounted
    results$cumulative_organ_damage_costs_discounted[(treatment - 1) * n + i] <- final_cumulative_organ_damage_costs_discounted
    results$cumulative_disease_activity_costs_discounted[(treatment - 1) * n + i] <- final_cumulative_disease_activity_costs_discounted
    results$cumulative_drug_costs_undiscounted[(treatment - 1) * n + i] <- final_cumulative_drug_costs_undiscounted
    results$cumulative_organ_damage_costs_undiscounted[(treatment - 1) * n + i] <- final_cumulative_organ_damage_costs_undiscounted
    results$cumulative_disease_activity_costs_undiscounted[(treatment - 1) * n + i] <- final_cumulative_disease_activity_costs_undiscounted
  }
}

  # Mean costs and QALYs for Treatment 1
  mean_costs_treatment_1_discounted <- mean(results$cumulative_costs_discounted[results$treatment == 1])
  mean_QALY_treatment_1_discounted <- mean(results$cumulative_QALY_discounted[results$treatment == 1])
  mean_costs_treatment_1_undiscounted <- mean(results$cumulative_costs_undiscounted[results$treatment == 1])
  mean_QALY_treatment_1_undiscounted <- mean(results$cumulative_QALY_undiscounted[results$treatment == 1])
  mean_life_years_treatment_1 <- mean(results$life_years[results$treatment == 1])
  mean_SLICC_score_1 <- mean(results$SLICC[results$treatment == 1])
  mean_AMS_score_1 <- mean(results$AMS[results$treatment == 1])
  
  # Mean cost components for Treatment 1 (Discounted)
  mean_drug_costs_treatment_1_discounted <- mean(results$cumulative_drug_costs_discounted[results$treatment == 1])
  mean_organ_damage_costs_treatment_1_discounted <- mean(results$cumulative_organ_damage_costs_discounted[results$treatment == 1])
  mean_disease_activity_costs_treatment_1_discounted <- mean(results$cumulative_disease_activity_costs_discounted[results$treatment == 1])
  
  # Mean cost components for Treatment 1 (Undiscounted)
  mean_drug_costs_treatment_1_undiscounted <- mean(results$cumulative_drug_costs_undiscounted[results$treatment == 1])
  mean_organ_damage_costs_treatment_1_undiscounted <- mean(results$cumulative_organ_damage_costs_undiscounted[results$treatment == 1])
  mean_disease_activity_costs_treatment_1_undiscounted <- mean(results$cumulative_disease_activity_costs_undiscounted[results$treatment == 1])
  
  # Mean costs and QALYs for Treatment 2
  mean_costs_treatment_2_discounted <- mean(results$cumulative_costs_discounted[results$treatment == 2])
  mean_QALY_treatment_2_discounted <- mean(results$cumulative_QALY_discounted[results$treatment == 2])
  mean_costs_treatment_2_undiscounted <- mean(results$cumulative_costs_undiscounted[results$treatment == 2])
  mean_QALY_treatment_2_undiscounted <- mean(results$cumulative_QALY_undiscounted[results$treatment == 2])
  mean_life_years_treatment_2 <- mean(results$life_years[results$treatment == 2])
  mean_SLICC_score_2 <- mean(results$SLICC[results$treatment == 2])
  mean_AMS_score_2 <- mean(results$AMS[results$treatment == 2])
  
  # Mean cost components for Treatment 2 (Discounted)
  mean_drug_costs_treatment_2_discounted <- mean(results$cumulative_drug_costs_discounted[results$treatment == 2])
  mean_organ_damage_costs_treatment_2_discounted <- mean(results$cumulative_organ_damage_costs_discounted[results$treatment == 2])
  mean_disease_activity_costs_treatment_2_discounted <- mean(results$cumulative_disease_activity_costs_discounted[results$treatment == 2])
  
  # Mean cost components for Treatment 2 (Undiscounted)
  mean_drug_costs_treatment_2_undiscounted <- mean(results$cumulative_drug_costs_undiscounted[results$treatment == 2])
  mean_organ_damage_costs_treatment_2_undiscounted <- mean(results$cumulative_organ_damage_costs_undiscounted[results$treatment == 2])
  mean_disease_activity_costs_treatment_2_undiscounted <- mean(results$cumulative_disease_activity_costs_undiscounted[results$treatment == 2])

# Calculate the ICER (for both discounted and undiscounted)
incremental_costs_discounted <- mean_costs_treatment_2_discounted - mean_costs_treatment_1_discounted
incremental_QALY_discounted <- mean_QALY_treatment_2_discounted - mean_QALY_treatment_1_discounted
ICER_discounted <- incremental_costs_discounted / incremental_QALY_discounted

incremental_costs_undiscounted <- mean_costs_treatment_2_undiscounted - mean_costs_treatment_1_undiscounted
incremental_QALY_undiscounted <- mean_QALY_treatment_2_undiscounted - mean_QALY_treatment_1_undiscounted
ICER_undiscounted <- incremental_costs_undiscounted / incremental_QALY_undiscounted

# Display the results
cat("Mean Costs for Treatment 1 (Belimumab, Discounted): £", mean_costs_treatment_1_discounted, "\n")
cat("Mean QALYs for Treatment 1 (Belimumab, Discounted): ", mean_QALY_treatment_1_discounted, "\n")
cat("Mean Life Years for Treatment 1: ", mean_life_years_treatment_1, "\n\n")

cat("Mean Costs for Treatment 1 (Belimumab, Undiscounted): £", mean_costs_treatment_1_undiscounted, "\n")
cat("Mean QALYs for Treatment 1 (Belimumab, Undiscounted): ", mean_QALY_treatment_1_undiscounted, "\n")
cat("Mean Life Years for Treatment 1: ", mean_life_years_treatment_1, "\n\n")

cat("Mean Costs for Treatment 2 (TAK279, Discounted): £", mean_costs_treatment_2_discounted, "\n")
cat("Mean QALYs for Treatment 2 (TAK279, Discounted): ", mean_QALY_treatment_2_discounted, "\n")
cat("Mean Life Years for Treatment 2: ", mean_life_years_treatment_2, "\n\n")

cat("Mean Costs for Treatment 2 (TAK279, Undiscounted): £", mean_costs_treatment_2_undiscounted, "\n")
cat("Mean QALYs for Treatment 2 (TAK279, Undiscounted): ", mean_QALY_treatment_2_undiscounted, "\n")
cat("Mean Life Years for Treatment 2: ", mean_life_years_treatment_2, "\n\n")

cat("Incremental Costs (Discounted): £", incremental_costs_discounted, "\n")
cat("Incremental QALYs (Discounted): ", incremental_QALY_discounted, "\n")
cat("ICER (Discounted): £", ICER_discounted, "per QALY\n\n")

cat("Incremental Costs (Undiscounted): £", incremental_costs_undiscounted, "\n")
cat("Incremental QALYs (Undiscounted): ", incremental_QALY_undiscounted, "\n")
cat("ICER (Undiscounted): £", ICER_undiscounted, "per QALY\n")

# write.table(results,'results_Sep8th_60k.csv')