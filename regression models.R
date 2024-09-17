################################################################################

###################  Define function for regression models  ####################

################################################################################

###########################################
##### section A: model coefficients #######
###########################################

######### defining the coefficients for the 1. disease activity, 2. mortality, 3. steroid use, 4.1-12. organ damage models (12), 5. clean utility regression

#### 1. Define the coefficients from the regression model for SLEDAI scores ####
coefficients_sle <- list(
  intercept = 1.4914,
  sle_lag = -0.4600,
  male_gender = -0.0798,
  log_age = -0.2414,
  renal_inv = -0.3014,
  african_american = 0.3829,
  increased_dna_binding = 0.2759,
  low_complement = 0.4838,
  haematological_involvement = 0.1043,
  anaemia = 0.1521
)

#### 2. define the coefficients from the mortality model ####
coefficients_mortality <- list(
  intercept = -8.8971,
  age_at_diagnosis = 1.0344,
  AMS = 1.2462,
  renal_inv = 1.7534,
  haematological_involvement = 2.3247,
  HCQ = 0.5721,
  cardiovascular_score = 1.402,
  renal_score = 1.7408,
  musculoskeletal_score = 1.3801,
  peripheral_vascular_score = 2.4659,
  gastrointestinal_score = 1.9663,
  malignancy = 2.7104,
  african_american = 2.1041,
  Cholesterol = 1.0041,
  infection = 2.7912,
  distribution_param = 1.8143
)

#### 3. Define the coefficients from the steroid model ####
coefficients_steroid <- list(
  current_sle = 0.777,
  intercept = 3.4754
)

#### 4. organ damage models  ####

# 4.1 Define the coefficients from the musculoskeletal model
coefficients_musculoskeletal <- list(
  intercept = -7.4545,
  log_age = 2.6446,
  cumulative_steroid = 1.0012,
  cytotoxic_Tx = 1.2955,
  SLICC = 1.1084,
  distribution_param = 0.8544
)

# 4.2 Define the coefficients from the ocular model
coefficients_ocular <- list(
  intercept = -12.814,
  log_age = 9.7554,
  neuroinv = 2.1807,
  cumulative_steroid = 1.0013,
  hypertension = 1.4793,
  distribution_param = 0.8027
)

# 4.3 Define the coefficients from the neuropsychiatric model
coefficients_neuro <- list(
  intercept = -8.0924,
  log_age = 2.4353,
  neuroinv = 7.9842,
  cumulative_steroid = 1.0004,
  cholesterol = 1.0034,
  hypertension = 1.8038,
  distribution_param = 0.8254
)

# 4.4 Define the coefficients from the peripheral vascular damage model
coefficients_peri_vascular <- list(
  intercept = -7.3112,
  smoker=1.9814,
  cholesterol = 1.0049,
  anticoagulant_positive = 2.2758
)

# 4.5 Define the coefficients from the gastrointestinal model
coefficients_gastro <- list(
  intercept = -5.0532,
  vasculitis=5.2815,
  cumulative_steroid = 1.0008
)

# 4.6 Define the coefficients from the gonadal failure model
coefficients_gonadal <- list(
  intercept = -7.2658,
  cumulative_steroid = 1.0017,
  cytotoxic_Tx = 2.3074,
  HCQ = 0.3896,
  cholesterol = 1.0044
)

# 4.7 Define the coefficients from the diabetes model
coefficients_diabetes <- list(
  intercept = -14.672,
  log_age = 9.8192,
  renal_inv = 2.4083,
  cumulative_steroid = 1.0012,
  african_american = 2.0703
)

# 4.8 Define the coefficients from the malignancy model
coefficients_malignancy <- list(
  intercept = -4.809,
  log_duration = 1.3486,
  age_at_diagnosis = 1.0259,
  SLICC = 1.2013,
  cholesterol = 0.9906
)

# 4.9 Define the coefficients from the cardiovascular model
coefficients_cardiovascular <- list(
  intercept = -8.8971,
  log_duration = 1.8532,
  age_at_diagnosis = 1.0479,
  AMS = 1.1585,
  serositisinv = 3.0638,
  cumulative_steroid = 1.0012,
  cholesterol = 1.0027,
  hypertension = 2.2692,
  distribution_param = -0.0445
)

# 4.10 Define the coefficients from the skin model
coefficients_skin <- list(
  intercept = 10.2613,
  skininv = 0.0866,
  HCQ = 4.4634,
  african_american = 0.0873,
  smoker = 0.1584,
  distribution_param = 1.5915
)

# 4.11 Define the coefficients from the renal model
coefficients_renal <- list(
  intercept = 10.2576,
  AMS = 0.819,
  renal_inv = 0.0257,
  cholesterol = 0.9922,
  distribution_param = 1.4942
)

# 4.12 Define the coefficients from the pulmonary model
coefficients_pulmonary <- list(
  intercept = 11.5193,
  log_age = 0.1824,
  renal_inv = 0.4234,
  serositisinv = 0.1062,
  dna_binding = 0.6176,
  SLICC = 0.8542,
  aCL = 0.2918,
  distribution_param = 1.1486
)

# 4.13 Define the coefficients for 'clean' utility
coefficients_clean_utility <- list(
  intercept = 1.275,
  log_age = -0.140,
  african_american = -0.036,
  current_sle = -0.009
)

#######################################
##### section B: model functions ######
#######################################
#### defining the model functions ####

#### 1. disease activity model ####
#### 1.1 SLEDAI score model (disease activity) 

# Define a function to predict SLEDAI scores
predict_sledai <- function(previous_sle, gender, age, renalinv, ethnicity, dna_binding, complement, haematologicalinv, anaemia) {
  log_age <- log(age)
  predicted_sle <- coefficients_sle$intercept +
    coefficients_sle$sle_lag * previous_sle +
    coefficients_sle$male_gender * gender +
    coefficients_sle$log_age * log_age +
    coefficients_sle$renal_inv * renalinv +
    coefficients_sle$african_american * ethnicity +
    coefficients_sle$increased_dna_binding * dna_binding +
    coefficients_sle$low_complement * complement +
    coefficients_sle$haematological_involvement * haematologicalinv +
    coefficients_sle$anaemia * anaemia
  current_SS <- previous_sle + predicted_sle # note the model is returning the CHANGE in SS score 
  return(current_SS)
}

#### 1.2. calculating AMS
# Define a function to calculate AMS
calculate_ams <- function(s, t) {
  if (length(s) != length(t)) {
    stop("The lengths of s and t must be the same")
  }
  
  # Ensure there are at least two points to calculate AMS
  if (length(s) < 2) {
    stop("Not enough points to calculate AMS")
  }
  
  n <- length(s)
  numerator <- sum((s[2:n] + s[1:(n-1)]) / 2 * (t[2:n] - t[1:(n-1)]))
  denominator <- sum(t[2:n] - t[1:(n-1)])
  AMS <- numerator / denominator
  return(AMS)
}

#### 2. calculating prob of death ####
# Define the combined mortality probability function with AMS
mortality_probability <- function(age_at_diagnosis,AMS,haematologicalinv, renalinv, HCQ, cardiovascular_s, renal_s, musculoskeletal_s,
                                  peripheral_vascular_s, gastrointestinal_s, malignancy, ethnicity, cholesterol, infection, l, u) {
  
  lam <- exp(coefficients_mortality$intercept +
               log(coefficients_mortality$age_at_diagnosis) * age_at_diagnosis +
               log(coefficients_mortality$AMS) * AMS +
               log(coefficients_mortality$haematological_involvement) * haematologicalinv +
               log(coefficients_mortality$renal_inv) * renalinv +
               log(coefficients_mortality$HCQ) * HCQ +
               log(coefficients_mortality$cardiovascular_score) * cardiovascular_s +
               log(coefficients_mortality$renal_score) * renal_s +
               log(coefficients_mortality$musculoskeletal_score) * musculoskeletal_s +
               log(coefficients_mortality$peripheral_vascular_score) * peripheral_vascular_s + 
               log(coefficients_mortality$gastrointestinal_score) * gastrointestinal_s + 
               log(coefficients_mortality$malignancy) * malignancy +
               log(coefficients_mortality$african_american) * ethnicity +
               log(coefficients_mortality$Cholesterol) * cholesterol +
               log(coefficients_mortality$infection) * infection)
  
  haz <- function(time) {
    lam * coefficients_mortality$distribution_param * time^(coefficients_mortality$distribution_param - 1)
  }
  
  ch <- integrate(haz, lower = l, upper = u)$value
  pr <- 1 - exp(-ch)
  
  return(pr)
}

#### 3. calculating the steroid use ####
predict_steroid <- function(previous_sle, gender, age, renalinv, ethnicity, dna_binding, complement, haematologicalinv, anaemia) {
  # Predict SLEDAI scores over the time interval
  current_sle <- predict_sledai (previous_sle, gender, age, renalinv, ethnicity, dna_binding, complement, haematologicalinv, anaemia)
  steroid_use <-
    (coefficients_steroid$intercept +
       coefficients_steroid$current_sle * current_sle)*365/12 # calculating here the monthly use of steroid (mg)
  return(steroid_use)
}

#### 4. calculating the prob of musculoskeletal damage ####

musculoskeletal_probability <- function(age,cumulative_steroid,cytotoxic_Tx,SLICC,l,u){
  # Define the hazard function using all covariates
  lam <- exp(coefficients_musculoskeletal$intercept +
               log(coefficients_musculoskeletal$log_age) * log(age) +
               log(coefficients_musculoskeletal$cumulative_steroid) * cumulative_steroid +
               log(coefficients_musculoskeletal$cytotoxic_Tx) * cytotoxic_Tx + 
               log(coefficients_musculoskeletal$SLICC) * SLICC)  
  
  haz <- function(time) {
    lam * coefficients_musculoskeletal$distribution_param * time^(coefficients_musculoskeletal$distribution_param - 1)
  } # note the musculoskeletal damage model is a Weibull model
  
  # Integrate the hazard function to get the cumulative hazard within [l, u]
  ch <- integrate(haz, lower = l, upper = u)$value
  
  # Calculate the probability of death based on the cumulative hazard
  pr <- 1 - exp(-ch)
  
  return(pr)
}      


#### 4. calculating the prob of ocular damage ####
ocular_probability <- function(age,neuroinv,cumulative_steroid,hypertension,l,u){
  # Define the hazard function using all covariates
  lam <- exp(coefficients_ocular$intercept +
               log(coefficients_ocular$log_age) * log(age) +
               log(coefficients_ocular$neuroinv) * neuroinv +
               log(coefficients_ocular$cumulative_steroid) * cumulative_steroid +
               log(coefficients_ocular$hypertension) * hypertension)  
  
  haz <- function(time) {
    lam * coefficients_ocular$distribution_param * time^(coefficients_ocular$distribution_param - 1)
  } # note the ocular damage model is a Weibull model
  
  # Integrate the hazard function to get the cumulative hazard within [l, u]
  ch <- integrate(haz, lower = l, upper = u)$value
  
  # Calculate the probability of death based on the cumulative hazard
  pr <- 1 - exp(-ch)
  
  return(pr)
}            

#### 5. calculating the prob of neuropsychiatric damage ####
neuro_probability <- function(age,neuroinv,cumulative_steroid,cholesterol,hypertension,l,u){
  # Define the hazard function using all covariates
  lam <- exp(coefficients_neuro$intercept +
               log(coefficients_neuro$log_age) * log(age) +
               log(coefficients_neuro$neuroinv) * neuroinv +
               log(coefficients_neuro$cumulative_steroid) * cumulative_steroid +
               log(coefficients_neuro$cholesterol) * cholesterol + 
               log(coefficients_neuro$hypertension) * hypertension)  
  
  haz <- function(time) {
    lam * coefficients_neuro$distribution_param * time^(coefficients_neuro$distribution_param - 1)
  } # note the neuropsychiatric damage model is a Weibull model
  
  # Integrate the hazard function to get the cumulative hazard within [l, u]
  ch <- integrate(haz, lower = l, upper = u)$value
  
  # Calculate the probability of death based on the cumulative hazard
  pr <- 1 - exp(-ch)
  
  return(pr)
}            

#### 6. calculating the prob of peripheral vascular damage ####
### !! this is an exponential model ###
probability_peri_vascular <- function(smoker,cholesterol,anticoagulant_positive,l,u){
  # Define the hazard function using all covariates
  lam <- exp(coefficients_peri_vascular$intercept +
               log(coefficients_peri_vascular$smoker) * smoker +
               log(coefficients_peri_vascular$cholesterol) * cholesterol + 
               log(coefficients_peri_vascular$anticoagulant_positive) * anticoagulant_positive)
  
  haz <- function(time) { lam*(time/time) } # since this is an exponential model, the baseline hazard rate is constant with time
  
  # Integrate the hazard function to get the cumulative hazard within [l, u]
  ch <- integrate(haz, lower = l, upper = u)$value
  
  # Calculate the probability of death based on the cumulative hazard
  pr <- 1 - exp(-ch)
  
  return(pr)
}   

#### 7. calculating the prob of gastrointestinal damage ####
### !! this is an exponential model ###
probability_gastro <- function(vasculitis,cumulative_steroid,l,u){
  # Define the hazard function using all covariates
  lam <- exp(coefficients_gastro$intercept +
               log(coefficients_gastro$vasculitis) * vasculitis +
               log(coefficients_gastro$cumulative_steroid) * cumulative_steroid)
  haz <- function(time) { lam*(time/time) } # since this is an exponential model, the baseline hazard rate is constant with time
  
  # Integrate the hazard function to get the cumulative hazard within [l, u]
  ch <- integrate(haz, lower = l, upper = u)$value
  
  # Calculate the probability of death based on the cumulative hazard
  pr <- 1 - exp(-ch)
  
  return(pr)
}         


#### 8. calculating the prob of gonadal damage ####
### !! this is an exponential model ###
probability_gonadal <- function(cumulative_steroid,cytotoxic_Tx,HCQ,cholesterol,l,u){
  # Define the hazard function using all covariates
  lam <- exp(coefficients_gonadal$intercept +
               log(coefficients_gonadal$cumulative_steroid) * cumulative_steroid +
               log(coefficients_gonadal$cytotoxic_Tx) * cytotoxic_Tx +
               log(coefficients_gonadal$HCQ) * HCQ + 
               log(coefficients_gonadal$cholesterol) * cholesterol)
  
  haz <- function(time) { lam*(time/time) } # since this is an exponential model, the baseline hazard rate is constant with time
  
  # Integrate the hazard function to get the cumulative hazard within [l, u]
  ch <- integrate(haz, lower = l, upper = u)$value
  
  # Calculate the probability of death based on the cumulative hazard
  pr <- 1 - exp(-ch)
  
  return(pr)
}                      

#### 9. calculating the prob of diabetes ####
### !! this is an exponential model ###
probability_diabetes <- function(age,renalinv,cumulative_steroid,ethnicity,l,u){
  # Define the hazard function using all covariates
  lam <- exp(coefficients_diabetes$intercept +
               log(coefficients_diabetes$log_age) * log(age) +
               log(coefficients_diabetes$renal_inv) * renalinv +
               log(coefficients_diabetes$cumulative_steroid) * cumulative_steroid + 
               log(coefficients_diabetes$african_american) * ethnicity)
  
  haz <- function(time) { lam*(time/time) } # since this is an exponential model, the baseline hazard rate is constant with time
  
  # Integrate the hazard function to get the cumulative hazard within [l, u]
  ch <- integrate(haz, lower = l, upper = u)$value
  
  # Calculate the probability of death based on the cumulative hazard
  pr <- 1 - exp(-ch)
  
  return(pr)
}                    

#### 10. calculating the prob of malignancy ####
### !! this is an exponential model ###
probability_malignancy <- function(duration,age_at_diagnosis,SLICC,cholesterol,l,u){
  # Define the hazard function using all covariates
  lam <- exp(coefficients_malignancy$intercept +
               log(coefficients_malignancy$log_duration) * log(duration) +
               log(coefficients_malignancy$age_at_diagnosis) * age_at_diagnosis +
               log(coefficients_malignancy$SLICC) * SLICC + 
               log(coefficients_malignancy$cholesterol) * cholesterol)
  
  haz <- function(time) { lam*(time/time) } # since this is an exponential model, the baseline hazard rate is constant with time
  
  # Integrate the hazard function to get the cumulative hazard within [l, u]
  ch <- integrate(haz, lower = l, upper = u)$value
  
  # Calculate the probability of death based on the cumulative hazard
  pr <- 1 - exp(-ch)
  
  return(pr)
}          

#### 11. calculating the prob of cardiovascular damage ####
# note this is a Gompertz model

probability_cardio <- function(duration, age_at_diagnosis, AMS, serositisinv, cumulative_steroid, cholesterol, hypertension, l, u) {
  
  lam <- exp(coefficients_cardiovascular$intercept +
               log(coefficients_cardiovascular$log_duration) * log(duration) +
               log(coefficients_cardiovascular$age_at_diagnosis) * age_at_diagnosis +
               log(coefficients_cardiovascular$AMS) * AMS +
               log(coefficients_cardiovascular$serositisinv) * serositisinv + 
               log(coefficients_cardiovascular$cumulative_steroid) * cumulative_steroid + 
               log(coefficients_cardiovascular$cholesterol) * cholesterol + 
               log(coefficients_cardiovascular$hypertension) * hypertension)
  
  haz <- function(time) { lam * exp(coefficients_cardiovascular$distribution_param * time) }
  
  ch <- integrate(haz, lower = l, upper = u)$value
  pr <- 1 - exp(-ch)
  
  return(pr)
}

#### 12. calculating the prob of skin damage ####
# note this is a log-logistic model
probability_skin <- function(skininv,HCQ, ethnicity, smoker,l,u){
  # Define the hazard function using all covariates
  lam <- 1/exp(coefficients_skin$intercept +
               log(coefficients_skin$skininv) * skininv +
               log(coefficients_skin$HCQ) * HCQ +
               log(coefficients_skin$african_american) * ethnicity +
               log(coefficients_skin$smoker) * smoker
  )
  
  haz <- function(time) { ((lam^(1/coefficients_skin$distribution_param))*(time^((1/coefficients_skin$distribution_param)-1)))/(coefficients_skin$distribution_param*(1+((lam*time)^(1/coefficients_skin$distribution_param)))) } # note this is a log-logistic model
  
  # Integrate the hazard function to get the cumulative hazard within [l, u]
  ch <- integrate(haz, lower = l, upper = u)$value
  
  # Calculate the probability of death based on the cumulative hazard
  pr <- 1 - exp(-ch)
  
  return(pr)
} 

#### 13. calculating the prob of renal damage ####
# note this is a log-logistic model & time ratio
probability_renal <- function(AMS, renalinv,cholesterol,l,u){
  
  # Define the hazard function using all covariates
  lam <- 1/exp(coefficients_renal$intercept +
               log(coefficients_renal$AMS) * AMS +
               log(coefficients_renal$renal_inv) * renalinv +
               log(coefficients_renal$cholesterol) * cholesterol
  )
  
  # Define the hazard function based on the log-logistic distribution
  time_ratio <- function(time) { 
    ((lam^(1/coefficients_renal$distribution_param)) * 
       (time^((1/coefficients_renal$distribution_param)-1))) / 
      (coefficients_renal$distribution_param * 
         (1 + ((lam*time)^(1/coefficients_renal$distribution_param)))) 
  }
  
  # Integrate the hazard function to get the cumulative hazard within [l, u]
  ch <- integrate(time_ratio, lower = l, upper = u, subdivisions = 2000)$value
  
  # Calculate the probability of death based on the cumulative hazard
  pr <- 1 - exp(-ch)
  
  return(pr)
}

#### 14. calculating the prob of pulmonary damage ####
# note this is a log-logistic model
probability_pulmonary <- function(age,renalinv,serositisinv,dna_binding,SLICC,aCL,l,u){
  
  # Define the hazard function using all covariates
  lam <- 1/exp(coefficients_pulmonary$intercept +
               log(coefficients_pulmonary$log_age) * log(age) +
               log(coefficients_pulmonary$renal_inv) * renalinv +
               log(coefficients_pulmonary$serositisinv) * serositisinv +
               log(coefficients_pulmonary$dna_binding) * dna_binding +
               log(coefficients_pulmonary$SLICC) * SLICC +
               log(coefficients_pulmonary$aCL) * aCL
  ) 
  
 # note this is a log-logistic model
  time_ratio <- function(time) { 
    (lam^(1/coefficients_pulmonary$distribution_param) * time^((1/coefficients_pulmonary$distribution_param) - 1)) /
      (coefficients_pulmonary$distribution_param * (1 + (lam * time)^(1/coefficients_pulmonary$distribution_param)))
  }
  
  
  # Integrate the hazard function to get the cumulative hazard within [l, u]
  ch <- integrate(time_ratio, lower = l, upper = u, subdivisions = 2000)$value
  
  # Calculate the probability of damage based on the cumulative hazard
  pr <- 1 - exp(-ch)
  
  return(pr)
} 


#### 15. calculating the clean utility ####
clean_u <- function(age,ethnicity,current_sle){
  clean_utility <- coefficients_clean_utility$intercept + 
    coefficients_clean_utility$log_age*log(age) + 
    coefficients_clean_utility$african_american*ethnicity + 
    coefficients_clean_utility$current_sle*current_sle
  return(clean_utility)
}
