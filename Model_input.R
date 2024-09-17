###################
#### SLE model ####
###################

### Model input ###

# Define the number of individuals
n <- 5000

# load dplyr package to use tibble function to define the data structure.
library(dplyr)
library(tidyr)
# install.packages('readxl')
library(readxl)

# Define the path to the Excel file
file_path <- "Life table.xlsx"  # Adjust the path as necessary

# Read the entire sheet without specifying a complex range
life_table_raw <- read_excel(path = file_path, sheet = 1, col_names = FALSE)

# Assign appropriate column names based on assumed structure
# Adjust column indices if your structure is different
colnames(life_table_raw) <- c("Age_Male", "qx_Male", "Age_Female", "qx_Female")

# Reshape the life table to long format for easier access
life_table <- life_table_raw %>%
  select(Age = Age_Male, qx_Male, qx_Female)

# Reshape to long format
life_table_long <- life_table %>%
  pivot_longer(
    cols = c(qx_Male, qx_Female),
    names_to = "Gender",
    names_prefix = "qx_",
    values_to = "qx"
  )

# Ensure no duplicates
life_table_long <- life_table_long %>%
  distinct(Age, Gender, .keep_all = TRUE)

# Reshape back to wide format (if needed)
life_table_wide <- life_table_long %>%
  pivot_wider(names_from = Gender, values_from = qx)

# Inspect the result
head(life_table_wide)

# Define SMR values directly as simple variables for the age ranges
# this would make it easier to do DSA and PSA later
smr_16_24 <- 19.2
smr_25_39 <- 8.0
smr_40_59 <- 3.7
smr_60_plus <- 1.4

# add SMR
# Add SMR values directly into the life_table_wide based on age ranges
life_table_wide <- life_table_wide %>%
  mutate(SMR = case_when(
    Age >= 16 & Age <= 24 ~ smr_16_24,
    Age >= 25 & Age <= 39 ~ smr_25_39,
    Age >= 40 & Age <= 59 ~ smr_40_59,
    Age >= 60 ~ smr_60_plus,
    TRUE ~ 1  # Default SMR of 1 for ages not covered explicitly
  ))

set.seed(123)

# Number of total population (e.g., 1000)
total_population <- 50000

# Probability distributions for each organ damage score
probabilities <- list(
  cardiovascular = c(0.94, 0.051, 0.009, 0.000, 0.000),
  diabetes = c(0.976, 0.024, 0.000, 0.000, 0.000),
  gastrointestinal = c(0.964, 0.034, 0.002, 0.000, 0.000),
  malignancy = c(0.996, 0.004, 0.000, 0.000, 0.000),
  musculoskeletal = c(0.874, 0.088, 0.032, 0.004, 0.002),
  neuropsychiatric = c(0.889, 0.092, 0.015, 0.004, 0.000),
  ocular = c(0.936, 0.062, 0.002, 0.000, 0.000),
  peripheral_vascular = c(0.944, 0.051, 0.004, 0.001, 0.000),
  premature_gonadal_failure = c(0.989, 0.011, 0.000, 0.000, 0.000),
  pulmonary = c(0.970, 0.028, 0.002, 0.000, 0.000),
  renal = c(0.974, 0.026, 0.000, 0.000, 0.000),
  skin = c(0.921, 0.071, 0.006, 0.002, 0.000)
)

# Create the large baseline population
total_population_pool <- tibble(
  id = 1:total_population, 
  age = pmax(pmin(rgamma(total_population, shape = 10.292, scale = 3.453), 100), 18),  # Gamma distribution for age
  gender = rbinom(total_population, 1, prob = 0.05),  # Bernoulli for gender
  ethnicity = rbinom(total_population, 1, prob = 0.04),  # Bernoulli for ethnicity
  disease_duration = NA,
  age_at_diagnosis = NA,
  SLEDAI_baseline = pmax(rgamma(total_population, shape = 6.717, scale = 1.454), 6),  # SLEDAI baseline
  dna_binding = rbinom(total_population, 1, prob = 0.74),  # DNA binding
  complement = rbinom(total_population, 1, prob = 0.59),
  vasculitis = rbinom(total_population, 1, prob = 0.08),
  neuropsychiatric_inv = rbinom(total_population, 1, prob = 0.02),
  renal_inv = rbinom(total_population, 1, prob = 0.20),
  serositis_inv = rbinom(total_population, 1, prob = 0.04),
  haematological_inv = rbinom(total_population, 1, prob = 0.07),
  skin_inv = rbinom(total_population, 1, prob = 0.82),
  cardiovascular_score = apply(rmultinom(total_population, 1, probabilities$cardiovascular), 2, function(x) which(x == 1) - 1),
  diabetes_score = apply(rmultinom(total_population, 1, probabilities$diabetes), 2, function(x) which(x == 1) - 1),
  gastrointestinal_score = apply(rmultinom(total_population, 1, probabilities$gastrointestinal), 2, function(x) which(x == 1) - 1),
  malignancy_score = apply(rmultinom(total_population, 1, probabilities$malignancy), 2, function(x) which(x == 1) - 1),
  musculoskeletal_score = apply(rmultinom(total_population, 1, probabilities$musculoskeletal), 2, function(x) which(x == 1) - 1),
  neuropsychiatric_score = apply(rmultinom(total_population, 1, probabilities$neuropsychiatric), 2, function(x) which(x == 1) - 1),
  ocular_score = apply(rmultinom(total_population, 1, probabilities$ocular), 2, function(x) which(x == 1) - 1),
  peripheral_vascular_score = apply(rmultinom(total_population, 1, probabilities$peripheral_vascular), 2, function(x) which(x == 1) - 1),
  premature_gonadal_failure_score = apply(rmultinom(total_population, 1, probabilities$premature_gonadal_failure), 2, function(x) which(x == 1) - 1),
  pulmonary_score = apply(rmultinom(total_population, 1, probabilities$pulmonary), 2, function(x) which(x == 1) - 1),
  renal_score = apply(rmultinom(total_population, 1, probabilities$renal), 2, function(x) which(x == 1) - 1),
  skin_score = apply(rmultinom(total_population, 1, probabilities$skin), 2, function(x) which(x == 1) - 1),
  steroid_baseline = pmin(rgamma(total_population, shape = 2.120, scale = 5.977), 30),  # Steroid use
  cytotoxic_Tx = rbinom(total_population, 1, prob = 0.42),  # Cytotoxic Tx
  HCQ = sample(c(0, 1), total_population, replace = TRUE),
  smoker = rbinom(total_population, 1, prob = 0.389),
  cholesterol = pmin(rgamma(total_population, shape = 13.63, scale = 13.615), 240),
  hypertension = rbinom(total_population, 1, prob = 0.158),
  aCL = rbinom(total_population, 1, prob = 0.9367),
  anticoagulant_positive = rbinom(total_population, 1, prob = 0.096),
  anaemia = rbinom(total_population, 1, prob = 0.5),
  infection = rbinom(total_population, 1, prob = 0.6634),
  weights = pmax(pmin(rgamma(total_population, shape = 1.96, scale = 35.897), 97.9), 39)
)

# Calculate disease_duration and other fields
total_population_pool <- total_population_pool %>%
  mutate(disease_duration = rgamma(total_population, shape = 0.981, scale = 5.404) %>% 
           pmin(age - 5),
         age_at_diagnosis = age - disease_duration,
         SLICC_baseline = cardiovascular_score + diabetes_score + gastrointestinal_score + malignancy_score +
           musculoskeletal_score + neuropsychiatric_score + ocular_score + peripheral_vascular_score +
           premature_gonadal_failure_score + pulmonary_score + renal_score + skin_score
  )

# Sample from the large pre-generated baseline population
baseline <- total_population_pool %>% slice(1:n)

# Print the sampled patients
print(baseline)

#### setting up the SLICC score per organ
SLICC_score <- list(
  Cardiovascular=1.42,
  Diabetes=1,
  Gastrointestinal=1.09,
  Malignancy=1,
  Musculoskeletal=1.41,
  Neuropsychiatric=1.37,
  Ocular=1.23,
  Peripheral_vascular=1.21,
  gonadal_failure=1,
  Pulmonary=1.31,
  Renal=1.83,
  Skin=1.14
)

#### setting up the cost parameters
# drug costs (annual)
belimumab_cost <- (10.125*baseline$weights+154)*13 # 10.125 is the per mg price, 154 is the admin price per infusion
TAK279_cost <- 8994.64

# disease activity related costs
disease_activity_costs <- list(
  SS_0 = 1487.70,
  SS_1 = 1659.96,
  SS_2 = 1832.20,
  SS_3 = 1954.25,
  SS_4 = 2026.02,
  SS_5 = 2097.79,
  SS_6 = 2169.57,
  SS_7 = 2241.34,
  SS_8 = 2313.10,
  SS_9 = 2396.87,
  SS_10 = 2492.57,
  SS_11 = 2588.26,
  SS_12 = 2683.96,
  SS_13_20 = 2779.69
)

# organ damage related costs
# Year 1
organ_damage_costs_y1 <- list(
  Cardio = 5392.14,                # Cardiovascular
  Diabetes = 3054.63,
  Gastro = 3559.13,                # Gastrointestinal
  Malignancy = 8154.87,
  Musculoskeletal = 8251.40,
  Neuro = 7838.83,                 # Neuropsychiatric
  Ocular = 2080.09,
  Peripheral_vascular = 3769.44,
  Gonadal = 0,                     # Gonadal failure
  Pulmonary = 17109.59,
  Renal = 2835.13,
  Skin = 0
)

# Year 2 and onwards
organ_damage_costs_y2 <- list(
  Cardio = 1490.54,                # Cardiovascular
  Diabetes = 3054.63,
  Gastro = 0.00,                   # Gastrointestinal
  Malignancy = 0.00,
  Musculoskeletal = 2742.04,
  Neuro = 3201.73,                 # Neuropsychiatric
  Ocular = 27.58,
  Peripheral_vascular = 814.80,
  Gonadal = 0.00,                  # Gonadal failure
  Pulmonary = 17165.90,
  Renal = 4184.31,
  Skin = 0.00
)

### setting up utility multiplier for organ damages
# utility multiplier for year 1
organ_damage_multiplier_y1 <- list(
  cardio = 0.779,
  diabetes = 0.91,
  gastro = 0.786,
  malignancy = 0.837,
  musculoskeletal = 0.655,
  neuro = 0.713,
  ocular = 0.974,
  peripheral_vascular = 0.863,
  gonadal = 1,
  pulmonary = 0.713,
  renal = 0.972,
  skin = 0.943
)

# utility multiplier for year 2 onwards
organ_damage_multiplier_y2 <- list(
  cardio = 0.806,
  diabetes = 0.91,
  gastro = 0.906,
  malignancy = 0.837,
  musculoskeletal = 0.729,
  neuro = 0.772,
  ocular = 0.992,
  peripheral_vascular = 0.873,
  gonadal = 1,
  pulmonary = 0.719,
  renal = 0.955,
  skin = 0.943
)


### year 1 disease activity
# this is used assuming patients who 1. are considered non-responders 2. are responders but discontinue in the first year
Y1_SS_score <- list(
  SOC = -0.349,
  Tx = -0.343-0.280,
  Tx_1 = -0.343,
  Tx_2 = -0.343,
  Tx_1_responder = -0.280,
  Tx_2_responder = -0.280
)


# responding rate for belimumab and TAK279
responder_belimumab <- 0.5240
# responder_TAK279 <- 0.5820
responder_TAK279 <- 0.5820

# discontinuation rate for belimumab and TAK279
### year 1
discontinue_belimumab_y1 <- 0.08
discontinue_TAK279_y1 <- 0.08 #0.088

# discontinuation rate for belimumab and TAK279
### year 2 on wards
discontinue_belimumab_y2 <- 0.117
discontinue_TAK279_y2 <- 0.117

# model calibration factors for belimumab and TAK279
calibration_belimumab <- 0.415
calibration_TAK279 <- 0.415
calibration_duration <- 6

# discount rates
qaly_d <- 0.035
cost_d <- 0.035
