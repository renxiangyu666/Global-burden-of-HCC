# =============================================================================================
# Population Attributable Fraction (PAF) Calculation in 2022 by 10, 15, 20 year latency period 
# =============================================================================================
# Purpose: Calculate PAFs and attributable HCC (AHCC) cases of nine risk factors for HCC globally
# Input: 
## Country-level age-standardized prevalence and region-level pooled RR data by  risk and sex classifications;
## Country-level HCC case in 2022 estimated by GLOBOCAN 2022 by sex;
# Output: Individual and combined PAFs results by subregion, region, and global levels
# ---------------------------------------------------------------------------------------------

# Load Libraries & Set Working Directory
library(psych)
library(readxl)
library(dplyr)
library(tidyr)
library(forecast)
library(purrr)
setwd('C:/Users/24399/Desktop/Core code for HCC and PAF estimation/PAF calculation/')

# ----------------------------------------------------
##PART 1: Merge Risk Factor Data and Standardize Nomenclature 
# ----------------------------------------------------

# 1. Load and Process age-standardized Alcohol & Smoking Dataset 
P1 <- read.csv("original files/alcohol and smoke.csv") %>% 
  # Standardize column names
  rename(cause_id = rei_id,
         cause_name = rei_name) %>% 
  # Convert percentage values to proportions (0-1 scale)
  mutate(across(c(val, upper, lower), ~ . * 0.01))

# Load and Process diabetes,HBV,HCV,NASH Dataset 
P2 <- read.csv("original files/diabetes,HBV,HCV,NASH.csv")

# 2. Combine Datasets 
combined_data <- bind_rows(P1, P2) %>% 
  # Recode categorical variables
  mutate(
    sex_name = recode(sex_name,
                      "Both" = "B",
                      "Female" = "F", 
                      "Male" = "M"),
    cause_name = recode(cause_name,
                        "Chronic hepatitis B including cirrhosis" = "HBV",
                        "Chronic hepatitis C including cirrhosis" = "HCV",
                        "Diabetes mellitus type 2" = "Diabetes",
                        "High alcohol use" = "Alcohol",
                        "Nonalcoholic fatty liver disease including cirrhosis" = "NAFLD/NASH",
                        "Smoking" = "Smoke")
  ) %>% 
  # Standardize prevalence column names
  rename(prevalence = val,
         prevalenceH = upper,
         prevalenceL = lower)

# 3. Export Merged Dataset 
write.csv(combined_data, 
          'prevalence predic 2022/GBD_combine.csv', 
          row.names = FALSE)

# -----------------------------------------------------------
##PART 2: Prevalence Forecasting for 2022 (Country Level) 
# -----------------------------------------------------------

# 1. Load preprocessed dataset from Part 1
prevalence_GBD <- read.csv('prevalence predic 2022/GBD_combine.csv') %>% 
  filter(sex_name %in% c("M", "F")) %>% 
  mutate(
    id1 = paste0(location_name, cause_name, sex_name),  # Create unique category ID
    year = as.integer(year)
  ) %>% 
  arrange(id1,year)  # Ensure temporal ordering

#### 2. prevalence prediction
ts_data <- ts(prevalence_GBD$prevalence,start=c(1990,1), frequency = 1)
categories <- unique(prevalence_GBD$id1)
forecasts <- list()
for (cat in categories) {
  cat_data <- subset(prevalence_GBD, id1 == cat)
  cat_ts <- ts(cat_data$prevalence, frequency = 1)
  holt_model <- holt(cat_ts)
  forecast_result <- forecast(holt_model, h = 1)
  forecasts[[cat]] <- forecast_result
}

P<-read.csv("prevalence predic 2022/forecast_AlcoholAfricaB.csv")

#Prediction of Prevalence for Various Categories in 2022
for (cat in categories) {
  forecast_df <- data.frame(
    prevalence = forecasts[[cat]]$mean,
    id1 = cat,
    year= 2022
  )
  P<- rbind(P, forecast_df)
}

#### 3. prevalence prediction:95% CI lower (PrevalenceL)
ts_data <- ts(prevalence_GBD$prevalenceL,start=c(1990,1), frequency = 1)
categories <- unique(prevalence_GBD$id1)
forecasts <- list()
for (cat in categories) {
  cat_data <- subset(prevalence_GBD, id1 == cat)
  cat_ts <- ts(cat_data$prevalenceL, frequency = 1)
  holt_model <- holt(cat_ts)
  forecast_result <- forecast(holt_model, h = 1)
  forecasts[[cat]] <- forecast_result
}

PL<-read.csv("prevalence predic 2022/forecast_AlcoholAfricaB_L region.csv")

#Prediction of PrevalenceL for Various Categories in 2022
for (cat in categories) {
  forecast_df <- data.frame(
    prevalenceL = forecasts[[cat]]$mean,
    id1 = cat,
    year= 2022
  )
  PL <- rbind(PL, forecast_df)
}

#### 4. prevalence prediction:95% CI upper (PrevalenceH)
ts_data <- ts(prevalence_GBD$prevalenceH,start=c(1990,1), frequency = 1)
categories <- unique(prevalence_GBD$id1)
forecasts <- list()
for (cat in categories) {
  cat_data <- subset(prevalence_GBD, id1 == cat)
  cat_ts <- ts(cat_data$prevalenceH, frequency = 1)
  holt_model <- holt(cat_ts)
  forecast_result <- forecast(holt_model, h = 1)
  forecasts[[cat]] <- forecast_result
}

PH<-read.csv("prevalence predic 2022/forecast_AlcoholAfricaB_U region.csv")

#Prediction of PrevalenceH for Various Categories in 2022
for (cat in categories) {
  forecast_df <- data.frame(
    prevalenceH = forecasts[[cat]]$mean,
    id1 = cat,
    year= 2022
  )
  PH <- rbind(PH, forecast_df)
}

## 5. Merge the P/L/H forecast data for 2022 with the original data from 1990-2021
P_2022<-cbind(P,PL[1],PH[1])
combine <- bind_rows(prevalence_GBD, P_2022[-1,])
filled_data <- combine %>%
  group_by(id1) %>% 
  mutate(across(
    .cols = -year,
    .fns = ~ {
      value_2021 <- .x[year == 2021] 
      ifelse(year == 2022 & is.na(.x), value_2021, .x)  
    }
  )) %>%
  ungroup()

write.csv(filled_data,  "prevalence predic 2022/pred_GBD1990-2022.csv")


# ----------------------------------------------------------------------
##PART 3: Country ISO Standardization & Overseas Territories Imputation 
# ----------------------------------------------------------------------
### 1. Process GBD Risk Factors ------------------------------------------------
# Load preprocessed prevalence data
all_data <- read.csv("prevalence predic 2022/pred_GBD1990-2022.csv") %>% 
  mutate(location_name = recode(location_name,
                                "Republic of Côte d'Ivoire" = "Ivoire",
                                "Côte d'Ivoire" = "Ivoire",
                                "French Republic" = "France"
  ))

# Load ISO code reference
ISO_ref <- readxl::read_xlsx("original files/COD3all.xlsx", sheet = "GBD2021")

# Create French overseas territories template
french_territories <- tribble(
  ~location_name,           ~ISO,
  "France, Guadeloupe",     "GLP",
  "France, La Réunion",     "REU",  # Fixed encoding issue
  "France, Martinique",     "MTQ",
  "French Guyana",          "GUF",
  "French Polynesia",       "PYF",
  "New Caledonia",          "NCL"
)

# Generate overseas territories data
french_data <- map_df(french_territories$location_name, ~ {
  all_data %>% 
    filter(location_name == "France") %>% 
    mutate(location_name = .x,
           ISO = french_territories$ISO[french_territories$location_name == .x])
})

# Combine original and imputed data
all_data_iso <- all_data %>% 
  left_join(ISO_ref, by = "location_name") %>% 
  bind_rows(french_data)

### 2. Process Obesity Data ---------------------------------------------------
obesity_data <- read.csv("original files/NCD_RisC_Lancet_2024_BMI_age_standardised_country.csv") %>% 
  mutate(
    location_name = recode(location_name,
                           "Cote d'Ivoire" = "Ivoire",
                           "Côte d'Ivoire" = "Ivoire"
    ),
    sex_name = recode(sex_name,
                      "Women" = "F",
                      "Men" = "M"
    )
  )

# Create special territory mappings
special_territories <- tribble(
  ~original,                  ~new_name,                      ~new_iso,
  "France",                   "France, Guadeloupe",          "GLP",
  "France",                   "France, La Réunion",          "REU",
  "France",                   "France, Martinique",          "MTQ",
  "France",                   "French Guyana",               "GUF",
  "France",                   "New Caledonia",               "NCL",
  "Zambia",                   "Tanzania, United Republic of","TZA",
  "Montenegro",               "Bosnia Herzegovina",          "BIH",
  "Puerto Rico",              "Guam",                        "GUM"
)

# Generate imputed obesity data
imputed_obesity <- map_df(1:nrow(special_territories), ~ {
  obesity_data %>% 
    filter(location_name == special_territories$original[.x]) %>% 
    mutate(location_name = special_territories$new_name[.x],
           ISO = special_territories$new_iso[.x])
})

obesity_iso <- bind_rows(obesity_data, imputed_obesity)

### 3. Data Integration ---------------------------------------------------
# Load additional datasets
clonorchiasis <- read.csv("original files/C. sinensis.csv")

aflatoxin <- read.csv("original files/B1_allcountry.csv") %>% 
  mutate(cause_name = "Aflatoxin B1",
         PAFU = as.numeric(PAFH)) %>% 
  select(-case)

# Final integration
combined_df <- bind_rows(
  all_data_iso,
  obesity_iso,
  clonorchiasis,
  aflatoxin
) %>% 
  select(-Region, -Continent) %>% 
  # Validate gender categories
  mutate(sex_name = if_else(sex_name %in% c("M", "F"), sex_name, NA_character_))

### Export Results -------------------------------------------------------------
writexl::write_xlsx(combined_df, 'prevalence1990-2022/original files/prevalence1990-2022.xlsx')


# ------------------------------------
## PART4: HCC Case Data Integration 
# ------------------------------------

### 1. Data Preparation 
# Load HCC case data from 1990 to 2022
mat <- read.csv("original files/HCC case 1990_2022.csv")

# Recode Sex values (1 = M, 2 = F, 3 = B) and aggregate cases by Time, Alpha.3.code, Sex, Region, and Continent
mat <- mat %>%
  mutate(Sex = recode(Sex, "1" = "M", "2" = "F", "3" = "B")) %>%
  group_by(Time, Alpha.3.code, Sex, Region, Continent) %>%
  summarise(case = round(sum(case),0), .groups = 'drop') %>%
  ungroup()

# 2. Extract 2022 data and create lagged datasets for 2012, 2007, and 2002
mat_2022 <- mat %>% filter(Time == "2022")
mat_2012 <- mat_2022 %>% mutate(Time = 2012)
mat_2007 <- mat_2022 %>% mutate(Time = 2007)
mat_2002 <- mat_2022 %>% mutate(Time = 2002)

# Combine lagged datasets into one
mat <- rbind(mat_2012, mat_2007, mat_2002)

# Create a unique ID for each row
mat$id <- paste0(mat$Alpha.3.code, mat$Time, mat$Sex)
combined_df$id <- paste0(combined_df$ISO, combined_df$year, combined_df$sex_name)

# 3. Join the combined dataset with the case data
all_datacase <- left_join(combined_df, mat %>% select(id, Region, Continent, case), by = "id")
merged_data <- all_datacase[!is.na(all_datacase$case), ]

# 4. Export final dataset
write.csv(merged_data, 'original files/HCC_case_integrated.csv', row.names = FALSE)


# ----------------------------------------------------------------------
# PART 5: PAF Calculation Pipeline
# Purpose: Calculate population attributable fraction (PAF) 
# ----------------------------------------------------------------------

# 1. Data Preparation & RR Matching 
# Load reference data
RR <- read_xlsx("original files/RR_Matching_File.xlsx", sheet = "RR")
merged_data <- read.csv("original files/HCC_case_integrated.csv")

# Separate aflatoxin (pre-calculated) from other factors
paf_aflatoxin <- merged_data %>% filter(cause_name == "Aflatoxin B1")
paf_other <- merged_data %>% filter(cause_name != "Aflatoxin B1")

# Create matching ID and merge with RR values
paf_other <- paf_other %>%
  mutate(id2 = paste0(cause_name, Continent, sex_name)) %>%
  left_join(RR, by = "id2")

# 2. PAF Calculation at country level
# Core PAF formula: PAF = [p(RR-1)] / [p(RR-1)+1]
calculate_paf <- function(prevalence, rr) {
  round(prevalence * (rr - 1) / (prevalence * (rr - 1) + 1), 3)
}

# Calculate estimates for all factors
paf_other <- paf_other %>%
  mutate(
    PAF = calculate_paf(prevalence, RR),
    PAFL = calculate_paf(prevalence, CIL),
    PAFU = calculate_paf(prevalence, CIH)
  )

# Combine all factors
PAF_combined <- bind_rows(paf_aflatoxin, paf_other)

# Calculate attributable cases
PAF_combined <- PAF_combined %>%
  mutate(
    AHCC = case * PAF,
    AHCCL = case * PAFL,
    AHCCU = case * PAFU
  )

# 3. Sex Aggregation -----------------------------------------------------
# Create 'Both Sexes' estimates
PAF_sex_aggregated <- PAF_combined %>%
  group_by(year, cause_name, Region, ISO, Continent) %>%
  summarise(
    case = sum(case),
    AHCC = sum(AHCC),
    AHCCL = sum(AHCCL),
    AHCCU = sum(AHCCU),
    .groups = 'drop'
  ) %>%
  mutate(
    sex_name = "B",
    PAF = AHCC / case,
    PAFL = AHCCL / case,
    PAFU = AHCCU / case
  )

# Final combined dataset
PAF_final <- bind_rows(PAF_combined, PAF_sex_aggregated)

## 4. Export Results
# Unformatted Version
writexl::write_xlsx(PAF_final, 'individual PAF results/Unformatted_PAF_Country_2002-2012.xlsx')

# Formatting & Export -------------------------------------------------
# Function for consistent formatting
format_paf_output <- function(df) {
  df %>%
    mutate(
      across(c(case, AHCC, AHCCL, AHCCU), ~round(.)),
      across(c(PAF, PAFL, PAFU), ~sprintf("%.1f", .x * 100)),
      across(c(case, AHCC, AHCCL, AHCCU), ~format(., big.mark = ",")),
      PAF_CI = paste0(PAF, " (", PAFL, "–", PAFU, ")")
    )
}

# Apply formatting
PAF_formatted <- PAF_final %>%
  format_paf_output() %>%
  # Order factors
  mutate(
    cause_name = factor(cause_name, 
                        levels = c("HBV", "HCV", "C. sinensis", "NAFLD/NASH", 
                                   "Obesity", "Diabetes", "Alcohol", "Smoke", "Aflatoxin B1")
    )
  ) %>%
  arrange(cause_name)

# Export results
writexl::write_xlsx(PAF_formatted, 'individual PAF results/Formatted_PAF_country_2002-2012.xlsx')


## 5. Single-Factor PAF Aggregation Pipeline (Subregion/Region/Global)
# Load unformatted country-level PAF results
PAF <- read_xlsx("individual PAF results/Unformatted_PAF_Country_2002-2012.xlsx")

# Aggregate by Subregion Level
PAF_subregion <- PAF %>%
  group_by(year, cause_name, Region, sex_name) %>%
  summarise(
    case = sum(case),
    AHCC = sum(AHCC),      # Attributable HCC cases
    AHCCL = sum(AHCCL),    # Lower confidence bound
    AHCCU = sum(AHCCU),    # Upper confidence bound
    .groups = 'drop'
  ) %>%
  mutate(CATregion = "subregion")

# Aggregate by Continent Level

PAF_region <- PAF %>%
  group_by(year, cause_name, Continent, sex_name) %>%
  summarise(
    case = sum(case),
    AHCC = sum(AHCC),
    AHCCL = sum(AHCCL),
    AHCCU = sum(AHCCU),
    .groups = 'drop'
  ) %>%
  mutate(
    CATregion = "region",
    Region = Continent  # Standardize column name
  )

# 6. Calculate Global Totals
PAF_globe <- PAF %>%
  group_by(year, cause_name, sex_name) %>%
  summarise(
    case = sum(case),
    AHCC = sum(AHCC),
    AHCCL = sum(AHCCL),
    AHCCU = sum(AHCCU),
    .groups = 'drop'
  ) %>%
  mutate(
    Region = "Globe",
    CATregion = "region"  # Maintain hierarchy consistency
  )

# 6. Combine All Aggregations
PAF_all <- bind_rows(PAF_subregion, PAF_region, PAF_globe)

# 7. Calculate PAF Metrics
PAF_all <- PAF_all %>%
  mutate(
    PAF = round(AHCC / case, 3),       # Population Attributable Fraction
    PAFL = round(AHCCL / case, 3),     # Lower bound
    PAFU = round(AHCCU / case, 3)      # Upper bound
  )


# 8. Export Unformatted Results
writexl::write_xlsx(PAF_all, 'individual PAF results/Unformatted_PAF_Region_2002-2012.xlsx')


# ----------------------------------------------------------------------
# PART 6: Combined PAF Calculation Pipeline at country level
# ----------------------------------------------------------------------

# 1. Load metabolic risk factor data, single-factor PAF results
metabolic <- read_xlsx("original files/combined metaboPAF.xlsx", sheet = 1)

PAF_Country <- read_xlsx("individual PAF results/Unformatted_PAF_Country_2002-2012.xlsx", sheet = 1)
PAF_Country<-PAF_Country %>%
  filter(!cause_name %in% c("NAFLD/NASH", "Obesity","Diabetes"),!sex_name %in% c("B"))

# 2. Data Filtering & Merging
# Merge datasets using common columns
common_cols <- intersect(colnames(PAF_Country), colnames(metabolic))
merged_df <- bind_rows(
  PAF_Country %>% select(any_of(common_cols)),
  metabolic %>% select(any_of(common_cols))
)

# Final filtering after merge
PAF_comb <- merged_df %>%
  filter(!cause_name %in% c("metabolic", "NAFLD/NASH", "Obesity"),
         !sex_name %in% c("B"))

# 3. Risk Category Mapping
# Load risk category mapping,Join category information
CATcause <- readxl::read_xlsx("original files/RR_Matching_File.xlsx", sheet = "CATcause")
PAF_comb <- left_join(PAF_comb, CATcause, by = "cause_name")

# 4. Combined PAF Calculation
# Calculate combined PAF at country level
# All Risk Factors
combined_paf_country <- PAF_comb  %>%
  group_by(ISO, Region, Continent, year, sex_name) %>%
  summarise(
    # Handle missing values 
    PAF = coalesce(PAF, 0),
    PAFL = coalesce(PAFL, 0), 
    PAFU = coalesce(PAFU, 0),
    case = case,
    # Cumulative product calculation
    cumprod_PAF = cumprod(1 - PAF),
    cumprod_PAFL = cumprod(1 - PAFL),
    cumprod_PAFU = cumprod(1 - PAFU)
  ) %>%
  mutate(
    # Combined PAF calculation
    combined_PAF = 1 - cumprod_PAF,
    combined_PAFL = 1 - cumprod_PAFL,
    combined_PAFU = 1 - cumprod_PAFU,
    CATcause = "All risk factors"
  ) %>%
  slice(n())  # Keep last observation per group

# Group Risk Factors (Infectious factors, Metabolic factors,Behavioral/toxic factors)
specific_paf_country <- PAF_comb %>%
  group_by(ISO, Region, Continent, year, sex_name, CATcause) %>%
  summarise(
    PAF = coalesce(PAF, 0),
    PAFL = coalesce(PAFL, 0),
    PAFU = coalesce(PAFU, 0),
    case = case,
    cumprod_PAF = cumprod(1 - PAF),
    cumprod_PAFL = cumprod(1 - PAFL),
    cumprod_PAFU = cumprod(1 - PAFU)
  ) %>%
  mutate(
    combined_PAF = 1 - cumprod_PAF,
    combined_PAFL = 1 - cumprod_PAFL,
    combined_PAFU = 1 - cumprod_PAFU
  ) %>%
  slice(n())

# Data Integration
final_paf_dataset <- bind_rows(combined_paf_country, specific_paf_country)%>%
  mutate(
    AHCC = round(case * round(combined_PAF, 3)),
    AHCCL = round(case * round(combined_PAFL, 3)),
    AHCCU = round(case * round(combined_PAFU, 3))
  )

# Both Sexes Calculation
both_sexes_dataset <- final_paf_dataset %>% 
  group_by(ISO,Region,Continent,year,CATcause) %>%
  summarise(AHCC = coalesce(AHCC, 0),
            AHCCL = coalesce(AHCCL, 0),
            AHCCU = coalesce(AHCCU, 0),
            case = sum(case),
            AHCC = sum(AHCC),
            AHCCL = sum(AHCCL),
            AHCCU = sum(AHCCU),
            sex_name = "B")%>%
  mutate(
    combined_PAF = round(AHCC / case, 3),
    combined_PAFL = round(AHCCL / case, 3),
    combined_PAFU = round(AHCCU / case, 3)
  ) %>%
  slice(n())

PAFcountry_all <- bind_rows(final_paf_dataset, both_sexes_dataset)

# 5. Output Formatting
# Unformatted output
writexl::write_xlsx(PAFcountry_all, 'combined PAF results/Unformatted_combinedPAF_Country_2002-2012.xlsx')


# ----------------------------------------------------------------------
# PART 7: Calculate combined PAF at subregion/region/globe level
# ----------------------------------------------------------------------

# 1. Environment Setup & Data Loading
PAF <- readxl::read_xlsx("combined PAF results/Unformatted_combinedPAF_Country_2002-2012.xlsx", sheet = 1)

# 2. Multi-level Aggregation
# Subregional aggregation
PAF_subregion <- PAF %>%
  group_by(year, CATcause, Region, sex_name) %>%
  summarise(
    case = sum(case),
    AHCC = sum(AHCC),       # Attributable cases
    AHCCL = sum(AHCCL),     # Lower bound
    AHCCU = sum(AHCCU),     # Upper bound
    .groups = 'drop'
  ) %>%
  mutate(CATregion = "subregion")

# Continental aggregation
PAF_region <- PAF %>%
  group_by(year, CATcause, Continent, sex_name) %>%
  summarise(
    case = sum(case),
    AHCC = sum(AHCC),
    AHCCL = sum(AHCCL),
    AHCCU = sum(AHCCU),
    .groups = 'drop'
  ) %>%
  mutate(
    CATregion = "region",
    Region = Continent  # Standardize column name
  )

# Global aggregation
PAF_globe <- PAF %>%
  group_by(year, CATcause, sex_name) %>%
  summarise(
    case = sum(case),
    AHCC = sum(AHCC),
    AHCCL = sum(AHCCL),
    AHCCU = sum(AHCCU),
    .groups = 'drop'
  ) %>%
  mutate(
    Region = "Globe",
    CATregion = "region"  # Maintain hierarchy
  )

# 3. Combine All Levels
PAF_all <- bind_rows(PAF_subregion, PAF_region, PAF_globe) %>%
  mutate(
    PAF = round(AHCC / case, 3),     # Calculate PAF percentage
    PAFL = round(AHCCL / case, 3),   # Lower bound
    PAFU = round(AHCCU / case, 3)    # Upper bound
  )

# 4. Export Raw Results
writexl::write_xlsx(PAF_all, 'combined PAF results/Unformatted_combinedPAF_region2002-2012.xlsx')


