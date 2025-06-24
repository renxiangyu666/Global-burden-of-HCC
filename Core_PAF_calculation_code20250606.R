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


# ------------------------------------
## PART1: HCC Case Data Integration 
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
# PART 2: PAF Calculation Pipeline
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

# Combine All Aggregations
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
# PART 3: Combined PAF Calculation Pipeline at country level
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
# PART 4: Calculate combined PAF at subregion/region/globe level
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


