############################HCC country
############################Purpose: calculate the HCC cases/deaths and Age standardized ratio(ASR) by country 
############################Input:GLOBOCAN 2022 data, population structure data, HCC ratio data, and standardized age structure
############################Output: formatted HCC cases/deaths and ASR result by country and sex

#####data for incidence and mortality was download from GLOBOCAN 2022
library(dplyr)
library(writexl)
library(readxl)
setwd("C:/Users/24399/Desktop/Core code for HCC and PAF estimation/HCC burden/")
######### PART1: Incidence
male_age0_4<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-males-age-0-4-in-2022-liver-and-intrahepatic-bile-ducts.csv")
male_age0_4 <-male_age0_4 %>%
  rename(N0_4 = Number)
male_age5_9<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-males-age-5-9-in-2022-liver-and-intrahepatic-bile-ducts.csv")
male_age5_9 <-male_age5_9 %>%
  rename(N5_9 = Number)
male_age10_14<-read_xlsx("GLOBOCAN 2022 data/countries/incidence/dataset-inc-males-age-10-14-in-2022-liver-and-intrahepatic-bile-ducts.xlsx")
male_age10_14 <-male_age10_14 %>%
  rename(N10_14 = Number)
male_age15_19<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-males-age-15-19-in-2022-liver-and-intrahepatic-bile-ducts.csv")
male_age15_19 <-male_age15_19 %>%
  rename(N15_19 = Number)
male_age20_24<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-males-age-20-24-in-2022-liver-and-intrahepatic-bile-ducts.csv")
male_age20_24 <-male_age20_24 %>%
  rename(N20_24 = Number)
male_age25_29<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-males-age-25-29-in-2022-liver-and-intrahepatic-bile-ducts.csv")
male_age25_29 <-male_age25_29 %>%
  rename(N25_29 = Number)
male_age30_34<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-males-age-30-34-in-2022-liver-and-intrahepatic-bile-ducts.csv")
male_age30_34 <-male_age30_34 %>%
  rename(N30_34 = Number)
male_age35_39<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-males-age-35-39-in-2022-liver-and-intrahepatic-bile-ducts.csv")
male_age35_39 <-male_age35_39 %>%
  rename(N35_39 = Number)
male_age40_44<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-males-age-40-44-in-2022-liver-and-intrahepatic-bile-ducts.csv")
male_age40_44 <-male_age40_44 %>%
  rename(N40_44 = Number)
male_age45_49<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-males-age-45-49-in-2022-liver-and-intrahepatic-bile-ducts.csv")
male_age45_49 <-male_age45_49 %>%
  rename(N45_49 = Number)
male_age50_54<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-males-age-50-54-in-2022-liver-and-intrahepatic-bile-ducts.csv")
male_age50_54 <-male_age50_54 %>%
  rename(N50_54 = Number)
male_age55_59<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-males-age-55-59-in-2022-liver-and-intrahepatic-bile-ducts.csv")
male_age55_59 <-male_age55_59 %>%
  rename(N55_59 = Number)
male_age60_64<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-males-age-60-64-in-2022-liver-and-intrahepatic-bile-ducts.csv")
male_age60_64 <-male_age60_64 %>%
  rename(N60_64 = Number)
male_age65_69<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-males-age-65-69-in-2022-liver-and-intrahepatic-bile-ducts.csv")
male_age65_69 <-male_age65_69 %>%
  rename(N65_69 = Number)
male_age70_74<-read_xlsx("GLOBOCAN 2022 data/countries/incidence/dataset-inc-males-age-70-74-in-2022-liver-and-intrahepatic-bile-ducts.xlsx")
male_age70_74 <-male_age70_74 %>%
  rename(N70_74 = Number)
male_age75_79<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-males-age-75-79-in-2022-liver-and-intrahepatic-bile-ducts.csv")
male_age75_79 <-male_age75_79 %>%
  rename(N75_79 = Number)
male_age80_84<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-males-age-80-84-in-2022-liver-and-intrahepatic-bile-ducts.csv")
male_age80_84 <-male_age80_84 %>%
  rename(N80_84 = Number)
male_age85_<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-males-age-85-85-in-2022-liver-and-intrahepatic-bile-ducts.csv")
male_age85_ <-male_age85_ %>%
  rename(N85. = Number)
data_frames_list <- list(male_age85_, male_age80_84, male_age75_79, male_age70_74, male_age65_69, male_age60_64,
                         male_age55_59, male_age50_54, male_age5_9, male_age45_49, male_age40_44, male_age35_39, 
                         male_age30_34, male_age25_29, male_age20_24, male_age15_19, male_age10_14, male_age0_4)
data_frames_lists <- lapply(data_frames_list, function(df) {
  n <- ncol(df)                     # get number of row
  df[, -c(1,2,3,(n-5):n)]     # delete rows
})
combined_male <- Reduce(function(x, y) merge(x, y, by = c("Sex",  "Country","Label"), all = TRUE), data_frames_lists)

female_age0_4<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-females-age-0-4-in-2022-liver-and-intrahepatic-bile-ducts.csv")
female_age0_4 <-female_age0_4 %>%
  rename(N0_4 = Number)
female_age5_9<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-females-age-5-9-in-2022-liver-and-intrahepatic-bile-ducts.csv")
female_age5_9 <-female_age5_9 %>%
  rename(N5_9 = Number)
female_age10_14<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-females-age-10-14-in-2022-liver-and-intrahepatic-bile-ducts.csv")
female_age10_14 <-female_age10_14 %>%
  rename(N10_14 = Number)
female_age15_19<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-females-age-15-19-in-2022-liver-and-intrahepatic-bile-ducts.csv")
female_age15_19 <-female_age15_19 %>%
  rename(N15_19 = Number)
female_age20_24<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-females-age-20-24-in-2022-liver-and-intrahepatic-bile-ducts.csv")
female_age20_24 <-female_age20_24 %>%
  rename(N20_24 = Number)
female_age25_29<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-females-age-25-29-in-2022-liver-and-intrahepatic-bile-ducts.csv")
female_age25_29 <-female_age25_29 %>%
  rename(N25_29 = Number)
female_age30_34<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-females-age-30-34-in-2022-liver-and-intrahepatic-bile-ducts.csv")
female_age30_34 <-female_age30_34 %>%
  rename(N30_34 = Number)
female_age35_39<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-females-age-35-39-in-2022-liver-and-intrahepatic-bile-ducts.csv")
female_age35_39 <-female_age35_39 %>%
  rename(N35_39 = Number)
female_age40_44<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-females-age-40-44-in-2022-liver-and-intrahepatic-bile-ducts.csv")
female_age40_44 <-female_age40_44 %>%
  rename(N40_44 = Number)
female_age45_49<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-females-age-45-49-in-2022-liver-and-intrahepatic-bile-ducts.csv")
female_age45_49 <-female_age45_49 %>%
  rename(N45_49 = Number)
female_age50_54<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-females-age-50-54-in-2022-liver-and-intrahepatic-bile-ducts.csv")
female_age50_54 <-female_age50_54 %>%
  rename(N50_54 = Number)
female_age55_59<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-females-age-55-59-in-2022-liver-and-intrahepatic-bile-ducts.csv")
female_age55_59 <-female_age55_59 %>%
  rename(N55_59 = Number)
female_age60_64<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-females-age-60-64-in-2022-liver-and-intrahepatic-bile-ducts.csv")
female_age60_64 <-female_age60_64 %>%
  rename(N60_64 = Number)
female_age65_69<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-females-age-65-69-in-2022-liver-and-intrahepatic-bile-ducts.csv")
female_age65_69 <-female_age65_69 %>%
  rename(N65_69 = Number)
female_age70_74<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-females-age-70-74-in-2022-liver-and-intrahepatic-bile-ducts.csv")
female_age70_74 <-female_age70_74 %>%
  rename(N70_74 = Number)
female_age75_79<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-females-age-75-79-in-2022-liver-and-intrahepatic-bile-ducts.csv")
female_age75_79 <-female_age75_79 %>%
  rename(N75_79 = Number)
female_age80_84<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-females-age-80-84-in-2022-liver-and-intrahepatic-bile-ducts.csv")
female_age80_84 <-female_age80_84 %>%
  rename(N80_84 = Number)
female_age85_<-read.csv("GLOBOCAN 2022 data/countries/incidence/dataset-inc-females-age-85-85-in-2022-liver-and-intrahepatic-bile-ducts.csv")
female_age85_ <-female_age85_ %>%
  rename(N85. = Number)
data_frames_list <- list(female_age85_, female_age80_84, female_age75_79, female_age70_74, female_age65_69, female_age60_64,
                         female_age55_59, female_age50_54, female_age5_9, female_age45_49, female_age40_44, female_age35_39, 
                         female_age30_34, female_age25_29, female_age20_24, female_age15_19, female_age10_14, female_age0_4)
data_frames_lists <- lapply(data_frames_list, function(df) {
  n <- ncol(df)                     # get row numbers
  df[, -c(1,2,3,(n-5):n)]     # delete rows
})
# combine data 
combined_female <- Reduce(function(x, y) merge(x, y, by = c("Sex",  "Country","Label"), all = TRUE), data_frames_lists)
combined_female$Sex=2
combined_male$Sex=1
combined_female<-combined_female[-nrow(combined_female),]
combined_male<-combined_male[-nrow(combined_male),]
combined_all<-merge(combined_female,combined_male,by=c("Sex","Country","Label","N85.","N80_84","N75_79","N70_74",
                                                       "N65_69","N60_64","N55_59","N50_54","N45_49","N40_44","N35_39",
                                                       "N30_34","N25_29","N20_24","N15_19","N10_14","N5_9","N0_4"),all=T)
cols_with_n <- grep("^N", names(combined_all), value = TRUE)
case<-combined_all
case$Label[which(case$Label=="Türkiye")]<-"Turkey"
ratio<-read_xlsx("HCC country ratio.xlsx")
case_all<-merge(case,ratio,by=c("Sex","Label"),all.x = T)

# get case columns
n_columns <- grep("^N", names(case_all), value = TRUE)

# calculate HCC cases
for (col in n_columns) {
  case_all[[col]] <- round(case_all[[col]] * case_all$Ratio)
}
alpha_code<-male_age0_4[,c(1,5)]
alpha_code$Label[which(alpha_code$Label=="Türkiye")]<-"Turkey"
case_all<-merge(alpha_code,case_all)
n<-ncol(case_all)
case_HCC<-case_all[,-n]
case<-case_HCC
case$Total <- rowSums(case[, grep("^N", names(case))], na.rm = TRUE)
#write_xlsx(case,"HCC_incidence_country.xlsx")

#######standardized age structure
age_st<-data.frame(
  AgeGrp = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84","85+"),
  pop = c(0.12, 0.10, 0.09, 0.09, 0.08, 0.08, 0.06, 0.06, 0.06, 0.06, 0.05, 0.04, 0.04, 0.03, 0.02, 0.01, 0.005, 0.005)
)
########population structure
population_structure<-read_xlsx("population structure.xlsx")
population_structure$AgeGrp[population_structure$AgeGrp=='45421']<-'5-9'
population_structure$AgeGrp[population_structure$AgeGrp=='45579']<-'10-14'
population_structure <- population_structure[which(population_structure$LocTypeName=="Country/Area"),]
population_structure<-population_structure[which(population_structure$Time==2022),]
table(population_structure$Location)
population_structure_world <- population_structure %>%
  dplyr::select(Time, AgeGrp, PopTotal,PopMale,PopFemale,Location)

table(population_structure$Location)
# Create age group 85+
population_structure_85plus <- population_structure_world %>%
  filter(AgeGrp %in% c("85-89", "90-94", "95-99", "100+")) %>%
  group_by(Time, Location) %>%
  summarise(
    PopTotal = sum(PopTotal, na.rm = TRUE),
    PopMale = sum(PopMale, na.rm = TRUE),
    PopFemale = sum(PopFemale, na.rm = TRUE),
    .groups = 'drop'  # This will prevent the grouping message and drop the grouping
  ) %>%
  mutate(AgeGrp = "85+")

# remove four age groups
population_structure_world <- population_structure_world %>%
  filter(!AgeGrp %in% c("85-89", "90-94", "95-99", "100+"))

population_structure_final <- population_structure_world %>%
  bind_rows(population_structure_85plus) %>%
  arrange(Time, AgeGrp,Location)
#######rename some country names to match the case data
population_structure_final$Location[population_structure_final$Location=="Bolivia (Plurinational State of)"]<-"Bolivia"
population_structure_final$Location[population_structure_final$Location=="Bosnia and Herzegovina"]<-"Bosnia Herzegovina"
population_structure_final$Location[population_structure_final$Location=="Cabo Verde"]<-"Cape Verde"
population_structure_final$Location[population_structure_final$Location=="Congo"]<-"Congo, Republic of"
population_structure_final$Location[population_structure_final$Location=="Democratic Republic of the Congo"]<-"Congo, Democratic Republic of"
population_structure_final$Location[population_structure_final$Location=="France"]<-"France (metropolitan)"
population_structure_final$Location[population_structure_final$Location=="French Guiana"]<-"French Guyana"
population_structure_final$Location[population_structure_final$Location=="Guadeloupe"]<-"France, Guadeloupe"
population_structure_final$Location[population_structure_final$Location=="Iran (Islamic Republic of)"]<-"Iran, Islamic Republic of"
population_structure_final$Location[population_structure_final$Location=="Republic of Korea"]<-"Korea, Republic of"
population_structure_final$Location[population_structure_final$Location=="Dem. People's Republic of Korea"]<-"Korea, Democratic People Republic of"
population_structure_final$Location[population_structure_final$Location=="Martinique"]<-"France, Martinique"
population_structure_final$Location[population_structure_final$Location=="Republic of Moldova"]<-"Moldova"
population_structure_final$Location[population_structure_final$Location=="Netherlands"]<-"The Netherlands"
population_structure_final$Location[population_structure_final$Location=="Réunion"]<-"France, La Réunion"
population_structure_final$Location[population_structure_final$Location=="United Republic of Tanzania"]<-"Tanzania, United Republic of"
population_structure_final$Location[population_structure_final$Location=="Venezuela (Bolivarian Republic of)"]<-"Venezuela"
population_structure_final$Location[population_structure_final$Location=="Gambia"]<-"The Republic of the Gambia"
population_structure_final$Location[population_structure_final$Location=="TÜRKIYE"]<-"Turkey"
population_structure_final$Location[population_structure_final$Location=="C么te d'Ivoire"]<-"Côte d'Ivoire"
population_structure_final$Location[population_structure_final$Location=="R茅union"]<-"France, La Réunion"
population_structure_final$Location[population_structure_final$Location=="State of Palestine"]<-"Gaza Strip and West Bank"
# Prepare vectors to store results
table(population_structure_final$Location)
library(dplyr)
library(tidyr)


###########calculate ASIR
long_data <- case %>%
  pivot_longer(
    cols = starts_with("N"),   # Selects all columns that start with 'N'
    names_to = "AgeGrp",       # New column for age groups
    names_prefix = "N",        # Remove 'N' prefix in age group names
    values_to = "Cases"        # New column for the case numbers
  ) %>%
  mutate(AgeGrp = replace(AgeGrp, AgeGrp == "85", "85+"))  # Adjusting age group naming for 85+
long_data$TOTAL<-NULL

summarized_data <- long_data %>%
  group_by(Country, Label, AgeGrp) %>%
  summarise(Cases = sum(Cases),    # Summing up the cases
            Sex = 0,               # Setting Sex to 0
            .groups = "drop")      # Drop the grouping
case_male<-case[which(case$Sex==1),]
continents<-case_male$Label
table(population_structure$Location)
asir<- vector("numeric", length(continents))
table(population_structure_final$Location)
# Loop through each country

for (i in seq_along(continents)) {
  continent <- continents[i]
  #continent<-"Turkey"
  # Filter data for the current continent
  pop_data <- population_structure_final[population_structure_final$Location == continent,]
  case_data <- long_data[long_data$Label == continent,]
  male_cases_data <- case_data[case_data$Sex == 1, ]
  female_cases_data <- case_data[case_data$Sex == 2, ]
  pop_data <- pop_data[order(pop_data$AgeGrp), ]
  male_cases_data <- male_cases_data[order(male_cases_data$AgeGrp), ]
  female_cases_data <- female_cases_data[order(female_cases_data$AgeGrp), ]
  age_st<-age_st[order(age_st$AgeGrp),]
  ASR<- ((male_cases_data$Cases+female_cases_data$Cases)/ pop_data$PopTotal) * 100000/1000
  Weighted_ASR <- ASR * age_st$pop
  
  # Sum weighted ASRs to get ASIR
  asir[i] <- sum(Weighted_ASR, na.rm = TRUE)
}
ASIR<-data.frame(Location = rep(NA, 185))
ASIR$Location<-continents
ASIR$ASIR<-asir
alpha<-case[,c(1,2)]
alpha$Location<-alpha$Label
alpha$Label<-NULL
alpha<-unique(alpha)
ASIR<-merge(ASIR,alpha,by="Location",all.x = T)
##write_xlsx(ASIR,"ASIR_country.xlsx")

###########calculate the case for both sex and combine case and ASIR
#####calculate the case for both sex by unregion
case_incidence<-case%>%
  group_by(Label,Alpha.3.code)%>%
  summarize(
    total=sum(Total)
  )
case_incidence$Location<-case_incidence$Label
case_incidence$Label<-NULL
######keep the ASR of both sex
######merge these two data frame 
incidence<-merge(ASIR,case_incidence,by=c("Location","Alpha.3.code"))



#############PART2: mortality
case<-read_xlsx("PLC_mortality_country.xlsx") ##### merged GLOBOCAN 2022 data 
ratio<-read_xlsx("HCC country ratio.xlsx")
case_all<-merge(case,ratio,by=c("Sex","Label"))

# get all case lines
n_columns <- grep("^N", names(case_all), value = TRUE)

# calculate HCC death cases 
for (col in n_columns) {
  case_all[[col]] <- round(case_all[[col]] * case_all$Ratio)
}



n<-ncol(case_all)
case_HCC<-case_all[,-n]
case<-case_HCC
case<-merge(case,alpha_code,by="Label")

#write_xlsx(case,"HCC_mortality_country.xlsx")
#case<-read_xlsx("HCC_mortality_country.xlsx")


#calculate ASMR
long_data <- case %>%
  pivot_longer(
    cols = starts_with("N"),   # Selects all columns that start with 'N'
    names_to = "AgeGrp",       # New column for age groups
    names_prefix = "N",        # Remove 'N' prefix in age group names
    values_to = "Cases"        # New column for the case numbers
  ) %>%
  mutate(AgeGrp = replace(AgeGrp, AgeGrp == "85", "85+"))  # Adjusting age group naming for 85+
long_data$TOTAL<-NULL



summarized_data <- long_data %>%
  group_by(Country, Label, AgeGrp) %>%
  summarise(Cases = sum(Cases),    # Summing up the cases
            Sex = 0,               # Setting Sex to 0
            .groups = "drop")      # Drop the grouping
case_male<-case[which(case$Sex==1),]
continents<-case_male$Label
table(population_structure$Location)
asmr<- vector("numeric", length(continents))
table(population_structure_final$Location)
# Loop through each continent

for (i in seq_along(continents)) {
  continent <- continents[i]
  #continent<-"Turkey"
  # Filter data for the current continent
  pop_data <- population_structure_final[population_structure_final$Location == continent,]
  case_data <- long_data[long_data$Label == continent,]
  male_cases_data <- case_data[case_data$Sex == 1, ]
  female_cases_data <- case_data[case_data$Sex == 2, ]
  pop_data <- pop_data[order(pop_data$AgeGrp), ]
  male_cases_data <- male_cases_data[order(male_cases_data$AgeGrp), ]
  female_cases_data <- female_cases_data[order(female_cases_data$AgeGrp), ]
  age_st<-age_st[order(age_st$AgeGrp),]
  ASR<- ((male_cases_data$Cases+female_cases_data$Cases)/ pop_data$PopTotal) * 100000/1000
  Weighted_ASR <- ASR * age_st$pop
  
  # Sum weighted ASRs to get ASIR
  asmr[i] <- sum(Weighted_ASR, na.rm = TRUE)
}
ASMR<-data.frame(Location = rep(NA, 185))
ASMR$Location<-continents
ASMR$ASMR<-asmr
alpha_code$Location<-alpha_code$Label
alpha_code$Label<-NULL

ASMR<-merge(ASMR,alpha_code,by="Location",all.x = T)
#write_xlsx(ASMR,"ASMR_country.xlsx")

###########calculate the deaths for both sex and combine death and ASMR
#####calculate the case for both sex by unregion
case$Total <- rowSums(case[, grep("^N", names(case))], na.rm = TRUE)
case_mortality<-case%>%
  group_by(Label,Alpha.3.code)%>%
  summarize(
    total=sum(Total)
  )
case_mortality$Location<-case_mortality$Label
case_mortality$Label<-NULL
######keep the ASR of both sex
######merge these two data frame 
mortality<-merge(ASMR,case_mortality,by=c("Location","Alpha.3.code"))


########combine the incidence and mortality
incidence$case<-incidence$total
incidence$total<-NULL
mortality$death<-mortality$total
mortality$total<-NULL
inc_mort<-merge(incidence,mortality,by=c("Location","Alpha.3.code"))
####write_xlsx(inc_mort,"HCC country.xlsx")