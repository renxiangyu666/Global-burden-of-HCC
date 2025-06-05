############################HCC UNregion
############################Purpose: Calculate HCC cases/deaths and Age standardized ratio(ASR)
############################Input: country-level HCC cases/deaths data, population structure data, and standardized age structure
############################Output: Formatted HCC cases/death and ASR result for subregions, regions and global levels

library(dplyr)
library(readxl)
library(writexl)
#############incidence
#########HCC case calculation
setwd("C:/Users/24399/Desktop/Core code for HCC and PAF estimation/HCC burden/")
case_country<-read_xlsx("HCC_incidence_country.xlsx")
n_columns <- grep("^N", names(case_country), value = TRUE)
case_UNregion <- case_country %>%
  group_by(Sex, Unregion,Region,Continent) %>%
  summarise(across(all_of(n_columns), sum, na.rm = TRUE))
#calculate the case for every continent 
case_continent <- case_country %>%
  group_by(Sex,Continent) %>%
  summarise(across(all_of(n_columns), sum, na.rm = TRUE))
 case_continent$Unregion<-case_continent$Continent
 case_continent$Region<-case_continent$Continent
#calculate the case for globe
 case_Globe <- case_country %>%
   group_by(Sex) %>%
   summarise(across(all_of(n_columns), sum, na.rm = TRUE))
 case_Globe$Continent<-'Globe'
 case_Globe$Unregion<-'Globe'
 case_Globe$Region<-'Globe'
#combine these three datasets
case<-rbind(case_continent,case_UNregion,case_Globe)
#########calculate the total case
case$Label<-case$Unregion
case$Total <- rowSums(case[, grep("^N", names(case))], na.rm = TRUE)
#write_xlsx(case,"HCC_incidence_UNregion.xlsx")



#########age standard data for ASR calculation
age_standard<-data.frame(
  age_group = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84","85+"),
  proportion = c(0.12, 0.10, 0.09, 0.09, 0.08, 0.08, 0.06, 0.06, 0.06, 0.06, 0.05, 0.04, 0.04, 0.03, 0.02, 0.01, 0.005, 0.005)
)


age_stand<-age_standard$proportion

age_stand
age_st<-data.frame(Pop = rep(NA, 18))
age_st$Pop<-age_stand


#########population data for ASR calculation
library(dplyr)
population_structure<-read_xlsx("population structure.xlsx")
population_structure$Location[population_structure$Location=="World"]<-"Globe"
population_structure$Location[population_structure$Location=="Eastern Europe"]<-"Central and Eastern Europe"
population_structure$AgeGrp[population_structure$AgeGrp=='45421']<-'5-9'
population_structure$AgeGrp[population_structure$AgeGrp=='45579']<-'10-14'
population_structure<-population_structure[which(population_structure$Time==2022),]
population_structure<-population_structure[-which(population_structure$LocID==927),]

population_structure_world <- population_structure %>%
  dplyr::select(Time, AgeGrp, PopTotal,PopMale,PopFemale,Location)

# create 85+ agegroup，by adding 85-89, 90-94, 95-99, 100+
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

# remove agegroup 
population_structure_world <- population_structure_world %>%
  filter(!AgeGrp %in% c("85-89", "90-94", "95-99", "100+"))

population_structure_final <- population_structure_world %>%
  bind_rows(population_structure_85plus) %>%
  arrange(Time, AgeGrp,Location)

# combine continents list
labels_to_combine <- c("Micronesia", "Polynesia","Melanesia")
pop<-population_structure_final[which(population_structure_final$Location=="Micronesia"|population_structure_final$Location=="Polynesia"|population_structure_final$Location=="Melanesia"),]
pop_combined <- pop %>%
  filter(Location %in% labels_to_combine) %>%
  group_by(Time, AgeGrp) %>%
  summarise(
    PopTotal = sum(PopTotal, na.rm = TRUE),
    PopMale = sum(PopMale, na.rm = TRUE),
    PopFemale = sum(PopFemale, na.rm = TRUE),
    Location = "Melanesia/Micronesia/Polynesia",
    .groups = 'drop'
  )
population_structure_final<-rbind(pop_combined,population_structure_final)

labels_to_combine <- c("Central Asia", "Southern Asia")
pop<-population_structure_final[which(population_structure_final$Location=="Central Asia"|population_structure_final$Location=="Southern Asia"),]
pop_combined <- pop %>%
  filter(Location %in% labels_to_combine) %>%
  group_by(Time, AgeGrp) %>%
  summarise(
    PopTotal = sum(PopTotal, na.rm = TRUE),
    PopMale = sum(PopMale, na.rm = TRUE),
    PopFemale = sum(PopFemale, na.rm = TRUE),
    Location = "South Central Asia",
    .groups = 'drop'
  )
population_structure_final<-rbind(pop_combined,population_structure_final)

labels_to_combine <- c("Northern America", "South America","Central America","Caribbean")
pop<-population_structure_final[which(population_structure_final$Location=="Northern America"|population_structure_final$Location=="South America"|
                                        population_structure_final$Location=="Central America"|population_structure_final$Location=="Caribbean"),]
pop_combined <- pop %>%
  filter(Location %in% labels_to_combine) %>%
  group_by(Time, AgeGrp) %>%
  summarise(
    PopTotal = sum(PopTotal, na.rm = TRUE),
    PopMale = sum(PopMale, na.rm = TRUE),
    PopFemale = sum(PopFemale, na.rm = TRUE),
    Location = "America",
    .groups = 'drop'
  )
population_structure_final<-rbind(pop_combined,population_structure_final)




#####calculate ASIRs
# Prepare vectors to store results

library(dplyr)
library(tidyr)
orders <- c("Globe","Melanesia/Micronesia/Polynesia","Australia/New Zealand","Oceania",
            "Western Europe","Southern Europe","Northern Europe","Central and Eastern Europe","Europe",
            "South America","Caribbean","Central America","Northern America","America",
            "Western Asia","South Central Asia","South-Eastern Asia","Eastern Asia","Asia",
            "Western Africa","Southern Africa","Northern Africa","Middle Africa", "Eastern Africa","Africa")
# calculate ASIR
long_data <- case %>%
  pivot_longer(
    cols = starts_with("N"),   # Selects all columns that start with 'N'
    names_to = "AgeGrp",       # New column for age groups
    names_prefix = "N",        # Remove 'N' prefix in age group names
    values_to = "Cases"        # New column for the case numbers
  ) %>%
  mutate(AgeGrp = replace(AgeGrp, AgeGrp == "85", "85+"))  # Adjusting age group naming for 85+
long_data$TOTAL<-NULL

asir_male <- vector("numeric", length(orders))
asir_female <- vector("numeric", length(orders))
asir_both<-vector("numeric", length(orders))
# Loop through each continent

for (i in seq_along(orders)) {
  continent <- orders[i]
  
  #continent<-"Southern Africa"
  # continent<-"Globe"
  # Filter data for the current continent
  pop_data <- population_structure_final[population_structure_final$Location == continent,]
  case_data <- long_data[long_data$Label == continent,]
  male_cases_data <- case_data[case_data$Sex == 1, ]
  female_cases_data <- case_data[case_data$Sex == 2, ]
  pop_data <- pop_data[order(pop_data$AgeGrp), ]
  male_cases_data <- male_cases_data[order(male_cases_data$AgeGrp), ]
  female_cases_data <- female_cases_data[order(female_cases_data$AgeGrp), ]
  ASR_Male <- (male_cases_data$Cases/ pop_data$PopMale) * 100000/1000
  ASR_Female <- (female_cases_data$Cases/ pop_data$PopFemale) * 100000/1000
  ASR_both<-((female_cases_data$Cases+male_cases_data$Cases)/ pop_data$PopTotal) * 100000/1000
  Weighted_ASR_Male <- ASR_Male * age_stand
  Weighted_ASR_Female <- ASR_Female * age_stand
  Weighted_ASR_both<-ASR_both*age_stand
  # Sum weighted ASRs to get ASIR
  asir_male[i] <- sum(Weighted_ASR_Male, na.rm = TRUE)
  asir_female[i] <- sum(Weighted_ASR_Female, na.rm = TRUE)
  asir_both[i] <-sum(Weighted_ASR_both,na.rm=T)
}

asir_female
asir_male
ASIR<-data.frame(Location = rep(NA, 25))
ASIR$Location<-orders
ASIR$male<-asir_male
ASIR$female<-asir_female
ASIR$Both<-asir_both
#write_xlsx(ASIR,"UNregion_ASIR.xlsx") ####this file will be used to create Figure 1A


###########calculate the case for both sex and combine case and ASIR
#####calculate the case for both sex by unregion
case_incidence<-case%>%
  group_by(Label,Continent)%>%
  summarize(
    total=sum(Total)
  )
case_incidence$Location<-case_incidence$Label
case_incidence$Label<-NULL
######keep the ASR of both sex
ASIR<-ASIR[,c(1,4)]
######merge these two data frame 
incidence<-merge(ASIR,case_incidence,by="Location")



################Mortality 
########calculate the HCC deaths and ASMR by Unregion 
#####calculate the HCC deaths for both sex 
###########HCC death number
case_country<-read_xlsx("HCC_mortality_country.xlsx")
n_columns <- grep("^N", names(case_country), value = TRUE)
case_UNregion <- case_country %>%
  group_by(Sex, Unregion,Region,Continent) %>%
  summarise(across(all_of(n_columns), sum, na.rm = TRUE))
case_continent <- case_country %>%
  group_by(Sex,Continent) %>%
  summarise(across(all_of(n_columns), sum, na.rm = TRUE))
case_continent$Unregion<-case_continent$Continent
case_continent$Region<-case_continent$Continent
case_Globe <- case_country %>%
  group_by(Sex) %>%
  summarise(across(all_of(n_columns), sum, na.rm = TRUE))
case_Globe$Continent<-'Globe'
case_Globe$Unregion<-'Globe'
case_Globe$Region<-'Globe'
case<-rbind(case_continent,case_UNregion,case_Globe)
case$Label<-case$Unregion
case$Total <- rowSums(case[, grep("^N", names(case))], na.rm = TRUE)
###write_xlsx(case,"HCC_mortality_UNregion.xlsx")




# Prepare vectors to store results

library(dplyr)
library(tidyr)
orders <- c("Globe","Melanesia/Micronesia/Polynesia","Australia/New Zealand","Oceania",
            "Western Europe","Southern Europe","Northern Europe","Central and Eastern Europe","Europe",
            "South America","Caribbean","Central America","Northern America","America",
            "Western Asia","South Central Asia","South-Eastern Asia","Eastern Asia","Asia",
            "Western Africa","Southern Africa","Northern Africa","Middle Africa", "Eastern Africa","Africa")

long_data <- case %>%
  pivot_longer(
    cols = starts_with("N"),   # Selects all columns that start with 'N'
    names_to = "AgeGrp",       # New column for age groups
    names_prefix = "N",        # Remove 'N' prefix in age group names
    values_to = "Cases"        # New column for the case numbers
  ) %>%
  mutate(AgeGrp = replace(AgeGrp, AgeGrp == "85", "85+"))  # Adjusting age group naming for 85+
long_data$TOTAL<-NULL

ASMR_male <- vector("numeric", length(orders))
ASMR_female <- vector("numeric", length(orders))
ASMR_both<-vector("numeric", length(orders))
# Loop through each continent

for (i in seq_along(orders)) {
  continent <- orders[i]
  
  #continent<-"Southern Africa"
  # continent<-"Globe"
  # Filter data for the current continent
  pop_data <- population_structure_final[population_structure_final$Location == continent,]
  case_data <- long_data[long_data$Label == continent,]
  male_cases_data <- case_data[case_data$Sex == 1, ]
  female_cases_data <- case_data[case_data$Sex == 2, ]
  pop_data <- pop_data[order(pop_data$AgeGrp), ]
  male_cases_data <- male_cases_data[order(male_cases_data$AgeGrp), ]
  female_cases_data <- female_cases_data[order(female_cases_data$AgeGrp), ]
  ASR_Male <- (male_cases_data$Cases/ pop_data$PopMale) * 100000/1000
  ASR_Female <- (female_cases_data$Cases/ pop_data$PopFemale) * 100000/1000
  ASR_both<-((female_cases_data$Cases+male_cases_data$Cases)/ pop_data$PopTotal) * 100000/1000
  Weighted_ASR_Male <- ASR_Male * age_stand
  Weighted_ASR_Female <- ASR_Female * age_stand
  Weighted_ASR_both<-ASR_both*age_stand
  # Sum weighted ASRs to get ASMR
  ASMR_male[i] <- sum(Weighted_ASR_Male, na.rm = TRUE)
  ASMR_female[i] <- sum(Weighted_ASR_Female, na.rm = TRUE)
  ASMR_both[i] <-sum(Weighted_ASR_both,na.rm=T)
}

ASMR_female
ASMR_male
ASMR<-data.frame(Location = rep(NA, 25))
ASMR$Location<-orders
ASMR$male<-ASMR_male
ASMR$female<-ASMR_female
ASMR$Both<-ASMR_both
#write_xlsx(ASMR,"UNregion_ASMR.xlsx") ###this file would be used to create Figure1A



###########calculate the deaths for both sex and combine death and ASMR
#####calculate the case for both sex by unregion
case_mortality<-case%>%
  group_by(Label,Continent)%>%
  summarize(
    total=sum(Total)
  )
case_mortality$Location<-case_mortality$Label
case_mortality$Label<-NULL
######keep the ASR of both sex
ASMR<-ASMR[,c(1,4)]
######merge these two data frame 
mortality<-merge(ASMR,case_mortality,by="Location")


########combine the incidence and mortality
incidence$case<-incidence$total
incidence$ASIR<-incidence$Both
incidence$Both<-NULL
incidence$total<-NULL
mortality$death<-mortality$total
mortality$ASMR<-mortality$Both
mortality$Both<-NULL
mortality$total<-NULL
inc_mort<-merge(incidence,mortality,by=c("Location","Continent"))
####write_xlsx(inc_mort,"HCC Region.xlsx")