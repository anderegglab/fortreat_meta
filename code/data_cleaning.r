##---------------------------------------------------------------
## Data Cleaning
##
##---------------------------------------------------------------

##---------------------------------------------------------------
## 0. Head
##---------------------------------------------------------------

## libraries

source("code/functions.r")

## read data
insect_data <- read.csv("data/raw_data/insect_analysis.csv")
drought_data <- read.csv("data/raw_data/drought_analysis.csv")
fire_data <- read.csv("data/raw_data/fire_analysis.csv")

## clean up idiosyncracies between files
insect_data <- insect_data[,!(colnames(insect_data) %in% c("PercentComposition", "thin_severity"))]

colnames(fire_data)[colnames(fire_data) == "siteName"] <- "SiteName"

drought_data <- drought_data[,colnames(drought_data) != "thin_severity"]
drought_data$SiteName <- NA
drought_data <- drought_data[,colnames(insect_data)]

insect_data$disturbance_type <- "insect"
drought_data$disturbance_type <- "drought"
fire_data$disturbance_type <- "fire"

## combine datasets
full_data <- rbind(insect_data, fire_data, drought_data)

## rename columns with egregiously long names
colnames(full_data)[colnames(full_data) == "X.BA.removed"] <- "ba_removed"
colnames(full_data)[colnames(full_data) == "MortalityAttributionUncertain..1.yes."] <- "mortality_attribution_uncertain"
colnames(full_data)[colnames(full_data) == "GroupingFlags..same.number...should.be.grouped.later..recycle.numbers.between.studies..If.nothing.should.be.grouped..write.NA."] <- "grouping_flags"
colnames(full_data)[colnames(full_data) == "ImputedSD...1.yes."] <- "impute_sd"
colnames(full_data)[colnames(full_data) == "ImputedBAremoved...1.yes."] <- "impute_ba"
colnames(full_data)[colnames(full_data) == "Carbon.vs.Mortality..1.Carbon..2.Mortality."] <- "carbon_vs_mortality"


##---------------------------------------------------------------
## 1. Flag studies that we still need to sort
##---------------------------------------------------------------

## need to remove some of these problem studies

## TODO
## Westling and Kerns 2021 (studyID = 284)
## reports mortality as dead trees per ha, no way to convert to % mortality
full_data <- full_data[full_data$studyID != 284,]

## TODO
## Muzika et al 2000 (studyID = 201)
## can't calculate SE in correct units, because transformation requires dividing by initial mortality
## for now, impute SE
full_data[full_data$studyID == 201, c("se_treatment", "se_control")] <- NA
full_data[full_data$studyID == 201, "impute_sd"] <- 1

## TODO
## Young et al. 2020 (studyID = 250)
## Need to re-extract using Figure 2a-2b. For now remove
full_data <- full_data[full_data$studyID != 250,]


## TODO
## Zhang et al. (studyID = 193)
## mean_control value = 1.3? Remove for now
full_data <- full_data[!(full_data$studyID == 193 & full_data$mean_control > 1), ]


##---------------------------------------------------------------
## 2. mortality vs. survivorship
##---------------------------------------------------------------

full_data[full_data$carbon_vs_mortality == 2, "mean_treatment"]
full_data[full_data$carbon_vs_mortality == 2, "mean_control"]
full_data[full_data$carbon_vs_mortality == 2, "se_control"]

hist(full_data[full_data$carbon_vs_mortality == 2, "mean_treatment"])


## remove NAs
full_data <- full_data[!is.na(full_data$mean_treatment),]

## 100% mortality
nrow(full_data[full_data$carbon_vs_mortality == 2 & full_data$mean_control == 1,])
nrow(full_data[full_data$carbon_vs_mortality == 2 & full_data$mean_treatment == 1,])

## 0% mortality
nrow(full_data[full_data$carbon_vs_mortality == 2 & full_data$mean_control == 0,])
nrow(full_data[full_data$carbon_vs_mortality == 2 & full_data$mean_treatment == 0,])

## looks like flipping to survivorship is the move
full_data[full_data$carbon_vs_mortality == 2, c("mean_treatment", "mean_control")] <-
  1 - full_data[full_data$carbon_vs_mortality == 2, c("mean_treatment", "mean_control")]

##---------------------------------------------------------------
## 2. Group data
##---------------------------------------------------------------

## run grouping script
full_data_grouped <- group_data(full_data)

##---------------------------------------------------------------
## 3. Calculate effect sizes and standard errors
##---------------------------------------------------------------
## for now we will add a small number when survivorship is 0:
full_data_grouped[full_data_grouped$carbon_vs_mortality == 2 & full_data_grouped$mean_treatment == 0, c("mean_treatment", "sd_treatment", "se_treatment")] <-
  full_data_grouped[full_data_grouped$carbon_vs_mortality == 2 & full_data_grouped$mean_treatment == 0, c("mean_treatment", "sd_treatment", "se_treatment")] + 0.01

full_data_grouped[full_data_grouped$carbon_vs_mortality == 2 & full_data_grouped$mean_control == 0, c("mean_control", "sd_control", "se_control")] <-
  full_data_grouped[full_data_grouped$carbon_vs_mortality == 2 & full_data_grouped$mean_control == 0, c("mean_control", "sd_control", "se_control")] + 0.01

## calculate log response ratio
full_data_grouped$lrr <- lrr(full_data_grouped$mean_treatment,
                             full_data_grouped$mean_control)

## calculate se
full_data_grouped$lrr_se <- lrr_se(full_data_grouped$mean_treatment,
                               full_data_grouped$mean_control,
                               full_data_grouped$se_treatment,
                               full_data_grouped$se_control)
## NAs remain NAs (missing values in og data)

full_data_grouped[3:8,]


##---------------------------------------------------------------
## 4. save cleaned data file
##---------------------------------------------------------------

write.csv(full_data_grouped, "data/processed_data/data_cleaned.csv")
