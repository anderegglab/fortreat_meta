##---------------------------------------------------------------
## Figure 1
##
##---------------------------------------------------------------

## load packages
library("metafor")
library("ggplot2")

## read in data
data <- read.csv("data/processed_data/data_cleaned.csv")

## models
freq_fit_fire <- readRDS("data/model_objects/freq_fit_fire.rds")
freq_fit_insect <- readRDS("data/model_objects/freq_fit_insect.rds")
freq_fit_drought <- readRDS("data/model_objects/freq_fit_drought.rds")

##

figure_1 <- ##
