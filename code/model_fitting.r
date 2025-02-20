
##---------------------------------------------------------------
## Model Fitting for Meta Analysis
##
load("/Users/evangeline/Desktop/Github/fortreat_meta/data/model_objects/model_fitting.RData")
save.image("~/Desktop/Github/fortreat_meta/data/model_objects/model_fitting.RData")
##---------------------------------------------------------------

##---------------------------------------------------------------
## 0. Head
##---------------------------------------------------------------

## libraries
library("metafor")
library("mice")
library("brms")
library("ggplot2")

setwd('/Users/evangeline/Desktop/Github/fortreat_meta/')
source("code/functions.r")

## read data
data <- read.csv("data/processed_data/data_cleaned.csv")

##---------------------------------------------------------------
## 1. Impute missing data
##---------------------------------------------------------------

## There are lots of missing standard errors, and
## some missing BA_removed. Because imputing BA removed requires
## info from other disturbance types, and some way of combining
## treatment names, I will just demo imputation with the standard
## errors here.

## for imputation we will use the 'mice' package

## MICE stands for multiple imputation with chained equations and
## is performed using the following workflow:

## 1. create N copies of the dataset
## 2. predict missing values in each copy ('mice' automatically
##    chooses a method depending on the variable type, for example
##    predictive mean matching (PMM) with continuous data)
## 3. fit N models on each of the N dataset copies
## 4. pool results using a method-dependent pooling algorithm (or
##    if using Bayes, simply append all the posteriors)

## see figure 1 in https://vbn.aau.dk/ws/portalfiles/portal/257318283/ejbrm_volume15_issue1_article450.pdf
## for a visual summary of this workflow

## here is another nice summary of how the predictive part (step 2) works:
## https://stats.stackexchange.com/questions/421545/multiple-imputation-by-chained-equations-mice-explained

## To run MICE, we need to choose the form of the imputation model,
## including which variables we will use as predictors, and any
## interaction or non-linear terms. However, if we choose to use
## random forest imputation, we do not have to specify interaction
## terms or nonliinear relationships manually.

## For this demo, let's just try it with a single predictor:
## 1) treatment class (thinning vs. rx fire vs. both)

## to do this we first need to categorize each treatment, as they dont have
## standardized names

cat_trt <- function(trt) {
  if (((grepl("burn", tolower(trt)) | grepl("fire", tolower(trt))) & !grepl("unburn", tolower(trt))) &
      (grepl("thin", tolower(trt)) | grepl("density", tolower(trt)) | grepl("umz", tolower(trt)))) return("both")
  else if((grepl("burn", tolower(trt)) | grepl("fire", tolower(trt))) & !grepl("unburn", tolower(trt))) return("rx_fire")
  else if(grepl("thin", tolower(trt)) | grepl("density", tolower(trt)) | grepl("umz", tolower(trt)))
    return("thinning")
  else return(NA)
}

## also lets just look at mortality for now, ignoring biomass/carbon data
## so again we need to clean up the responseVariable column
## TODO deal with all the various non-mortality variables

data$burn <- "no"
data$thin <- "no"
for (i in 1:nrow(data)) {
  print(i)
  data[i,"trt_class"] <- cat_trt(data[i,"treatment"])
  if (!is.na(data[i, "trt_class"])) {
    if (data[i,"trt_class"] == "thinning") {
      data[i,"thin"] <- "yes"
    } else if (data[i,"trt_class"] == "rx_fire") {
      data[i,"burn"] <- "yes"
    } else if (data[i,"trt_class"] == "both") {
      data[i,"thin"] <- "yes"
      data[i,"burn"] <- "yes"
    }
  }
}
data$trt_class ## looks correct

## remove na trt_class for now
## TODO go back to studies and figure out the uncertain trt_classes
## specifically, studies 72 and 419 have confusing treatment descriptions
data <- data[!is.na(data$trt_class),]

data$trt_class <- factor(data$trt_class, levels = c("rx_fire", "thinning", "both"))

## TODO move to functions file, keep here now for readability
impute_data <- function(data, vars = c("lrr", "lrr_se", "trt_class")) {

  mice_data <- data[, vars]

  ## make predictor matrix
  predictor_matrix <- make.predictorMatrix(mice_data)
  predictor_matrix ## looks good

  impute_method <- make.method(mice_data)
  impute_method ## no method specified for complete variables

  ## impute data
  imputed_data <- mice(mice_data,
                       m = 20,
                       predictorMatrix = predictor_matrix,
                       method = impute_method,
                       seed = 1)

  return(imputed_data)

}

# mortality/survivorship
fire_imputed <- impute_data(data[data$disturbance_type == "fire" & data$carbon_vs_mortality == 2,])
insect_imputed <- impute_data(data[data$disturbance_type == "insect" & data$carbon_vs_mortality == 2,])
drought_imputed <- impute_data(data[data$disturbance_type == "drought" & data$carbon_vs_mortality == 2,])

# carbon
fire_imputed_C <- impute_data(data[data$disturbance_type == "fire" & data$carbon_vs_mortality == 1,])
insect_imputed_C <- impute_data(data[data$disturbance_type == "insect" & data$carbon_vs_mortality == 1,])
drought_imputed_C <- impute_data(data[data$disturbance_type == "drought" & data$carbon_vs_mortality == 1,])

str(fire_imputed) ## looks good
complete(fire_imputed_C,2)

saveRDS(fire_imputed, "data/model_objects/imputed_fire.rds")
saveRDS(insect_imputed, "data/model_objects/imputed_insect.rds")
saveRDS(drought_imputed, "data/model_objects/imputed_drought.rds")
saveRDS(fire_imputed_C, "data/model_objects/imputed_fire_C.rds")
saveRDS(insect_imputed_C, "data/model_objects/imputed_insect_C.rds")
saveRDS(drought_imputed_C, "data/model_objects/imputed_drought_C.rds")

##---------------------------------------------------------------
## 4. Run models
##---------------------------------------------------------------

## we are going to try two different methods, frequentist and Bayes, just for fun

## frequentist first, using 'metafor'
# mortality/survivorship
freq_fit_fire <- with(fire_imputed,
                      rma(yi = lrr,
                          sei = lrr_se,
                          mods = ~ trt_class))
saveRDS(freq_fit_fire, "data/model_objects/freq_fit_fire.rds")

freq_fit_insect <- with(insect_imputed,
                      rma(yi = lrr,
                          sei = lrr_se,
                          mods = ~ trt_class))
saveRDS(freq_fit_insect, "data/model_objects/freq_fit_insect.rds")

freq_fit_drought <- with(drought_imputed,
                         rma(yi = lrr,
                             sei = lrr_se,
                             mods = ~ trt_class))
saveRDS(freq_fit_drought, "data/model_objects/freq_fit_drought.rds")

# carbon
freq_fit_fire_C <- with(fire_imputed_C,
                      rma(yi = lrr,
                          sei = lrr_se,
                          mods = ~ trt_class))
saveRDS(freq_fit_fire_C, "data/model_objects/freq_fit_fire_C.rds")

freq_fit_insect_C <- with(insect_imputed_C,
                        rma(yi = lrr,
                            sei = lrr_se,
                            mods = ~ trt_class))
saveRDS(freq_fit_insect_C, "data/model_objects/freq_fit_insect_C.rds")

freq_fit_drought_C <- with(drought_imputed_C,
                         rma(yi = lrr,
                             sei = lrr_se,
                             mods = ~ trt_class))
saveRDS(freq_fit_drought_C, "data/model_objects/freq_fit_drought_C.rds")

## now bayes
## Commenting this out for now -- lets keep it simple.
## if you do run it, it takes a minute (or 5)

## bayes_fit_fire <- brm_multiple(lrr | se(lrr_se) ~ trt_class,
##                                data = fire_imputed,
##                                chains = 4, cores = 4,
##                                silent = 2, refresh = 0,
##                                open_progress = FALSE)
## saveRDS(bayes_fit_fire, "data/model_objects/bayes_fit_fire.rds")

## bayes_fit_insect <- brm_multiple(lrr | se(lrr_se) ~ trt_class,
##                                data = insect_imputed,
##                                chains = 4, cores = 4,
##                                silent = 2, refresh = 0,
##                                open_progress = FALSE)
## saveRDS(bayes_fit_insect, "data/model_objects/bayes_fit_insect.rds")


## bayes_fit_drought <- brm_multiple(lrr | se(lrr_se) ~ trt_class,
##                                data = drought_imputed,
##                                chains = 4, cores = 4,
##                                silent = 2, refresh = 0,
##                                open_progress = FALSE)
## saveRDS(bayes_fit_drought, "data/model_objects/bayes_fit_drought.rds")

##---------------------------------------------------------------
## 5. Results
##---------------------------------------------------------------

## first frequentist
pool_fire <- summary(pool(freq_fit_fire))
pool_fire[-1] <- round(pool_fire[-1], digits = 3)
pool_fire

pool_insect <- summary(pool(freq_fit_insect))
pool_insect[-1] <- round(pool_insect[-1], digits = 3)
pool_insect

pool_drought <- summary(pool(freq_fit_drought))
pool_drought[-1] <- round(pool_drought[-1], digits = 3)
pool_drought


## ## then bayes
## summary(bayes_fit_fire)
## plot(bayes_fit_fire)

## summary(bayes_fit_insect)
## plot(bayes_fit_insect)

## summary(bayes_fit_drought)
## plot(bayes_fit_drought)

