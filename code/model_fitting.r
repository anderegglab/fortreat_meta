
##---------------------------------------------------------------
## Model Fitting for Meta Analysis
## Demo with Insect Data
##
## Jacob Levine; 1/20/25
##---------------------------------------------------------------

##---------------------------------------------------------------
## 0. Head
##---------------------------------------------------------------

## libraries
library("metafor")
library("mice")
library("brms")
library("ggplot2")

source("functions.r")

## read data
insect_data <- read.csv("data/insect_data.csv")

## rename columns with egregiously long names
colnames(insect_data)[33] <- "grouping_flags"
colnames(insect_data)[28] <- "carbon_or_mortality"

##---------------------------------------------------------------
## 1. Group data
##---------------------------------------------------------------

## TODO Move functions to separate file, leave here now for readability
## run grouping script
insect_data_grouped <- group_data(insect_data)

##---------------------------------------------------------------
## 2. Calculate effect sizes and standard errors
##---------------------------------------------------------------
## calculate log response ratio
insect_data_grouped$lrr <- lrr(insect_data_grouped$mean_treatment,
                               insect_data_grouped$mean_control)
## produces NaNs and -Inf when treatment mean is 0. Need to decide what
## to do here, but for now I will toss them.

## TODO try flipping to survivorship, and do sensitivity analysis when
## adding small numbers. If models are sensitive, use mean dif for
## mortality

## There is also a weird NA for treatment_mean, will toss that row too.

insect_data_grouped <- insect_data_grouped[!is.nan(insect_data_grouped$lrr) &
                                           !is.na(insect_data_grouped$lrr) &
                                           insect_data_grouped$lrr != -Inf,]

## calculate se
insect_data_grouped$lrr_se <- lrr_se(insect_data_grouped$mean_treatment,
                               insect_data_grouped$mean_control,
                               insect_data_grouped$se_treatment,
                               insect_data_grouped$se_control)
## NAs remain NAs (missing values in og data)

##---------------------------------------------------------------
## 3. Impute missing data
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

## Another wrinkle: "every relationship included in the analysis model
## needs to be included in the imputation model regardless of whether
## they contain missing values or not." This means that we cannot
## actually perform the imputation until we know what predictors we would
## like to include in the analysis model. It also creates an issue
## with imputing BA removed. If we want to include additional data from
## other disturbance types when imputing BA removed, that data needs to be
## used for predicting any other missing values, AND it needs to be included
## in the model we eventually fit to analyze the data. Because we were
## planning to use separate models to analyze each disturbance type, this
## poses a problem. Two solutions I see are that are definitely kosher are: 1)
## fit everything in a single modeling framework, where there is an interaction
## between each dependent variable and a categorical variable describing the
## disturbance type; or 2) only predict the BA removed in each disturbance type
## using the data from that disturbance type. A third option which may or may not
## be okay is to perform MICE on N copies of the full dataset (all disturbances)
## and then split these by disturbance and run the analysis models, pooling
## results by disturbance. I *feel* like this is okay, but not 100% sure.

## For this demo, let's just try it with a single predictor:
## 1) treatment class (thinning vs. rx fire vs. both)

## to do this we first need to categorize each treatment, as they dont have
## standardized names

cat_trt <- function(trt) {
  if ((grepl("burn", tolower(trt)) | grepl("fire", tolower(trt))) &
      (grepl("thin", tolower(trt)) | grepl("density", tolower(trt)))) return("both")
  else if(grepl("burn", tolower(trt)) | grepl("fire", tolower(trt))) return("rx_fire")
  else if (grepl("thin", tolower(trt)) | grepl("density", tolower(trt))) return("thinning")
  else return(NA)
}
## TODO check what "Fifty (50) UMZ (low density stand)" means,
## assume thin for now

## also lets just look at mortality for now, ignoring biomass/carbon data
## so again we need to clean up the responseVariable column
## TODO deal with all the various non-mortality variables,
## short on atm so going to leave as is

insect_data_grouped$burn <- "no"
insect_data_grouped$thin <- "no"
for (i in 1:nrow(insect_data_grouped)) {
  insect_data_grouped[i,"trt_class"] <- cat_trt(insect_data_grouped[i,"treatment"])
  if (insect_data_grouped[i,"trt_class"] == "thinning") {
    insect_data_grouped[i,"thin"] <- "yes"
  } else if (insect_data_grouped[i,"trt_class"] == "rx_fire") {
    insect_data_grouped[i,"burn"] <- "yes"
  } else if (insect_data_grouped[i,"trt_class"] == "both") {
    insect_data_grouped[i,"thin"] <- "yes"
    insect_data_grouped[i,"burn"] <- "yes"
  }
}
insect_data_grouped$trt_class ## looks correct

## remove all extra columns, so we have clean dataset:
mice_data <- insect_data_grouped[insect_data_grouped$carbon_or_mortality == 2,
                                 c("lrr", "lrr_se", "trt_class")]
mice_data$trt_class <- as.factor(mice_data$trt_class) ## ensure factor encoded correctly

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

str(imputed_data) ## looks good

##---------------------------------------------------------------
## 4. Run models
##---------------------------------------------------------------

## we are going to try two different methods, frequentist and Bayes, just for fun

## frequentist first, using 'metafor'
freq_fit <- with(imputed_data,
                 rma(yi = lrr,
                     sei = lrr_se,
                     mods = ~ trt_class))

## now bayes
## takes a minute (or 5)
bayes_fit <- brm_multiple(lrr | se(lrr_se) ~ trt_class,
                          data = imputed_data,
                          chains = 4, cores = 4,
                          silent = 2, refresh = 0,
                          open_progress = FALSE)

##---------------------------------------------------------------
## 5. Results
##---------------------------------------------------------------

## first frequentist
## no significant results
pool <- summary(pool(freq_fit))
pool[-1] <- round(pool[-1], digits = 3)
pool

## then bayes
## intestingly, the HMC fit shows clear positive effect of rx fire,
## a clearly negative intercept (rx fire + thin), and a negative but unclear
## effect of thinning alone.
summary(bayes_fit)
plot(bayes_fit)

conditional_effects(bayes_fit, points = TRUE)
