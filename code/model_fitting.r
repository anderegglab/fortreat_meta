
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

## read data
insect_data <- read.csv("data/insect_data.csv")

## rename columns with egregiously long names
colnames(insect_data)[33] <- "grouping_flags"

##---------------------------------------------------------------
## 1. Group data
##---------------------------------------------------------------

## TODO Move functions to separate file, leave here now for readability

## define functions
overall_mean <- function(mean, n) {
  return(sum(mean * n) / sum(n))
}

overall_sd <- function(mean, n, sd) {
  return({
    sqrt((1 / (sum(n) - 1)) *
         (sum((n - 1) * sd^2) + sum(n * (mean - overall_mean(mean, n)))))
    })
}

combine_se  <- function(mean, n, se, N) {
  return({
    sqrt(((N[1] * se[1]^2 + N[2] * se[2]^2) / (N[1] + N[2])) +
         ((n[1] * n[1] * (mean[1] - mean[2])^2) / ((n[1] + n[2]) * (N[1] + N[2]))))
  })
}

overall_se <- function(mean, n, se) {
  while(TRUE) {
    N <- n^2 - n
    ov_se <- combine_se(mean[1:2], n[1:2], se[1:2], N[1:2])
    if (length(mean) == 2) {
      break
    }
    mean <- c(mean(mean[1:2]), mean[3:length(mean)])
    n <- c(sum(n[1:2]), n[3:length(n)])
    se <- c(ov_se, se[3:length(se)])
  }
  return(ov_se)
}

group_data <- function(data, remove_zeros = FALSE, add_constant = FALSE) {

  if (remove_zeros) {
    data <- data[data$treatment_se != 0 & data$control_se != 0,]
  } else if(add_constant) {
    data[data$se_treatment == 0 & !is.na(data$se_treatment), c("se_treatment", "sd_treatment")] <- 1e-10
    data[data$se_control == 0 & !is.na(data$se_control), c("se_control", "sd_control")] <- 1e-10
  }

  data_grouped <- data[0,]

  ## loop through studies
  for (id in unique(data$studyID)) {

    sdata <- data[data$studyID == id,] ## pull all data from ID

    if (!is.na(sdata[1,"grouping_flags"])) {

      for (grp in unique(sdata$grouping_flags)) {

        gdata <- sdata[sdata$grouping_flags == grp,] ## pull data with same group
        ndata <- sdata[1,]

        ndata[,c("mean_treatment", "n_treatment", "sd_treatment", "se_treatment")] <-
         c(overall_mean(gdata$mean_treatment, gdata$n_treatment),
           sum(gdata$n_treatment),
           overall_sd(gdata$mean_treatment, gdata$n_treatment, gdata$sd_treatment),
           overall_se(gdata$mean_treatment, gdata$n_treatment, gdata$se_treatment))

        ndata[,c("mean_control", "n_control", "sd_control", "se_control")] <-
          c(overall_mean(gdata$mean_control, gdata$n_control),
            sum(gdata$n_control),
            overall_sd(gdata$mean_control, gdata$n_control, gdata$sd_control),
            overall_se(gdata$mean_control, gdata$n_control, gdata$se_control))

        data_grouped <- rbind(data_grouped, ndata)
      }

    } else {

      data_grouped <- rbind(data_grouped, sdata)

    }
  }
  return(data_grouped)
}

## TODO CHECK STUDY 306:
## one of the controls has 0 sd, but positive mean with n = 4. Seems unlikely?

## run grouping script
insect_data_grouped <- group_data(insect_data)

## This includes some standard errors in the grouping which were entered as 0
## We can either A: leave as is, or

## B. remove all observations with SE=0
insect_data_grouped_nozeros <- group_data(insect_data, remove_zeros = TRUE)

## C. add very small number to se and sd that equal zero
insect_data_grouped_wconst <- group_data(insect_data, add_constant = TRUE)

##---------------------------------------------------------------
## 2. Calculate effect sizes and standard errors
##---------------------------------------------------------------

## fn to calculate log response ratio
lrr <- function(mean_1, mean_2) {
  return(log(mean_1 / mean_2))
}

## fn to calculate standard error of log response ratio
lrr_se <- function(mean_1, mean_2, se_1, se_2) {
  return({
    sqrt((se_1^2 / mean_1^2) + (se_2^2  / mean_2^2))
  })
}

## calculate log response ratio
insect_data_grouped$lrr <- lrr(insect_data_grouped$mean_treatment,
                               insect_data_grouped$mean_control)
## produces NaNs and -Inf when treatment mean is 0. Need to decide what
## to do here, but for now I will toss them.

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
cat_response <- function(response) {
  if (grepl("mort", tolower(response))) return("mortality")
  else return(NA)
}

for (i in 1:nrow(insect_data_grouped)) {
  insect_data_grouped[i,"trt_class"] <- cat_trt(insect_data_grouped[i,"treatment"])
  insect_data_grouped[i,"response_class"] <- cat_response(insect_data_grouped[i,"responseVariable"])
}
insect_data_grouped$trt_class ## looks correct

## remove all extra columns, so we have clean dataset:
mice_data <- insect_data_grouped[!is.na(insect_data_grouped$response_class),
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
bayes_fit <- brm_multiple(lrr | se(lrr_se) ~ trt_class, data = imputed_data,
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
