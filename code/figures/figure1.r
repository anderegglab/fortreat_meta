##---------------------------------------------------------------
## Figure 1
##
load("/Users/evangeline/Desktop/Github/fortreat_meta/data/model_objects/figure1.RData")
save.image("~/Desktop/Github/fortreat_meta/data/model_objects/figure1.RData")
##---------------------------------------------------------------

## load packages
library("metafor")
library("ggplot2")
setwd('/Users/evangeline/Desktop/Github/fortreat_meta/')
cat_trt <- function(trt) {
  if (((grepl("burn", tolower(trt)) | grepl("fire", tolower(trt))) & !grepl("unburn", tolower(trt))) &
      (grepl("thin", tolower(trt)) | grepl("density", tolower(trt)) | grepl("umz", tolower(trt)))) return("both")
  else if((grepl("burn", tolower(trt)) | grepl("fire", tolower(trt))) & !grepl("unburn", tolower(trt))) return("rx_fire")
  else if(grepl("thin", tolower(trt)) | grepl("density", tolower(trt)) | grepl("umz", tolower(trt)))
    return("thinning")
  else return(NA)
}

#------- read in data --------
data <- read.csv("data/processed_data/data_cleaned.csv")
for (i in 1:nrow(data)) {
  data[i,"trt_class"] <- cat_trt(data[i,"treatment"])
}

## models outputs and imputed dataset
freq_fit_fire <- readRDS("data/model_objects/freq_fit_fire.rds")
freq_fit_insect <- readRDS("data/model_objects/freq_fit_insect.rds")
freq_fit_drought <- readRDS("data/model_objects/freq_fit_drought.rds")
imputed_fire <- readRDS("data/model_objects/imputed_fire.rds")
imputed_insect <- readRDS("data/model_objects/imputed_insect.rds")
imputed_drought <- readRDS("data/model_objects/imputed_drought.rds")

freq_fit_fire_C <- readRDS("data/model_objects/freq_fit_fire_C.rds")
freq_fit_insect_C <- readRDS("data/model_objects/freq_fit_insect_C.rds")
freq_fit_drought_C <- readRDS("data/model_objects/freq_fit_drought_C.rds")
imputed_fire_C <- readRDS("data/model_objects/imputed_fire_C.rds")
imputed_insect_C <- readRDS("data/model_objects/imputed_insect_C.rds")
imputed_drought_C <- readRDS("data/model_objects/imputed_drought_C.rds")

#
###-------- mortality/survivorship ------------
#------- insect --------
## first frequentist
## no significant results
summary(pool(freq_fit_insect))
pool_insect <- pool(freq_fit_insect)
pool_insect

# Extract coefficients and confidence intervals
coefs_freq_insect <- pool_insect$pooled$estimate  # Extract estimates
names(coefs_freq_insect) <- pool_insect$pooled$term  # Assign variable names
print(coefs_freq_insect)
se_freq_insect <- sqrt(pool_insect$pooled$t)  # Standard errors
cilower_freq_insect <- coefs_freq_insect - 1.96 * se_freq_insect
ciupper_freq_insect <- coefs_freq_insect + 1.96 * se_freq_insect
# extract sample size
ss_insect <- data.frame(table(complete(imputed_insect,1)$trt_class))

# Convert to a data frame
imputed_insect$data$trt_class <- as.factor(imputed_insect$data$trt_class)
plotdata_insect <- data.frame(Treatment = levels(imputed_insect$data$trt_class),
                            EffectSize = coefs_freq_insect, 
                            LowerCI = cilower_freq_insect,
                            UpperCI = ciupper_freq_insect,
                            SS=ss_insect[,2])

# Plot the estimated effects per treatment category
theme_set(theme_bw())
ggplot(plotdata_insect, aes(x = Treatment, y = EffectSize)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.05) +
  labs(title = "Insect: Effect Sizes by Treatment", x = "Treatment",y = "Effect Size")+
  geom_hline(yintercept = 0, color='gray') + 
  theme(panel.grid=element_blank(), axis.title = element_text(size=12), axis.text = element_text(size=12)) + geom_text(aes(label = SS, y = UpperCI + 0.05), size = 4)

#
#------- fire --------
## first frequentist
## no significant results
summary(pool(freq_fit_fire))
pool_fire <- pool(freq_fit_fire)
pool_fire

# Extract coefficients and confidence intervals
coefs_freq_fire <- pool_fire$pooled$estimate  # Extract estimates
names(coefs_freq_fire) <- pool_fire$pooled$term  # Assign variable names
print(coefs_freq_fire)
se_freq_fire <- sqrt(pool_fire$pooled$t)  # Standard errors
cilower_freq_fire <- coefs_freq_fire - 1.96 * se_freq_fire
ciupper_freq_fire <- coefs_freq_fire + 1.96 * se_freq_fire
# extract sample size
ss_fire <- data.frame(table(complete(imputed_fire,1)$trt_class))

# Convert to a data frame
imputed_fire$data$trt_class <- as.factor(imputed_fire$data$trt_class)
plotdata_fire <- data.frame(Treatment = levels(imputed_fire$data$trt_class),
                              EffectSize = coefs_freq_fire, 
                              LowerCI = cilower_freq_fire,
                              UpperCI = ciupper_freq_fire,
                              SS=ss_fire[,2])

# Plot the estimated effects per treatment category
theme_set(theme_bw())
ggplot(plotdata_fire, aes(x = Treatment, y = EffectSize)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.05) +
  labs(title = "Fire: Effect Sizes by Treatment", x = "Treatment",y = "Effect Size")+
  geom_hline(yintercept = 0, color='gray') + 
  theme(panel.grid=element_blank(), axis.title = element_text(size=12), axis.text = element_text(size=12)) + geom_text(aes(label = SS, y = UpperCI + 0.1), size = 4)

#
#------- drought --------
## first frequentist
## no significant results
summary(pool(freq_fit_drought))
pool_drought <- pool(freq_fit_drought)
pool_drought

# Extract coefficients and confidence intervals
coefs_freq_drought <- pool_drought$pooled$estimate  # Extract estimates
names(coefs_freq_drought) <- pool_drought$pooled$term  # Assign variable names
print(coefs_freq_drought)
se_freq_drought <- sqrt(pool_drought$pooled$t)  # Standard errors
cilower_freq_drought <- coefs_freq_drought - 1.96 * se_freq_drought
ciupper_freq_drought <- coefs_freq_drought + 1.96 * se_freq_drought
# extract sample size
ss_drought <- data.frame(table(complete(imputed_drought,1)$trt_class))

# Convert to a data frame
imputed_drought$data$trt_class <- as.factor(imputed_drought$data$trt_class)
plotdata_drought <- data.frame(Treatment = levels(imputed_drought$data$trt_class),
                            EffectSize = coefs_freq_drought, 
                            LowerCI = cilower_freq_drought,
                            UpperCI = ciupper_freq_drought,
                            SS=ss_drought[,2])

# Plot the estimated effects per treatment category
ggplot(plotdata_drought, aes(x = Treatment, y = EffectSize)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.05) +
  labs(title = "Drought: Effect Sizes by Treatment", x = "Treatment",y = "Effect Size")+
  geom_hline(yintercept = 0, color='gray') + 
  theme(panel.grid=element_blank(), axis.title = element_text(size=12), axis.text = element_text(size=12)) + geom_text(aes(label = SS, y = UpperCI + 0.05), size = 4)

#
###-------------------------------
###---------- carbon -------------
#------- insect --------
## first frequentist
## no significant results
summary(pool(freq_fit_insect_C))
pool_insect_C <- pool(freq_fit_insect_C)
pool_insect_C

# Extract coefficients and confidence intervals
coefs_freq_insect_C <- pool_insect_C$pooled$estimate  # Extract estimates
names(coefs_freq_insect_C) <- pool_insect_C$pooled$term  # Assign variable names
print(coefs_freq_insect_C)
se_freq_insect_C <- sqrt(pool_insect_C$pooled$t)  # Standard errors
cilower_freq_insect_C <- coefs_freq_insect_C - 1.96 * se_freq_insect_C
ciupper_freq_insect_C <- coefs_freq_insect_C + 1.96 * se_freq_insect_C
# extract sample size
ss_insect_C <- data.frame(table(complete(imputed_insect_C,1)$trt_class))

# Convert to a data frame
imputed_insect_C$data$trt_class <- as.factor(imputed_insect_C$data$trt_class)
plotdata_insect_C <- data.frame(Treatment = levels(imputed_insect_C$data$trt_class),
                              EffectSize = coefs_freq_insect_C, 
                              LowerCI = cilower_freq_insect_C,
                              UpperCI = ciupper_freq_insect_C,
                              SS=ss_insect_C[,2])

# Plot the estimated effects per treatment category
ggplot(plotdata_insect_C, aes(x = Treatment, y = EffectSize)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.05) +
  labs(title = "Insect: Effect Sizes by Treatment", x = "Treatment",y = "Effect Size")+
  geom_hline(yintercept = 0, color='gray') + 
  theme(panel.grid=element_blank(), axis.title = element_text(size=12), axis.text = element_text(size=12)) + geom_text(aes(label = SS, y = UpperCI + 0.1), size = 4)

#
#------- fire --------
## first frequentist
## no significant results
summary(pool(freq_fit_fire_C))
pool_fire_C <- pool(freq_fit_fire_C)
pool_fire_C

# Extract coefficients and confidence intervals
coefs_freq_fire_C <- pool_fire_C$pooled$estimate  # Extract estimates
names(coefs_freq_fire_C) <- pool_fire_C$pooled$term  # Assign variable names
print(coefs_freq_fire_C)
se_freq_fire_C <- sqrt(pool_fire_C$pooled$t)  # Standard errors
cilower_freq_fire_C <- coefs_freq_fire_C - 1.96 * se_freq_fire_C
ciupper_freq_fire_C <- coefs_freq_fire_C + 1.96 * se_freq_fire_C
# extract sample size
ss_fire_C <- data.frame(table(complete(imputed_fire_C,1)$trt_class))

# Convert to a data frame
imputed_fire_C$data$trt_class <- as.factor(imputed_fire_C$data$trt_class)
plotdata_fire_C <- data.frame(Treatment = levels(imputed_fire_C$data$trt_class),
                            EffectSize = coefs_freq_fire_C, 
                            LowerCI = cilower_freq_fire_C,
                            UpperCI = ciupper_freq_fire_C,
                            SS=ss_fire_C[,2])

# Plot the estimated effects per treatment category
ggplot(plotdata_fire_C, aes(x = Treatment, y = EffectSize)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.05) +
  labs(title = "Fire: Effect Sizes by Treatment", x = "Treatment",y = "Effect Size")+
  geom_hline(yintercept = 0, color='gray') + 
  theme(panel.grid=element_blank(), axis.title = element_text(size=12), axis.text = element_text(size=12)) + geom_text(aes(label = SS, y = UpperCI + 0.1), size = 4)

#
#------- drought --------
## first frequentist
## no significant results
summary(pool(freq_fit_drought_C))
pool_drought_C <- pool(freq_fit_drought_C)
pool_drought_C

# Extract coefficients and confidence intervals
coefs_freq_drought_C <- pool_drought_C$pooled$estimate  # Extract estimates
names(coefs_freq_drought_C) <- pool_drought_C$pooled$term  # Assign variable names
print(coefs_freq_drought_C)
se_freq_drought_C <- sqrt(pool_drought_C$pooled$t)  # Standard errors
cilower_freq_drought_C <- coefs_freq_drought_C - 1.96 * se_freq_drought_C
ciupper_freq_drought_C <- coefs_freq_drought_C + 1.96 * se_freq_drought_C
# extract sample size
ss_drought_C <- data.frame(table(complete(imputed_drought_C,1)$trt_class))

# Convert to a data frame
imputed_drought_C$data$trt_class <- as.factor(imputed_drought_C$data$trt_class)
plotdata_drought_C <- data.frame(Treatment = levels(imputed_drought_C$data$trt_class),
                               EffectSize = coefs_freq_drought_C, 
                               LowerCI = cilower_freq_drought_C,
                               UpperCI = ciupper_freq_drought_C,
                               SS=ss_drought_C[,2])

# Plot the estimated effects per treatment category
ggplot(plotdata_drought_C, aes(x = Treatment, y = EffectSize)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.05) +
  labs(title = "Drought: Effect Sizes by Treatment", x = "Treatment",y = "Effect Size")+
  geom_hline(yintercept = 0, color='gray') + 
  theme(panel.grid=element_blank(), axis.title = element_text(size=12), axis.text = element_text(size=12)) + geom_text(aes(label = SS, y = UpperCI + 0.1), size = 4)

#