
##---------------------------------------------------------------
## Figure 2 -- overall effect of treatments on survivorship
## Sensitivity analysis -- Dropping all uncertain studies (N = 2, no true control and uncertain mortality attribution)
##
## author: Cedric Zahnd
##         Jacob Levine; contact: jacob.levine@utah.edu
##---------------------------------------------------------------

##---------------------------------------------------------------
## 0. Head
##---------------------------------------------------------------

## libraries
library("metafor")
library("mice")
library("brms")
library("ggplot2")
library("patchwork")
library("cowplot")

source("code/00_functions.r")

## read data
data <- read.csv("data/processed_data/data_cleaned.csv")

## Drop uncertain studies:
data <- data[data$Nis2 == 0 & data$NoTrueControl == 0 & data$mortality_attribution_uncertain == 0 & data$zeroSurv == 0,]

## TODO For now I'm not saving any tables or Rdata files, can change later if needed.

##---------------------------------------------------------------
## 1 Overall
##---------------------------------------------------------------

data$disturbance_type <- factor(data$disturbance_type, levels = c("fire", "drought", "insect"))

impute_data <- function(data, vars = c("lrr", "lrr_se", "disturbance_type"), m = 20) {
  
  mice_data <- data[,vars]
  
  ## make predictor matrix
  predictor_matrix <- make.predictorMatrix(mice_data)
  predictor_matrix ## looks good
  
  impute_method <- make.method(mice_data)
  impute_method ## no method specified for complete variables
  
  imputed_data <- mice(mice_data, method = impute_method, predictorMatrix = predictor_matrix,
                       maxit = 40, seed = 1, m = m)
  
  return(imputed_data)
  
}

mort_imputed <- impute_data(data[data$carbon_vs_mortality == 2,], m = 100)
plot(mort_imputed)

mort_fit <- with(mort_imputed,
                 rma(yi = lrr,
                     sei = lrr_se,
                     mods = ~ 0 + disturbance_type))
#saveRDS(mort_fit, "data/model_objects/mort_fit_overall.rds")

pool_mort <- summary(pool(mort_fit))
pool_mort[-1] <- round(pool_mort[-1], digits = 3)
pool_mort

#table_gen(pool_mort, "overall_trt_mort.csv")

## overall survivorship plot
pdata <- data.frame(disturbance_type = factor(c("fire", "drought", "insect"), levels = c("fire", "drought", "insect")))
pdata$mean <- pool_mort$estimate
pdata$se <- pool_mort$std.error
pdata$lower <- pdata$mean - 1.97*pdata$se
pdata$upper <- pdata$mean + 1.97*pdata$se

A <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8) +
  geom_jitter(data = data[data$carbon_vs_mortality == 2,], aes(x = lrr, y = disturbance_type, color = disturbance_type), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata, aes(x = mean, y = disturbance_type, color = disturbance_type), size = 8) +
  geom_linerange(data = pdata, aes(y = disturbance_type, xmin = lower, xmax = upper, color = disturbance_type), size = 3) +
  scale_color_manual(values = c(red, blue, yellow)) +
  scale_x_continuous(limits = c(-0.5, 5.5)) +
  ggtitle("Treatment Effects: Survivorship") +
  annotate("text", y = "fire", x = 5, label = "*", size = 15) +
  xlab("Log Response Ratio") +
  theme_bw() +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14, face = "italic"), plot.title = element_text(face = "bold", size = 18))
A

##---------------------------------------------------------------
## 2 Broken out by treatment type
##---------------------------------------------------------------

impute_data <- function(data, vars = c("lrr", "lrr_se", "burn_bin", "thin_bin"), m = 20) {
  
  mice_data <- data[,vars]
  mice_data <- cbind(mice_data, thin.burn = 0)
  
  ## make predictor matrix
  predictor_matrix <- make.predictorMatrix(mice_data)
  predictor_matrix ## looks good
  
  impute_method <- make.method(mice_data)
  impute_method ## no method specified for complete variables
  
  x <- mice(mice_data, max = 0)
  meth <- x$meth
  meth["thin.burn"] <- "~I(thin_bin*burn_bin)"
  
  pred <- x$pred
  pred[c("burn_bin", "thin_bin"), c("thin.burn")] <- 0
  imputed_data <- mice(mice_data, meth = meth, pred = pred, maxit = 40, seed = 1, m = m)
  
  return(imputed_data)
  
}

fire_mort_imputed <- impute_data(data[data$disturbance_type == "fire" & data$carbon_vs_mortality == 2,], m = 100)
drought_mort_imputed <- impute_data(data[data$disturbance_type == "drought" & data$carbon_vs_mortality == 2,], m = 100)
insect_mort_imputed <- impute_data(data[data$disturbance_type == "insect" & data$carbon_vs_mortality == 2,], m = 100)

fire_mort_fit <- with(fire_mort_imputed,
                      rma(yi = lrr,
                          sei = lrr_se,
                          mods = ~ 0 + thin_bin + thin_bin:burn_bin))
#saveRDS(fire_mort_fit, "data/model_objects/fire_mort_fit.rds")

drought_mort_fit <- with(drought_mort_imputed,
                         rma(yi = lrr,
                             sei = lrr_se,
                             mods = ~ 0 + thin_bin * burn_bin))
#saveRDS(drought_mort_fit, "data/model_objects/drought_mort_fit.rds")

insect_mort_fit <- with(insect_mort_imputed,
                        rma(yi = lrr,
                            sei = lrr_se,
                            mods = ~ 0 + thin_bin * burn_bin))
#saveRDS(insect_mort_fit, "data/model_objects/insect_mort_fit.rds")

pool_fire_mort <- summary(pool(fire_mort_fit))
pool_fire_mort[-1] <- round(pool_fire_mort[-1], digits = 3)
pool_fire_mort

#table_gen(pool_fire_mort, "trtclass_fire_mort.csv")

pool_drought_mort <- summary(pool(drought_mort_fit))
pool_drought_mort[-1] <- round(pool_drought_mort[-1], digits = 3)
pool_drought_mort

#table_gen(pool_drought_mort, "trtclass_drought_mort.csv")

pool_insect_mort <- summary(pool(insect_mort_fit))
pool_insect_mort[-1] <- round(pool_insect_mort[-1], digits = 3)
pool_insect_mort

#table_gen(pool_insect_mort, "trtclass_insect_mort.csv")

pdata <- data.frame(trt_class = factor(c("thinning", "both"), levels = c("thinning", "both")))
pdata$disturbance_type = "fire"
pdata$mean[1] <- pool_fire_mort$estimate[1]
pdata$mean[2] <- pool_fire_mort$estimate[1] + pool_fire_mort$estimate[2]
pdata$se <- pool_fire_mort$std.error
pdata$lower <- pdata$mean - 1.97*pdata$se
pdata$upper <- pdata$mean + 1.97*pdata$se

B <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8) +
  geom_jitter(data = data[data$disturbance_type == "fire" & data$carbon_vs_mortality == 2,], aes(x = lrr, y = trt_class, color = disturbance_type), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata, aes(x = mean, y = trt_class, color = disturbance_type), size = 6) +
  geom_linerange(data = pdata, aes(y = trt_class, xmin = lower, xmax = upper, color = disturbance_type), size = 3) +
  scale_color_manual(values = c(red)) +
  annotate("text", y = "both", x = 4.8, label = "*", size = 15) +
  annotate("text", y = "thinning", x = 4.8, label = "*", size = 15) +
  xlab("Log Response Ratio") +
  xlim(-0.5, 5) +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.y = element_text(face = "italic", size = 12), axis.title.x = element_text(size = 16),)
B

pdata <- data.frame(trt_class = factor(c("thinning", "rx_fire", "both"), levels = c("thinning", "rx_fire", "both")))
pdata$disturbance_type = "drought"
pdata$mean[1] <- pool_drought_mort$estimate[1]
pdata$mean[2] <- pool_drought_mort$estimate[2]
pdata$mean[3] <- pool_drought_mort$estimate[1] + pool_drought_mort$estimate[2] + pool_drought_mort$estimate[3]
pdata$se <- pool_drought_mort$std.error
pdata$lower <- pdata$mean - 1.97*pdata$se
pdata$upper <- pdata$mean + 1.97*pdata$se

C <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8) +
  geom_jitter(data = data[data$disturbance_type == "drought" & data$carbon_vs_mortality == 2,], aes(x = lrr, y = trt_class, color = trt_class), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata, aes(x = mean, y = trt_class, color = trt_class), size = 6) +
  geom_linerange(data = pdata, aes(y = trt_class, xmin = lower, xmax = upper, color = trt_class), size = 3) +
  annotate("text", y = "thinning", x = 4.8, label = "*", size = 15, col = "white") +
  scale_color_manual(values = c(blue, "grey", "grey")) +
  xlim(-0.5, 5) +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.y = element_text(face = "italic", size = 12))
C

pdata <- data.frame(trt_class = factor(c("thinning", "rx_fire", "both"), levels = c("thinning", "rx_fire", "both")))
pdata$disturbance_type = "insect"
pdata$mean[1] <- pool_insect_mort$estimate[1]
pdata$mean[2] <- pool_insect_mort$estimate[2]
pdata$mean[3] <- pool_insect_mort$estimate[1] + pool_insect_mort$estimate[2] + pool_insect_mort$estimate[3]
pdata$se <- pool_insect_mort$std.error
pdata$lower <- pdata$mean - 1.97*pdata$se
pdata$upper <- pdata$mean + 1.97*pdata$se

D <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8) +
  geom_jitter(data = data[data$disturbance_type == "insect" & data$carbon_vs_mortality == 2,], aes(x = lrr, y = trt_class, color = trt_class), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata, aes(x = mean, y = trt_class, color = trt_class), size = 6) +
  geom_linerange(data = pdata, aes(y = trt_class, xmin = lower, xmax = upper, color = trt_class), size = 3) +
  annotate("text", y = "thinning", x = 4.8, label = "*", size = 15) +
  scale_color_manual(values = c(yellow, yellow, "grey")) +
  xlim(-0.5, 5) +
  theme_bw() +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_text(face = "italic", size = 12))
D


A + D + C + B +
  plot_layout(
    design = "
AB
AC
AE")

ggsave("figures/sensitivity/sensitivity_figure2.png", width = 8, height = 8)
ggsave("figures/illustrator/sensitivity_figure2.pdf", width = 8, height = 8)


## --------------------------------------------------------------
## Figure out why mortality - drought - thinning changes:
## --------------------------------------------------------------

impute_data <- function(data, vars = c("lrr", "lrr_se", "burn_bin", "thin_bin"), m = 20) {
  
  mice_data <- data[,vars]
  mice_data <- cbind(mice_data, thin.burn = 0)
  
  ## make predictor matrix
  predictor_matrix <- make.predictorMatrix(mice_data)
  predictor_matrix ## looks good
  
  impute_method <- make.method(mice_data)
  impute_method ## no method specified for complete variables
  
  x <- mice(mice_data, max = 0)
  meth <- x$meth
  meth["thin.burn"] <- "~I(thin_bin*burn_bin)"
  
  pred <- x$pred
  pred[c("burn_bin", "thin_bin"), c("thin.burn")] <- 0
  imputed_data <- mice(mice_data, meth = meth, pred = pred, maxit = 40, seed = 1, m = m)
  
  return(imputed_data)
  
}

## read data
data_full <- read.csv("data/processed_data/data_cleaned.csv")

## Drop Nis2: - Not significant (dropped 1 / 4)
data <- data_full[data_full$Nis2 == 0,]
## Drop NoTrueC: - Not significant (dropped 1 / 4)
data <- data_full[data_full$NoTrueControl == 0,]
## Drop MortAttr: - Significant (dropped 1 / 3)
data <- data_full[data_full$mortality_attribution_uncertain == 0,]

##
## Both Nis2 and NoTrueC drop study 436 Cochran and Seidl 1999
## That is therefore the study driving that effect.
## Remaining: 7 studies, 15 observations
##

drought_mort_imputed <- impute_data(data[data$disturbance_type == "drought" & data$carbon_vs_mortality == 2,], m = 100)

drought_mort_fit <- with(drought_mort_imputed,
                         rma(yi = lrr,
                             sei = lrr_se,
                             mods = ~ 0 + thin_bin * burn_bin))
#saveRDS(drought_mort_fit, "data/model_objects/drought_mort_fit.rds")

#table_gen(pool_fire_mort, "trtclass_fire_mort.csv")

pool_drought_mort <- summary(pool(drought_mort_fit))
pool_drought_mort[-1] <- round(pool_drought_mort[-1], digits = 3)
pool_drought_mort

#table_gen(pool_drought_mort, "trtclass_drought_mort.csv")


pdata <- data.frame(trt_class = factor(c("thinning", "rx_fire", "both"), levels = c("thinning", "rx_fire", "both")))
pdata$disturbance_type = "drought"
pdata$mean[1] <- pool_drought_mort$estimate[1]
pdata$mean[2] <- pool_drought_mort$estimate[2]
pdata$mean[3] <- pool_drought_mort$estimate[1] + pool_drought_mort$estimate[2] + pool_drought_mort$estimate[3]
pdata$se <- pool_drought_mort$std.error
pdata$lower <- pdata$mean - 1.97*pdata$se
pdata$upper <- pdata$mean + 1.97*pdata$se

C <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8) +
  geom_jitter(data = data[data$disturbance_type == "drought" & data$carbon_vs_mortality == 2,], aes(x = lrr, y = trt_class, color = trt_class), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata, aes(x = mean, y = trt_class, color = trt_class), size = 8) +
  geom_linerange(data = pdata, aes(y = trt_class, xmin = lower, xmax = upper, color = trt_class), size = 3) +
  annotate("text", y = "both", x = 4.8, label = "", size = 15, col = "grey") +
  scale_color_manual(values = c("#377eb8", "grey", "grey")) +
  xlim(-0.5, 5) +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank(), axis.text.y = element_text(face = "italic", size = 12))
C






