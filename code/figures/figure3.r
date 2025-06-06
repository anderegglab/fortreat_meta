
##---------------------------------------------------------------
## Figure 3 -- overall effect of treatments on carbon
##
## author: Jacob Levine; contact: jacob.levine@utah.edu
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

carbon_imputed <- impute_data(data[data$carbon_vs_mortality == 1,], m = 100)

carbon_fit <- with(carbon_imputed,
                 rma(yi = lrr,
                     sei = lrr_se,
                     mods = ~ 0 + disturbance_type))
saveRDS(carbon_fit, "data/model_objects/carbon_fit_overall.rds")

pool_carbon <- summary(pool(carbon_fit))
pool_carbon[-1] <- round(pool_carbon[-1], digits = 3)
pool_carbon

table_gen(pool_carbon, "overall_trt_carbon.csv")

pdata <- data.frame(disturbance_type = factor(c("fire", "drought", "insect"), levels = c("fire", "drought", "insect")))
pdata$mean <- pool_carbon$estimate
pdata$se <- pool_carbon$std.error
pdata$lower <- pdata$mean - 1.97*pdata$se
pdata$upper <- pdata$mean + 1.97*pdata$se

A <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8) +
  geom_jitter(data = data[data$carbon_vs_mortality == 1,], aes(x = lrr, y = disturbance_type, color = disturbance_type), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata, aes(x = mean, y = disturbance_type, color = disturbance_type), size = 8) +
  geom_linerange(data = pdata, aes(y = disturbance_type, xmin = lower, xmax = upper, color = disturbance_type), size = 3) +
  scale_color_manual(values = c(red, blue, yellow)) +
  annotate("text", y = "drought", x = 1.5, label = "*", size = 15) +
  xlab("Log Response Ratio") +
  ggtitle("Treatment Effects: Carbon") +
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

fire_carbon_imputed <- impute_data(data[data$disturbance_type == "fire" & data$carbon_vs_mortality == 1,], m = 100)
insect_carbon_imputed <- impute_data(data[data$disturbance_type == "insect" & data$carbon_vs_mortality == 1,], m = 100)
drought_carbon_imputed <- impute_data(data[data$disturbance_type == "drought" & data$carbon_vs_mortality == 1,], m = 100)

fire_carbon_fit <- with(fire_carbon_imputed,
                      rma(yi = lrr,
                          sei = lrr_se,
                          mods = ~ 0 + thin_bin * burn_bin))
saveRDS(fire_carbon_fit, "data/model_objects/fire_carbon_fit.rds")

drought_carbon_fit <- with(drought_carbon_imputed,
                      rma(yi = lrr,
                          sei = lrr_se,
                          mods = ~ 0 + thin_bin * burn_bin))
saveRDS(drought_carbon_fit, "data/model_objects/drought_carbon_fit.rds")

insect_carbon_fit <- with(insect_carbon_imputed,
                      rma(yi = lrr,
                          sei = lrr_se,
                          mods = ~ 0 + thin_bin * burn_bin))
saveRDS(insect_carbon_fit, "data/model_objects/insect_carbon_fit.rds")

pool_fire_carbon <- summary(pool(fire_carbon_fit))
pool_fire_carbon[-1] <- round(pool_fire_carbon[-1], digits = 3)
pool_fire_carbon

table_gen(pool_fire_carbon, "trtclass_fire_carbon.csv")

pool_drought_carbon <- summary(pool(drought_carbon_fit))
pool_drought_carbon[-1] <- round(pool_drought_carbon[-1], digits = 3)
pool_drought_carbon

table_gen(pool_drought_carbon, "trtclass_drought_carbon.csv")

pool_insect_carbon <- summary(pool(insect_carbon_fit))
pool_insect_carbon[-1] <- round(pool_insect_carbon[-1], digits = 3)
pool_insect_carbon

table_gen(pool_insect_carbon, "trtclass_insect_carbon.csv")

pdata <- data.frame(trt_class = factor(c("thinning", "rx_fire", "both"), levels = c("thinning", "rx_fire", "both")))
pdata$disturbance_type = "fire"
pdata$mean[1] <- pool_fire_carbon$estimate[1]
pdata$mean[2] <- pool_fire_carbon$estimate[2]
pdata$mean[3] <- pool_fire_carbon$estimate[1] + pool_fire_carbon$estimate[2] + pool_fire_carbon$estimate[3]
pdata$se <- pool_fire_carbon$std.error
pdata$lower <- pdata$mean - 1.97*pdata$se
pdata$upper <- pdata$mean + 1.97*pdata$se

B <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8) +
  geom_jitter(data = data[data$disturbance_type == "fire" & data$carbon_vs_mortality == 1,], aes(x = lrr, y = trt_class, color = disturbance_type), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata, aes(x = mean, y = trt_class, color = disturbance_type), size = 6) +
  geom_linerange(data = pdata, aes(y = trt_class, xmin = lower, xmax = upper, color = disturbance_type), size = 3) +
  annotate("text", y = "thinning", x = 4.6, label = "*", size = 15) +
  scale_color_manual(values = c(red)) +
  xlab("Log Response Ratio") +
  xlim(-1.5, 5) +
  theme_bw() +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(face = "italic", size = 12))
B


pdata <- data.frame(trt_class = factor(c("thinning", "rx_fire", "both"), levels = c("thinning", "rx_fire", "both")))
pdata$disturbance_type = "drought"
pdata$mean[1] <- pool_drought_carbon$estimate[1]
pdata$mean[2] <- pool_drought_carbon$estimate[2]
pdata$mean[3] <- pool_drought_carbon$estimate[1] + pool_drought_carbon$estimate[2] + pool_drought_carbon$estimate[3]
pdata$se <- pool_drought_carbon$std.error
pdata$lower <- pdata$mean - 1.97*pdata$se
pdata$upper <- pdata$mean + 1.97*pdata$se

C <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8) +
  geom_jitter(data = data[data$disturbance_type == "drought" & data$carbon_vs_mortality == 1,], aes(x = lrr, y = trt_class, color = disturbance_type), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata, aes(x = mean, y = trt_class, color = disturbance_type), size = 6) +
  geom_linerange(data = pdata, aes(y = trt_class, xmin = lower, xmax = upper, color = disturbance_type), size = 3) +
  annotate("text", y = "both", x = 4.6, label = "*", size = 15) +
  annotate("text", y = "thinning", x = 4.6, label = "*", size = 15) +
  scale_color_manual(values = c(blue)) +
  xlim(-1.5, 5) +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank(),
        axis.text.y = element_text(face = "italic", size = 12))
C

pdata <- data.frame(trt_class = factor(c("thinning", "rx_fire", "both"), levels = c("thinning", "rx_fire", "both")))
pdata$disturbance_type = "insect"
pdata$mean[1] <- pool_insect_carbon$estimate[1]
pdata$mean[2] <- pool_insect_carbon$estimate[2]
pdata$mean[3] <- pool_insect_carbon$estimate[1] + pool_insect_carbon$estimate[2] + pool_insect_carbon$estimate[3]
pdata$se <- pool_insect_carbon$std.error
pdata$lower <- pdata$mean - 1.97*pdata$se
pdata$upper <- pdata$mean + 1.97*pdata$se

D <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8) +
  geom_jitter(data = data[data$disturbance_type == "insect" & data$carbon_vs_mortality == 1,], aes(x = lrr, y = trt_class, color = disturbance_type), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata, aes(x = mean, y = trt_class, color = disturbance_type), size = 6) +
  geom_linerange(data = pdata, aes(y = trt_class, xmin = lower, xmax = upper, color = disturbance_type), size = 3) +
  annotate("text", y = "both", x = 4.6, label = "", size = 15) + # Weird but makes the xaxis be in the right order
  scale_color_manual(values = c(yellow)) +
  xlim(-1.5, 5) +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank(),
        axis.text.y = element_text(face = "italic", size = 12))
D


A + D + C + B +
  plot_layout(
  design = "
AB
AC
AE")

ggsave("figures/figure3.png", width = 8, height = 8)
ggsave("figures/illustrator/figure3.pdf", width = 8, height = 8)
