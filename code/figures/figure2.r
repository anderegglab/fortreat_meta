
##---------------------------------------------------------------
## Figure 2 -- overall effect of treatments on survivorship
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

source("code/functions.r")

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

mort_imputed <- impute_data(data[data$carbon_vs_mortality == 2,], m = 100)
plot(mort_imputed)

mort_fit <- with(mort_imputed,
                 rma(yi = lrr,
                     sei = lrr_se,
                     mods = ~ 0 + disturbance_type))
saveRDS(mort_fit, "data/model_objects/mort_fit_overall.rds")

pool_mort <- summary(pool(mort_fit))
pool_mort[-1] <- round(pool_mort[-1], digits = 3)
pool_mort

## overall survivorship plot
pdata <- data.frame(disturbance_type = factor(c("fire", "drought", "insect"), levels = c("fire", "drought", "insect")))
pdata$mean <- pool_mort$estimate
pdata$se <- pool_mort$std.error
pdata$lower <- pdata$mean - 1.97*pdata$se
pdata$upper <- pdata$mean + 1.97*pdata$se

A <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8) +
  geom_jitter(data = data[data$carbon_vs_mortality == 2,], aes(x = lrr, y = disturbance_type, color = disturbance_type), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata, aes(x = mean, y = disturbance_type, color = disturbance_type), size = 12) +
  geom_linerange(data = pdata, aes(y = disturbance_type, xmin = lower, xmax = upper, color = disturbance_type), size = 3) +
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  scale_x_continuous(limits = c(-0.5, 5.5)) +
  ggtitle("Treatment Effects: Survivorship") +
  annotate("text", y = "fire", x = 5, label = "***", size = 15) +
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
saveRDS(fire_mort_fit, "data/model_objects/fire_mort_fit.rds")

drought_mort_fit <- with(drought_mort_imputed,
                      rma(yi = lrr,
                          sei = lrr_se,
                          mods = ~ 0 + thin_bin * burn_bin))
saveRDS(drought_mort_fit, "data/model_objects/drought_mort_fit.rds")

insect_mort_fit <- with(insect_mort_imputed,
                      rma(yi = lrr,
                          sei = lrr_se,
                          mods = ~ 0 + thin_bin * burn_bin))
saveRDS(insect_mort_fit, "data/model_objects/insect_mort_fit.rds")

pool_fire_mort <- summary(pool(fire_mort_fit))
pool_fire_mort[-1] <- round(pool_fire_mort[-1], digits = 3)
pool_fire_mort

pool_drought_mort <- summary(pool(drought_mort_fit))
pool_drought_mort[-1] <- round(pool_drought_mort[-1], digits = 3)
pool_drought_mort

pool_insect_mort <- summary(pool(insect_mort_fit))
pool_insect_mort[-1] <- round(pool_insect_mort[-1], digits = 3)
pool_insect_mort

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
  geom_point(data = pdata, aes(x = mean, y = trt_class, color = disturbance_type), size = 8) +
  geom_linerange(data = pdata, aes(y = trt_class, xmin = lower, xmax = upper, color = disturbance_type), size = 3) +
  scale_color_manual(values = c("#e41a1c")) +
  annotate("text", y = "both", x = 4.8, label = "***", size = 15) +
  annotate("text", y = "thinning", x = 4.8, label = "***", size = 15) +
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
  geom_jitter(data = data[data$disturbance_type == "drought" & data$carbon_vs_mortality == 2,], aes(x = lrr, y = trt_class, color = disturbance_type), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata, aes(x = mean, y = trt_class, color = disturbance_type), size = 8) +
  geom_linerange(data = pdata, aes(y = trt_class, xmin = lower, xmax = upper, color = disturbance_type), size = 3) +
  annotate("text", y = "both", x = 4.8, label = "***", size = 15) +
  scale_color_manual(values = c("#377eb8")) +
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
  geom_jitter(data = data[data$disturbance_type == "insect" & data$carbon_vs_mortality == 2,], aes(x = lrr, y = trt_class, color = disturbance_type), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata, aes(x = mean, y = trt_class, color = disturbance_type), size = 8) +
  geom_linerange(data = pdata, aes(y = trt_class, xmin = lower, xmax = upper, color = disturbance_type), size = 3) +
  annotate("text", y = "thinning", x = 4.8, label = "***", size = 15) +
  scale_color_manual(values = c("#4daf4a")) +
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

ggsave("figures/figure2.png", width = 14, height = 12)
