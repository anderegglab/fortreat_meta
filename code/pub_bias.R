##---------------------------------------------------------------
## Tests for Publication Bias
## Survivorship
##---------------------------------------------------------------

## libraries
library("dplyr")
library("metafor")
library("mice")
library("brms")
library("ggplot2")
library("patchwork")
library("cowplot")

source("code/00_functions.r")

## read data
data <- read.csv("data/processed_data/data_cleaned.csv")
data <- data[data$carbon_vs_mortality == 2, c("studyID", "disturbance_type", "trt_class","thin_bin", "burn_bin", "lrr", "lrr_se")]

## add colors and pchs:
data$cls <- "#ff3319"
data$cls[data$disturbance_type == "drought"] <- "#0057ba"
data$cls[data$disturbance_type == "insect"] <-  "#ffab00"
data$pchs <- 16
data$pchs[is.na(data$lrr_se)] <- 17

## read in overall mortality (survivorship) model
mort_fit <- readRDS("data/model_objects/mort_fit_overall.rds")

pool_mort <- summary(pool(mort_fit))
pool_mort[-1] <- round(pool_mort[-1], digits = 3)
pool_mort

## assess publication bias for each MICE thing separately
bias_tests <- lapply(mort_fit$analyses, function(mod) {
  regtest(mod, model = "rma")
})

pvals <- sapply(bias_tests, function(x) if (!is.null(x)) x$pval else NA)
hist(pvals, breaks = 100) ; abline(v = 0.05, col = "red", lwd = 4)
summary(pvals, na.rm = TRUE)
sum(pvals < 0.05) / length(pvals) # 75%

## extract standard errors:
temp <- data.frame(vi1 = mort_fit$analyses[[1]]$vi)
for(i in 1:100){
  temp[,i] <- mort_fit$analyses[[i]]$vi
}
temp$mvi <- rowMeans(temp)
data$se <- sqrt(temp$mvi)

## extract residuals:
temp <- data.frame(res1 = rstandard(mort_fit$analyses[[1]])[[1]])
for(i in 1:100){
  temp[,i] <- rstandard(mort_fit$analyses[[i]])[[1]]
}
temp$mres <- rowMeans(temp)
data$resid <- temp$mres

## Funnel plot
cx <- 1

png("figures/funnel_mortality.png", width = 7, height = 6, units = "in", res = 72)
funnel(data$lrr, sei = data$se, col = data$cls, pch = data$pchs, cex.lab = cx, cex = cx, cex.axis = cx, xlab = "LRR", main = "Survivorship")
legend("topleft", col = c("#ff3319", "#0057ba", "#ffab00", "black", "black"), pch = c(rep(16,4),17), cex  = cx, legend = c("Fire", "Drought", "Insects", "Measured SE", "Imputed SE"))
text(2, 0.1, "Egger's test for asymmetry:", adj = c(0,0), cex = cx)
text(2, 0.2, "75% p < 0.05", adj = c(0,0), cex = cx)
box()
dev.off()

## Test average residuals instead:
regtest(data$resid, sei = data$se, model = "rma")

png("figures/funnel_mortality_resid.png", width = 7, height = 6, units = "in", res = 72)
funnel(data$resid, sei = data$se, col = data$cls, pch = data$pchs, cex.lab = cx, cex = cx, cex.axis = cx, xlab = "Residuals", main = "Survivorship")
legend("topleft", col = c("#ff3319", "#0057ba", "#ffab00", "black", "black"), pch = c(rep(16,4),17), cex  = cx, legend = c("Fire", "Drought", "Insects", "Measured SE", "Imputed SE"))
text(1, 0.1, "Egger's test for asymmetry:", adj = c(0,0), cex = cx)
text(1, 0.2, "z = 3.94, p < 0.001", adj = c(0,0), cex = cx)
box()
dev.off()

##---------------------------------------------------------------
## Reanalyse without North & Hurteau (studyID = 7)
##---------------------------------------------------------------

data <- data[data$studyID != 7, c(1:9)]

##---------------------------------------------------------------
## Overall
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

mort_imputed <- impute_data(data, m = 100)
plot(mort_imputed)

mort_fit <- with(mort_imputed,
                 rma(yi = lrr,
                     sei = lrr_se,
                     mods = ~ 0 + disturbance_type))

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
  geom_jitter(data = data[data$disturbance_type == "fire",], aes(x = lrr, y = disturbance_type, color = disturbance_type), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata[pdata$disturbance_type == "fire",], aes(x = mean, y = disturbance_type, color = disturbance_type), size = 8) +
  geom_linerange(data = pdata[pdata$disturbance_type == "fire",], aes(y = disturbance_type, xmin = lower, xmax = upper, color = disturbance_type), size = 3) +
  scale_color_manual(values = c(red)) +
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

fire_mort_imputed <- impute_data(data[data$disturbance_type == "fire",], m = 100)

fire_mort_fit <- with(fire_mort_imputed,
                      rma(yi = lrr,
                          sei = lrr_se,
                          mods = ~ 0 + thin_bin + thin_bin:burn_bin))

pool_fire_mort <- summary(pool(fire_mort_fit))
pool_fire_mort[-1] <- round(pool_fire_mort[-1], digits = 3)
pool_fire_mort

pdata <- data.frame(trt_class = factor(c("thinning", "both"), levels = c("thinning", "both")))
pdata$disturbance_type = "fire"
pdata$mean[1] <- pool_fire_mort$estimate[1]
pdata$mean[2] <- pool_fire_mort$estimate[1] + pool_fire_mort$estimate[2]
pdata$se <- pool_fire_mort$std.error
pdata$lower <- pdata$mean - 1.97*pdata$se
pdata$upper <- pdata$mean + 1.97*pdata$se

B <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8) +
  geom_jitter(data = data[data$disturbance_type == "fire",], aes(x = lrr, y = trt_class, color = disturbance_type), height = 0.2, size = 3, alpha = 0.5) +
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

A + B +
  plot_layout(
    design = "
AB")

ggsave("figures/sensitivity/sensitivity_bias_figure2.png", width = 8, height = 4)
ggsave("figures/illustrator/sensitivity_bias_figure2.pdf", width = 8, height = 4)

##---------------------------------------------------------------
## Redo publication bias analysis
##---------------------------------------------------------------

## assess publication bias for each MICE thing separately
bias_tests <- lapply(mort_fit$analyses, function(mod) {
  regtest(mod, model = "rma")
})

pvals <- sapply(bias_tests, function(x) if (!is.null(x)) x$pval else NA)
hist(pvals, breaks = 100) ; abline(v = 0.05, col = "red", lwd = 4)
summary(pvals, na.rm = TRUE)
sum(pvals < 0.05) / length(pvals) # 84%

## extract standard errors:
temp <- data.frame(vi1 = mort_fit$analyses[[1]]$vi)
for(i in 1:100){
  temp[,i] <- mort_fit$analyses[[i]]$vi
}
temp$mvi <- rowMeans(temp)
data$se <- sqrt(temp$mvi)

## extract residuals:
temp <- data.frame(res1 = rstandard(mort_fit$analyses[[1]])[[1]])
for(i in 1:100){
  temp[,i] <- rstandard(mort_fit$analyses[[i]])[[1]]
}
temp$mres <- rowMeans(temp)
data$resid <- temp$mres

## Funnel plot
cx <- 1

png("figures/funnel_mortality_without7.png", width = 7, height = 6, units = "in", res = 72)
funnel(data$lrr, sei = data$se, col = data$cls, pch = data$pchs, cex.lab = cx, cex = cx, cex.axis = cx, xlab = "LRR", main = "Survivorship")
#legend("topleft", col = c("#ff3319", "#0057ba", "#ffab00", "black", "black"), pch = c(rep(16,4),17), cex  = cx, legend = c("Fire", "Drought", "Insects", "Measured SE", "Imputed SE"))
text(1, 0.1, "Egger's test for asymmetry:", adj = c(0,0), cex = cx)
text(1, 0.2, "84% p < 0.05", adj = c(0,0), cex = cx)
box()
dev.off()

#data[data$lrr > 2.5,]

## Test average residuals instead:
regtest(data$resid, sei = data$se, model = "rma")

png("figures/funnel_mortality_resid_without7.png", width = 7, height = 6, units = "in", res = 72)
funnel(data$resid, sei = data$se, col = data$cls, pch = data$pchs, cex.lab = cx, cex = cx, cex.axis = cx, xlab = "Residuals", main = "Survivorship")
#legend("topleft", col = c("#ff3319", "#0057ba", "#ffab00", "black", "black"), pch = c(rep(16,4),17), cex  = cx, legend = c("Fire", "Drought", "Insects", "Measured SE", "Imputed SE"))
text(1, 0.1, "Egger's test for asymmetry:", adj = c(0,0), cex = cx)
text(1, 0.2, "z = 2.53, p = 0.011", adj = c(0,0), cex = cx)
box()
dev.off()

#regtest(rstandard(mort_fit$analyses[[1]])[[1]], sei = sqrt(mort_fit$analyses[[1]]$vi), model = "rma")
#funnel(rstandard(mort_fit$analyses[[1]])[[1]], sei = sqrt(mort_fit$analyses[[1]]$vi))

##---------------------------------------------------------------
## Carbon
##---------------------------------------------------------------

## read data
data <- read.csv("data/processed_data/data_cleaned.csv")
data <- data[data$carbon_vs_mortality == 1, c("studyID", "disturbance_type", "trt_class","thin_bin", "burn_bin", "lrr", "lrr_se")]

## add colors and pchs:
data$cls <- "#ff3319"
data$cls[data$disturbance_type == "drought"] <- "#0057ba"
data$cls[data$disturbance_type == "insect"] <-  "#ffab00"
data$pchs <- 16
data$pchs[is.na(data$lrr_se)] <- 17

## read in overall mortality (survivorship) model
carbon_fit <- readRDS("data/model_objects/carbon_fit_overall.rds")

pool_carbon <- summary(pool(carbon_fit))
pool_carbon[-1] <- round(pool_carbon[-1], digits = 3)
pool_carbon

## assess publication bias for each MICE thing separately
bias_tests <- lapply(carbon_fit$analyses, function(mod) {
  regtest(mod, model = "rma")
})

pvals <- sapply(bias_tests, function(x) if (!is.null(x)) x$pval else NA)
hist(pvals, breaks = 100) ; abline(v = 0.05, col = "red", lwd = 4)
summary(pvals, na.rm = TRUE)
sum(pvals < 0.05) / length(pvals) # 8%

## extract standard errors:
temp <- data.frame(vi1 = carbon_fit$analyses[[1]]$vi)
for(i in 1:100){
  temp[,i] <- carbon_fit$analyses[[i]]$vi
}
temp$mvi <- rowMeans(temp)
data$se <- sqrt(temp$mvi)

## extract residuals:
temp <- data.frame(res1 = rstandard(carbon_fit$analyses[[1]])[[1]])
for(i in 1:100){
  temp[,i] <- rstandard(carbon_fit$analyses[[i]])[[1]]
}
temp$mres <- rowMeans(temp)
data$resid <- temp$mres

## Funnel plot
cx <- 1

png("figures/funnel_carbon.png", width = 7, height = 6, units = "in", res = 72)
funnel(data$lrr, sei = data$se, col = data$cls, pch = data$pchs, cex.lab = cx, cex = cx, cex.axis = cx, xlab = "LRR", main = "Carbon stocks")
#legend("topleft", col = c("#ff3319", "#0057ba", "#ffab00", "black", "black"), pch = c(rep(16,4),17), cex  = cx, legend = c("Fire", "Drought", "Insects", "Measured SE", "Imputed SE"))
text(-3, 0.1, "Egger's test for asymmetry:", adj = c(0,0), cex = cx)
text(-3, 0.2, "8% p < 0.05", adj = c(0,0), cex = cx)
box()
dev.off()

## Test average residuals instead:
regtest(data$resid, sei = data$se, model = "rma")

png("figures/funnel_carbon_resid.png", width = 7, height = 6, units = "in", res = 72)
funnel(data$resid, sei = data$se, col = data$cls, pch = data$pchs, cex.lab = cx, cex = cx, cex.axis = cx, xlab = "Residuals", main = "Carbon stocks")
#legend("topleft", col = c("#ff3319", "#0057ba", "#ffab00", "black", "black"), pch = c(rep(16,4),17), cex  = cx, legend = c("Fire", "Drought", "Insects", "Measured SE", "Imputed SE"))
text(1, 0.1, "Egger's test for asymmetry:", adj = c(0,0), cex = cx)
text(1, 0.2, "z = 1.21, p = 0.227", adj = c(0,0), cex = cx)
box()
dev.off()
