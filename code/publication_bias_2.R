##---------------------------------------------------------------
## Tests for Publication Bias
##
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

##---------------------------------------------------------------
## Survivorship
##---------------------------------------------------------------

## read in overall mortality (survivorship) model
mort_fit <- readRDS("data/model_objects/mort_fit_overall.rds")

pool_mort <- summary(pool(mort_fit))
pool_mort[-1] <- round(pool_mort[-1], digits = 3)
pool_mort

## assess publication bias for each MICE thing separately
bias_tests <- lapply(mort_fit$analyses, function(mod) {
  regtest(mod, model = "lm")
})

pvals <- sapply(bias_tests, function(x) if (!is.null(x)) x$pval else NA)
hist(pvals, breaks = 100)
abline(v = 0.05, col = "red", lwd = 4)
summary(pvals, na.rm = TRUE)
sum(pvals < 0.05) / length(pvals)

ggplot(data.frame(pvals = pvals), aes(x = pvals)) +
  geom_histogram(binwidth = 0.05, color = "white", fill = "black") +
  geom_vline(xintercept = 0.05, color = "red", size = 2) +
  annotate("text", x = 0.2, y = 55, label = "66% < 0.05", size = 12) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 60)) +
  ggtitle("Eggers Test: Survivorship") +
  xlab("P Values") +
  ylab("Count") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), plot.title = element_text(face = "bold", size = 18))

ggsave("figures/publication_bias_mortality.png", width = 12, height = 8)

## assess average publication bias:

out <- list()
for(i in 1:length(mort_fit$analyses)){
  temp <- rstandard(mort_fit$analyses[[i]])
  temp <- data.frame(resid = temp$resid, se = temp$se)
  temp$rownr <- 1:nrow(temp)
  out[[i]] <- temp ; rm(temp)
}

mort_res <- do.call(rbind, out)
mort_avg <- mort_res %>% group_by(rownr) %>% summarise(mres = mean(resid), mse = mean(se))

data <- read.csv("data/processed_data/data_cleaned.csv")
mort_data <- data[data$carbon_vs_mortality == 2,c("studyID", "disturbance_type", "trt_class", "impute_sd", "lrr_se")]
mort_data$resid <- mort_avg$mres ; mort_data$se <- mort_avg$mse
mort_data$cls <- "blue"
mort_data$cls[is.na(mort_data$lrr_se)] <- "red"

regtest(mort_data$resid, sei = mort_data$se, model = "lm")

cx <- 1.5
png("figures/publication_bias_funnel_mortality.png", width = 12, height = 8, units = "in", res = 72)
funnel(mort_data$resid, sei = mort_data$se, col = mort_data$cls, xlab = "Residual value", cex.lab = cx, cex = cx, cex.axis = cx)
legend("topright", col = c("blue", "red"), legend = c("Measured SE", "Imputed SE"), pch = 16, cex = cx)
text(-4, 0.1, "Egger's test for asymmetry:", adj = c(0,0), cex = cx)
text(-4, 0.25, "t = 5.73, df = 126, p < 0.001", adj = c(0,0), cex = cx)
dev.off()

##---------------------------------------------------------------
## For Carbon
##---------------------------------------------------------------

carbon_fit <- readRDS("data/model_objects/carbon_fit_overall.rds")

pool_carbon <- summary(pool(carbon_fit))
pool_carbon[-1] <- round(pool_carbon[-1], digits = 3)
pool_carbon

## assess publication bias for each MICE thing separately
bias_tests <- lapply(carbon_fit$analyses, function(mod) {
  regtest(mod, model = "lm")
})

pvals <- sapply(bias_tests, function(x) if (!is.null(x)) x$pval else NA)
hist(pvals, breaks = 100)
abline(v = 0.05, col = "red", lwd = 4)
summary(pvals, na.rm = TRUE)
sum(pvals < 0.05) / length(pvals)

ggplot(data.frame(pvals = pvals), aes(x = pvals)) +
  geom_histogram(binwidth = 0.05, fill = "black", color = "white") +
  geom_vline(xintercept = 0.05, color = "red", size = 2) +
  annotate("text", x = 0.2, y = 9.5, label = "0% < 0.05", size = 12) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 10)) +
  ggtitle("Eggers Test: Carbon") +
  xlab("P Values") +
  ylab("Count") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), plot.title = element_text(face = "bold", size = 18))

ggsave("figures/publication_bias_carbon.png", width = 12, height = 8)

## assess average publication bias:

out <- list()
for(i in 1:length(carbon_fit$analyses)){
  temp <- rstandard(carbon_fit$analyses[[i]])
  temp <- data.frame(resid = temp$resid, se = temp$se)
  temp$rownr <- 1:nrow(temp)
  out[[i]] <- temp ; rm(temp)
}

carbon_res <- do.call(rbind, out)
carbon_avg <- carbon_res %>% group_by(rownr) %>% summarise(mres = mean(resid), mse = mean(se))

data <- read.csv("data/processed_data/data_cleaned.csv")
carbon_data <- data[data$carbon_vs_mortality == 1,c("studyID", "disturbance_type", "trt_class", "impute_sd", "lrr_se")]
carbon_data$resid <- carbon_avg$mres ; carbon_data$se <- carbon_avg$mse
carbon_data$cls <- "blue"
carbon_data$cls[is.na(carbon_data$lrr_se)] <- "red"

regtest(carbon_data$resid, sei = carbon_data$se, model = "lm")

cx <- 1.5
png("figures/publication_bias_funnel_carbon.png", width = 12, height = 8, units = "in", res = 72)
funnel(carbon_data$resid, sei = carbon_data$se, col = carbon_data$cls, xlab = "Residual value", cex.lab = cx, cex = cx, cex.axis = cx)
legend("topright", col = c("blue", "red"), legend = c("Measured SE", "Imputed SE"), pch = 16, cex = cx)
text(-3, 0.1, "Egger's test for asymmetry:", adj = c(0,0), cex = cx)
text(-3, 0.25, "t = 0.58, df = 89, p < 0.562", adj = c(0,0), cex = cx)
dev.off()

