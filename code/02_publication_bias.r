##---------------------------------------------------------------
## Tests for Publication Bias
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

## Egger's tests: ------------------------------------------------

## Mortality / Survivorship: -------------------------------------

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

bias_tests[[1]]$fit$beta

## combine results using Rubin's rules
M <- length(bias_tests)   # list of 100 regtest objects

slopes <- numeric(M)
vars   <- numeric(M)

for (m in 1:M) {
  obj <- bias_tests[[m]]
  slopes[m] <- obj$fit$beta["sei", 1]
  vars[m]   <- obj$fit$vb["sei", "sei"]
}

## Rubin's rules
beta_bar <- mean(slopes)
U_bar    <- mean(vars)
B        <- var(slopes)
T_var    <- U_bar + (1 + 1/M) * B
SE       <- sqrt(T_var)
t_stat   <- beta_bar / SE

## Degrees of freedom (Barnard–Rubin)
nu <- (M - 1) * (1 + U_bar / ((1 + 1/M) * B))^2

## p-values
p_val_t <- 2 * pt(abs(t_stat), df = nu, lower.tail = FALSE)
p_val_t # 0.036
beta_bar

#pvals <- sapply(bias_tests, function(x) if (!is.null(x)) x$pval else NA)
#hist(pvals, breaks = 100) ; abline(v = 0.05, col = "red", lwd = 4)
#summary(pvals, na.rm = TRUE)
#sum(pvals < 0.05) / length(pvals) # 75%

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

pdf("figures/funnel_mortality.pdf", width = 7, height = 6)
funnel(data$lrr, sei = data$se, col = data$cls, pch = data$pchs, cex.lab = cx, cex = cx, cex.axis = cx, xlab = "LRR", main = "Survivorship")
legend("topleft", col = c("#ff3319", "#0057ba", "#ffab00", "black", "black"), pch = c(rep(16,4),17), cex  = cx, legend = c("Fire", "Drought", "Insects", "Measured SE", "Imputed SE"))
text(2, 0.1, "Egger's test for asymmetry:", adj = c(0,0), cex = cx)
text(2, 0.2, "p = 0.036", adj = c(0,0), cex = cx)
box()
dev.off()

## Carbon: ------------------------------------------------------

## read data
data <- read.csv("data/processed_data/data_cleaned.csv")
data <- data[data$carbon_vs_mortality == 1, c("studyID", "disturbance_type", "trt_class","thin_bin", "burn_bin", "lrr", "lrr_se")]

## add colors and pchs:
data$cls <- "#ff3319"
data$cls[data$disturbance_type == "drought"] <- "#0057ba"
data$cls[data$disturbance_type == "insect"] <-  "#ffab00"
data$pchs <- 16
data$pchs[is.na(data$lrr_se)] <- 17

## read in overall carbon model
carb_fit <- readRDS("data/model_objects/carbon_fit_overall.rds")

pool_carb <- summary(pool(carb_fit))
pool_carb[-1] <- round(pool_carb[-1], digits = 3)
pool_carb

## assess publication bias for each MICE thing separately
bias_tests <- lapply(carb_fit$analyses, function(mod) {
  regtest(mod, model = "rma")
})

bias_tests[[1]]$fit$beta


## combine results using Rubin's rules
M <- length(bias_tests)   # list of 100 regtest objects

slopes <- numeric(M)
vars   <- numeric(M)

for (m in 1:M) {
  obj <- bias_tests[[m]]
  slopes[m] <- obj$fit$beta["sei", 1]
  vars[m]   <- obj$fit$vb["sei", "sei"]
}

## Rubin's rules
beta_bar <- mean(slopes)
U_bar    <- mean(vars)
B        <- var(slopes)
T_var    <- U_bar + (1 + 1/M) * B
SE       <- sqrt(T_var)
t_stat   <- beta_bar / SE

## Degrees of freedom (Barnard–Rubin)
nu <- (M - 1) * (1 + U_bar / ((1 + 1/M) * B))^2

## p-values
p_val_t <- 2 * pt(abs(t_stat), df = nu, lower.tail = FALSE)
p_val_t # 0.139
beta_bar

#pvals <- sapply(bias_tests, function(x) if (!is.null(x)) x$pval else NA)
#hist(pvals, breaks = 100) ; abline(v = 0.05, col = "red", lwd = 4)
#summary(pvals, na.rm = TRUE)
#sum(pvals < 0.05) / length(pvals) # 75%

## extract standard errors:
temp <- data.frame(vi1 = carb_fit$analyses[[1]]$vi)
for(i in 1:100){
  temp[,i] <- carb_fit$analyses[[i]]$vi
}
temp$mvi <- rowMeans(temp)
data$se <- sqrt(temp$mvi)

## extract residuals:
temp <- data.frame(res1 = rstandard(carb_fit$analyses[[1]])[[1]])
for(i in 1:100){
  temp[,i] <- rstandard(carb_fit$analyses[[i]])[[1]]
}
temp$mres <- rowMeans(temp)
data$resid <- temp$mres

## Funnel plot
cx <- 1

pdf("figures/funnel_carbon.pdf", width = 7, height = 6)
funnel(data$lrr, sei = data$se, col = data$cls, pch = data$pchs, cex.lab = cx, cex = cx, cex.axis = cx, xlab = "LRR", main = "Carbon")
legend("topleft", col = c("#ff3319", "#0057ba", "#ffab00", "black", "black"), pch = c(rep(16,4),17), cex  = cx, legend = c("Fire", "Drought", "Insects", "Measured SE", "Imputed SE"))
text(2, 0.1, "Egger's test for asymmetry:", adj = c(0,0), cex = cx)
text(2, 0.2, "p = 0.139", adj = c(0,0), cex = cx)
box()
dev.off()


##--------------------------------------------------------------
## Macaskill's tests (more robust to small sample size and heterogeneity)                                                  ##--------------------------------------------------------------

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

## for mortality
data <- read.csv("data/processed_data/data_cleaned.csv")

mort_imputed <- impute_data(data[data$carbon_vs_mortality == 2,], m = 100)

n <- data[data$carbon_vs_mortality == 2,"n_control"] + data[data$carbon_vs_mortality == 2,"n_treatment"]
yi <-data[data$carbon_vs_mortality == 2,"lrr"]
vi <- data[data$carbon_vs_mortality == 2,"lrr_se"]^2

ind <- which(is.na(vi))

effect_size <- numeric(100)
p <- numeric(100)
vars <- numeric(100)
intercept <- numeric(100)
for (i in 1:100) {
    vi[ind] <- mort_imputed$imp$lrr_se[[i]]^2
    s <- summary(lm(yi ~ n, weights = 1/vi))
    intercept[i] <- s$coefficients[1,1]
    effect_size[i] <- s$coefficients[2,1]
    vars[i] <- s$coefficients[2,2]^2
    p[i] <- s$coefficients[2,4]
}

M <- length(yi)

beta_bar <- mean(effect_size)
intercept_bar <- mean(intercept)

U_bar <- mean(vars)
B <- var(effect_size)
T_var <- U_bar + (1 + 1/M) * B
SE <- sqrt(T_var)
t_stat <- beta_bar / SE
nu <- (M - 1) * (1 + U_bar / ((1 + 1/M) * B))^2
p_val <- 2 * pt(abs(t_stat), df = nu, lower.tail = FALSE)
p_val ## overall p value, 0.267
beta_bar ## overall effect size

## plot output

pdf("figures/macaskill_survival.pdf", width = 7, height = 6)
plot(n, yi, xlab = "n", ylab = "LRR", pch = 16, cex = 1.5, main = "Survivorship", cex.lab = 1.5, cex.axis = 1.5)
abline(intercept_bar, beta_bar, col = "red", lwd = 2)
text(275, 4.3, "Macaskill's test:", adj = c(0,0), cex = 1)
text(275, 4, "p = 0.267", adj = c(0,0), cex = 1)
box()
dev.off()


## for carbon
mort_imputed <- impute_data(data[data$carbon_vs_mortality == 1,], m = 100)

n <- data[data$carbon_vs_mortality == 1,"n_control"] + data[data$carbon_vs_mortality == 1,"n_treatment"]
yi <-data[data$carbon_vs_mortality == 1,"lrr"]
vi <- data[data$carbon_vs_mortality == 1,"lrr_se"]^2

ind <- which(is.na(vi))

effect_size <- numeric(100)
p <- numeric(100)
vars <- numeric(100)
for (i in 1:100) {
    vi[ind] <- mort_imputed$imp$lrr_se[[i]]^2
    s <- summary(lm(yi ~ n, weights = 1/vi))
    effect_size[i] <- s$coefficients[2,1]
    vars[i] <- s$coefficients[2,2]^2
    p[i] <- s$coefficients[2,4]
}

M <- length(yi)

beta_bar <- mean(effect_size)

U_bar <- mean(vars)
B <- var(effect_size)
T_var <- U_bar + (1 + 1/M) * B
SE <- sqrt(T_var)
t_stat <- beta_bar / SE
nu <- (M - 1) * (1 + U_bar / ((1 + 1/M) * B))^2
p_val <- 2 * pt(abs(t_stat), df = nu, lower.tail = FALSE)
p_val ## overall p value, 0.603
beta_bar ## overall effect size


pdf("figures/macaskill_carbon.pdf", width = 7, height = 6)
plot(n, yi, xlab = "n", ylab = "LRR", pch = 16, cex = 1.5, main = "Carbon", cex.lab = 1.5, cex.axis = 1.5)
abline(intercept_bar, beta_bar, col = "red", lwd = 2)
text(60, 1.2, "Macaskill's test:", adj = c(0,0), cex = 1)
text(60, 1.05, "p = 0.603", adj = c(0,0), cex = 1)
box()
dev.off()
