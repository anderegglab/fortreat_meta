##---------------------------------------------------------------
## Figure 5 -- Treatment effects by Biome (survivorship)
##
## author: Jacob Levine; contact: jacob.levine@utah.edu
##---------------------------------------------------------------

## libraries
library("metafor")
library("mice")
library("brms")
library("ggplot2")
library("cowplot")
library("patchwork")

source("code/00_functions.r")

## read data
data <- read.csv("data/processed_data/data_cleaned.csv")

##---------------------------------------------------------------
## 1. Impute missing data and fit models
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

mort_imputed_tcon <- impute_data(data[data$Olson_Biome == "Temperate Conifer Forests" & data$carbon_vs_mortality == 2,], m = 100)

mort_fit_tcon <- with(mort_imputed_tcon,
                 rma(yi = lrr,
                     sei = lrr_se,
                     mods = ~ 0 + disturbance_type))

pool_mort_tcon <- summary(pool(mort_fit_tcon))
pool_mort_tcon[-1] <- round(pool_mort_tcon[-1], digits = 3)
pool_mort_tcon

table_gen(pool_mort_tcon, "biome_tcon_mort.csv")

data[data$Olson_Biome == "Temperate Broadleaf and Mixed Forests" & data$carbon_vs_mortality == 2,]
mort_imputed_tbro <- impute_data(data[data$Olson_Biome == "Temperate Broadleaf and Mixed Forests" & data$carbon_vs_mortality == 2,], m = 100)

mort_fit_tbro <- with(mort_imputed_tbro,
                 rma(yi = lrr,
                     sei = lrr_se,
                     mods = ~ 0 + disturbance_type))

pool_mort_tbro <- summary(pool(mort_fit_tbro))
pool_mort_tbro[-1] <- round(pool_mort_tbro[-1], digits = 3)
pool_mort_tbro

table_gen(pool_mort_tbro, "biome_tbro_mort.csv")

mort_imputed_med <- impute_data(data[data$Olson_Biome == "Mediterranean Forests, Woodlands and Scrub" & data$carbon_vs_mortality == 2,], m = 100)

mort_fit_med <- with(mort_imputed_med,
                 rma(yi = lrr,
                     sei = lrr_se,
                     mods = ~ 0 + disturbance_type))

pool_mort_med <- summary(pool(mort_fit_med))
pool_mort_med[-1] <- round(pool_mort_med[-1], digits = 3)
pool_mort_med

table_gen(pool_mort_med, "biome_med_mort.csv")


##---------------------------------------------------------------
## 2. BarPlot
##---------------------------------------------------------------

data$disturbance_type <- factor(data$disturbance_type, levels = c("fire", "drought", "insect"))

A <- ggplot(data[data$carbon_vs_mortality == 2, ]) +
  geom_bar(aes(x = Olson_Biome, fill = disturbance_type)) +
  labs(x = "Biome", y = "Count", fill = "Disturbance") +
  annotate("text", x = "Deserts and Xeric Shrublands", y = 10, label = "Deserts and \n Xeric Shrublands", size = 4) +
  annotate("text", x = "Mediterranean Forests, Woodlands and Scrub", y = 13, label = "Mediterranean Forests, \n Woodlands and Scrub", size = 4) +
  annotate("text", x = "Temperate Broadleaf and Mixed Forests", y = 17, label = "Temperate Broadleaf \n and Mixed Forests", size = 4) +
  annotate("text", x = "Temperate Conifer Forests", y = 105, label = "Temperate Conifer \n Forests", size = 4) +
  ggtitle("Biome Effects: Survivorship") +
  scale_fill_manual(values = c(red, yellow, blue)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 110)) +
  theme_bw() +
  theme(legend.position = "inside", legend.position.inside = c(0.1, 0.9), axis.text.x = element_blank(),
        plot.title = element_text(face = "bold", size = 18))
A

##---------------------------------------------------------------
## 3. Effect Plots
##---------------------------------------------------------------

## Temperate Conifer Forests
pd <- data.frame(disturbance_typefire = c(1, 0, 0),
           disturbance_typedrought = c(0, 1, 0),
           disturbance_typeinsect = c(0, 0, 1))

out <- list(mean = matrix(nrow = 3, ncol = 100), ci.lb = matrix(nrow = 3, ncol = 100), ci.ub = matrix(nrow = 3, ncol = 100))

for (i in 1:100) {

  p <- predict.rma(mort_fit_tcon$analyses[[i]], newmods = as.matrix(pd))
  out[[1]][,i] <- p[[1]]
  out[[2]][,i] <- p[[3]]
  out[[3]][,i] <- p[[4]]

}

pdata <- data.frame(disturbance_type = factor(c("fire", "drought", "insect"), levels = c("fire", "drought", "insect")))

pdata$mean <- rowMeans(out[[1]])
pdata$lower <- rowMeans(out[[2]])
pdata$upper <- rowMeans(out[[3]])

B <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8) +
  geom_jitter(data = data[data$Olson_Biome == "Temperate Conifer Forests" & data$carbon_vs_mortality == 2,], aes(x = lrr, y = disturbance_type, color = disturbance_type), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata, aes(x = mean, y = disturbance_type, color = disturbance_type), size = 8) +
  geom_linerange(data = pdata, aes(y = disturbance_type, xmin = lower, xmax = upper, color = disturbance_type), size = 3) +
  annotate("text", x = 4.9, y = "fire", label = "*", size = 14) +
  scale_color_manual(values = c(red, blue, yellow)) +
  ggtitle("Temperate Conifer Forests") +
  xlim(-1.2, 5) +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank())
B

## Temperate Broadleaf and Mixed Forests
pd <- data.frame(disturbance_typedrought = c(1, 0),
                 disturbance_typeinsect = c(0, 1))

out <- list(mean = matrix(nrow = 2, ncol = 100), ci.lb = matrix(nrow = 2, ncol = 100), ci.ub = matrix(nrow = 2, ncol = 100))

for (i in 1:100) {

  p <- predict.rma(mort_fit_tbro$analyses[[i]], newmods = as.matrix(pd))
  out[[1]][,i] <- p[[1]]
  out[[2]][,i] <- p[[3]]
  out[[3]][,i] <- p[[4]]

}

pdata <- data.frame(disturbance_type = factor(c("drought", "insect"), levels = c("fire", "drought", "insect")))

pdata$mean <- rowMeans(out[[1]])
pdata$lower <- rowMeans(out[[2]])
pdata$upper <- rowMeans(out[[3]])

C <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8) +
  geom_jitter(data = data[data$Olson_Biome == "Temperate Broadleaf and Mixed Forests" & data$carbon_vs_mortality == 2,], aes(x = lrr, y = disturbance_type, color = disturbance_type), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata, aes(x = mean, y = disturbance_type, color = disturbance_type), size = 8) +
  geom_linerange(data = pdata, aes(y = disturbance_type, xmin = lower, xmax = upper, color = disturbance_type), size = 3) +
  #annotate("text", x = 4.9, y = "insect", label = "***", size = 14) +
  #annotate("text", x = 4.9, y = "drought", label = "***", size = 14) +
  scale_color_manual(values = c(blue, yellow)) +
  ggtitle("Temperate Broadleaf and Mixed Forests") +
  xlim(-1.2, 5) +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank())
C

## Mediterranean Forests, Woodlands and Scrub
pd <- data.frame(disturbance_typedrought = c(1, 0),
                 disturbance_typefire = c(0, 1))

out <- list(mean = matrix(nrow = 2, ncol = 100), ci.lb = matrix(nrow = 2, ncol = 100), ci.ub = matrix(nrow = 2, ncol = 100))

for (i in 1:100) {

  p <- predict.rma(mort_fit_med$analyses[[i]], newmods = as.matrix(pd))
  out[[1]][,i] <- p[[1]]
  out[[2]][,i] <- p[[3]]
  out[[3]][,i] <- p[[4]]

}

pdata <- data.frame(disturbance_type = factor(c("drought", "fire"), levels = c("fire", "drought", "insect")))

pdata$mean <- rowMeans(out[[1]])
pdata$lower <- rowMeans(out[[2]])
pdata$upper <- rowMeans(out[[3]])

D <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8) +
  geom_jitter(data = data[data$Olson_Biome == "Mediterranean Forests, Woodlands and Scrub" & data$carbon_vs_mortality == 2,], aes(x = lrr, y = disturbance_type, color = disturbance_type), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata, aes(x = mean, y = disturbance_type, color = disturbance_type), size = 8) +
  geom_linerange(data = pdata, aes(y = disturbance_type, xmin = lower, xmax = upper, color = disturbance_type), size = 3) +
  annotate("text", x = 4.9, y = "fire", label = "*", size = 14) +
  scale_color_manual(values = c(red, blue)) +
  xlim(-1.2, 5) +
  ggtitle("Mediterranean Forests, Woodlands and Scrub") +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank())
D


A + D + B + C +
  plot_layout(
  design = "
AB
AC
AE")

ggsave("figures/figure5.png", width = 13, height = 10)
ggsave("figures/illustrator/figure5.pdf", width = 13, height = 10)
