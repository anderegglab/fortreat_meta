##---------------------------------------------------------------
## Figure 6 -- Treatment effects by Biome (carbon)
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

source("code/functions.r")

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

carbon_imputed_tcon <- impute_data(data[data$Olson_Biome == "Temperate Conifer Forests" & data$carbon_vs_mortality == 1,], m = 100)

carbon_fit_tcon <- with(carbon_imputed_tcon,
                 rma(yi = lrr,
                     sei = lrr_se,
                     mods = ~ 0 + disturbance_type))

pool_carbon_tcon <- summary(pool(carbon_fit_tcon))
pool_carbon_tcon[-1] <- round(pool_carbon_tcon[-1], digits = 3)
pool_carbon_tcon

carbon_imputed_tbro <- impute_data(data[data$Olson_Biome == "Temperate Broadleaf and Mixed Forests" & data$carbon_vs_mortality == 1,], m = 100)

carbon_fit_tbro <- with(carbon_imputed_tbro,
                 rma(yi = lrr,
                     sei = lrr_se,
                     mods = ~ 0 + disturbance_type))

pool_carbon_tbro <- summary(pool(carbon_fit_tbro))
pool_carbon_tbro[-1] <- round(pool_carbon_tbro[-1], digits = 3)
pool_carbon_tbro

carbon_imputed_med <- impute_data(data[data$Olson_Biome == "Mediterranean Forests, Woodlands and Scrub" & data$carbon_vs_mortality == 1,], m = 100)

carbon_fit_med <- with(carbon_imputed_med,
                 rma(yi = lrr,
                     sei = lrr_se))

pool_carbon_med <- summary(pool(carbon_fit_med))
pool_carbon_med[-1] <- round(pool_carbon_med[-1], digits = 3)
pool_carbon_med


carbon_imputed_des <- impute_data(data[data$Olson_Biome == "Deserts and Xeric Shrublands" & data$carbon_vs_mortality == 1,], m = 100)

carbon_fit_des <- with(carbon_imputed_des,
                 rma(yi = lrr,
                     sei = lrr_se))

pool_carbon_des <- summary(pool(carbon_fit_des))
pool_carbon_des[-1] <- round(pool_carbon_des[-1], digits = 3)
pool_carbon_des


##---------------------------------------------------------------
## 2. Barplot
##---------------------------------------------------------------

data$disturbance_type <- factor(data$disturbance_type, levels = c("fire", "drought", "insect"))

A <- ggplot(data[data$carbon_vs_mortality == 1, ]) +
  geom_bar(aes(x = Olson_Biome, fill = disturbance_type)) +
  labs(x = "Biome", y = "Count", fill = "Disturbance") +
  annotate("text", x = "Deserts and Xeric Shrublands", y = 12, label = "Deserts and \n Xeric Shrublands", size = 5) +
  annotate("text", x = "Mediterranean Forests, Woodlands and Scrub", y = 32, label = "Mediterranean Forests, \n Woodlands and Scrub", size = 5) +
  annotate("text", x = "Temperate Broadleaf and Mixed Forests", y = 18, label = "Temperate Broadleaf \n and Mixed Forests", size = 5) +
  annotate("text", x = "Temperate Conifer Forests", y = 45, label = "Temperate Conifer \n Forests", size = 5) +
  ggtitle("Biome Effects: Carbon") +
  scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 100)) +
  theme_bw() +
  theme(legend.position = "inside", legend.position.inside = c(0.1, 0.9), axis.text.x = element_blank(),
        plot.title = element_text(face = "bold", size = 18))
A

## mortality
pd <- data.frame(disturbance_typefire = c(1, 0, 0),
           disturbance_typedrought = c(0, 1, 0),
           disturbance_typeinsect = c(0, 0, 1))

out <- list(mean = matrix(nrow = 3, ncol = 100), ci.lb = matrix(nrow = 3, ncol = 100), ci.ub = matrix(nrow = 3, ncol = 100))

for (i in 1:100) {

  p <- predict.rma(carbon_fit_tcon$analyses[[i]], newmods = as.matrix(pd))
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
  geom_jitter(data = data[data$Olson_Biome == "Temperate Conifer Forests" & data$carbon_vs_mortality == 1,], aes(x = lrr, y = disturbance_type, color = disturbance_type), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata, aes(x = mean, y = disturbance_type, color = disturbance_type), size = 8) +
  geom_linerange(data = pdata, aes(y = disturbance_type, xmin = lower, xmax = upper, color = disturbance_type), size = 3) +
  annotate("text", x = 4.9, y = "fire", label = "***", size = 14) +
  annotate("text", x = 4.9, y = "drought", label = "***", size = 14) +
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  ggtitle("Temperate Conifer Forests") +
  xlim(-1.2, 5) +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank())
B

pd <- data.frame(disturbance_typefire = c(1, 0, 0),
           disturbance_typedrought = c(0, 1, 0),
           disturbance_typeinsect = c(0, 0, 1))

out <- list(mean = matrix(nrow = 3, ncol = 100), ci.lb = matrix(nrow = 3, ncol = 100), ci.ub = matrix(nrow = 3, ncol = 100))

for (i in 1:100) {

  p <- predict.rma(carbon_fit_tbro$analyses[[i]], newmods = as.matrix(pd))
  out[[1]][,i] <- p[[1]]
  out[[2]][,i] <- p[[3]]
  out[[3]][,i] <- p[[4]]

}

pdata <- data.frame(disturbance_type = factor(c("fire", "drought", "insect"), levels = c("fire", "drought", "insect")))

pdata$mean <- rowMeans(out[[1]])
pdata$lower <- rowMeans(out[[2]])
pdata$upper <- rowMeans(out[[3]])

C <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8) +
  geom_jitter(data = data[data$Olson_Biome == "Temperate Broadleaf and Mixed Forests" & data$carbon_vs_mortality == 1,], aes(x = lrr, y = disturbance_type, color = disturbance_type), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata, aes(x = mean, y = disturbance_type, color = disturbance_type), size = 8) +
  geom_linerange(data = pdata, aes(y = disturbance_type, xmin = lower, xmax = upper, color = disturbance_type), size = 3) +
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  ggtitle("Temperate Conifer Forests") +
  xlim(-1.2, 5) +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank())
C


## mortality
pd <- data.frame(disturbance_typedrought = c(1))

out <- list(mean = matrix(nrow = 1, ncol = 100), ci.lb = matrix(nrow = 1, ncol = 100), ci.ub = matrix(nrow = 1, ncol = 100))

for (i in 1:100) {

  p <- predict.rma(carbon_fit_med$analyses[[i]])
  out[[1]][,i] <- p[[1]]
  out[[2]][,i] <- p[[3]]
  out[[3]][,i] <- p[[4]]

}

pdata <- data.frame(disturbance_type = factor(c("drought"), levels = c("fire", "drought", "insect")))

pdata$mean <- rowMeans(out[[1]])
pdata$lower <- rowMeans(out[[2]])
pdata$upper <- rowMeans(out[[3]])

D <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8) +
  geom_jitter(data = data[data$Olson_Biome == "Mediterranean Forests, Woodlands and Scrub" & data$carbon_vs_mortality == 1,], aes(x = lrr, y = disturbance_type, color = disturbance_type), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata, aes(x = mean, y = disturbance_type, color = disturbance_type), size = 8) +
  geom_linerange(data = pdata, aes(y = disturbance_type, xmin = lower, xmax = upper, color = disturbance_type), size = 3) +
  annotate("text", x = 4.9, y = "drought", label = "***", size = 14) +
  scale_color_manual(values = c("#377eb8")) +
  xlim(-1.2, 5) +
  ggtitle("Mediterranean Forests, Woodlands and Scrub") +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank())
D


## mortality
pd <- data.frame(disturbance_typefire = c(1))

out <- list(mean = matrix(nrow = 1, ncol = 100), ci.lb = matrix(nrow = 1, ncol = 100), ci.ub = matrix(nrow = 1, ncol = 100))

for (i in 1:100) {

  p <- predict.rma(carbon_fit_des$analyses[[i]])
  out[[1]][,i] <- p[[1]]
  out[[2]][,i] <- p[[3]]
  out[[3]][,i] <- p[[4]]

}

pdata <- data.frame(disturbance_type = factor(c("fire"), levels = c("fire", "drought", "insect")))

pdata$mean <- rowMeans(out[[1]])
pdata$lower <- rowMeans(out[[2]])
pdata$upper <- rowMeans(out[[3]])

E <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8) +
  geom_jitter(data = data[data$Olson_Biome == "Deserts and Xeric Shrublands" & data$carbon_vs_mortality == 1,], aes(x = lrr, y = disturbance_type, color = disturbance_type), height = 0.2, size = 3, alpha = 0.5) +
  geom_point(data = pdata, aes(x = mean, y = disturbance_type, color = disturbance_type), size = 8) +
  geom_linerange(data = pdata, aes(y = disturbance_type, xmin = lower, xmax = upper, color = disturbance_type), size = 3) +
  annotate("text", x = 4.9, y = "fire", label = "***", size = 14) +
  scale_color_manual(values = c("#e41a1c")) +
  xlim(-1.2, 5) +
  ggtitle("Deserts and Xeric Shrublands") +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank())
E

A + E + D + C + B +
  plot_layout(
  design = "
AB
AC
AE
AF")

ggsave("figures/figure6.png", width = 17, height = 10)
