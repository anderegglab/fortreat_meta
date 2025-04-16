
##---------------------------------------------------------------
## Model Fitting for Meta Analysis
##
##---------------------------------------------------------------

##---------------------------------------------------------------
## 0. Head
##---------------------------------------------------------------

## libraries
library("metafor")
library("mice")
library("brms")
library("ggplot2")
library("broom")
library("patchwork")

source("code/00_functions.r")

## read data
data <- read.csv("data/processed_data/data_cleaned.csv")

##---------------------------------------------------------------
## 1. Overall effect
##---------------------------------------------------------------

data$ba_removed <- as.numeric(gsub("< ", "", data$ba_removed))

impute_data <- function(data, vars = c("lrr", "lrr_se", "ba_removed", "disturbance_type"), m = 20) {

  #mice_data <- data[data$disturbance_type == "fire" & data$carbon_vs_mortality == 2,][, vars]
  mice_data <- data[,vars]
  #mice_data[mice_data$burn_bin == 1, "ba_removed"] <- 0
  mice_data$ba_removed <- as.numeric(scale(mice_data$ba_removed))
  mice_data <- cbind(mice_data, ba.dist = 0)

  ## make predictor matrix
  predictor_matrix <- make.predictorMatrix(mice_data)
  predictor_matrix ## looks good

  impute_method <- make.method(mice_data)
  impute_method ## no method specified for complete variables

  x <- mice(mice_data, max = 0)
  meth <- x$meth
  meth["ba.dist"] <- "~I(ba_removed*disturbance_type)"

  pred <- x$pred
  pred[c("ba_removed", "disturbance_type"), c("ba.dist")] <- 0
  imputed_data <- mice(mice_data, meth = meth, pred = pred, maxit = 40, seed = 1, m = 100)

  return(imputed_data)

}

#sub <- data[data$carbon_vs_mortality == 2,]
#sub$ba_removed_scaled <- scale(sub$ba_removed)

#sub[sub$disturbance_type == "fire", "ba_removed_scaled"][order(sub[sub$disturbance_type == "fire", "ba_removed_scaled"])]

mort_sev_imputed <- impute_data(data[data$carbon_vs_mortality == 2,], m = 100)

#mort_sev_imputed <- impute_data(sub[sub$carbon_vs_mortality == 2 & sub$ba_removed_scaled != -1.1546552 & sub$ba_removed_scaled != -0.9037813,], m = 100)
mort_fit <- with(mort_sev_imputed,
                 rma(yi = lrr,
                     sei = lrr_se,
                     mods = ~ 0 + ba_removed * disturbance_type))

pool_mort <- summary(pool(mort_fit))
pool_mort[-1] <- round(pool_mort[-1], digits = 4)
pool_mort

table_gen(pool_mort, "baremoved_burn.csv")

## mortality
pd <- data.frame(disturbance_typefire = c(rep(1, times = 20), rep(0, times = 40)),
           disturbance_typedrought = c(rep(0, times = 20), rep(1, times = 20), rep(0, times = 20)),
           disturbance_typeinsect = c(rep(0, times = 40), rep(1, times = 20)),
           ba_removed = c(rep(0, times = 20), seq(-2, 2.5, length.out = 20), rep(0, times = 20)),
           ba_removeddisturbance_typefire = c(seq(-2, 2.5, length.out = 20), rep(0, times = 40)),
           ba_removeddisturbance_typeinsect = c(rep(0, times = 40), seq(-2, 2.5, length.out = 20)))

colnames(pd)[5] <- "ba_removed:disturbance_typefire"
colnames(pd)[6] <- "ba_removed:disturbance_typeinsect"

out <- list(mean = matrix(nrow = 60, ncol = 100), ci.lb = matrix(nrow = 60, ncol = 100), ci.ub = matrix(nrow = 60, ncol = 100))

for (i in 1:100) {

  p <- predict.rma(mort_fit$analyses[[i]], newmods = as.matrix(pd))
  out[[1]][,i] <- p[[1]]
  out[[2]][,i] <- p[[3]]
  out[[3]][,i] <- p[[4]]

}

pdata <- data.frame(disturbance_type = factor(rep(c("fire", "drought", "insect"), each = 20), levels = c("fire", "drought", "insect")),
                    ba_removed = rep(seq(-2, 2.5, length.out = 20), times = 3))

pdata$mean <- rowMeans(out[[1]])
pdata$lower <- rowMeans(out[[2]])
pdata$upper <- rowMeans(out[[3]])


p <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata, aes(x = ba_removed, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata, aes(x = ba_removed, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$carbon_vs_mortality == 2,], aes(x = scale(ba_removed), y = lrr, color = disturbance_type), size = 5) +
  scale_x_continuous(limits = c(-2, 2.5), expand = c(0,0)) +
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank())
p


p1 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata[pdata$disturbance_type == "fire",], aes(x = ba_removed, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata[pdata$disturbance_type == "fire",], aes(x = ba_removed, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$disturbance_type == "fire" & data$carbon_vs_mortality == 2,], aes(x = scale(ba_removed), y = lrr, color = disturbance_type), size = 5) +
  scale_x_continuous(limits = c(-2, 2.5), expand = c(0,0)) +
  scale_y_continuous(limits = c(-1, 5), expand = c(0,0)) +
  scale_color_manual(values = c("#e41a1c")) +
  scale_fill_manual(values = c("#e41a1c")) +
  annotate("text", y = 4.5, x = -1.5, label = "***", size = 15) +
  ylab("Log Response Ratio") +
  ggtitle("Fire") +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 18))
p1

p2 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata[pdata$disturbance_type == "drought",], aes(x = ba_removed, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata[pdata$disturbance_type == "drought",], aes(x = ba_removed, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$disturbance_type == "drought" & data$carbon_vs_mortality == 2,], aes(x = scale(ba_removed), y = lrr, color = disturbance_type), size = 5) +
  scale_x_continuous(limits = c(-2, 2.5), expand = c(0,0)) +
  scale_y_continuous(limits = c(-1, 5), expand = c(0,0)) +
  scale_color_manual(values = c("#377eb8")) +
  scale_fill_manual(values = c("#377eb8")) +
  xlab("BA Removed (scaled)") +
  ggtitle("Drought") +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_text(size = 16), axis.title.y = element_blank(),
        plot.title = element_text(face = "bold", size = 18))
p2

p3 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata[pdata$disturbance_type == "insect",], aes(x = ba_removed, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata[pdata$disturbance_type == "insect",], aes(x = ba_removed, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$disturbance_type == "insect" & data$carbon_vs_mortality == 2,], aes(x = scale(ba_removed), y = lrr, color = disturbance_type), size = 5) +
  scale_x_continuous(limits = c(-2, 2.5), expand = c(0,0)) +
  scale_y_continuous(limits = c(-1, 5), expand = c(0,0)) +
  scale_color_manual(values = c("#4daf4a")) +
  scale_fill_manual(values = c("#4daf4a")) +
  ggtitle("Insect") +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank(),
        plot.title = element_text(face = "bold", size = 18))
p3

mort_plot <- p1 + p2 + p3 + plot_layout(
               design = "
ABC
")
mort_plot

ggsave("figures/figure4A.png", width = 20, height = 8)



carbon_sev_imputed <- impute_data(data[data$carbon_vs_mortality == 1,], m = 100)

carbon_fit <- with(carbon_sev_imputed,
                 rma(yi = lrr,
                     sei = lrr_se,
                     mods = ~ 0 + ba_removed * disturbance_type))

pool_carbon <- summary(pool(carbon_fit))
pool_carbon[-1] <- round(pool_carbon[-1], digits = 4)
pool_carbon

table_gen(pool_carbon, "baremoved_carbon.csv")

## mortality
pd <- data.frame(disturbance_typefire = c(rep(1, times = 20), rep(0, times = 40)),
           disturbance_typedrought = c(rep(0, times = 20), rep(1, times = 20), rep(0, times = 20)),
           disturbance_typeinsect = c(rep(0, times = 40), rep(1, times = 20)),
           ba_removed = c(rep(0, times = 20), seq(-2, 2.5, length.out = 20), rep(0, times = 20)),
           ba_removeddisturbance_typefire = c(seq(-2, 2.5, length.out = 20), rep(0, times = 40)),
           ba_removeddisturbance_typeinsect = c(rep(0, times = 40), seq(-2, 2.5, length.out = 20)))

colnames(pd)[5] <- "ba_removed:disturbance_typefire"
colnames(pd)[6] <- "ba_removed:disturbance_typeinsect"

out <- list(mean = matrix(nrow = 60, ncol = 100), ci.lb = matrix(nrow = 60, ncol = 100), ci.ub = matrix(nrow = 60, ncol = 100))

for (i in 1:100) {

  p <- predict.rma(carbon_fit$analyses[[i]], newmods = as.matrix(pd))
  out[[1]][,i] <- p[[1]]
  out[[2]][,i] <- p[[3]]
  out[[3]][,i] <- p[[4]]

}

pdata <- data.frame(disturbance_type = factor(rep(c("fire", "drought", "insect"), each = 20), levels = c("fire", "drought", "insect")),
                    ba_removed = rep(seq(-2, 2.5, length.out = 20), times = 3))

pdata$mean <- rowMeans(out[[1]])
pdata$lower <- rowMeans(out[[2]])
pdata$upper <- rowMeans(out[[3]])

p1 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata[pdata$disturbance_type == "fire",], aes(x = ba_removed, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata[pdata$disturbance_type == "fire",], aes(x = ba_removed, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$disturbance_type == "fire" & data$carbon_vs_mortality == 1,], aes(x = scale(ba_removed), y = lrr, color = disturbance_type), size = 5) +
  scale_x_continuous(limits = c(-2, 2.5), expand = c(0,0)) +
  scale_y_continuous(limits = c(-2, 5), expand = c(0,0)) +
  scale_color_manual(values = c("#e41a1c")) +
  scale_fill_manual(values = c("#e41a1c")) +
  ylab("Log Response Ratio") +
  ggtitle("Fire") +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 18))
p1

p2 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata[pdata$disturbance_type == "drought",], aes(x = ba_removed, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata[pdata$disturbance_type == "drought",], aes(x = ba_removed, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$disturbance_type == "drought" & data$carbon_vs_mortality == 1,], aes(x = scale(ba_removed), y = lrr, color = disturbance_type), size = 5) +
  scale_x_continuous(limits = c(-2, 2.5), expand = c(0,0)) +
  scale_y_continuous(limits = c(-2, 5), expand = c(0,0)) +
  scale_color_manual(values = c("#377eb8")) +
  scale_fill_manual(values = c("#377eb8")) +
  xlab("BA Removed (scaled)") +
  ggtitle("Drought") +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_text(size = 16), axis.title.y = element_blank(),
        plot.title = element_text(face = "bold", size = 18))
p2

p3 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata[pdata$disturbance_type == "insect",], aes(x = ba_removed, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata[pdata$disturbance_type == "insect",], aes(x = ba_removed, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$disturbance_type == "insect" & data$carbon_vs_mortality == 1,], aes(x = scale(ba_removed), y = lrr, color = disturbance_type), size = 5) +
  scale_x_continuous(limits = c(-2, 2.5), expand = c(0,0)) +
  scale_y_continuous(limits = c(-2, 5), expand = c(0,0)) +
  scale_color_manual(values = c("#4daf4a")) +
  scale_fill_manual(values = c("#4daf4a")) +
  ggtitle("Insect") +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank(),
        plot.title = element_text(face = "bold", size = 18))
p3

carbon_plot <- p1 + p2 + p3 + plot_layout(
               design = "
ABC
")
carbon_plot

ggsave("figures/figure4B.png", width = 20, height = 8)

