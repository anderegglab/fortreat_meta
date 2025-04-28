
##---------------------------------------------------------------
## Figure 7 -- climate impacts
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

source("code/00_functions.r")

## read data
data <- read.csv("data/processed_data/data_cleaned.csv")

##---------------------------------------------------------------
## 1. Overall effect
##---------------------------------------------------------------

impute_data <- function(data, vars = c("lrr", "lrr_se", "cwd", "disturbance_type"), m = 20) {

  mice_data <- data[,vars]
  mice_data$cwd <- as.numeric(scale(mice_data$cwd))
  mice_data <- cbind(mice_data, cwd.dist = 0)

  ## make predictor matrix
  predictor_matrix <- make.predictorMatrix(mice_data)
  predictor_matrix ## looks good

  impute_method <- make.method(mice_data)
  impute_method ## no method specified for complete variables

  x <- mice(mice_data, max = 0)
  meth <- x$meth
  meth["cwd.dist"] <- "~I(cwd*disturbance_type)"

  pred <- x$pred
  pred[c("cwd", "disturbance_type"), c("cwd.dist")] <- 0
  imputed_data <- mice(mice_data, meth = meth, pred = pred, maxit = 40, seed = 1, m = 100)

  return(imputed_data)

}

mort_cwd_imputed <- impute_data(data[data$carbon_vs_mortality == 2,], m = 100)

mort_fit <- with(mort_cwd_imputed,
                 rma(yi = lrr,
                     sei = lrr_se,
                     mods = ~ 0 + cwd * disturbance_type))

pool_mort <- summary(pool(mort_fit))
pool_mort[-1] <- round(pool_mort[-1], digits = 4)
pool_mort

table_gen(pool_mort, "cwd_mort.csv")

## mortality
pd <- data.frame(disturbance_typefire = c(rep(1, times = 20), rep(0, times = 40)),
           disturbance_typedrought = c(rep(0, times = 20), rep(1, times = 20), rep(0, times = 20)),
           disturbance_typeinsect = c(rep(0, times = 40), rep(1, times = 20)),
           cwd = c(rep(0, times = 20), seq(-2, 3, length.out = 20), rep(0, times = 20)),
           cwddisturbance_typefire = c(seq(-2, 3, length.out = 20), rep(0, times = 40)),
           cwddisturbance_typeinsect = c(rep(0, times = 40), seq(-2, 3, length.out = 20)))

colnames(pd)[5] <- "cwd:disturbance_typefire"
colnames(pd)[6] <- "cwd:disturbance_typeinsect"

out <- list(mean = matrix(nrow = 60, ncol = 100), ci.lb = matrix(nrow = 60, ncol = 100), ci.ub = matrix(nrow = 60, ncol = 100))

for (i in 1:100) {

  p <- predict.rma(mort_fit$analyses[[i]], newmods = as.matrix(pd))
  out[[1]][,i] <- p[[1]]
  out[[2]][,i] <- p[[3]]
  out[[3]][,i] <- p[[4]]

}

pdata <- data.frame(disturbance_type = factor(rep(c("fire", "drought", "insect"), each = 20), levels = c("fire", "drought", "insect")),
                    cwd = rep(seq(-2, 3, length.out = 20), times = 3))

pdata$mean <- rowMeans(out[[1]])
pdata$lower <- rowMeans(out[[2]])
pdata$upper <- rowMeans(out[[3]])


p <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata, aes(x = cwd, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata, aes(x = cwd, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$carbon_vs_mortality == 2,], aes(x = scale(cwd), y = lrr, color = disturbance_type), size = 5) +
  scale_x_continuous(limits = c(-2, 3), expand = c(0,0)) +
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank())
p


p1 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata[pdata$disturbance_type == "fire",], aes(x = cwd, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata[pdata$disturbance_type == "fire",], aes(x = cwd, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$disturbance_type == "fire" & data$carbon_vs_mortality == 2,], aes(x = scale(cwd), y = lrr, color = disturbance_type), size = 5) +
  scale_x_continuous(limits = c(-2, 3), expand = c(0,0)) +
  scale_y_continuous(limits = c(-1, 5), expand = c(0,0)) +
  ylab("Log Response Ratio") +
  ggtitle("Fire") +
  scale_color_manual(values = c("#e41a1c")) +
  scale_fill_manual(values = c("#e41a1c")) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 18))
p1

p2 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata[pdata$disturbance_type == "drought",], aes(x = cwd, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata[pdata$disturbance_type == "drought",], aes(x = cwd, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$disturbance_type == "drought" & data$carbon_vs_mortality == 2,], aes(x = scale(cwd), y = lrr, color = disturbance_type), size = 5) +
  ggtitle("Drought") +
  xlab("Climatic Water Deficit") +
  scale_x_continuous(limits = c(-2, 3), expand = c(0,0)) +
  scale_y_continuous(limits = c(-1, 5), expand = c(0,0)) +
  scale_color_manual(values = c("#377eb8")) +
  scale_fill_manual(values = c("#377eb8")) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_text(size = 16), axis.title.y = element_blank(),
        plot.title = element_text(face = "bold", size = 18))
p2

p3 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata[pdata$disturbance_type == "insect",], aes(x = cwd, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata[pdata$disturbance_type == "insect",], aes(x = cwd, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$disturbance_type == "insect" & data$carbon_vs_mortality == 2,], aes(x = scale(cwd), y = lrr, color = disturbance_type), size = 5) +
  ggtitle("Insect") +
  scale_x_continuous(limits = c(-2, 3), expand = c(0,0)) +
  scale_y_continuous(limits = c(-1, 5), expand = c(0,0)) +
  scale_color_manual(values = c("#4daf4a")) +
  scale_fill_manual(values = c("#4daf4a")) +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank(),
        plot.title = element_text(face = "bold", size = 18))
p3


p1 + p2 + p3 +
  plot_layout(
  design = "ABC")

ggsave("figures/figure7A.png", width = 20, height = 8)


carbon_sev_imputed <- impute_data(data[data$carbon_vs_mortality == 1,], m = 100)

carbon_fit <- with(carbon_sev_imputed,
                 rma(yi = lrr,
                     sei = lrr_se,
                     mods = ~ 0 + cwd * disturbance_type))

pool_carbon <- summary(pool(carbon_fit))
pool_carbon[-1] <- round(pool_carbon[-1], digits = 4)
pool_carbon

table_gen(pool_carbon, "cwd_carbon.csv")

## mortality
pd <- data.frame(disturbance_typefire = c(rep(1, times = 20), rep(0, times = 40)),
           disturbance_typedrought = c(rep(0, times = 20), rep(1, times = 20), rep(0, times = 20)),
           disturbance_typeinsect = c(rep(0, times = 40), rep(1, times = 20)),
           cwd = c(rep(0, times = 20), seq(-2, 3, length.out = 20), rep(0, times = 20)),
           cwddisturbance_typefire = c(seq(-2, 3, length.out = 20), rep(0, times = 40)),
           cwddisturbance_typeinsect = c(rep(0, times = 40), seq(-2, 3, length.out = 20)))

colnames(pd)[5] <- "cwd:disturbance_typefire"
colnames(pd)[6] <- "cwd:disturbance_typeinsect"

out <- list(mean = matrix(nrow = 60, ncol = 100), ci.lb = matrix(nrow = 60, ncol = 100), ci.ub = matrix(nrow = 60, ncol = 100))

for (i in 1:100) {

  p <- predict.rma(carbon_fit$analyses[[i]], newmods = as.matrix(pd))
  out[[1]][,i] <- p[[1]]
  out[[2]][,i] <- p[[3]]
  out[[3]][,i] <- p[[4]]

}

pdata <- data.frame(disturbance_type = factor(rep(c("fire", "drought", "insect"), each = 20), levels = c("fire", "drought", "insect")),
                    cwd = rep(seq(-2, 3, length.out = 20), times = 3))

pdata$mean <- rowMeans(out[[1]])
pdata$lower <- rowMeans(out[[2]])
pdata$upper <- rowMeans(out[[3]])

p1 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata[pdata$disturbance_type == "fire",], aes(x = cwd, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata[pdata$disturbance_type == "fire",], aes(x = cwd, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$disturbance_type == "fire" & data$carbon_vs_mortality == 1,], aes(x = scale(cwd), y = lrr, color = disturbance_type), size = 5) +
  scale_x_continuous(limits = c(-2, 3), expand = c(0,0)) +
  scale_y_continuous(limits = c(-2,3), expand = c(0,0)) +
  ylab("Log Response Ratio") +
  ggtitle("Fire") +
  scale_color_manual(values = c("#e41a1c")) +
  scale_fill_manual(values = c("#e41a1c")) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 18))
p1

p2 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata[pdata$disturbance_type == "drought",], aes(x = cwd, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata[pdata$disturbance_type == "drought",], aes(x = cwd, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$disturbance_type == "drought" & data$carbon_vs_mortality == 1,], aes(x = scale(cwd), y = lrr, color = disturbance_type), size = 5) +
  ggtitle("Drought") +
  xlab("Climatic Water Deficit") +
  scale_x_continuous(limits = c(-2, 3), expand = c(0,0)) +
  scale_y_continuous(limits = c(-2, 3), expand = c(0,0)) +
  scale_color_manual(values = c("#377eb8")) +
  scale_fill_manual(values = c("#377eb8")) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_text(size = 16), axis.title.y = element_blank(),
        plot.title = element_text(face = "bold", size = 18))
p2

p3 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata[pdata$disturbance_type == "insect",], aes(x = cwd, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata[pdata$disturbance_type == "insect",], aes(x = cwd, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$disturbance_type == "insect" & data$carbon_vs_mortality == 1,], aes(x = scale(cwd), y = lrr, color = disturbance_type), size = 5) +
  ggtitle("Insect") +
  scale_x_continuous(limits = c(-2, 3), expand = c(0,0)) +
  scale_y_continuous(limits = c(-2, 3), expand = c(0,0)) +
  scale_color_manual(values = c("#4daf4a")) +
  scale_fill_manual(values = c("#4daf4a")) +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank(),
        plot.title = element_text(face = "bold", size = 18))
p3


p1 + p2 + p3 +
  plot_layout(
  design = "ABC")

ggsave("figures/figure7B.png", width = 20, height = 8)

##---------------------------------------------------------------
## TEMP
##---------------------------------------------------------------

impute_data <- function(data, vars = c("lrr", "lrr_se", "mat", "disturbance_type"), m = 20) {

  mice_data <- data[,vars]
  mice_data$mat <- as.numeric(scale(mice_data$mat))
  mice_data <- cbind(mice_data, mat.dist = 0)

  ## make predictor matrix
  predictor_matrix <- make.predictorMatrix(mice_data)
  predictor_matrix ## looks good

  impute_method <- make.method(mice_data)
  impute_method ## no method specified for complete variables

  x <- mice(mice_data, max = 0)
  meth <- x$meth
  meth["mat.dist"] <- "~I(mat*disturbance_type)"

  pred <- x$pred
  pred[c("mat", "disturbance_type"), c("mat.dist")] <- 0
  imputed_data <- mice(mice_data, meth = meth, pred = pred, maxit = 40, seed = 1, m = 100)

  return(imputed_data)

}

mort_mat_imputed <- impute_data(data[data$carbon_vs_mortality == 2,], m = 100)

mort_fit <- with(mort_mat_imputed,
                 rma(yi = lrr,
                     sei = lrr_se,
                     mods = ~ 0 + mat * disturbance_type))

pool_mort <- summary(pool(mort_fit))
pool_mort[-1] <- round(pool_mort[-1], digits = 4)
pool_mort

table_gen(pool_mort, "mat_mort.csv")

## mortality
pd <- data.frame(disturbance_typefire = c(rep(1, times = 20), rep(0, times = 40)),
           disturbance_typedrought = c(rep(0, times = 20), rep(1, times = 20), rep(0, times = 20)),
           disturbance_typeinsect = c(rep(0, times = 40), rep(1, times = 20)),
           mat = c(rep(0, times = 20), seq(-2, 3, length.out = 20), rep(0, times = 20)),
           matdisturbance_typefire = c(seq(-2, 3, length.out = 20), rep(0, times = 40)),
           matdisturbance_typeinsect = c(rep(0, times = 40), seq(-2, 3, length.out = 20)))

colnames(pd)[5] <- "mat:disturbance_typefire"
colnames(pd)[6] <- "mat:disturbance_typeinsect"

out <- list(mean = matrix(nrow = 60, ncol = 100), ci.lb = matrix(nrow = 60, ncol = 100), ci.ub = matrix(nrow = 60, ncol = 100))

for (i in 1:100) {

  p <- predict.rma(mort_fit$analyses[[i]], newmods = as.matrix(pd))
  out[[1]][,i] <- p[[1]]
  out[[2]][,i] <- p[[3]]
  out[[3]][,i] <- p[[4]]

}

pdata <- data.frame(disturbance_type = factor(rep(c("fire", "drought", "insect"), each = 20), levels = c("fire", "drought", "insect")),
                    mat = rep(seq(-2, 3, length.out = 20), times = 3))

pdata$mean <- rowMeans(out[[1]])
pdata$lower <- rowMeans(out[[2]])
pdata$upper <- rowMeans(out[[3]])


p <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata, aes(x = mat, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata, aes(x = mat, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$carbon_vs_mortality == 2,], aes(x = scale(mat), y = lrr, color = disturbance_type), size = 5) +
  scale_x_continuous(limits = c(-2, 3), expand = c(0,0)) +
  scale_y_continuous(limits = c(-1, 5), expand = c(0,0)) +
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank())
p

p1 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata[pdata$disturbance_type == "fire",], aes(x = mat, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata[pdata$disturbance_type == "fire",], aes(x = mat, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$disturbance_type == "fire" & data$carbon_vs_mortality == 2,], aes(x = scale(mat), y = lrr, color = disturbance_type), size = 5) +
  scale_x_continuous(limits = c(-2, 3), expand = c(0,0)) +
  scale_y_continuous(limits = c(-1,5), expand = c(0,0)) +
  ylab("Log Response Ratio") +
  ggtitle("Fire") +
  scale_color_manual(values = c("#e41a1c")) +
  scale_fill_manual(values = c("#e41a1c")) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 18))
p1

p2 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata[pdata$disturbance_type == "drought",], aes(x = mat, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata[pdata$disturbance_type == "drought",], aes(x = mat, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$disturbance_type == "drought" & data$carbon_vs_mortality == 2,], aes(x = scale(mat), y = lrr, color = disturbance_type), size = 5) +
  ggtitle("Drought") +
  xlab("Mean Annual Temperature") +
  scale_x_continuous(limits = c(-2, 3), expand = c(0,0)) +
  scale_y_continuous(limits = c(-1, 5), expand = c(0,0)) +
  scale_color_manual(values = c("#377eb8")) +
  scale_fill_manual(values = c("#377eb8")) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_text(size = 16), axis.title.y = element_blank(),
        plot.title = element_text(face = "bold", size = 18))
p2

p3 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata[pdata$disturbance_type == "insect",], aes(x = mat, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata[pdata$disturbance_type == "insect",], aes(x = mat, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$disturbance_type == "insect" & data$carbon_vs_mortality == 2,], aes(x = scale(mat), y = lrr, color = disturbance_type), size = 5) +
  ggtitle("Insect") +
  scale_x_continuous(limits = c(-2, 3), expand = c(0,0)) +
  scale_y_continuous(limits = c(-1, 5), expand = c(0,0)) +
  scale_color_manual(values = c("#4daf4a")) +
  scale_fill_manual(values = c("#4daf4a")) +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank(),
        plot.title = element_text(face = "bold", size = 18))
p3


p1 + p2 + p3 +
  plot_layout(
  design = "ABC")

ggsave("figures/figure7C.png", width = 20, height = 8)


carbon_sev_imputed <- impute_data(data[data$carbon_vs_mortality == 1,], m = 100)

carbon_fit <- with(carbon_sev_imputed,
                 rma(yi = lrr,
                     sei = lrr_se,
                     mods = ~ 0 + mat * disturbance_type))

pool_carbon <- summary(pool(carbon_fit))
pool_carbon[-1] <- round(pool_carbon[-1], digits = 4)
pool_carbon

table_gen(pool_carbon, "mat_carbon.csv")

## mortality
pd <- data.frame(disturbance_typefire = c(rep(1, times = 20), rep(0, times = 40)),
           disturbance_typedrought = c(rep(0, times = 20), rep(1, times = 20), rep(0, times = 20)),
           disturbance_typeinsect = c(rep(0, times = 40), rep(1, times = 20)),
           mat = c(rep(0, times = 20), seq(-2, 3, length.out = 20), rep(0, times = 20)),
           matdisturbance_typefire = c(seq(-2, 3, length.out = 20), rep(0, times = 40)),
           matdisturbance_typeinsect = c(rep(0, times = 40), seq(-2, 3, length.out = 20)))

colnames(pd)[5] <- "mat:disturbance_typefire"
colnames(pd)[6] <- "mat:disturbance_typeinsect"

out <- list(mean = matrix(nrow = 60, ncol = 100), ci.lb = matrix(nrow = 60, ncol = 100), ci.ub = matrix(nrow = 60, ncol = 100))

for (i in 1:100) {

  p <- predict.rma(carbon_fit$analyses[[i]], newmods = as.matrix(pd))
  out[[1]][,i] <- p[[1]]
  out[[2]][,i] <- p[[3]]
  out[[3]][,i] <- p[[4]]

}

pdata <- data.frame(disturbance_type = factor(rep(c("fire", "drought", "insect"), each = 20), levels = c("fire", "drought", "insect")),
                    mat = rep(seq(-2, 3, length.out = 20), times = 3))

pdata$mean <- rowMeans(out[[1]])
pdata$lower <- rowMeans(out[[2]])
pdata$upper <- rowMeans(out[[3]])

p1 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata[pdata$disturbance_type == "fire",], aes(x = mat, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata[pdata$disturbance_type == "fire",], aes(x = mat, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$disturbance_type == "fire" & data$carbon_vs_mortality == 1,], aes(x = scale(mat), y = lrr, color = disturbance_type), size = 5) +
  scale_x_continuous(limits = c(-2, 3), expand = c(0,0)) +
  scale_y_continuous(limits = c(-5,3), expand = c(0,0)) +
  annotate("text", x = -1, y = 2.5, label = "***", size = 16) +
  ylab("Log Response Ratio") +
  ggtitle("Fire") +
  scale_color_manual(values = c("#e41a1c")) +
  scale_fill_manual(values = c("#e41a1c")) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 18))
p1

p2 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata[pdata$disturbance_type == "drought",], aes(x = mat, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata[pdata$disturbance_type == "drought",], aes(x = mat, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$disturbance_type == "drought" & data$carbon_vs_mortality == 1,], aes(x = scale(mat), y = lrr, color = disturbance_type), size = 5) +
  ggtitle("Drought") +
  xlab("Mean Annual Temperature") +
  scale_x_continuous(limits = c(-2, 3), expand = c(0,0)) +
  scale_y_continuous(limits = c(-5, 3), expand = c(0,0)) +
  scale_color_manual(values = c("#377eb8")) +
  scale_fill_manual(values = c("#377eb8")) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_text(size = 16), axis.title.y = element_blank(),
        plot.title = element_text(face = "bold", size = 18))
p2

p3 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 2) +
  geom_ribbon(data = pdata[pdata$disturbance_type == "insect",], aes(x = mat, ymin = lower, ymax = upper, fill = disturbance_type), alpha = 0.3) +
  geom_line(data = pdata[pdata$disturbance_type == "insect",], aes(x = mat, y = mean, color = disturbance_type), size = 3) +
  geom_point(data = data[data$disturbance_type == "insect" & data$carbon_vs_mortality == 1,], aes(x = scale(mat), y = lrr, color = disturbance_type), size = 5) +
  ggtitle("Insect") +
  scale_x_continuous(limits = c(-2, 3), expand = c(0,0)) +
  scale_y_continuous(limits = c(-5, 3), expand = c(0,0)) +
  scale_color_manual(values = c("#4daf4a")) +
  scale_fill_manual(values = c("#4daf4a")) +
  theme_bw() +
  theme(legend.position = "none", axis.title = element_blank(),
        plot.title = element_text(face = "bold", size = 18))
p3

p1 + p2 + p3 +
  plot_layout(
  design = "ABC")

ggsave("figures/figure7D.png", width = 20, height = 8)
