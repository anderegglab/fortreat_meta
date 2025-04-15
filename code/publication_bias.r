##---------------------------------------------------------------
## Tests for Publication Bias
##
##---------------------------------------------------------------

## libraries
library("metafor")
library("mice")
library("brms")
library("ggplot2")
library("patchwork")
library("cowplot")

source("code/functions.r")

##---------------------------------------------------------------
## Survivorship
##---------------------------------------------------------------

## read in overall mortality (survivorship) model
mort_fit <- readRDS("data/model_objects/mort_fit_overall.rds")

pool_mort <- summary(pool(mort_fit))
pool_mort[-1] <- round(pool_mort[-1], digits = 3)
pool_mort

## assess publication bias
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
  annotate("text", x = 0.2, y = 55, label = "69% < 0.05", size = 12) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 60)) +
  ggtitle("Eggers Test: Survivorship") +
  xlab("P Values") +
  ylab("Count") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), plot.title = element_text(face = "bold", size = 18))

funnel(mort_fit$analyses[[1]])

ggsave("figures/publication_bias_mortality.png", width = 12, height = 8)

##---------------------------------------------------------------
## For Carbon
##---------------------------------------------------------------

carbon_fit <- readRDS("data/model_objects/carbon_fit_overall.rds")

pool_carbon <- summary(pool(carbon_fit))
pool_carbon[-1] <- round(pool_carbon[-1], digits = 3)
pool_carbon

funnel(carbon_fit$analyses[[4]])

## assess publication bias
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
  annotate("text", x = 0.2, y = 9.5, label = "1% < 0.05", size = 12) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 10)) +
  ggtitle("Eggers Test: Carbon") +
  xlab("P Values") +
  ylab("Count") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), plot.title = element_text(face = "bold", size = 18))

ggsave("figures/publication_bias_carbon.png", width = 12, height = 8)
