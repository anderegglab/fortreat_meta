

data <- read.csv("~/Downloads/individual_tree_data.csv")
plot_data <- read.csv("~/Downloads/plot_data.csv")

length(unique(data$plot))

ndata <- data.frame(plot = unique(data$plot))

for (p in unique(data$plot)) {

  subdata <- data[data$plot == p,]
  ndata[ndata$plot == p, "mortality"] <- nrow(subdata[subdata$rdead == 1,]) / nrow(subdata)

}


ndata$burned <- 0
for (i in 1:nrow(ndata)) {
  if (plot_data[plot_data$plot == ndata[i, "plot"], "burn_status"] == "burned") {
    ndata[i, "burned"] <- 1
  }
}


mean_burned <- mean(ndata[ndata$burned == 1, "mortality"])
mean_unburned <- mean(ndata[ndata$burned == 0, "mortality"])

sd_burned <- sd(ndata[ndata$burned == 1, "mortality"])
sd_unburned <- sd(ndata[ndata$burned == 0, "mortality"])

n_burned <- nrow(ndata[ndata$burned == 1,])
n_unburned <- nrow(ndata[ndata$burned == 0,])

sd_burned / sqrt(n_burned)
sd_unburned / sqrt(n_unburned)
