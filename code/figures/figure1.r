##########
# Create a world and western US map as Figure 1. Add them together later in Illustrator
# Cedric Zahnd
##########
# Load packages
library(maps)

# Load data
data <- read.csv("data/processed_data/data_cleaned.csv")

# Subset to keep it tidy
data <- data[,c("studyID", "trt_class", "disturbance_type", "carbon_vs_mortality", "latitude.coordinate", "longitude.coordinate")]

# Colors and symbols
clfire <- "#ff3319" ; cldrou <- "#0057ba" ; clinse <- "#ffab00"
alph <- 50
clfire2 <- rgb(255,51,25,alph,maxColorValue = 255) # fire
cldrou2 <- rgb(0,87,186,alph,maxColorValue = 255) # drought
clinse2 <- rgb(255,171,0,alph,maxColorValue = 255) # insect

data$cl1 <- clfire ; data$cl2 <- clfire2
data$cl1[data$disturbance_type == "drought"] <- cldrou ; data$cl2[data$disturbance_type == "drought"] <- cldrou2  
data$cl1[data$disturbance_type == "insect"] <- clinse ; data$cl2[data$disturbance_type == "insect"] <- clinse2

data$pchs <- 21
data$pchs[data$trt_class == "rx_fire"] <- 24
data$pchs[data$trt_class == "both"] <- 23

#
cx1 <- 1.4

pdf("figures/worldmap.pdf", width = 9.2, height = 4)
map("world", fill = T, col = "lightgrey", border = "white", ylim = c(-60,90), xlim = c(-175,190), myborder = 0, mar = c(0.25,0.1,0.5,0.1),lforce="e")
with(data, points(longitude.coordinate, latitude.coordinate, pch = pchs, col = cl1, bg = cl2, cex = cx1))
legend("bottom", inset = 0.02, bg = "white", col = c("black", "black", clfire, "black", cldrou, "black", clinse, "black"), pch  = c(26,26,21,16,21,24,21,23), pt.bg = c(NA, NA, clfire2, "black", cldrou2, "black", clinse2, "black"),
       legend = c("Disturbance:", "Treatment:", "Fire", "Thinning", "Drought", "Rx Fire", "Insect", "Both"), ncol = 4)
box()
dev.off()

pdf("figures/usamap.pdf", width = 3, height = 3)
map("state", fill = T, col = "lightgrey", border = "white", xlim = c(-130,-100), ylim = c(25,52), myborder = 0, mar = c(0.25,0.1,0.5,0.1),lforce="e")
with(data, points(longitude.coordinate, latitude.coordinate, pch = pchs, col = cl1, bg = cl2, cex = cx1-0.2))
legend("top", legend = "Western USA", bty = "n")
box()
dev.off()


