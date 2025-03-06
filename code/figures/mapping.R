##########
# Creates 3 figures:
# - World map with each study site by response variable, treatment and disturbance (map_world.pdf)
# - USA map with each study site by response variable, treatment and disturbance (map_us.pdf)
# - Legend for these two maps
# I suggest to combine the maps and legend in illustrator once we're ready to submit.
#
# Notes:
# Categorisation of treatments is copied from Jacobs scripts
# Here, I still add coordinates "manually" where missing, but these have been updated in the analysis sheets...
# Maps are a bit a mess right now, evtl worth recoding in ggplot or so.
# 
# TODO: 
# -Integrate categorisation into the data cleaning script(s)
# -Remove manual change of coordinates
# -Unify disturbance naming (there are lots of names for wildfire right now)
##########

# Load packages
library(maps)

# Load data
data <- read.csv("data/processed_data/data_cleaned.csv")

# Subset to keep it tidy
data <- data[,c("studyID", "treatment", "treatment.notes", "disturbance", "carbon_vs_mortality", "SiteName", "latitude.coordinate", "longitude.coordinate")]

# Categorize treatments:
cat_trt <- function(trt) {
  if (((grepl("burn", tolower(trt)) | grepl("fire", tolower(trt))) & !grepl("unburn", tolower(trt))) &
      (grepl("thin", tolower(trt)) | grepl("density", tolower(trt)) | grepl("umz", tolower(trt)))) return("both")
  else if((grepl("burn", tolower(trt)) | grepl("fire", tolower(trt))) & !grepl("unburn", tolower(trt))) return("rx_fire")
  else if(grepl("thin", tolower(trt)) | grepl("density", tolower(trt)) | grepl("umz", tolower(trt)))
    return("thinning")
  else return(NA)
}

data$burn <- "no"
data$thin <- "no"
for (i in 1:nrow(data)) {
  print(i)
  data[i,"trt_class"] <- cat_trt(data[i,"treatment"])
  if (!is.na(data[i, "trt_class"])) {
    if (data[i,"trt_class"] == "thinning") {
      data[i,"thin"] <- "yes"
    } else if (data[i,"trt_class"] == "rx_fire") {
      data[i,"burn"] <- "yes"
    } else if (data[i,"trt_class"] == "both") {
      data[i,"thin"] <- "yes"
      data[i,"burn"] <- "yes"
    }
  }
}
data$trt_class[is.na(data$trt_class)] <- "thinning"
data$trt_class ## looks correct

# Add coordinates where missing:
data$latitude.coordinate[data$studyID == 33] <- 41.449763 ; data$longitude.coordinate[data$studyID == 33] <- -120.251172 
data$latitude.coordinate[data$studyID == 59] <- 43.952778 ; data$longitude.coordinate[data$studyID == 59] <- -103.681242 
data$latitude.coordinate[data$studyID == 181] <- 47.7 ; data$longitude.coordinate[data$studyID == 181] <- -91.935
data$latitude.coordinate <- as.numeric(data$latitude.coordinate)

# Name all non-insect or drought disturbances "wildfire":
data$disturbance[!data$disturbance %in% c("insects","drought&insect","drought")] <- "wildfire"

# Add color and point symbols depending on treatment, disturbance and response type
data$pchs <- NA ; data$clrs <- NA

alph <- 100
cl1 <- rgb(255,165,0,alph,maxColorValue = 255) # insects
cl2 <- rgb(160,32,240,alph,maxColorValue = 255) # drought and insects
cl3 <- rgb(0,0,255,alph,maxColorValue = 255) # drought
cl4 <- rgb(255,0,0,alph,maxColorValue = 255) # wildfire

for(i in 1:nrow(data)){
  if(data$disturbance[i] == "insects"){
    data$clrs[i] <- cl1
  }else if(data$disturbance[i] == "drought&insect"){
    data$clrs[i] <- cl2
  }else if(data$disturbance[i] == "drought"){
    data$clrs[i] <- cl3
  }else if(data$disturbance[i] == "wildfire"){
    data$clrs[i] <- cl4
  }
  
  if(data$carbon_vs_mortality[i] == 1){
    if(data$trt_class[i] == "thinning"){
      data$pchs[i] <- 15
    }else if(data$trt_class[i] == "rx_fire"){
      data$pchs[i] <- 17
    }else if(data$trt_class[i] == "both"){
      data$pchs[i] <- 19
    }
  }else if(data$carbon_vs_mortality[i] == 2){
    if(data$trt_class[i] == "thinning"){
      data$pchs[i] <- 0
    }else if(data$trt_class[i] == "rx_fire"){
      data$pchs[i] <- 2
    }else if(data$trt_class[i] == "both"){
      data$pchs[i] <- 1
    }
  }
}

# ---------------------------------------
# Make the maps:
# CAREFUL: This will overwrite previous pdf versions!

# World map:
pdf("figures/map_word.pdf", 12, 6)
map("world", fill = T, col = "lightgrey")
with(data, points(longitude.coordinate, latitude.coordinate, pch = pchs, col = clrs, cex = 1))
dev.off()

# US map:
pdf("figures/map_us.pdf", 8, 6)
map("world", fill = T, col = "lightgrey", xlim = c(-130,-60), ylim = c(22,55))
with(data, points(longitude.coordinate, latitude.coordinate, pch = pchs, col = clrs, cex = 1))
dev.off()

# Legend
pdf("figures/map_legend.pdf", 5, 5)
plot(0:1,0:1,type = "n", xaxt = "n", yaxt = "n", xlab = NA, ylab = NA, frame.plot = F)
legend("center", legend = c("Thinning", "Fire", "Both", "", "Wildfire", "Drought", "Insects", "Insects & Drought", "", "Carbon", "Mortality"),
       col = c("black","black","black","black","red","blue","orange","purple","black","black", "black"),
       pch = c(15,17,19,NA,18,18,18,18,NA,18,5))
dev.off()


