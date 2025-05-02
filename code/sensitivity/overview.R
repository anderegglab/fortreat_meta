
data <- read.csv("data/processed_data/data_cleaned.csv")
head(data)

n2 <- data[data$Nis2 == 1,]
noc <- data[data$NoTrueControl == 1,]
mor <- data[data$mortality_attribution_uncertain == 1,]
zer <- data[data$zeroSurv == 1,]

n2[,c("studyID", "disturbance_type", "trt_class", "carbon_vs_mortality")]
unique(n2$studyID)
### N = 2:
# carbon: fire & drought
# Surviv: drought

# carbon - fire - rx_fire (1 / 1)
# carbon - drought - thinning (2 / 7)
# surviv - drought - thinning (1 / 4)

noc[,c("studyID", "disturbance_type", "trt_class", "carbon_vs_mortality")]
unique(noc$studyID)
### No true control:
# carbon: drought
# surviv: drought & insects

# carbon - drought - thinning (1 / 4)
# surviv - drought - thinning (1 / 4)
# surviv - insects - thinning (1 / 14)

mor[,c("studyID", "disturbance_type", "trt_class", "carbon_vs_mortality")]
unique(mor$studyID)
### Mortality uncertain:
# carbon: drought
# surviv: drought & insects

# carbon - drought - thinning (1 / 3)
# surviv - drought - thinning (1 / 3)
# surviv - insects - rx_fire (3 / 9)

zer[,c("studyID", "disturbance_type", "trt_class", "carbon_vs_mortality")]
unique(zer$studyID)
### Zero mortality:
# surviv: drought & insects

# surviv - fire - thinning (2 / 6)
# surviv - fire - both (2 / 2)


# Total:
# carbon: fire & drought
# surviv: drought & insects & fire

# carbon - fire - rx_fire
# carbon - drought - thinning
# surviv - drought - thinning
# surviv - insects - thinning
# surviv - insects - rx_fire
# surviv - fire - thinning
# surviv - fire - both

