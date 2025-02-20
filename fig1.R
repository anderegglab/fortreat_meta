library(tidyverse)

palette(c("black", "red"))

data <- read.csv("/Users/acpostleon/Desktop/Lab/GroupMetaAnalysis/Data/data_cleaned.csv")



data <- data %>%
  mutate(treatment = str_to_lower(treatment)) %>%
  mutate(treatment.notes = str_to_lower(treatment.notes)) %>%
  mutate(burn = case_when(str_detect(treatment, "burn")== TRUE | str_detect(treatment, "fire")== TRUE | str_detect(treatment.notes, "fire")==TRUE | str_detect(treatment.notes, "burn")==TRUE  ~ "burned",
                str_detect(treatment, "burn")== FALSE & str_detect(treatment, "fire")== FALSE & str_detect(treatment.notes, "fire")==FALSE & str_detect(treatment.notes, "burn")==FALSE  ~ "unburned")) %>%
  mutate(thin = case_when(str_detect(treatment, "thin")== TRUE | str_detect(treatment, "harvest")== TRUE | str_detect(treatment.notes, "thin")==TRUE | str_detect(treatment.notes, "harvest")==TRUE  ~ "thin",
                          str_detect(treatment, "thin")== FALSE & str_detect(treatment, "harvest")== FALSE & str_detect(treatment.notes, "thin")==FALSE & str_detect(treatment.notes, "harvest")==FALSE  ~ "not thinned")) %>%
  mutate(treatment_type = case_when(burn=="unburned" & thin =="thin" ~ "thin",
                                    burn=="burned" & thin =="not thinned" ~ "burn",
                                    burn=="burned" & thin == "thin" ~ "thin + burn")) %>%
  mutate(treatment_number = case_when(treatment_type == "thin" ~ 3,
                                      treatment_type == "burn" ~ 2,
                                      treatment_type == "thin + burn" ~ 1)) %>%
  mutate(treat_carbon_number = case_when(carbon_vs_mortality == 1 ~ treatment_number+ 0.1,
                                         carbon_vs_mortality == 2 ~ treatment_number - 0.1))


png("/Users/acpostleon/Desktop/Lab/GroupMetaAnalysis/Graphs/disturbance_response.png", width=12, height=6, units="in", res=300)
par(mfrow=c(1,3))
for (dist in unique(data$disturbance_type)){
  data2 <- data %>%
    filter(disturbance_type == dist)
  plot(data2$lrr, data2$treat_carbon_number, col=alpha(data2$carbon_vs_mortality, 0.6), pch=19,cex=1.3,
       ylim=c(0.5, 3.2), yaxt="n", ylab="", xlab="Change in disturbance response\n(log response ratio)", main=dist)
  points()
  axis(side=2, at=c(1,2,3), labels =c("thin +\nburn  ", "burn", "thin"), las=2)
  legend("bottomleft", col = alpha(c(1,2),0.6), legend = c("carbon", "mortality"), pch=19, cex=1.1)
  
}
dev.off()



