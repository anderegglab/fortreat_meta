##---------------------------------------------------------------
## Figure 1
##
##---------------------------------------------------------------

## load packages
library("metafor")
library("ggplot2")
library("maps")
library("sf")
library("MetBrewer")
library("patchwork")

## read in data
data <- read.csv("data/processed_data/data_cleaned.csv")

sf_data <- st_as_sf(data, coords = c("longitude.coordinate", "latitude.coordinate"), crs = 4326)
sf_data$disturbance_type <- factor(data$disturbance_type, levels = c("fire", "drought", "insect"))

world <- st_read("../../Princeton/precipitation_intensification/data/World_Countries_Generalized.shp")
world <- st_transform(world, crs = 4326)

p1 <- ggplot(world) +
  geom_sf(fill = "gray", color = "white") +
  geom_sf(data = sf_data, aes(color = disturbance_type), size = 5, alpha = 0.9) +
  theme_bw() +
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  labs(color = "Disturbance") +
  #scale_color_gradient2(low = "#f7fbff", mid = "#6baed6", high = "#084594", midpoint = 800) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Disturbance") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1.1, "cm"),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "italic", size = 12),
        legend.position = "inside",
        legend.position.inside = c(0.05, 0.5))
p1
ggsave("figures/figure1A.png", width = 14, height = 10)

p2 <- ggplot(world) +
  geom_sf(fill = "gray", color = "white") +
  geom_sf(data = sf_data, aes(color = Olson_Biome), size = 5, alpha = 0.9) +
  theme_bw() +
  #scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  labs(color = "Olson Biome") +
  scale_color_manual(values=met.brewer("Lakota", 4)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Biome") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1.1, "cm"),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "italic", size = 12),
        legend.position = "inside",
        legend.position.inside = c(0.15, 0.3))
p2
ggsave("figures/figure1B.png", width = 14, height = 10)


p3 <- ggplot(world) +
  geom_sf(fill = "gray", color = "white") +
  geom_sf(data = sf_data, aes(color = trt_class), size = 5, alpha = 0.9) +
  theme_bw() +
  #scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  scale_color_manual(values=met.brewer("Klimt", 3)) +
  labs(color = "Olson Biome") +
  #scale_color_gradient2(low = "#f7fbff", mid = "#6baed6", high = "#084594", midpoint = 800) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Treatment") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1.1, "cm"),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "italic", size = 12),
        legend.position = "inside",
        legend.position.inside = c(0.15, 0.3))
p3
ggsave("figures/figure1C.png", width = 14, height = 10)


library(usmap)

states <- us_map()
states <- st_transform(states, crs = 4326)
states <- states[states$abbr != "AK" & states$abbr != "HI",]


p4 <- ggplot(states) +
  geom_sf(fill = "gray", color = "white") +
  geom_sf(data = sf_data, aes(color = disturbance_type), size = 5, alpha = 0.9) +
  theme_bw() +
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  labs(color = "Disturbance") +
  xlim(c(-130,-60)) +
  ylim(c(25,50)) +
  #scale_color_gradient2(low = "#f7fbff", mid = "#6baed6", high = "#084594", midpoint = 800) +
  #scale_x_continuous(expand = c(0,0)) +
  #scale_y_continuous(expand = c(0,0)) +
  ggtitle("Disturbance (US)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1.1, "cm"),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "italic", size = 12),
        legend.position = "inside",
        legend.position.inside = c(0.05, 0.5))
p4
ggsave("figures/figure1D.png", width = 14, height = 10)

p5 <- ggplot(states) +
  geom_sf(fill = "gray", color = "white") +
  geom_sf(data = sf_data, aes(color = Olson_Biome), size = 5, alpha = 0.9) +
  theme_bw() +
  xlim(c(-130,-60)) +
  ylim(c(25,50)) +
  #scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  labs(color = "Olson Biome") +
  scale_color_manual(values=met.brewer("Lakota", 4)) +
  #scale_x_continuous(expand = c(0,0)) +
  #scale_y_continuous(expand = c(0,0)) +
  ggtitle("Biome (US)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 20),
        legend.key.size = unit(0.9, "cm"),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "italic", size = 12),
        legend.position = "inside",
        legend.position.inside = c(0.15, 0.17))
p5
ggsave("figures/figure1E.png", width = 14, height = 10)

p6 <- ggplot(states) +
  geom_sf(fill = "gray", color = "white") +
  geom_sf(data = sf_data, aes(color = trt_class), size = 5, alpha = 0.9) +
  theme_bw() +
  xlim(c(-130,-60)) +
  ylim(c(25,50)) +
  #scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a")) +
  scale_color_manual(values=met.brewer("Klimt", 3)) +
  labs(color = "Treatment") +
  #scale_color_gradient2(low = "#f7fbff", mid = "#6baed6", high = "#084594", midpoint = 800) +
  #scale_x_continuous(expand = c(0,0)) +
  #scale_y_continuous(expand = c(0,0)) +
  ggtitle("Treatment (US)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 20),
        legend.key.size = unit(1.1, "cm"),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "italic", size = 12),
        legend.position = "inside",
        legend.position.inside = c(0.05, 0.2))
p6
ggsave("figures/figure1F.png", width = 14, height = 10)
