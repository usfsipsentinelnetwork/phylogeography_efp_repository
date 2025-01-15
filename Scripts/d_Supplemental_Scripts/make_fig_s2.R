library(ape)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(dendextend)
library(magrittr)
library(ggplot2)
library(patchwork)
library(devEMF)

load("Figures_Sept2024/n19clusters.RData")
load("Figures_Sept2024/n19clusters_worldmap.RData")

layout1 <- c(
  area(t = 1, l = 1, b = 2, r = 1),
  area(t = 3, l = 1, b = 5, r =1))

windows();gg.dendro.countries19 + WorldMap2 +  plot_layout(
  design = layout1,
  guides = 'collect') &
  theme(plot.tag = element_text(size=12, face='bold'),
        legend.position = 'bottom',
        legend.text = element_text(size=9))

emf(file = "Final_figures/S2.emf", 8,8, emfPlus=T)
(gg.dendro.countries19) %>%+ WorldMap2 +  plot_layout(
  design = layout1,
  guides = 'collect') &
  theme(plot.tag = element_text(size=12, face='bold'),
        legend.position = 'bottom',
        legend.text = element_text(size=9))
dev.off()

svg("Final_figures/S2.svg", 8,8, family="Cambria Math")
(gg.dendro.countries19)%>%+ WorldMap2 +  plot_layout(
  design = layout1,
  guides = 'collect') &
  theme(plot.tag = element_text(size=12, face='bold'),
        legend.position = 'bottom',
        legend.text = element_text(size=9))
dev.off()