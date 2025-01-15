library(ggplot2)
library(patchwork)
library(dplyr)
library(viridis)
library(devEMF)

#setwd("C:/Users/GeoffreyWilliams/OneDrive - USDA/Project_folders_Personal/Manuscripts/Phylogeographic_Patterns/v2_July7_2024/Repository")
load('Covariates_input/world_usca_admin2.RData')
load('Covariates_input/zdf.RData')

zdf1 <- zdf
zdf <- zdf1

# Color palette for climate classification
climate.colors=c("#960000", "#FF0000", "#FF6E6E", "#FFCCCC", "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64", "#007800", "#005000", "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00", "#550055", "#820082", "#C800C8", "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF", "#64FFFF", "#F5FFFF")

zdf$main <- as.factor(zdf$main)
levels(zdf$main) <- c("", "Equatorial", "Arid", "Warm Temperate", "Snow", "Polar")
climate.colors.main   <- climate.colors[c(1,5,9,18,30)] %>% c('white',.)

climate.main<-
  ggplot() + geom_tile(data= zdf, aes(y = coords.x2, x = coords.x1, fill = main)) +
  scale_fill_manual(values = climate.colors.main) +
  geom_sf(data = world_usca_admin2, aes(), alpha=0,
          colour='black', linewidth=.5) +
  guides(fill = guide_legend(nrow=1)) +
  ggthemes::theme_map() +
  theme(legend.key.size = unit(.022, 'npc'),
        legend.position = 'bottom',
        legend.title = element_blank(),
        plot.margin = margin(0,0,0,0))

zdf$prec <- as.factor(zdf$prec)
levels(zdf$prec) <- c("", "Fully humid", "Monsoon", "Dry summer", "Steppe", "Dry winter", "Desert")
zdf$prec <- factor(zdf$prec, levels=
                     c("Fully humid", "Monsoon", "Dry summer", "Steppe", "Dry winter", "Desert", ""))

climate.colors.precip <- climate.colors[c(6,7,10,12,15,2)]%>% c(.,'white')

climate.precip<-
  ggplot() + geom_tile(data= zdf, aes(y = coords.x2, x = coords.x1, fill = prec)) +
  scale_fill_manual(values = climate.colors.precip) +
  geom_sf(data = world_usca_admin2, aes(), alpha=0,
          colour='black', linewidth=.5) +
  guides(fill = guide_legend(nrow=2, byrow=T)) +
  ggthemes::theme_map() +
  theme(legend.key.size = unit(.022, 'npc'),
        legend.position = 'bottom',
        legend.title = element_blank(),
        plot.margin = margin(0,0,0,0))

zdf$temp <- as.factor(zdf$temp)
levels(zdf$temp) <- c("", "Hot summer", "Warm summer", "Cool summer", "Extremely continental",
                      "Polar frost", "Hot arid", "Cold arid", "Polar tundra")
zdf$temp <- factor(zdf$temp, levels=
                     c("Hot summer", "Warm summer", "Cool summer", "Extremely continental",
                       "Polar frost", "Hot arid", "Cold arid", "Polar tundra",""))

climate.colors.temp   <- climate.colors[c(7,6,22:25,30:31)]%>% c(.,'white')

climate.temp<-
  ggplot() + geom_tile(data= zdf, aes(y = coords.x2, x = coords.x1, fill = temp)) +
  scale_fill_manual(values = climate.colors.temp)+#, na.translate = F) +
  geom_sf(data = world_usca_admin2, aes(), alpha=0,
          colour='black', linewidth=.5) +
  guides(fill = guide_legend(nrow=3, byrow=T)) +
  ggthemes::theme_map() +
  theme(legend.key.size = unit(.022, 'npc'),
        legend.position = 'bottom',
        legend.title = element_blank(),
        plot.margin = margin(0,0,0,0))

# for safekeeping
save(
  climate.main,
  climate.precip,
  climate.temp,
  file = 'Figures_Sept2024/climate.maps2.RData'
)

load('Figures_Sept2024/climate.maps2.RData')
load("Figures_Sept2024/ordination_heat_panels.RData")

layout_final <- '
AB
CD
EF
'

leg_heatmap <- ggpubr::as_ggplot(ggpubr::get_legend(heatmap_climate +
                                                      scale_fill_viridis(name = "Distance",
                                                                         option = 'magma',
                                                                         limits = c(0,.75),
                                                                         breaks = c(0,0.25,0.5,.75)) +
                                                      theme(legend.position.inside = c(.5,.5))))

#windows(1,1.25);leg_heatmap
#
#windows();((heatmap_climate + theme(legend.position = 'none'))/
#  (heatmap_climate + theme(legend.position = 'none')) )+ inset_element(leg_heatmap,
#                                                                  .7,.85,.8,1.15)

composite.figure.2.distances <-
  (
    ((climate.main + map_theme_custom + theme(legend.justification='top')) +
    (climate.precip + map_theme_custom + theme(legend.justification='top')) +
    (climate.temp + map_theme_custom + theme(legend.justification='top')) +
    (world_map2 + map_theme_custom + theme(legend.justification='top')) +
    free(heatmap_climate + coord_fixed(.4) + theme(#axis.title = 'element_blank',
                                                   legend.position = 'none')) +
    free(gg1.final + map_theme_custom)) +
    plot_layout(design = layout_final)+
    plot_annotation(tag_levels = list(c('A','B','C','D','E','F',''))) &
    theme(plot.tag = element_text(size=12, face='bold'))
  ) + inset_element(leg_heatmap,
                    -.3,
                    .9,
                    -.15,
                    1.3)

save(composite.figure.2.distances, file='Figures_Sept2024/Fig3_composite_ABCDEF.RData')

windows(
  8,
  8*1.1
);composite.figure.2.distances


###########################################################

emf(file = "Final_figures/Fig3_composite_ABCDEF.emf", 8,
    8.6, emfPlus=T)
composite.figure.2.distances
dev.off()

#svg("Final_figures/Fig3_composite_ABCDEF.svg", 8,
#    8.6, family="Cambria Math")
#composite.figure.2.distances
#dev.off()
