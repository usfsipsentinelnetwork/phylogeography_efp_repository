#library(devEMF)
#library(semPlot)
#library(semptools)
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)

# load necessary data

load("Covariates_input/countries_KG_community_matrix.RData")
load("Host_list/Output/ordinated_unifrac_data.Rdata")
load("Figures_Sept2024/Map_panels.RData")
load("Figures_Sept2024/S1_ordination_of_cluster.RData")

#    windows(); gg1#+theme(legend.position='none')

# we want to add numbers to clusters and legend ?



world_map <- WorldMap2 +
  #scale_x_discrete(expand=c(0,0)) +
  theme(plot.margin = margin(0,0,0,0))

#map_width <- layer_scales(world_map)$x$range$range %>% diff
#map_height <- layer_scales(world_map)$y$range$range %>% diff

#map_width/map_height

#(ps <- plotsize(world_map, width=1, height=1))

#(psratio1 <- ps$width/ps$height)

world_map$data$biogeographic_region[world_map$data$biogeographic_region=="Carribbean"]<-"Caribbean"

dendro_order_source<-
  c("SE North America",
    "NE North America",
    "Eurasian Palearctic",
    "Iberia N Africa",
    "Middle East",
    "W North America",
    "N North America",
    "SE Asia",
    "S Cone Pacific",
    "N Pacific",
    "Caribbean",
    "C America",
    "Amazonia",
    "Subsaharan Africa",
    "Arabia Sahara",
    "Hawaii",
    "Australia") %>% gsub("\\ ","\\.",.)

newdf <-
  data.frame(names = 
               factor(world_map$data$biogeographic_region %>% c("Australia"),
                      levels=dendro_order_source)) %>%
  arrange(as.numeric(names)) %>%
  mutate(extrasp = c(rep(" ", 9), rep("", 8))) %>%
  mutate(numbers = 1:17) %>%
  mutate(annotated = paste(extrasp, numbers, extrasp, "   ", names, sep = ""))

#dendro_order_source

world_map$data$biogeographic_region <- with(newdf[-17,]%>%arrange(as.character(names)),
                                            factor(annotated, levels=c(annotated, "17   Australia")))

c16_b <- c16#[-17]
names(c16_b) <- levels(world_map$data$biogeographic_region)

centroids <- read.csv('Covariates_input/landarea_and_centroids.csv') %>%
  dplyr::rename(biogeographic_region=X) %>%
  dplyr::mutate(biogeographic_region=sub("Carribbean","Caribbean",biogeographic_region)) %>%
  left_join(newdf, by= join_by('biogeographic_region'=='names'))

centroids[which(centroids$biogeographic_region=='Hawaii'),'lon'] <-
  centroids[which(centroids$biogeographic_region=='Hawaii'),'lon'] + 10

centroids[which(centroids$biogeographic_region=='Caribbean'),'lon'] <-
  centroids[which(centroids$biogeographic_region=='Caribbean'),'lon'] + 10

centroids[which(centroids$biogeographic_region=='N.Pacific'),'lon'] <-
  centroids[which(centroids$biogeographic_region=='N.Pacific'),'lon'] + 15

centroids[which(centroids$biogeographic_region=='NE.North.America'),'lon'] <-
  centroids[which(centroids$biogeographic_region=='NE.North.America'),'lon'] - 20
centroids[which(centroids$biogeographic_region=='NE.North.America'),'lat'] <-
  centroids[which(centroids$biogeographic_region=='NE.North.America'),'lat'] - 16

c16_c <- c16_b
names(c16_c) <- gsub("\\.",
                     " ",
                     names(c16_c))%>%
  gsub("North America", "CONUS/Canada", .)

# add centroids to map
world_map2 <- world_map + 
  geom_text(data=centroids %>% filter(!grepl("Hawaii", annotated)), mapping=aes(x=lon, y=lat, label=numbers)) +
  scale_fill_manual(values=c16_c, limits = c(
    names(c16_c) %>%
      grep("Hawaii",.,invert=T,value=T) %>%
      grep("Australia",.,invert=T,value=T) %>%
      sort
  ))+
  guides(fill = guide_legend(nrow=6, byrow=T)) +
  ggthemes::theme_map()+
  theme(legend.title = element_blank(),
    legend.position = 'bottom',
    legend.key.size = unit(.022, 'npc'), 
    legend.text = element_text(#size=7.5,
                               margin = margin(0,0,0,-9.25)),
    plot.margin = margin(0,0,0,0))

world_map2$data$biogeographic_region <-
  factor(world_map2$data$biogeographic_region,
         levels=newdf$annotated)

newdf <-
  left_join(newdf,
            (dat2.summary %>%
               mutate(
                 biogeographic_region =
                   sub('Carribbean',
                       'Caribbean',
                       biogeographic_region))),
            by=join_by('names' == 'biogeographic_region'))

#world_map2$data$biogeographic_region <-
#  factor(world_map2$data$biogeographic_region,
#         levels=newdf$annotated)

levels(world_map2$data$biogeographic_region) <- gsub("\\.",
                                                     " ",
                                                     levels(world_map2$data$biogeographic_region))

levels(world_map2$data$biogeographic_region) <-
  gsub("North America", "CONUS/Canada", levels(world_map2$data$biogeographic_region))

#windows();world_map2

#((world_map2 + world_map2) / (world_map2 + world_map2) ) + plot_layout()

#      windows(6,4); world_map2
      
#      windows(); gg1+theme(legend.position = 'none')+
#                    geom_text(
#                      data=newdf[-17,],
#                      mapping=aes(
#                        x= MDS1.mean,
#                        y= MDS2.mean,
#                        label = numbers))
      
#      composite.panel.1<-(world_map2 / (gg1+theme(legend.position = 'none')+
#                                          geom_text(
#                                            data=newdf[-17,],
#                                            mapping=aes(
#                                              x= MDS1.mean,
#                                              y= MDS2.mean,
#                                              label = numbers
#                                            )
#                                          )#+
#                                          #theme(legend.position='none')
#                                        )) + plot_layout(guides='collect')
#      windows();composite.panel.1

gg1.final <- gg1+theme(legend.position = 'none')+
                      geom_text(
                        data=newdf[-17,],
                        mapping=aes(
                          x= MDS1.mean,
                          y= MDS2.mean,
                          label = numbers))



load("Covariates_input/climate.df.RData")

#level_orders <-
#  c(
#    "Hawaii",
#    "Australia",
#    "Eurasian.Palearctic",
#    "W.North.America",
#    "N.North.America",
#    "NE.North.America",
#    "SE.North.America",
#    "Iberia.N.Africa",
#    "Middle.East",
#    "SE.Asia",
#    "S.Cone.Pacific",
#    "N.Pacific",
#    "Caribbean",
#    "C.America",
#    "Amazonia",
#    "Subsaharan.Afica",
#    "Arabia.Sahara"
#  )

climate.df$Var1 <- factor(climate.df$Var1 %>% gsub("Carribbean", "Caribbean", .), levels = newdf$names[c(16,6:7,2:1,3,17,4:5,8:15)])
climate.df$Var2 <- factor(climate.df$Var2 %>% gsub("Carribbean", "Caribbean", .), levels = newdf$names[c(16,6:7,2:1,3,17,4:5,8:15)])

climate.df <- arrange(climate.df, Var1, Var2)

levels(climate.df$Var2) <-levels(climate.df$Var2) %>% gsub("\\."," ",.)
levels(climate.df$Var1) <-levels(climate.df$Var1) %>% gsub("\\."," ",.)

#index_value <- NULL
#for(i in 1:17) {
#  index_value <- c(index_value, seq(1+(i-1)*17, length.out=i))
#}

climate.df2 <-
  climate.df %>% #[-index_value,] %>%
  #filter(Var2 != "Australia" & Var2 != "Hawaii") %>%
  #filter(Var1 != Var2) %>%
  filter(Var1 %in% c("Australia", "Eurasian Palearctic", "NE North America",
                     "SE North America", "N North America", "W North America",
                     "Hawaii"))

climate.df2$Var2 <- factor(climate.df2$Var2, levels = rev(levels(climate.df2$Var2)))

levels(climate.df2$Var1) <-
  gsub("North America", "CONUS/Canada", levels(climate.df2$Var1))
levels(climate.df2$Var2) <-
  gsub("North America", "CONUS/Canada", levels(climate.df2$Var2))

#    windows();ggplot(data=climate.df2, aes(Var1, Var2, fill=value)) +
#      geom_tile() +
#      guides(fill = guide_legend(title = "Dissimilarity"))+
#      scale_fill_distiller(palette = "Greys", direction=1)+
#      #xlab("Sink")+
#      #ylab("Source") +
#      theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#                              #axis.text.y = element_blank(),
#                              axis.title = element_blank())

heatmap_climate <-
  ggplot(data=climate.df2 %>%
           
           filter(!(
             (Var2 == "Eurasian Palearctic" & Var1 %in% c(unique(as.character(climate.df2$Var2)) %>% setdiff("Australia"))) |
               (Var2 == "Hawaii" & Var1 %in% c(unique(as.character(climate.df2$Var2)) %>% setdiff(c("Australia","Eurasian Palearctic")))) |
               (Var2 == "N CONUS/Canada" & Var1 %in% c(unique(as.character(climate.df2$Var2)) %>% setdiff(c("Australia","Eurasian Palearctic","Hawaii")))) |
               (Var2 == "NE CONUS/Canada" & Var1 %in% c(unique(as.character(climate.df2$Var2)) %>% setdiff(c("Australia","Eurasian Palearctic","Hawaii","N CONUS/Canada")))) |
               (Var2 == "SE CONUS/Canada" & Var1 %in% c(unique(as.character(climate.df2$Var2)) %>% setdiff(c("Australia","Eurasian Palearctic","Hawaii","N CONUS/Canada","NE CONUS/Canada")))) |
               (Var2 == "W CONUS/Canada" & Var1 %in% c(unique(as.character(climate.df2$Var2)) %>% setdiff(c("Australia","Eurasian Palearctic","Hawaii","N CONUS/Canada","NE CONUS/Canada","SE CONUS/Canada"))))
           )) %>%
           
#           filter(!(
#             grepl("CONUS/Canada", Var1) &
#               grepl("CONUS/Canada", Var2))) %>%
           filter(Var1 != Var2),
         aes(Var1, Var2, fill=value)) +
  geom_tile() +
  scale_x_discrete(
    limits = #sort(unique(climate.df2$Var1))
      c(
        "Australia",
        "Eurasian Palearctic",
        "Hawaii",
        "N CONUS/Canada",
        "NE CONUS/Canada",
        "SE CONUS/Canada",
        "W CONUS/Canada"
  )) + scale_y_discrete(
    limits = #sort(unique(as.character(climate.df2$Var2)))

      rev(c("Australia",
            "Eurasian Palearctic",
            "Hawaii",
            "N CONUS/Canada",
        "NE CONUS/Canada",
        "SE CONUS/Canada",
        "W CONUS/Canada"
        ) %>%
      (function(x) c(
        x,
        sort(setdiff(unique(as.character(climate.df2$Var2)),x))))) %>%
      
      setdiff(c("Australia","Hawaii"))
      
#      rev(c("N CONUS/Canada",
#        "NE CONUS/Canada",
#        "SE CONUS/Canada",
#        "W CONUS/Canada",
#        "Australia",
#        "Hawaii",
#        "Eurasian Palearctic"
#      ) %>% (function(x) c(x,sort(setdiff(climate.df2$Var2,x)))))
  ) + scale_fill_viridis(name = "Dissimilarity",
    option = 'magma',
    limits = c(0,.75),
    breaks = c(0,0.25,0.5,.75)) +
  xlab("\nSink")+ylab("Source\n")+
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, #vjust = 0.5,
                                                     hjust=1),
                          #axis.text.y = element_blank(),
                          #axis.title = element_blank(),
                          legend.position.inside = c(0.88,0.9),
                          legend.position='inside',
                          legend.key.width = unit(.015, 'npc'),
                          legend.key.height = unit(0.012,'npc')
                          )

windows(); heatmap_climate + coord_fixed(.4) + guides()

########################

# put everything together

map_theme_custom <-
  theme(legend.key.size = unit(.0175, 'npc'),
        legend.key.spacing = unit(.0015, 'npc'))

save(world_map2,
     heatmap_climate,
     gg1.final,
     map_theme_custom,
     file = "Figures_Sept2024/ordination_heat_panels.RData")
