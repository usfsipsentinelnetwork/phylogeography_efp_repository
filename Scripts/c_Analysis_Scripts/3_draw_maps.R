# making the maps

library(countrycode)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(rworldmap)
library(ggthemes)
library(doBy)
library(data.table)
library(magrittr)
library(scatterpie)
library(cowplot)

#load("Host_list/Output/ordinated_unifrac_data.Rdata")

load("Host_list/Output/unifrac_clusters_n_19.Rdata")

bioregions_vector %>% filter(grepl("Quebec", country_name))

bioregions_vector<-rbind(bioregions_vector, data.frame(country_name="Western Sahara", biogeographic_region="Arabia.Sahara"))
bioregions_vector<-rbind(bioregions_vector, data.frame(country_name="Greenland", biogeographic_region="NE.North.America"))

# dont run for n=19
bioregions_vector<-rbind(bioregions_vector, data.frame(country_name="Iceland", biogeographic_region="Eurasian.Palearctic"))
bioregions_vector[which(bioregions_vector$country_name=="Quebec"),"country_name"] <-"Québec"
bioregions_vector[which(bioregions_vector$country_name=="Newfoundland"),"country_name"] <-"NewfoundlandandLabrador"

countries <- bioregions_vector %>%
  mutate(country_name = replace(country_name, which(country_name=="Saint Martin"), "Sint Maarten")) %>%
  mutate(country_name_code = with(., countrycode(sourcevar=country_name, origin='country.name', destination='iso3c'))) %>%
  filter(!is.na(country_name_code))

c16 <- #pals::polychrome(29)
  c(
    #"dodgerblue2"
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    #"black",
    "gold1",
    #"skyblue2",
    "#FB9A99",
    "blue1", # lt pink
    "palegreen2",
    #"#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    
    "darkturquoise",# "green1", "yellow4", #"yellow3",
    "darkorange4", #"brown"#,
    #"#822E1C",
    "#BDCDFF","#FA0087","#2ED9FF","#E31A1C"#,  ## red

    #uncomment the comma above and these for n = 19
      
    #"maroon", "orchid1", "steelblue4", "black"
  )
country_map <- joinCountryData2Map(countries, joinCode = "ISO3", nameJoinColumn = "country_name_code")
mapCountryData(country_map, nameColumnToPlot='biogeographic_region', colourPalette=c16, addLegend=F, mapTitle="")

# ggplot

# first get states so can filter out names
state_prov <- ne_states(c("united states of america", "canada"), returnclass = 'sf')
state_prov$name2 <- gsub(" ", "", state_prov$name)

# start with countries
world <- ne_countries(scale = 'medium', returnclass = 'sf')

#plot(world %>% filter(admin=="France"))

################

## Fix French Guiana - it is currently displayed as part of mainland France. Needs splitting.
FRA <- world %>% filter(admin=="France") # Isolate France row
split <- st_cast(FRA, "POLYGON") # Split multipolygon into many polygons
#
#
#plot(split[1,][1]) # Corsica
#plot(split[2,][1]) # mainland
#plot(split[c(4:5),][1])
#plot(split[3,][1])
#plot(split[6:9,][1]) # Carribbean
#plot(split[10,][1]) # Guyana

split$part <- NA

split[c(1:2),]$part <- 'main'
split[3,]$part <- 'part2'
split[c(4:5),]$part <- 'part3'
split[6:9,]$part<- 'Carribbean'
split[10,]$part<- 'French Guiana'

#plot(split['part'])

#
guy <- split %>% filter(part == "French Guiana") # Isolate polygon for French Guiana

guy$iso_a3_eh <- "GUF"
guy$iso_a3 <- "GUF"
guy$name <- "French Guiana"
guy$admin<- "French Guiana"

plot(guy[1])
FRA_main <- split %>% filter(part == "main") # Isolate polygon for mainland France
plot(FRA_main[1])

FRA_carib<- split %>% filter(part == "Carribbean")
FRA_poly <- split %>% filter(part %in% c("part2", "part3")) # Isolate polygon for mainland France

## We can now replace all French polygons on main map with these 2 separate polygons. Other territories (e.g. French Polynesia) are too small to be viewed so can be ignored (depending on your map's purpose).
baguetteless <- world %>% filter(sovereignt != "France")

#plot(baguetteless[1])


plot(baguetteless[1])
plot(FRA_main[1])
plot(guy[1])

world <- rbind(baguetteless, FRA_main %>% select(-part), guy %>% select(-part))
plot(world[1])

# did it work?
#plot((world %>% filter(name=="France"))[1])
#plot(baguetteless[1])
#plot((FRA_main %>% filter(name=="France"))[1])
#plot((guy %>% filter(name=="France"))[1])


#
# merge countries with bioregion data
countries2 <- rename(countries, iso_a3_eh = country_name_code) %>% filter(!(country_name %in% c(state_prov$name2, "Sint Maarten")))

# the ones that will be dropped for not matching (nothing too concerning)
countries2$country_name[countries2$iso_a3_eh %in% (countries2$iso_a3_eh %>% setdiff(world$iso_a3_eh))]

# French Guiana and Norway are missing from map

# its dropping French Guiana

# THIS IS WHERE WE MERGE THE BIOGEOGRAPHIC INFORMATION
world2<-right_join(world, countries2, by='iso_a3_eh') %>% filter(!is.na(scalerank))

#plot(world2[1])

# did it work?
#plot((world2 %>% filter(name=="France"))[1])

# by now its already wrong
#world2[,c('name','biogeographic_region')] %>% arrange(name) %>% print(n=300)

#### MAP IS MISSING A BUNCH OF CANADA (QUEBEC, ETC)

# next
# merge states with bioregion data
state_prov_bioregions <- bioregions_vector %>% filter(country_name %in% c(state_prov$name2))
state_prov2<- right_join(state_prov, state_prov_bioregions, by=join_by('name2'=='country_name'))

# rbind with states/provinces
incommon <- names(world2) %>% intersect(names(state_prov2))
world_usca_admin <- rbind(world2[,incommon], state_prov2[,incommon])
world_usca_admin[,c('admin','name','biogeographic_region')] %>% arrange(name) %>% print(n=300)

plot((world_usca_admin %>% filter(admin %in% c("France", "French Guiana")))['biogeographic_region'])

# by now its already wrong
#world2[,c('name','biogeographic_region')] %>% arrange(name) %>% print(n=300)
#world_usca_admin[,c('name','biogeographic_region')] %>% arrange(name) %>% print(n=300)
#world_usca_admin %>% filter(admin=="France")
#plot((world_usca_admin %>% filter(admin=="France"))[1])

# remove the non contiguous regions (it seems)
regions <- levels(world_usca_admin$biogeographic_region)
sf::sf_use_s2(FALSE)
sf::sf_use_s2(FALSE)

########!!!!!!!!!!!!!!!!!!
# gonna need to make world_usca_admin3 with W Europe (+ Spain) and Australia

world_usca_admin2 <- NULL;
for (x in regions) {
  print(x)
  y <- filter(world_usca_admin, biogeographic_region == x)
  world_usca_admin2 <-
    rbind(
      world_usca_admin2,
      st_sf(
        biogeographic_region=x,
        geometry=
          nngeo::st_remove_holes(
            st_make_valid(
              st_union(
                st_make_valid(y))))))}

save(world_usca_admin2, file='world_usca_admin2.RData')

world_usca_admin3 <- NULL;
for (x in c('Eurasian.Palearctic','Iberia.N.Africa')) {
  print(x)
  y <- filter(world_usca_admin, biogeographic_region == x)
  
  if (x == 'Eurasian.Palearctic')
    y <- filter(y, !(name %in% c('Russia','Mongolia')))
  else
    y<- filter(y, name %in% c('Portugal','Spain'))
  print(y)
  world_usca_admin3 <-
    rbind(
      world_usca_admin3,
      st_sf(
        biogeographic_region=x,
        geometry=
          nngeo::st_remove_holes(
            st_make_valid(
              st_union(
                st_make_valid(y))))))}

# make plot - just of world
WorldMap2 <- ggplot(data = world_usca_admin2) +#ggplot(data = world_usca_admin) + # ggplot(data = world2) +
  geom_sf(aes(fill = biogeographic_region)) +
  #  theme_wsj(color = "white", base_family = "mono") +
  theme_map() +
  theme(legend.position="none") +
  scale_fill_manual(values=c16)+
  guides(fill=guide_legend(title=NULL))+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))

windows(15,10); WorldMap2   

# only for n = 19 clusters
#save(WorldMap2, file = "Figures_Sept2024/n19clusters_worldmap.RData")

c16.named <- c16
names(c16.named) <- world_usca_admin2$biogeographic_region

extract_region <- function(region_to_extract = "Hawaii",
                           x_factor_min=0,
                           x_factor_max=0,
                           y_factor_min=0,
                           y_factor_max=0,
                           scale=c16.named,
                           world_map_layer = world_usca_admin2) {
  extracted <-
    world_map_layer %>% filter(biogeographic_region == region_to_extract) %>%
    ggplot()+
    geom_sf(aes(fill = biogeographic_region)) +
    theme_map() +
    theme(legend.position="none") +
    scale_fill_manual(values=c16.named)+
    guides(fill=guide_legend(title=NULL))
  scale.x <- layer_scales(extracted)$x$range$range
  scale.y <- layer_scales(extracted)$y$range$range
  scale.x.new <- c(scale.x[1] + x_factor_min * diff(scale.x),
                   scale.x[2] + x_factor_max * diff(scale.x))
  scale.y.new <- c(scale.y[1] + y_factor_min * diff(scale.y),
                   scale.y[2] + y_factor_max * diff(scale.y))
  extracted +
    scale_x_continuous(limits=scale.x.new, expand=c(0,0))+
    scale_y_continuous(limits=scale.y.new, expand=c(0,0))
}

sink_map_HI <- extract_region("Hawaii", 0.77, 0, 0, -0.63)
sink_map_WNA <- extract_region("W.North.America")
sink_map_NNA <- extract_region("N.North.America", 0, -0.65)
sink_map_NEN <- extract_region("NE.North.America")
sink_map_SEN <- extract_region("SE.North.America")

sink_map_EU <- world_usca_admin3 %>% ggplot() +  # this one needs to be modified to just western Europe
  geom_sf(aes(fill = biogeographic_region))+
  theme_map() +
  theme(legend.position="none") +
  scale_fill_manual(values=c16.named)+
  guides(fill=guide_legend(title=NULL))+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))

scale.x <- layer_scales(sink_map_EU)$x$range$range
scale.y <- layer_scales(sink_map_EU)$y$range$range
scale.x.new <- c(scale.x[1] + 0.49* diff(scale.x),
                 scale.x[2])
scale.y.new <- c(scale.y[1] + 0.3 * diff(scale.y),
                 scale.y[2] + -0.1 * diff(scale.y))
sink_map_EU<-
  sink_map_EU +
    scale_x_continuous(limits=scale.x.new, expand=c(0,0))+
    scale_y_continuous(limits=scale.y.new, expand=c(0,0))

sink_map_AU <-  # this one needs to be modified to just Australia
  filter(world_usca_admin, name == 'Australia') %>%
  mutate(biogeographic_region = 'Australia') %>%
  ggplot() +
  geom_sf(aes(fill = biogeographic_region))+
  theme_map() +
  theme(legend.position="none") +
  scale_fill_manual(values=c16.named %>% c(Australia = "yellow"))+
  guides(fill=guide_legend(title=NULL))+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))
scale.x <- layer_scales(sink_map_AU)$x$range$range
scale.y <- layer_scales(sink_map_AU)$y$range$range
scale.x.new <- c(scale.x[1],
                 scale.x[2] - 0.1* diff(scale.x))
scale.y.new <- c(scale.y[1] + 0.2 * diff(scale.y),
                 scale.y[2])
sink_map_AU<-
  sink_map_AU +
  scale_x_continuous(limits=scale.x.new, expand=c(0,0))+
  scale_y_continuous(limits=scale.y.new, expand=c(0,0))

windows(2,2);sink_map_HI
windows(2,2);sink_map_WNA
windows(2,2);sink_map_NNA
windows(2,2);sink_map_NEN
windows(2,2);sink_map_SEN
windows(2,2);sink_map_EU
windows(2,2);sink_map_AU

save(WorldMap2, "Figures_July2024/1a_Map_final_panel1.RData")
save(WorldMap2,
     sink_map_HI,
     sink_map_WNA,
     sink_map_NNA,
     sink_map_NEN,
     sink_map_SEN,
     sink_map_EU,
     sink_map_AU,
     file="Figures_July2024/Map_panels.RData")

