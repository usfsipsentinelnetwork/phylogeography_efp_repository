###########################################################################################
##
## R source code to read and visualize Köppen-Geiger fields (Version of 19 August 2022)                                                                                    
##
## Climate classification after Kottek et al. (2006), downscaling after Rubel et al. (2017)
##
## Kottek, M., J. Grieser, C. Beck, B. Rudolf, and F. Rubel, 2006: World Map of the  
## Köppen-Geiger climate classification updated. Meteorol. Z., 15, 259-263.
##
## Rubel, F., K. Brugger, K. Haslinger, and I. Auer, 2017: The climate of the 
## European Alps: Shift of very high resolution Köppen-Geiger climate zones 1800-2100. 
## Meteorol. Z., DOI 10.1127/metz/2016/0816.
##
## (C) Climate Change & Infectious Diseases Group, University of Veterinary Medicine Vienna
##     Vetmeduni Vienna, Austria
##
###########################################################################################

#setwd("C:/Users/GeoffreyWilliams/OneDrive - USDA/Project_folders_Personal/Projects/Sentinels/Retrospective/Retrospective-May2023-Homework/Map_KG-Global")

load(".RData")

# required packages 
library(raster)
library(dplyr)
library(tidyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(maptools)
library(doBy)

library(rasterVis)
library(rworldxtra)
data(countriesHigh)
library(latticeExtra)
library(sf)
library(ggplot2)

# Read raster files
period='1986-2010'
r <- raster(paste('KG_', period, '.grd', sep=''))

# Color palette for climate classification
climate.colors=c("#960000", "#FF0000", "#FF6E6E", "#FFCCCC", "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64", "#007800", "#005000", "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00", "#550055", "#820082", "#C800C8", "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF", "#64FFFF", "#F5FFFF")

# Legend must correspond to all climate classes, insert placeholders
r0 <- r[1:32]; r[1:32] <- seq(1,32,1)

# Converts raster field to categorical data
r <- ratify(r); rat <- levels(r)[[1]]



# Legend is always drawn in alphabetic order
rat$climate <- c('Af', 'Am', 'As', 'Aw',
                 'BSh', 'BSk', 'BWh', 'BWk',
                 'Cfa', 'Cfb','Cfc', 'Csa', 'Csb', 'Csc', 'Cwa','Cwb', 'Cwc',
                 'Dfa', 'Dfb', 'Dfc','Dfd', 'Dsa', 'Dsb', 'Dsc', 'Dsd','Dwa', 'Dwb', 'Dwc', 'Dwd', 'EF','ET', 'Ocean')
for(i in c("A","B","C","D","E","W","S","f","s","w","m","h","k","a","b","c","d","F","T")) {
  rat[,i] <- grepl(i, rat$climate) %>% as.numeric
}

# Remove the placeholders
r[1:32] <- r0; levels(r) <- rat

# Select region 
  x1=-180; x2=180; y1=-56; y2=90; xat=20; yat=10

r <- crop(r, extent(x1, x2, y1, y2))

# Visualization		

#as.data.frame(r, xy=TRUE)

#if(.Platform$OS.type=="windows") {quartz<-function(...) windows(...)}
#quartz(width=13, height=10)#, dpi=100)

#print(levelplot(r, col.regions=climate.colors, xlab="", ylab="", 
#                scales=list(x=list(limits=c(xmin(r), xmax(r)), at=seq(xmin(r), xmax(r), xat)), 
#                            y=list(limits=c(ymin(r), ymax(r)), at=seq(ymin(r), ymax(r), yat))), colorkey=list(space="top", tck=0, maxpixels=ncell(r)))
#      +latticeExtra::layer(sp.polygons(countriesHigh, lwd=0.25)))

#print(levelplot(r, col.regions=climate.colors, xlab="", ylab="", 
#                scales=list(x=list(limits=c(xmin(r), xmax(r)), at=seq(xmin(r), xmax(r), xat)), 
#                            y=list(limits=c(ymin(r), ymax(r)), at=seq(ymin(r), ymax(r), yat))), colorkey=list(space="top", tck=0, maxpixels=ncell(r))))

#levelplot(r, col.regions=climate.colors, xlab="", ylab="", 
#          scales=list(x=list(limits=c(xmin(r), xmax(r)), at=seq(xmin(r), xmax(r), xat)), 
#                      y=list(limits=c(ymin(r), ymax(r)), at=seq(ymin(r), ymax(r), yat))), colorkey=list(space="top", tck=0, maxpixels=ncell(r)))+
#  theme(axis.line = element_blank())

#ggplot() + geom(data = z, aes(x = lon, y = lat, fill = as.factor(KG))) +
#  scale_fill_manual(values = climate.colors)

#test_spdf <- as(r, "SpatialPixelsDataFrame")
#test_df <- as.data.frame(test_spdf)
#colnames(test_df) <- c("value", "x", "y")

#ggplot(test_df) + geom_tile(aes(y = lat, x = lon, fill = as.factor(KG))) +
#  scale_fill_manual(values = climate.colors) +
#  coord_equal()

#plot(r)

z <- rasterToPoints(r, spatial=T); z <- spTransform(z, CRS=projection(r))

zdf <- data.frame(z, xy=T)

save(zdf, file = "zdf.RData")

#head(zdf)

#zdf$layer %>% factor %>% summary

#zdf$layer <- factor(zdf$layer)

#levels(zdf$layer) <- rat$climate


# Legend is always drawn in alphabetic order
rat$climate <- c('Af', 'Am', 'As', 'Aw',
                 'BSh', 'BSk', 'BWh', 'BWk',
                 'Cfa', 'Cfb','Cfc', 'Csa', 'Csb', 'Csc', 'Cwa','Cwb', 'Cwc',
                 'Dfa', 'Dfb', 'Dfc','Dfd', 'Dsa', 'Dsb', 'Dsc', 'Dsd','Dwa', 'Dwb', 'Dwc', 'Dwd', 'EF','ET', 'Ocean')
for(i in c("A","B","C","D","E","W","S","f","s","w","m","h","k","a","b","c","d","F","T")) {
  rat[,i] <- grepl(i, rat$climate) %>% as.numeric
}
rat[32,3:21] <- 0

zdf <-
  zdf %>%
  dplyr::select(-c(A,B,C,D,E,W,S,f,s,w,m,h,k,a,b,c,d,'F','T',main,prec,temp,climate)) %>%
  #dplyr::mutate(layer = factor(layer))
  dplyr::left_join(
    rat %>%
      dplyr::mutate(
        main = ifelse(A==1,'A',ifelse(B==1,'B',ifelse(C==1,'C',ifelse(D==1,'D',ifelse(E==1,'E',''))))),
        prec = ifelse(W==1,'W',ifelse(S==1,'S',ifelse(f==1,'f',ifelse(s==1,'s',ifelse(w==1,'w',ifelse(m==1,'m','')))))),
        temp = ifelse(h==1,'h',ifelse(k==1,'k',ifelse(a==1,'a',ifelse(b==1,'b',ifelse(c==1,'c',ifelse(d==1,'d',ifelse(F==1,'F',ifelse(T==1,'T',''))))))))
      ),
    by= join_by('layer'=='ID'))

head(zdf)

head(zdf$climate)

save(zdf, file = "zdf.RData")

climate.colors.main   <- climate.colors[c(1,5,9,18,30)] %>% c('white',.)
climate.colors.precip <- climate.colors[c(6,7,10,12,15,2)]%>% c('white',.)

setwd("C:/Users/GeoffreyWilliams/OneDrive - USDA/Project_folders_Personal/Manuscripts/Phylogeographic_Patterns/v2_July7_2024/Repository")
load('world_usca_admin2.RData')

load("zdf.RData")

climate.all<-
  ggplot() + geom_tile(data= zdf, aes(y = y, x = x, fill = as.factor(climate))) +
  scale_fill_manual(values = climate.colors) +
  geom_sf(data = world_usca_admin2, aes(), alpha=0,
          colour='black', linewidth=20) +
 # coord_equal() + 
  ggthemes::theme_map() +
  theme(legend.position = 'none')

zdf$main <- as.factor(zdf$main)
levels(zdf$main) <- c("", "Equatorial", "Arid", "Warm Temperate", "Snow", "Polar")

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

climate.precip<-
  ggplot() + geom_tile(data= zdf, aes(y = coords.x2, x = coords.x1, fill = as.factor(prec))) +
  scale_fill_manual(values = climate.colors.precip) +
  geom_sf(data = world_usca_admin2, aes(), alpha=0,
          colour='black', linewidth=.5) +
  guides(fill = guide_legend(nrow=1)) +
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
  ggplot() + geom_tile(data= zdf, aes(y = coords.x2, x = coords.x1, fill = as.factor(temp))) +
  scale_fill_manual(values = climate.colors.temp)+#, na.translate = F) +
  geom_sf(data = world_usca_admin2, aes(), alpha=0,
          colour='black', linewidth=.5) +
  guides(fill = guide_legend(nrow=2, byrow=T)) +
  ggthemes::theme_map() +
  theme(legend.key.size = unit(.022, 'npc'),
        legend.position = 'bottom',
        legend.title = element_blank(),
        plot.margin = margin(0,0,0,0))

windows(6,4);climate.main
windows(6,4);climate.precip
windows(6,4);climate.temp

save(climate.all, climate.main, climate.precip, climate.temp, file="Figures_July2024/climate.maps.RData")

#######
z <- rasterToPoints(r, spatial=T); z <- spTransform(z, CRS=projection(r))

# first get states so can filter out names
state_prov <- ne_states(c("united states of america", "canada"), returnclass = 'sp')
names(state_prov@data)[which(names(state_prov@data)=='name_en')]<-'name'

# start with countries
world <- ne_countries(scale = 'medium', returnclass = 'sp')
world<- world[-which(world@data$name %in% c("United States", "Canada")),]
countries_states <- spRbind(state_prov[,'name'], world[,'name'])

res <- over(z, countries_states)

KG_country<- data.frame(country=res$name, KG=z@data$layer) %>% filter(!is.na(country))

KG_expanded <- left_join(rename(KG_country, ID=KG), rat, by="ID")

countries_KG_community_matrix <- summaryBy(A+B+C+D+E+W+S+f+s+w+m+h+k+a+b+c+d+F+T~country, KG_expanded, FUN=sum, keep.names=T)
save(countries_KG_community_matrix, file="countries_KG_community_matrix.RData")
#############################

load("countries_KG_community_matrix.RData")

load('../Host_list/Output/ordinated_unifrac_data.Rdata')

bioregions_vector<-
bioregions_vector %>%
  mutate(country_name = gsub("([A-Z]{1}[a-z]*)([A-Z]{1}[a-z]*)([A-Z]{1}[a-z]*)", "\\1 \\2 \\3", country_name))%>%
  mutate(country_name = gsub("([A-Z]{1}[a-z]*)([A-Z]{1}[a-z]*)", "\\1 \\2", country_name))%>%
  mutate(country_name = sub("^Newfoundland$", "Newfoundland and Labrador", country_name))%>%
  mutate(country_name = sub("Bosnia and Herzegovina", "Bosnia and Herz.", country_name))%>%
  mutate(country_name = sub("Brunei Darussalam", "Brunei", country_name))%>%
  mutate(country_name = sub("CÃ´te d'Ivoire", "Côte d'Ivoire", country_name))%>%
  mutate(country_name = sub("Central African Republic", "Central African Rep.", country_name))%>%
  mutate(country_name = sub("Congo, The Democratic Republic of the", "Dem. Rep. Congo", country_name))%>%
  mutate(country_name = sub("Deleware", "Delaware", country_name))%>%
  mutate(country_name = sub("Dominican Republic", "Dominican Rep.", country_name))%>%
  mutate(country_name = sub("CuraÃ§ao", "Curaçao", country_name))%>%
  mutate(country_name = sub("Czechia", "Czech Rep.", country_name))%>%
  mutate(country_name = sub("Equatorial Guinea", "Eq. Guinea", country_name))%>%
  mutate(country_name = sub("French Polynesia", "Fr. Polynesia", country_name))%>%
  mutate(country_name = sub("Georgia Country", "Georgia", country_name))%>%
  mutate(country_name = sub("^Labrador$", "Newfoundland and Labrador", country_name))%>%
  mutate(country_name = sub("Lao People's Democratic Republic", "Lao PDR", country_name))%>%
  mutate(country_name = sub("North Korea", "Korea", country_name))%>%  
  mutate(country_name = sub("Eswatini", "Swaziland", country_name))%>%
  #mutate(country_name = sub("French Guiana", "", country_name))%>%
  mutate(country_name = sub("Quebec", "Québec", country_name))%>%
  mutate(country_name = sub("Russian Federation", "Russia", country_name))%>%
  mutate(country_name = sub("South Korea", "Korea", country_name))%>%
  mutate(country_name = sub("South Sudan", "S. Sudan", country_name))%>%
  mutate(country_name = sub("Syrian Arab Republic", "Syria", country_name))%>%
  mutate(country_name = sub("Viet Nam", "Vietnam", country_name))%>%
  mutate(country_name = sub("North Macedonia", "Macedonia", country_name))

# note that the following join does not include Greenland or Antarctica

bioregions_vector %>%
  rename(country=country_name) %>%
  left_join(countries_KG_community_matrix) %>% filter(grepl("^G", country))

bioregions_vector %>%
  rename(country=country_name) %>%
  right_join(countries_KG_community_matrix) %>% filter(grepl("^G", country))

# these countries are not bound to the climate frame, but are all city states, very small countries, or Greenland or Antarctica

(bioregions_vector %>%
  rename(country=country_name) %>%
  right_join(countries_KG_community_matrix)) %>%
  
  (function (x) x[which(is.na(x$biogeographic_region)),'country'])

merged_climate_zones<- bioregions_vector %>%
  rename(country=country_name) %>%
  left_join(countries_KG_community_matrix)

merged_climate_zones[which(is.na(merged_climate_zones$A)),"country"]

merged_climate_zones_aggregated<-
  aggregate(.~biogeographic_region, merged_climate_zones[,-1], sum)

merged_climate_zones_aggregated_standardized<-
vegan::decostand(merged_climate_zones_aggregated[,-1],"total")

row.names(merged_climate_zones_aggregated_standardized)<-merged_climate_zones_aggregated$biogeographic_region

climate.df<-
vegan::vegdist(merged_climate_zones_aggregated_standardized, method="euclid")%>%
  as.matrix%>%
 (reshape2::melt)#%>%
#  filter(Var1 %in% c("N.North.America","NC.North.America","W.North.America","NE.North.America","SE.North.America","Hawaii")) %>%
#  filter(!(
#    Var2 %in% c("N.North.America","NC.North.America","W.North.America","NE.North.America","SE.North.America") &
#      Var1 %in% c("N.North.America","NC.North.America","W.North.America","NE.North.America","SE.North.America")
#    ))

save(climate.df, file="climate.df.RData")

