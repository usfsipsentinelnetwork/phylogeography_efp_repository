
library(dplyr)
library(doBy)
library(reshape2)
library(tidyverse)
library(sjPlot)
library(cowplot)

load("Host_list/Output/ordination_of_cluster.RData")
load("Host_list/Output/ordinated_unifrac_data.Rdata")
load("Host_list/Output/countries_unifrac_3.RData")
load("Host_list/Output/countries_unifrac_clusters1.RData")
load("Host_list/Output/countries_unifrac_corrected_transformed.RData")



#gg1
##### analysis

####### !!!!!!!!! issues caught 2/12/25 !!!!!!!!!
######
# Duplicated origins are kept
# Also some duplicate destinations
# Duplicate of 'Melampsora laricis-populina' and 'Melampsora larici-populina'
#
# DO NOT use Expanded_data_EU_AU_NA_HI_noduplicateorigins_UpdatedFeb2025.csv or ...Sept2024_noduplicateorigins.csv
# Divide by number of origins (eg, in N America) for duplicated origins (Still not implemented as of 2/12/25)
#
# Need to add DISTINCT to deal with duplicate destinations
# 
#
#######
#######

pathogens <-
  read.csv("Pathogens/Input/Expanded_data_EU_AU_NA_HI_Updated_Feb2025.csv")%>%
#  read.csv("Pathogens/Input/Expanded_data_EU_AU_NA_HI_Sept2024.csv")%>%
  dplyr::select(c("Pathogen.Species", "Geographic.Origin.Cluster", "Invaded.Range.Cluster")) %>% distinct #%>%
  #filter(Geographic.Origin.Cluster != "Unknown")

pathogens$Invaded.Range.Cluster <- replace(pathogens$Invaded.Range.Cluster, pathogens$Invaded.Range.Cluster=="S.Cone.Pacific", "Australia")
  #  read.csv("Pathogens/Input/Table1_Pathogens_Established_NorthAmerica_analysis.2.7.24.csv")%>%
  #  dplyr::select(c("Pathogen.Species", "Geographic.Origin.Cluster", "Invaded.Range.in.North.America.Cluster")) %>%
 # filter(Geographic.Origin.Cluster != "") %>%
 # filter(Invaded.Range.in.North.America.Cluster != "")
  
pathogens <- pathogens %>%
  filter(Geographic.Origin.Cluster != Invaded.Range.Cluster)
  
# added 3/1/2025
pathogens <- pathogens %>%
  filter(!(grepl("North.America",Geographic.Origin.Cluster) & grepl("North.America",Invaded.Range.Cluster)))

pathogens$Geographic.Origin.Cluster<-as.factor(pathogens$Geographic.Origin.Cluster)

#path.sum <- summaryBy(Pathogen.Species ~ Geographic.Origin.Cluster + Invaded.Range.in.North.America.Cluster, data=pathogens, FUN=length)
#path.sum.wide <- dcast(path.sum, Invaded.Range.in.North.America.Cluster ~ Geographic.Origin.Cluster)


path.sum <- summaryBy(Pathogen.Species ~ Geographic.Origin.Cluster + Invaded.Range.Cluster, data=pathogens, FUN=length)

# added/modified 3/1/2025 to account for double origins
path.sum.factor <- summaryBy(Geographic.Origin.Cluster ~ Pathogen.Species,
                             data=dplyr::select(pathogens,Pathogen.Species,Geographic.Origin.Cluster)%>%
                               distinct, FUN=length)
path.sum <-
  pathogens %>% left_join(path.sum.factor %>% mutate(origin.factor = 1/Geographic.Origin.Cluster.length)) %>%
  summaryBy(origin.factor ~ Geographic.Origin.Cluster + Invaded.Range.Cluster, data=., FUN=sum) %>%
  right_join(path.sum)

path.sum.wide <- dcast(path.sum, Invaded.Range.Cluster ~ Geographic.Origin.Cluster, value.var='origin.factor.sum')


path.sum.wide[is.na(path.sum.wide)]<-0



rownames(path.sum.wide) <- path.sum.wide$Invaded.Range.Cluster

# NOT NEEDED FOR EXPANDED VERSION
#path.sum.wide <- path.sum.wide[,-which(grepl("North\\.America", colnames(path.sum.wide)))]%>%
#  cbind(path.sum.wide %>%
#              filter(Invaded.Range.in.North.America.Cluster=="Hawaii") %>%
#              dplyr::select(grep("North\\.America", colnames(path.sum.wide), value=T)) %>% .[,-1] %>%
#          cbind(data.frame(W.North.America=0, N.North.America=0))%>%
#          bind_rows(
#            data.frame(row.names=grep("North\\.America", rownames(path.sum.wide), value=T)) #%>% row.names
#                    ))

path.sum.wide <- path.sum.wide %>% dplyr::select(-Invaded.Range.Cluster)
path.sum.wide %>%cbind(Invaded.Range.Cluster=row.names(path.sum.wide),.)%>%
  write.csv("Pathogens/Output/path.sum.expanded-March4-2025.csv", row.names=F) # same as march 4

#path.sum.wide <- path.sum.wide[,-1]

# go back to script and modify this as well - March 1, 2025
# can do this here
#
abcd <- function (a,b,c,d,e)
  {cbind(a[which(c != d)],b[which(c != d)]) %>%
  dist %>%
  as.matrix %>%
  mean}
# now you are counting australia twice? - actually no, should work fine except when you are looking at S.Cone.Pacific? - should be rest of it
dats2 <- dats2 %>% mutate(biogeographic_region = sub("Carribbean","Caribbean",biogeographic_region))

dats3 <- dats2 %>%
  rbind(dats2[rownames(dats2)=="Australia",] %>% mutate(biogeographic_region = "Australia"))
#
all_pairs <-
  expand_grid(
    region1 = unique(dats3$biogeographic_region),
    region2 = unique(dats3$biogeographic_region)) %>%
  mutate(combo = paste(region1, region2, sep=' - ')) %>%
  rbind(
    expand_grid(
      region1 = unique(dats3$biogeographic_region),
      region2 = unique(dats3$biogeographic_region)) %>%
      mutate(combo = paste(region2, region1, sep=' - ')) 
  )
# added March 3, 2025
cdef <- function (c,d,e,f) {
  #f[attr(f,"Labels") %in% c[d==e], attr(f,"Labels") %in% c[d!=e]] %>% dist %>% mean
  
  if(!(c == "Australia" | d == "Australia")) {
    as.matrix(f)[attr(f,"Labels") %in% e$country[which(e$biogeographic_region == c)],
                 attr(f,"Labels") %in% e$country[which(e$biogeographic_region == d)]] %>% mean
  } else if (c == "Australia") {
    as.matrix(f)[attr(f,"Labels") == "Australia",
                 attr(f,"Labels") %in% setdiff(e$country[which(e$biogeographic_region == d)],"Australia")] %>% mean
  } else if (d == "Australia") {
    as.matrix(f)[attr(f,"Labels") %in% setdiff(e$country[which(e$biogeographic_region == c)],"Australia"),
                 attr(f,"Labels") == "Australia"] %>% mean
  }
  
 #%>% (function(x) log(1-x))
  
}

e <- dats3 %>%
  dplyr::select(biogeographic_region) %>%
  mutate(country = row.names(dats3))

as.matrix(uf.countries.2)[attr(uf.countries.2,"Labels") %in% e$country[which(e$biogeographic_region == "SE.North.America")],
             attr(uf.countries.2,"Labels") %in% e$country[which(e$biogeographic_region == "Middle.East")]]
as.matrix(uf.countries.2)[attr(uf.countries.2,"Labels") %in% e$country[which(e$biogeographic_region == "SE.North.America")],
                        attr(uf.countries.2,"Labels") %in% e$country[which(e$biogeographic_region == "Middle.East")]] %>% mean

cdef("Australia",
     "S.Cone.Pacific",
     dats3 %>%
       dplyr::select(biogeographic_region) %>%
       mutate(country = row.names(dats3)),
     uf.countries.2)

dat2.summary.march3.2025 <-
  
  expand_grid(region1 = unique(dats3$biogeographic_region),
              region2 = unique(dats3$biogeographic_region)) %>%
  
  filter(region1 != region2) %>%
  
  group_by(region1, region2) %>%
  summarise(
    unifrac_group_distance =
      cdef(region1,
           region2,
         dats3 %>%
         dplyr::select(biogeographic_region) %>%
         mutate(country = row.names(dats3)),
         uf.countries.2))
         
#         uf.countries)) %>%
#  
#  mutate(unifrac_group_distance_log = log(1-(unifrac_group_distance-.00001)))

# right now doesn't have Australia  

#dat2.summary.march3.2025 %>% na.omit %>% with(hist(log(1-unifrac_group_distance)))

distinct(dat2.summary.march3.2025)

hist(dat2.summary.march3.2025$unifrac_group_distance)#_log)

  #dats3 %>%
  #select(biogeographic_region) %>%
  #mutate(country = row.names(dats3)) %>%
  #expand_grid()
  
#  dats3 %>% mutate(country_name = rownames(dats3)) %>%
#  left_join(all_pairs, by = join_by(biogoegraphic_region = region1)) %>%
#  dplyr::filter(!(grepl("S.Cone.Pacific",combo) & grepl("Australia",combo) & country=="Australia")) %>%
#  dplyr::group_by(combo) %>%
#  dplyr::summarise(
#    unifrac_group_distance
#  )
#
dat2.summary.march1.2025 <-
  dats3 %>%
  dplyr::mutate(country = row.names(dats3)) %>%
  left_join(all_pairs, by = join_by(biogeographic_region == region1)) %>%
  
  # deal with fact we are counting Australia twice when comparing to S.Cone.Pacific
  dplyr::filter(!(grepl("S.Cone.Pacific",combo) & grepl("Australia",combo) & country=="Australia")) %>%
  dplyr::group_by(combo) %>%
  dplyr::summarise(
    nmds_group_distance = abcd(MDS1,MDS2,biogeographic_region,region2)
    ) %>%
  filter(!is.na(nmds_group_distance))
#
# end added March 1 2025
hist(dat2.summary.march1.2025$nmds_group_distance)

missing <-
  dat2.summary$biogeographic_region %>%
  setdiff(c("SE.North.America","NE.North.America","W.North.America","N.North.America","Hawaii","Australia")) %>%
  sub("Carribbean","Caribbean",.) %>%
  setdiff(colnames(path.sum.wide))

for (i in missing) {
  path.sum.wide[,i]<-0
}

# back to long format
path.sum <- path.sum.wide %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname) %>%
  dplyr::rename(Invaded.Range.Cluster=rowname) %>%
  dplyr::rename(Geographic.Origin.Cluster=key) %>%
  dplyr::rename(Pathogen.Species.length=value) %>% 
  filter(!(grepl("North\\.America",Invaded.Range.Cluster) & grepl("North\\.America",Geographic.Origin.Cluster))) %>%
  filter(Invaded.Range.Cluster != Geographic.Origin.Cluster) %>%
  na.omit


# new version march 1 2025

# are we considering movement within north america? i didn't think so...

# add distances
# dat2.summary <- dat2.summary %>% mutate(biogeographic_region = sub("Carribbean","Caribbean",biogeographic_region))


#dist_to_Nearctic2<-(dat2.summary[2:3] %>%
#                      dist(upper=T) %>%
#                      as.matrix) %>% as.data.frame
#row.names(dist_to_Nearctic2) <- colnames(dist_to_Nearctic2) <- dat2.summary$biogeographic_region
#dist_to_Nearctic2$Invaded.Range.Cluster <- row.names(dist_to_Nearctic2)
#dist_to_Nearctic2_melt <- melt(as.data.frame(dist_to_Nearctic2), variable.name="Geographic.Origin.Cluster", value.name="NMDS.Distance")
#dist_to_Nearctic2_melt$Geographic.Origin.Cluster <- as.character(dist_to_Nearctic2_melt$Geographic.Origin.Cluster)

dist_to_Nearctic2_melt <-
  
  # march 3
  dat2.summary.march3.2025 %>%
  rename(Invaded.Range.Cluster=region1, Geographic.Origin.Cluster=region2) %>%
  
  #dat2.summary.march1.2025 %>%
  #separate(combo, into=c("Invaded.Range.Cluster", "Geographic.Origin.Cluster"), sep = ' - ') %>%

  mutate(NMDS.Distance =
        #   nmds_group_distance
          
          # march 3
          unifrac_group_distance#_log
          
         ) %>%
  #select(-nmds_group_distance)
  dplyr::select(-unifrac_group_distance)

to_analyze <- path.sum %>%
  left_join(
    dist_to_Nearctic2_melt,
    by = c("Invaded.Range.Cluster", "Geographic.Origin.Cluster")) %>%
  dplyr::rename(Number.of.pests=Pathogen.Species.length)

summary(glm(Number.of.pests ~ NMDS.Distance, data=to_analyze %>% filter(!is.na(NMDS.Distance)), family="poisson"))

#with(to_analyze, plot(NMDS.Distance, Number.of.pests+0.001))

#with(to_analyze, plot(NMDS.Distance, Number.of.pests+0.001))

# removed spatauto part of script 3/1/25

load("Covariates_Input/climate.df.RData") # updated with Hawaii - NEED TO ADD EUROPE/AUSTRALIA
#load('Covariates_Input/trade_aggregated.RData')
#load('Covariates_Input/HI_trade_aggregated.RData')
#hawaiiNAm <- read.csv('Covariates_Input/HI_NorthAmerica_aggregated_1989-present-inflation-adjusted.csv') # - ADD EUROPE/AUSTRALIA

# add climate


to_analyze2 <- climate.df %>% mutate(Var2 = sub("Carribbean","Caribbean",Var2)) %>%
  dplyr::rename(Invaded.Range.Cluster=Var1) %>%
  dplyr::rename(Geographic.Origin.Cluster = Var2) %>%
  dplyr::rename(Climate.Distance = value) %>%
  left_join(to_analyze, .)


#plot(lm4)

##
library(PerformanceAnalytics)

windows(6,6)
chart.Correlation(to_analyze2[,3:5], histogram=TRUE, pch=19)

# add trade

trade_aggregated <- read.csv('Covariates_input/trade_aggregated_PPP_wide.csv', check.names=F)


trade_aggregated$biogeographic_region <-
  gsub("Carribbean","Caribbean",trade_aggregated$biogeographic_region)
trade_aggregated$`Reporter Name` <-
  gsub("Europe","Eurasian.Palearctic",trade_aggregated$`Reporter Name`)

#trade_aggregated$cumulativeT <-
#  replace(trade_aggregated$cumulativeT,
#          trade_aggregated$biogeographic_region)
  

## add trade
to_analyze3<-
  trade_aggregated %>% dplyr::select(biogeographic_region, "Reporter Name", cumulativeT) %>%
  dplyr::rename(Geographic.Origin.Cluster = biogeographic_region) %>%
  dplyr::rename(Invaded.Range.Cluster = "Reporter Name") %>%
  dplyr::rename(Cumulative.Import.Value.T = cumulativeT) %>%
  #dplyr::select(c(Geographic.Origin.Cluster,Cumulative.Import.Value.T)) %>%
  left_join(to_analyze2  %>%
              mutate(InvRange2 = 
                       ifelse(grepl("North.America", Invaded.Range.Cluster),
                              "North America", Invaded.Range.Cluster)) %>%
              mutate(Origin2 = 
                       ifelse(grepl("North.America", Geographic.Origin.Cluster),
                              "North America", Geographic.Origin.Cluster)),
            ., by=join_by(InvRange2==Invaded.Range.Cluster,
                          Origin2==Geographic.Origin.Cluster)) %>%
  filter(InvRange2 != Origin2) %>%
  mutate(NorthAm  = InvRange2=="North America") %>%
  dplyr::select(-c(InvRange2, Origin2))
  
summary(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + Cumulative.Import.Value.T, data=to_analyze3, family="poisson"))
car::Anova(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + Cumulative.Import.Value.T, data=to_analyze3, family="poisson"), type="III")

# INFLOR database

# add trees# - ADD EUROPE/AUSTRALIA
exotic_trees <- read.csv("Covariates_Input/exotic_trees.csv")[,-1]
exotic_trees_hi<-read.csv("Covariates_Input/HI_exotic_trees.csv")[,-1]
exotic_trees_au_eu<-read.csv("Covariates_Input/AU_EU_exotic_trees.csv")[,-1] %>%
  mutate(Invaded.Range.Cluster = gsub("Europe","Eurasian.Palearctic",Invaded.Range.Cluster))

to_analyze4 <- to_analyze3 %>%
  left_join(
    rbind(exotic_trees,exotic_trees_hi) %>%
      dplyr::rename(Invaded.Range.Cluster=Invaded.Range.in.North.America.Cluster) %>%
      rbind(exotic_trees_au_eu)%>%
      mutate(Geographic.Origin.Cluster=sub("Carribbean","Caribbean",Geographic.Origin.Cluster)))

summary(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + Cumulative.Import.Value.T +number_of_exotic_trees, data=to_analyze4, family="poisson"))
car::Anova(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + Cumulative.Import.Value.T+number_of_exotic_trees, data=to_analyze4, family="poisson"), type="III")

windows(6,6)
chart.Correlation(to_analyze4[,3:8], histogram=TRUE, pch=19)

####### with transformations
#######

summary(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + log(Cumulative.Import.Value.T) +log(number_of_exotic_trees+.01), data=to_analyze4, family="poisson"))
car::Anova(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + log(Cumulative.Import.Value.T) +log(number_of_exotic_trees+.01), data=to_analyze4, family="poisson"), type="II")
summary(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + log(Cumulative.Import.Value.T) +log(number_of_exotic_trees+.01)+
              NMDS.Distance:log(number_of_exotic_trees+.01), data=to_analyze4, family="poisson"), type="II")

to_analyze <- to_analyze4 #%>%
#  left_join(extrameta %>% select(X,area), by=join_by("Geographic.Origin.Cluster"=="X")) %>% rename(area_source=area) %>%
#  left_join(extrameta %>% select(X,area), by=join_by("Invaded.Range.in.North.America.Cluster"=="X")) %>% rename(area_sink=area)

##
write.csv(to_analyze, "Final_output/Final_analysis_data_pathogens-March4-2025-pnas.csv", row.names=F)
save(to_analyze, file="Final_output/Final_analysis_data_pathogens-March4-2025-pnas.RData")

############################################################
############################################################
