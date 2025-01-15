
library(dplyr)
library(doBy)
library(reshape2)
library(tidyverse)
library(sjPlot)
library(cowplot)

load("Host_list/Output/ordination_of_cluster.RData")
load("Host_list/Output/ordinated_unifrac_data.Rdata")
#gg1
##### analysis

pathogens <-
  read.csv("Pathogens/Input/Expanded_data_EU_AU_NA_HI_Sept2024.csv")%>%
  dplyr::select(c("Pathogen.Species", "Geographic.Origin.Cluster", "Invaded.Range.Cluster")) #%>%
  #filter(Geographic.Origin.Cluster != "Unknown")

pathogens$Invaded.Range.Cluster <- replace(pathogens$Invaded.Range.Cluster, pathogens$Invaded.Range.Cluster=="S.Cone.Pacific", "Australia")
  #  read.csv("Pathogens/Input/Table1_Pathogens_Established_NorthAmerica_analysis.2.7.24.csv")%>%
  #  dplyr::select(c("Pathogen.Species", "Geographic.Origin.Cluster", "Invaded.Range.in.North.America.Cluster")) %>%
 # filter(Geographic.Origin.Cluster != "") %>%
 # filter(Invaded.Range.in.North.America.Cluster != "")
  
pathogens <- pathogens %>%
  filter(Geographic.Origin.Cluster != Invaded.Range.Cluster)
  

pathogens$Geographic.Origin.Cluster<-as.factor(pathogens$Geographic.Origin.Cluster)

#path.sum <- summaryBy(Pathogen.Species ~ Geographic.Origin.Cluster + Invaded.Range.in.North.America.Cluster, data=pathogens, FUN=length)
#path.sum.wide <- dcast(path.sum, Invaded.Range.in.North.America.Cluster ~ Geographic.Origin.Cluster)


path.sum <- summaryBy(Pathogen.Species ~ Geographic.Origin.Cluster + Invaded.Range.Cluster, data=pathogens, FUN=length)
path.sum.wide <- dcast(path.sum, Invaded.Range.Cluster ~ Geographic.Origin.Cluster)


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

path.sum.wide <- path.sum.wide %>% select(-Invaded.Range.Cluster)
path.sum.wide %>%cbind(Invaded.Range.Cluster=row.names(path.sum.wide),.)%>%
  write.csv("Pathogens/Output/path.sum.expanded.csv", row.names=F)

#path.sum.wide <- path.sum.wide[,-1]

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
  dplyr::rename(Pathogen.Species.length=value) %>% na.omit

# are we considering movement within north america? i didn't think so...

# add distances
dat2.summary <- dat2.summary %>% mutate(biogeographic_region = sub("Carribbean","Caribbean",biogeographic_region))

dist_to_Nearctic2<-(dat2.summary[2:3] %>%
                      dist(upper=T) %>%
                      as.matrix) %>% as.data.frame

row.names(dist_to_Nearctic2) <- colnames(dist_to_Nearctic2) <- dat2.summary$biogeographic_region

dist_to_Nearctic2$Invaded.Range.Cluster <- row.names(dist_to_Nearctic2)

dist_to_Nearctic2_melt <- melt(as.data.frame(dist_to_Nearctic2), variable.name="Geographic.Origin.Cluster", value.name="NMDS.Distance")
dist_to_Nearctic2_melt$Geographic.Origin.Cluster <- as.character(dist_to_Nearctic2_melt$Geographic.Origin.Cluster)
to_analyze <- path.sum %>%
  left_join(
    dist_to_Nearctic2_melt,
    by = c("Invaded.Range.Cluster", "Geographic.Origin.Cluster")) %>%
  dplyr::rename(Number.of.pests=Pathogen.Species.length)

summary(glm(Number.of.pests ~ NMDS.Distance, data=to_analyze, family="poisson"))

#with(to_analyze, plot(NMDS.Distance, Number.of.pests+0.001))

#with(to_analyze, plot(NMDS.Distance, Number.of.pests+0.001))

library(progress)

pb <- progress_bar$new(total=dim(to_analyze)[1])
# add spatial autocorrelation covariate
# NOTE made a correction here that makes it mathematically correct!

# need to correct for other continents
# for example, iberia/n.africa, se.asia, and middle east are adjacent/border eurasian.palearctic

to_analyze$spatauto <- 0
for(i in 1:(dim(to_analyze)[1])) {
  
  #i<-1
  # make sure to skip Hawaii
  if (to_analyze[i,"Invaded.Range.Cluster"] %in% c("Hawaii","Australia","Eurasian.Palearctic")) {
    pb$tick()
    to_analyze[i,"spatauto"] <- 0
    next
  }
  
  other_NorthAmerica_cells <-
    setdiff(rownames(path.sum.wide) %>% setdiff(c("Hawaii","Australia","Eurasian.Palearctic")), # make sure to skip Hawaii
            to_analyze[i,"Invaded.Range.Cluster"])
  
  ## THIS USES NUMBER OF INVASIVE SPECIES IN AT LEAST ONE OTHER "NEIGHBORING" CELL - FROM THAT REGION
  
  # pathogens that are in this cell - current cluster - pathogen species in all other entries that are in the current cluster
  pathogens_from_same_origin_region <-
    pathogens %>%                      # search pathogen location combinations
    filter(                              # filter the ones from the same origin cluster
      as.character(Geographic.Origin.Cluster) ==
        as.character(to_analyze[i,"Geographic.Origin.Cluster"] %>% setdiff("Hawaii"))) %>%
    dplyr::select("Pathogen.Species") %>% unlist %>% unique ## and get those pathogen species - added unique!
  
  # unique(pathogens$Pathogen.Species) # now change to any pathogen species!
  
  # it should actually be the number of pathogens from that same origin region- that are present in other cells
  
  # now find out how many of those pathogens are present in at least one neighboring cell
  
  value <- 0
  for (p in pathogens_from_same_origin_region) {
    
    for (o in other_NorthAmerica_cells) {
      
      # which cells that contain the pathogen (p) are from other parts of north america (o)
      y <- which(pathogens$Pathogen.Species==p & pathogens$Invaded.Range.Cluster==o)
      
      if (length(y)) {
        
        # if at least one of the other parts of north america have the pathogen, add one to the value and go to the next!
        value <- value+1
        break
      }
      
    }
    
    # augment with Mexico/Central America - only gets here if haven't found one yet
    # not gonna add this right now
    #mex <- which(pathogens$Pathogen.Species==p & pathogens$Geographic.Origin.Cluster=="C.America")
    
  }
  pb$tick()
  to_analyze[i,"spatauto"] <- value
}
str(to_analyze)
# plot and analysis




## add other variables
## need to rederive

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
chart.Correlation(to_analyze2[,3:6], histogram=TRUE, pch=19)

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
  trade_aggregated %>% select(biogeographic_region, "Reporter Name", cumulativeT) %>%
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
  select(-c(InvRange2, Origin2))
  
## add hawaii
  
#  rbind(
#    trade_aggregated_hi %>% mutate(biogeographic_region = sub("Carribbean","Caribbean",biogeographic_region))%>% 
#      rbind(hawaiiNAm %>%
#              mutate(Adjusted=Total.1989.Present) %>%
#              dplyr::select(-Country) %>%
#              cbind(biogeographic_region=c("NE.North.America","W.North.America","SE.North.America","N.North.America"))) %>%
#    rename(Geographic.Origin.Cluster = biogeographic_region) %>%
#    mutate(Cumulative.Import.Value.T = Total.1989.Present/4/1000000000000) %>% # divide equally by four since there are four regions?
#    dplyr::select(c(Geographic.Origin.Cluster,Cumulative.Import.Value.T)) %>%
#    left_join(to_analyze2 %>% filter(Invaded.Range.in.North.America.Cluster == "Hawaii"), .) 
#  )

summary(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + Cumulative.Import.Value.T, data=to_analyze3, family="poisson"))
car::Anova(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + Cumulative.Import.Value.T, data=to_analyze3, family="poisson"), type="III")

summary(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + Cumulative.Import.Value.T + spatauto, data=to_analyze3, family="poisson"))
car::Anova(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + Cumulative.Import.Value.T + spatauto, data=to_analyze3, family="poisson"), type="III")

summary(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + Cumulative.Import.Value.T + spatauto:NorthAm, data=to_analyze3, family="poisson"))
car::Anova(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + Cumulative.Import.Value.T + spatauto:NorthAm, data=to_analyze3, family="poisson"), type="II")

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

summary(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + Cumulative.Import.Value.T +number_of_exotic_trees+ spatauto, data=to_analyze4, family="poisson"))
car::Anova(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + Cumulative.Import.Value.T +number_of_exotic_trees+ spatauto, data=to_analyze4, family="poisson"), type="III")

summary(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + Cumulative.Import.Value.T +number_of_exotic_trees+ spatauto:NorthAm, data=to_analyze4, family="poisson"))
car::Anova(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + Cumulative.Import.Value.T +number_of_exotic_trees+ spatauto:NorthAm, data=to_analyze4, family="poisson"), type="II")

summary(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + Cumulative.Import.Value.T +log(number_of_exotic_trees+.01)+ spatauto:NorthAm, data=to_analyze4, family="poisson"))
car::Anova(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + Cumulative.Import.Value.T +log(number_of_exotic_trees+.01)+ spatauto:NorthAm, data=to_analyze4, family="poisson"), type="II")

windows(6,6)
chart.Correlation(to_analyze4[,3:9], histogram=TRUE, pch=19)

#######
summary(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + log(Cumulative.Import.Value.T) +log(number_of_exotic_trees+.01)+ log(spatauto+.01):NorthAm, data=to_analyze4, family="poisson"))
car::Anova(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + log(Cumulative.Import.Value.T) +log(number_of_exotic_trees+.01)+ log(spatauto+.01):NorthAm, data=to_analyze4, family="poisson"), type="II")

summary(glm(Number.of.pests ~  Climate.Distance + log(Cumulative.Import.Value.T) +NMDS.Distance*log(number_of_exotic_trees+.01)+ log(spatauto+.01):NorthAm, data=to_analyze4, family="poisson"))
car::Anova(glm(Number.of.pests ~ Climate.Distance + log(Cumulative.Import.Value.T) +NMDS.Distance*log(number_of_exotic_trees+.01)+ log(spatauto+.01):NorthAm, data=to_analyze4, family="poisson"), type="II")
#######

summary(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + log(Cumulative.Import.Value.T) +log(number_of_exotic_trees+.01), data=to_analyze4, family="poisson"))
car::Anova(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance + log(Cumulative.Import.Value.T) +log(number_of_exotic_trees+.01), data=to_analyze4, family="poisson"), type="II")

summary(glm(Number.of.pests ~ NMDS.Distance +log(number_of_exotic_trees+.01), data=to_analyze4, family="poisson"))
car::Anova(glm(Number.of.pests ~ NMDS.Distance  +log(number_of_exotic_trees+.01), data=to_analyze4, family="poisson"), type="II")

summary(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance+log(number_of_exotic_trees+.01), data=to_analyze4, family="poisson"))
car::Anova(glm(Number.of.pests ~ NMDS.Distance  + Climate.Distance+log(number_of_exotic_trees+.01), data=to_analyze4, family="poisson"), type="II")

summary(glm(Number.of.pests ~ NMDS.Distance + log(Cumulative.Import.Value.T)+ Climate.Distance+log(number_of_exotic_trees+.01), data=to_analyze4, family="poisson"))
car::Anova(glm(Number.of.pests ~ NMDS.Distance  + log(Cumulative.Import.Value.T)+ Climate.Distance+log(number_of_exotic_trees+.01), data=to_analyze4, family="poisson"), type="II")

summary(glm(Number.of.pests ~ NMDS.Distance + Climate.Distance+log(number_of_exotic_trees+.01)+ log(spatauto+.01):NorthAm, data=to_analyze4, family="poisson"))
car::Anova(glm(Number.of.pests ~ NMDS.Distance  + Climate.Distance+log(number_of_exotic_trees+.01)+ log(spatauto+.01):NorthAm, data=to_analyze4, family="poisson"), type="II")


# add country area# - ADD EUROPE/AUSTRALIA

#extrameta <- read.csv('Covariates_Input/landarea_and_centroids.csv')%>% 
#  mutate(X=sub("Carribbean","Caribbean",X))

to_analyze <- to_analyze4 #%>%
#  left_join(extrameta %>% select(X,area), by=join_by("Geographic.Origin.Cluster"=="X")) %>% rename(area_source=area) %>%
#  left_join(extrameta %>% select(X,area), by=join_by("Invaded.Range.in.North.America.Cluster"=="X")) %>% rename(area_sink=area)

##
write.csv(to_analyze, "Final_output/Final_analysis_data_pathogens-Sept-2024.csv", row.names=F)
save(to_analyze, file="Final_output/Final_analysis_data_pathogens-Sept-2024.RData")

############################################################
############################################################
