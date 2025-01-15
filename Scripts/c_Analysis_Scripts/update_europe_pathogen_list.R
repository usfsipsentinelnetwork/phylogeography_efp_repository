# update europe
library('dplyr')
source("Scripts/Script1_Sentinels_functions_2023.1.R")

setwd('Pathogens/Input')

working_list <- read.csv("Expanded_data_EU_AU_NA_HI_7_8_2024_corrections.csv")

working_list_eu <- working_list %>% filter(Invaded.Range.Cluster == "Eurasian.Palearctic")

new_list <- read.csv("Invasive pathogens Europe 2024.csv") %>%
  
  dplyr::mutate(
    Species_short = 
      
      gsub(
        "^([A-Z][a-z]+ [a-z\\-]+)( [a-z.]+ [a-z\\-]+)?.*$",#( [\\(A-Z].*)( ssp\\. [a-z]+)?$",
        "\\1\\2",
        Species
      )
      
  )

new_list$Species_short[grep("Abad", new_list$Species_short)] <- "Phytophthora niederhauseri"

new_list$Species_short

# diff

shared <- intersect(new_list$Species_short, working_list_eu$Pathogen.Species)
additionals <- setdiff(new_list$Species_short, working_list_eu$Pathogen.Species)

shared

additionals

new_data <- filter(new_list, Species_short %in% additionals) %>%
  filter(Native.range != "Unknown")

filter(new_data, Species_short %in% working_list$Pathogen.Species)

##

new_data_to_add <-
  select(new_data, Species_short, Native.range) %>%
  dplyr::rename(Pathogen.Species = Species_short)

new_data_to_add$Invaded.Range.Cluster <- "Eurasian.Palearctic"

## already have data

already_in_list <- new_data_to_add$Pathogen.Species[new_data_to_add$Pathogen.Species %in% working_list$Pathogen.Species]
already_in_list_origins <- working_list %>% dplyr::filter(Pathogen.Species %in% already_in_list) %>% dplyr::select(-Invaded.Range.Cluster) %>% distinct

new_data_to_add <- 
  new_data_to_add %>% left_join(already_in_list_origins)

##

new_data_to_add$Native.range %>% unique
working_list$Geographic.Origin.Cluster %>% unique

## Europe

new_data_to_add$Geographic.Origin.Cluster[which(new_data_to_add$Native.range == "Europe")] <- "Eurasian.Palearctic"

## Central America

new_data_to_add$Geographic.Origin.Cluster[which(new_data_to_add$Native.range == "Central America")] <- "C.America"

## Australasia

new_data_to_add[which(new_data_to_add$Native.range == "Australasia"),]
new_data_to_add$Geographic.Origin.Cluster[which(new_data_to_add$Native.range == "Australasia")] <- "S.Cone.Pacific"

## Africa

new_data_to_add[which(new_data_to_add$Native.range == "Africa"),]
new_data_to_add$Geographic.Origin.Cluster[which(new_data_to_add$Native.range == "Africa")] <- "Subsaharan.Africa"

## South America
new_data_to_add$Geographic.Origin.Cluster[which(new_data_to_add$Native.range == "South America")] <- "Amazonia"

## Asia + North America

new_data_to_add<-

new_data_to_add %>%
  filter(!grepl("Ophiostoma novo-ulmi",Pathogen.Species)) %>%
  bind_rows(data.frame(Pathogen.Species = c("Ophiostoma ulmi","Ophiostoma novo-ulmi"),
                       Invaded.Range.Cluster = "Eurasian.Palearctic",
                       Geographic.Origin.Cluster = "SE.Asia"))

#new_data_to_add[which(new_data_to_add$Native.range == "Asia"),] 
## North America
#new_data_to_add[which(new_data_to_add$Native.range == "North America"),] 
#write.csv(
#  new_data_to_add[which(new_data_to_add$Native.range %in% c("Asia","North America")),],
#  "origins_to_determine_8-27-24.csv",
#  row.names=F
#)

# need to find which ones are also in north . america

new_data_to_add <- filter(new_data_to_add, !is.na(Geographic.Origin.Cluster))

#### CHECK ALREADY DOWNLOADED GBIF DATA

annotated.EU.Asia <- read.csv("origins_to_determine_8-27-24_filled.csv")

list.to.check.EU <-
  annotated.EU.Asia %>%
  filter(!grepl("North\\.America", Geographic.Origin.Cluster))%>%
  select(Pathogen.Species) %>% unlist %>% as.vector%>% unique 

list.to.check.EU.Asia <-
  annotated.EU.Asia %>%
  select(Pathogen.Species) %>% unlist %>% as.vector%>% unique 

FungalTraits.Pathogens <-
  read.csv("C:/Users/GeoffreyWilliams/OneDrive - USDA/Project_folders_Personal/Projects/Sentinels/Large_Files/FungalTraits 1.2_ver_16Dec_2020 - V.1.2.csv") %>%
  filter(grepl("[Pp]athogen", primary_lifestyle) |grepl("[Pp]athogen", Secondary_lifestyle) ) %>%
  select(GENUS) %>% unlist %>% as.vector %>% unique

ARS.Pathogens <-
  read.csv("C:/Users/GeoffreyWilliams/OneDrive - USDA/Project_folders_Personal/Projects/Sentinels/Large_Files/Fungus-Host_20211105.csv") %>%
  tidyr::separate(sciName, into=c("p1","p2","p3","p4","p5","p6","p7","p8"), sep = " ", remove=F) %>%
  filter(p1 %in% FungalTraits.Pathogens)

head(ARS.Pathogens)

GBIF.list.to.check.inNA.fromEU.Asia <- read.csv("GBIF_host_association.csv") %>% filter(countryCode %in% c('CA','US'))
GBIF.list.to.check.inHI.fromEU.Asia <- read.csv("GBIF_host_association.csv") %>% filter(stateProvince == "HI")
GBIF.list.to.check.inAU.fromEU.Asia <- read.csv("GBIF_host_association.csv") %>% filter(countryCode == "AU")

new.NA.invasives <- list.to.check.EU[(list.to.check.EU %in% unique(GBIF.list.to.check.inNA.fromEU.Asia$species) |
                                             list.to.check.EU %in% unique(filter(ARS.Pathogens, country %in% c("United States", "Canada"))$sciName)) &
                                            !(list.to.check.EU %in% filter(working_list, grepl("North\\.America", Invaded.Range.Cluster))$Pathogen.Species)]

new.HI.invasives <- list.to.check.EU.Asia[(list.to.check.EU.Asia %in% unique(GBIF.list.to.check.inHI.fromEU.Asia$species) |
                                             list.to.check.EU.Asia %in% unique(filter(ARS.Pathogens, stateProvince == "Hawaii")$sciName)) &
                                            !(list.to.check.EU.Asia %in% filter(working_list, grepl("Hawaii", Invaded.Range.Cluster))$Pathogen.Species)]

new.AU.invasives <- list.to.check.EU.Asia[(list.to.check.EU.Asia %in% unique(GBIF.list.to.check.inAU.fromEU.Asia$species) |
                                             list.to.check.EU.Asia %in% unique(filter(ARS.Pathogens, country == "Australia")$sciName)) &
                                            !(list.to.check.EU.Asia %in% filter(working_list, grepl("Australia", Invaded.Range.Cluster))$Pathogen.Species)]

list.to.check.other <- new_data_to_add %>% filter(!is.na(Geographic.Origin.Cluster))
#write.csv(list.to.check.other, "Input/Additional_Pathogens_to_GBIF.csv", row.names=F)

GBIF.list.to.check.inNA.fromEU.Asia2 <- read.csv("GBIF_host_association_additional_data_sep4.csv") %>% filter(countryCode %in% c('CA','US'))
GBIF.list.to.check.inHI.fromEU.Asia2 <- read.csv("GBIF_host_association_additional_data_sep4.csv") %>% filter(stateProvince == "HI")
GBIF.list.to.check.inAU.fromEU.Asia2 <- read.csv("GBIF_host_association_additional_data_sep4.csv") %>% filter(countryCode == "AU")

new.NA.invasives <- unique(c(new.NA.invasives,
                      list.to.check.other$Pathogen.Species[( list.to.check.other$Pathogen.Species %in% unique(GBIF.list.to.check.inNA.fromEU.Asia2) |
                        list.to.check.other$Pathogen.Species %in% unique(filter(ARS.Pathogens, country %in% c("United States", "Canada"))$sciName)) &
                          !(list.to.check.other$Pathogen.Species %in% filter(working_list, grepl("North\\.America", Invaded.Range.Cluster))$Pathogen.Species)]))

new.HI.invasives <- unique(c(new.HI.invasives,
                             list.to.check.other$Pathogen.Species[( list.to.check.other$Pathogen.Species %in% unique(GBIF.list.to.check.inHI.fromEU.Asia2) |
                               list.to.check.other$Pathogen.Species %in% unique(filter(ARS.Pathogens, stateProvince == "Hawaii")$sciName)) &
                                 !(list.to.check.other$Pathogen.Species %in% filter(working_list, grepl("Hawaii", Invaded.Range.Cluster))$Pathogen.Species)]))

new.AU.invasives <- unique(c(new.AU.invasives,
                             list.to.check.other$Pathogen.Species[( list.to.check.other$Pathogen.Species %in% unique(GBIF.list.to.check.inAU.fromEU.Asia2) |
                               list.to.check.other$Pathogen.Species %in% unique(filter(ARS.Pathogens, country == "Australia")$sciName)) &
                                 !(list.to.check.other$Pathogen.Species %in% filter(working_list, grepl("Australia", Invaded.Range.Cluster))$Pathogen.Species)]))

# now have to figure out which regions in North America pathogens have invaded

allgbif <-
  bind_rows(
    read.csv("GBIF_host_association.csv"),
    read.csv("GBIF_host_association_additional_data_sep4.csv")
  )

# N
new.N.NA.invasives <-
  new.NA.invasives[
    new.NA.invasives %in% unique((filter(ARS.Pathogens, stateProvince %in% c(
      "Kansas",          "Canada, Saskatchewan",           "Alaska",
      "Nebraska",        "Canada, Alberta",                "Central States",
      "South Dakota",    "Canada, Nunavut",                "Central and western states",
      "North Dakota",    "Canada, Northwest Territories",  "Eastern and central states",
      "Manitoba",        "Canada, Yukon Territory",        "Northern Great Plains states",
      "Northern states"
    ))$sciName %>% c(filter(allgbif, stateProvince %in% c(
      "Kansas",          "Saskatchewan",           "Alaska",
      "Nebraska",        "Alberta",
      "South Dakota",    "Nunavut",
      "North Dakota",    "Northwest Territories",
      "Manitoba",        "Yukon Territory", "Yukon"
    ))$species)))]
# NE
new.NE.NA.invasives <-
  new.NA.invasives[
    new.NA.invasives %in% (filter(ARS.Pathogens, stateProvince %in% c(
      "Iowa",            "Canada, Quebec",       "West Virginia",              "Canada, Prince Edward Island",
      "Minnesota",       "Canada, Ontario",      "Virginia",                   "Canada, Nova Scotia",
      "Illinois",        "Canada, New Brunsick", "Maryland",                   "Great Lake states", 
      "Indiana",         "Canada, Newfoundland", "Eastern and central states", "Connecticut",
      "Michigan",        "Delaware",             "New Jersey",                 "Rhode Island",
      "Northern states", "Massachusetts",        "Vermont",                    "New Hampshire",
      "Maine",           "New York",             "Ohio",                       "Pennsylvania",
      "Wisconsin",       "Northeastern states",  "Eastern states",             "Central states",
      "District of Columbia"
    ))$sciName %>% c(filter(allgbif, stateProvince %in% c(
      "Iowa",            "Quebec",       "West Virginia", "Prince Edward Island", "Minnesota", "Ontario",      "Virginia",
      "Nova Scotia",     "Illinois",     "New Brunsick",  "Maryland",             "Indiana",   "Newfoundland", "Connecticut",
      "Michigan",        "Delaware",     "New Jersey",    "Rhode Island",         "Northern states", "Massachusetts", "Vermont",
      "New Hampshire",   "Maine",         "New York",     "Ohio",                 "Pennsylvania",    "Wisconsin",     "District of Columbia"
    ))$species))]
# SE
new.SE.NA.invasives <-
  new.NA.invasives[
    new.NA.invasives %in% (filter(ARS.Pathogens, stateProvince %in% c(
      "Kentucky",       "Alabama",     "Louisiana", "Southern states",
      "Tennessee",      "Mississippi", "Georgia",   "Southeastern states",
      "North Carolina", "Missouri",    "Florida",   "South central states",
      "South Carolina", "Arkansas",    "Oklahoma",  "Texas"
    ))$sciName %>% c(filter(allgbif, stateProvince %in% c(
      "Kentucky",       "Alabama",     "Louisiana",
      "Tennessee",      "Mississippi", "Georgia",
      "North Carolina", "Missouri",    "Florida",
      "South Carolina", "Arkansas",    "Oklahoma",  "Texas"
    ))$species))]
# W
new.SW.NA.invasives <-
  new.NA.invasives[
    new.NA.invasives %in% (filter(ARS.Pathogens, stateProvince %in% c(
      "Oregon",       "New Mexico", "Montana", "Southern states",
      "Washington",   "Arizona",    "Idaho",   "Southwestern states",
      "California",   "Colorado",   "Utah",    "Northwestern states",
      "Nevada",       "Wyoming",    "Western states",
      "Central and western states", "Canada, British Columbia"
    ))$sciName %>% c(filter(allgbif, stateProvince %in% c(
      "Oregon",       "New Mexico", "Montana",
      "Washington",   "Arizona",    "Idaho",
      "California",   "Colorado",   "Utah", 
      "Nevada",       "Wyoming",
      "British Columbia"
    ))$species))]

#### PUT IT ALL TOGETHER

updated.list.combined <-
  bind_rows(
    working_list,
    new_data_to_add[,c(1,3,4)] %>% filter(!is.na(Geographic.Origin.Cluster)),
    annotated.EU.Asia[,c(1,3,4)] %>% filter(!is.na(Geographic.Origin.Cluster))
  ) %>% #filter(Invaded.Range.Cluster != Geographic.Origin.Cluster) %>%
  filter(!is.na(Geographic.Origin.Cluster)) %>%
  distinct

#working_list
#new_data_to_add %>% filter(!is.na(Geographic.Origin.Cluster))
#annotated.EU.Asia[,1:4] %>% filter(!is.na(Geographic.Origin.Cluster))

updated.list.combined <-
  bind_rows(
    updated.list.combined,
    data.frame(Pathogen.Species = new.N.NA.invasives) %>% left_join(updated.list.combined) %>% mutate(Invaded.Range.Cluster="N.North.America"),
    data.frame(Pathogen.Species = new.NE.NA.invasives) %>% left_join(updated.list.combined)%>% mutate(Invaded.Range.Cluster="NE.North.America"),
    data.frame(Pathogen.Species = new.SE.NA.invasives) %>% left_join(updated.list.combined)%>% mutate(Invaded.Range.Cluster="SE.North.America"),
    data.frame(Pathogen.Species = new.W.NA.invasives) %>% left_join(updated.list.combined)%>% mutate(Invaded.Range.Cluster="W.North.America"),
    data.frame(Pathogen.Species = new.HI.invasives) %>% left_join(updated.list.combined)%>% mutate(Invaded.Range.Cluster="Hawaii"),
    data.frame(Pathogen.Species = new.AU.invasives) %>% left_join(updated.list.combined)%>% mutate(Invaded.Range.Cluster="Australia"),
  ) %>% filter(Invaded.Range.Cluster != Geographic.Origin.Cluster) %>%
  filter(!is.na(Geographic.Origin.Cluster)) %>%
  distinct

# quality control - which ones have multiple origins?

distinct(select(updated.list.combined, Pathogen.Species, Geographic.Origin.Cluster)) %>%
  (function (x) {
    x[x$Pathogen.Species %in% names(table(x$Pathogen.Species)[table(x$Pathogen.Species) > 1]),]
  })

duplicate_origins <-
  
  distinct(select(updated.list.combined, Pathogen.Species, Geographic.Origin.Cluster)) %>%
  (function (x) {
    x[x$Pathogen.Species %in% names(table(x$Pathogen.Species)[table(x$Pathogen.Species) > 1]),]
  }) %>%
  select(Pathogen.Species) %>% unlist %>% as.vector %>% unique

# quality control - potential synonyms

updated.list.combined %>% 
  select(Pathogen.Species) %>%
  distinct %>%
  tidyr::separate(Pathogen.Species, c('g','s'), " ", F) %>%
  (function (x) {
    x[x$s %in% names(table(x$s)[table(x$s) > 1]),]
  }) %>%
  arrange(s,g,Pathogen.Species)

# rename one thing


updated.list.combined[which(updated.list.combined$Pathogen.Species=="Stigmina thujina"),'Pathogen.Species'] <- "Pseudocercospora thujina"
updated.list.combined[which(updated.list.combined$Invaded.Range.Cluster=="S.Cone.Pacific"),'Invaded.Range.Cluster'] <- "Australia"

#updated.list.combined <- read.csv("Pathogens/Input/Expanded_data_EU_AU_NA_HI_Sept2024.csv")
#updated.list.combined <- updated.list.combined %>% filter(Geographic.Origin.Cluster != "Unknown")
#updated.list.combined[which(updated.list.combined$Geographic.Origin.Cluster=="N.Pacfic"),'Geographic.Origin.Cluster'] <- "N.Pacific"
#updated.list.combined[which(updated.list.combined$Invaded.Range.Cluster=="SW.North.America"),'Invaded.Range.Cluster'] <- "W.North.America"


#Phaeoseptoria eucalypti!=Teratosphaeria eucalypti
#Eutypella parasitica!=Cryphonectria parasitica

####
#### dont forget to filter out pathogens where origins and destination are the same

write.csv(dplyr::arrange(updated.list.combined), "Expanded_data_EU_AU_NA_HI_Sept2024.csv", row.names=F)
#write.csv(dplyr::arrange(updated.list.combined), "Pathogens/Input/Expanded_data_EU_AU_NA_HI_Sept2024.csv", row.names=F)

# NOTE THIS HAS NOT BEEN CORRECTED/UPDATED
write.csv(dplyr::arrange(updated.list.combined %>% filter(!(Pathogen.Species %in% duplicate_origins))), "Expanded_data_EU_AU_NA_HI_Sept2024_noduplicateorigins.csv", row.names=F)
