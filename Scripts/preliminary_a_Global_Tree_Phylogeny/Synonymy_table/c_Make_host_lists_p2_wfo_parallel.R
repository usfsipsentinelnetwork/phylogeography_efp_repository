setwd("/home/ec2-user/Retrospective/")

#source("Scripts/Script1_Sentinels_functions_2023.R")

library(doParallel)
library(dplyr)
library(tidyr)
#library(V.PhyloMaker2)
#library(phyloseq)
library(WorldFlora)
WFO.remember("/home/ec2-user/Retrospective/WFO_Backbone/classification.csv")

#setwd("Host_list/USDA_DB/")

#na_plantlist <- read.csv( "host_trees_byregion_notsynonymized.csv") %>% dplyr::rename(host_full=host) %>% dplyr::rename(host=host_name)

## = commented out 3/27/23
##
#> names(us_plantlist)
#[1] "host"                 "host_name"            "host_genus"           "host_authority"
##
##setwd("../Other_DB/")
##
##### all names used in datasets need to be standardize to the same taxonomy
##### Extract names, add other relevant host names,
##
# 1 tree species names from global plant list (BGCI)
#> names(list1_bgci)
#[1] "country_name" "species"
##list1_bgci <- read.csv("Sarina_BGCI_species_country_df.csv") %>% dplyr::rename(host=species)
##
# 2 species names for V.PhyloPlot2
#> names(list2_vphyloplot)
#[1] "group"   "species" "genus"   "family" 
#data("tips.info.LCVP")
#data("tips.info.TPL")
##data("tips.info.WP")
##
#### INCLUDES NON-TREES - add later?
##list2_vphyloplot <-tips.info.WP %>% filter(group=="seedplant") %>% dplyr::mutate(host= sub("_", " ", species))
##
# 3 tree species names from **Ossala (GUTI) * has a year associated
#> names(list3_guti)
#[1] "GUTI.binomial"                "Country" ...
##list3_guti <- read.csv("Ossala_urban_forest_presence_absence_bycountry.csv") %>% dplyr::rename(host=GUTI.binomial)
## 
# 4 global alien flora database (GLONAF) * has years associated
## this one has both species names and species names with authorities
###> names(list4_glonaf)
###[1] "TaxonName"        "scientificName"   "Family" ... "LifeForm" ...
##
#### INCLUDES NON-TREES - add later?
##list4_glonaf <-
##  read.csv("GlobalAlienSpeciesFirstRecordDatabase_v2_reduced.csv") %>%
##  filter(LifeForm=="Vascular plants") %>%
##  dplyr::rename(host=TaxonName)
##
# 5 cabi * (some) years associated
##
#setdiff(names(read.csv("DataExtract20211203_Final_Host_Distributions_CABI.csv")),
#              names(read.csv("USForestryServiceDataExtract20221129_Final.csv")))
##
#setdiff(names(read.csv("USForestryServiceDataExtract20221129_Final.csv")),
#              names(read.csv("DataExtract20211203_Final_Host_Distributions_CABI.csv")))
##
##list5_cabi_a <-
##  read.csv("DataExtract20211203_Final_Host_Distributions_CABI.csv") %>%
##  dplyr::rename(species=sp) %>%
##  select(-c(Family, Genus, Sect, speth, authority)) %>%
##  rbind(
##    read.csv("USForestryServiceDataExtract20221129_Final.csv") %>%
##      dplyr::rename(species=ScientificName) %>%
##      dplyr::rename(LocationGeonameID = LocationGeonameId) %>%
##      dplyr::rename(MainlandStatus = MainlandStatusName)) %>%
##  dplyr::rename(host=species)
##
# 6 Spaulding -- ** incomplete ** (some) years associated
##
##list6_spaulding <-
##  read.csv("Spaulding_GBIF_othersources_host_chronology_incomplete.csv") %>%
##  filter(grepl("Spaulding", citation))
##  
# 7 tree species names used in ARS database
### also getting distributions from fungal-host associtions ** some years associated **
##
#### INCLUDES NON-TREES - add later?
##list7_ars <- read.csv("Fungus-Host_20211105.csv")
##
##### all unique names used
# just including known tree names for now keeps list length down from 190K++ to ~57.5K
##
##master_names_list_hosts_all_novplot <-
##  unique(c(na_plantlist$host_name,
##    list1_bgci$host,
#    list2_vphyloplot$host,
##    list3_guti$host,
##    list4_glonaf$host,
##    list5_cabi_a$host,
##    list6_spaulding$host,
##    list7_ars$host
##  )) %>% sort
##
##master_names_list_hosts_all <-
##  unique(c(na_plantlist$host_name,
##           list1_bgci$host,
##           list2_vphyloplot$host,
##           list3_guti$host,
##           list4_glonaf$host,
##           list5_cabi_a$host,
##           list6_spaulding$host,
##           list7_ars$host
##  )) %>% sort
##
##master_names_list_hosts_all_novplot[!(master_names_list_hosts_all_novplot %in% list2_vphyloplot$host)]
##
##setwd("../..")
##
##master_names_list_hosts_all %>% write.csv("Host_list/Output/names_used_all_sources.csv", row.names=F)
##
##master_names_list_hosts_all_novplot %>% write.csv("Host_list/Output/names_used_all_sources_novphyloplot.csv", row.names=F)
##
##
###### 
##### make a translation table (to be used for searching)
###working_trans_table <- NULL
###write.csv(read.csv("Host_list/Output/host_translation_table_all_new_wfo.csv"), "Host_list/Output/host_translation_table_all_working_wfo.csv", row.names=F)
##
#### correct earlier table
###working_trans_table <- read.csv("Host_list/Output/host_translation_table_all_working_wfo.csv")
###names_new <- master_names_list_hosts_all
###table1 <- data.frame()
###chunksize <- 500
##
###working_trans_table$host_sp_name_used <- ""
###length(names_new)
##
###for (i in 0:(dim(working_trans_table)[1] %/% chunksize)) {
###  start <- i*chunksize + 1
###  if (i*chunksize ==dim(working_trans_table)[1]) break
###  if (i == dim(working_trans_table)[1] %/% chunksize) {
###    stop <- start + dim(working_trans_table)[1] %% chunksize
###  } else {
###    stop <- start + chunksize-1
###  }
###  working_trans_table$host_sp_name_used[start:stop] <- names_new[start:stop]
###}
##
###with(working_trans_table, sum(spec.name.ORIG!=host_sp_name_used))
##
###with(working_trans_table[-which(is.na(working_trans_table$scientificName)),], sum(scientificName!=host_sp_name_used))
##
###

master_names_list_hosts_all <-
  read.csv("Host_list/Output/names_used_all_sources.csv")$x

working_trans_table <-
  read.csv("Host_list/Output/host_translation_table_all_working_wfo.csv")

###names_new <- master_names_list_hosts_all_novplot[which(!(master_names_list_hosts_all_novplot %in% unique(working_trans_table$spec.name.ORIG)))]

names_new <-
  master_names_list_hosts_all[which(!(master_names_list_hosts_all %in% unique(working_trans_table$spec.name.ORIG)))]

registerDoParallel(cores=63)

table1 <- working_trans_table #data.frame()
chunksize <- 1000

h <- 0:(length(names_new) %/% chunksize)

lists <- foreach(a=h, .combine=rbind, .packages=c('tidyr', 'WorldFlora')) %dopar% {
  start <- a*chunksize + 1
  if (a*chunksize == length(names_new)) break
  if (a == length(names_new) %/% chunksize) {stop <- length(names_new) } else {stop <- start + chunksize - 1}
  xx<-WFO.match(spec.data=names_new[start:stop], WFO.data=WFO.data)
  yy<-WFO.one(xx)
  filename <- paste("Host_list/Output/TranslationTable_WFO/", "names_", start, stop, sep="")
  write.csv(yy, filename, row.names=F)
  filename
}

setwd("Host_list/Output/TranslationTable_WFO/")

#lists <- dir()

for (x in lists) {
  table1 <- rbind(table1, read.csv(x))
}

setwd("../../..")

filenames_all <- data.frame()

for (y in h) {
  start <- y*chunksize + 1
  if (y == length(names_new) %/% chunksize) {stop <- length(names_new) } else {stop <- start + chunksize - 1}
  filenames_all <- rbind(filenames_all, data.frame(name=paste("Host_list/Output/TranslationTable_WFO/", "names_", start, stop, sep="")))
}

missing <- setdiff(filenames_all, lists)

#table1 <- rbind(table1, table_3)

table1 <- arrange(table1, 'spec.name', 'spec.name.ORIG')

write.csv(table1, "Host_list/Output/host_translation_table_all_new_wfo.csv", row.names=F)

print(paste("missing:", missing))

###

# incorporate annotations - if needed

###

# produce translated versions #