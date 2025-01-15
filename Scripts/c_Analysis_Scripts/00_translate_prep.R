# make phyloseq objects to calculate phylogenetic distances between regions and places, etc.

library(ape)
library(tidyr)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(reshape2)
library(foreach)
library(doParallel)
load("Host_list/Output/global_tree_phylogeny.RData")

melt.community <- function (x, idvars, speciesvar) {
  y <- NULL
  for (j in idvars) {
    y <- rbind(
      y,
      x %>%
        na.omit %>%
        dplyr::rename(k=j) %>%
        dplyr::rename(scientificName=speciesvar) %>%
        select(c(scientificName, k)) %>%
        filter(k) %>%
        dplyr::mutate(country_name=j) %>%
        select(-k) %>%
        as.data.frame
    )
  }
  y
}


#############################
# countries - for clustering

trees.countries <- read.csv("Host_list/Other_DB/Sarina_BGCI_species_country_df.csv")%>%distinct
tree_translation<-
  read.csv("Host_list/Output/host_translation_table_all_new_wfo.csv")%>%
  rbind(read.csv("Host_list/Output/missing_spp_5-12-2023.csv")[,-1]) %>%          #### added 5-16-2023
  select(spec.name,spec.name.ORIG,scientificName)%>%
  distinct

trees.countries.translated <- left_join(
  trees.countries,
  tree_translation %>% select(spec.name, spec.name.ORIG,scientificName),
  by = join_by("species" == "spec.name.ORIG"))

#checking
#dim(trees.countries.translated %>% distinct)
#dim(trees.countries)

# now modify trees.countries to separate Hawai'i and USDA Plants list for North America

usda.untranslated <- read.csv("Host_list/Output/host_trees_byregion_notsynonymized.csv")%>% distinct
#usda.translation.table <- read.csv("Host_list/Output/ars_north_america_search_synonym_table.csv")%>% distinct

usda.translated <-
  left_join(
    usda.untranslated ,#%>% select(host_name),
    tree_translation,
    by = join_by("host_name" == "spec.name.ORIG")
  )

#dim(usda.untranslated)
#dim(usda.translated)

hawaii.untranslated <- read.csv("Host_list/Other_DB/Tropicos_Hawaii_Names.csv")%>%distinct

hawaii.translated <-
  left_join(
    hawaii.untranslated,
    tree_translation,
    by = join_by("species" == "spec.name.ORIG")
  )

hawaii.translated <- na.omit(hawaii.translated)

hawaii.translated <-
  hawaii.translated %>%
  filter(
    scientificName %in%
      (trees.countries.translated %>%
         filter(country_name=="United States"))$scientificName) %>%
  select(scientificName) %>%
  mutate(country_name="Hawaii")

# need to remove United states
#
# and add the regional lists
#
# plus Hawai'i, separately as Hawai'i
#######

#na.omit(usda.translated)[-which(names(usda.translated)%in%c('host','host_name','host_genus','host_authority','spec.name'))]

usda.translated.melted <-
  usda.translated %>%
  na.omit %>%
  select(-c(host,host_name,host_genus,host_authority,spec.name))%>%
  (function(x)
    melt.community(
      x,
      speciesvar="scientificName",
      idvars=names(x)[which(!(names(x) == "scientificName"))]))

usda.stateprovince.untranslated <- read.csv("Host_list/Output/host_trees_bystateprovince_notsynonymized.csv")

usda.stateprovince.translated <-
  left_join(
    usda.stateprovince.untranslated,
    tree_translation,
    by = join_by("host_name" == "spec.name.ORIG")
  ) %>%
  select(-c(host,host_name,host_genus,host_authority, spec.name)) %>%
  na.omit %>%
  (function(x)
    melt.community(
      x,
      speciesvar="scientificName",
      idvars=names(x)[which(!(names(x) == "scientificName"))]))

# can try this either with province/state- or region-level for North America

trees.countries.translated <-
  trees.countries.translated %>%
  filter(!(country_name %in% c("Canada","United States"))) %>%
  select(scientificName, country_name) %>%
  #  rbind(usda.translated.melted, hawaii.translated) # region level
  rbind(hawaii.translated, usda.stateprovince.translated)  %>%
  distinct# state/province-level

#write.csv(trees.countries.translated, "Host_list/Output/trees.countries.all.translated.csv", row.names=F)
