setwd("/home/ec2-user/Retrospective/")

load("Host_list/Output/global_tree_phylogeny_pruned.RData")

load("Host_list/Output/nodes_final.Rdata")

list1_bgci <- read.csv("Host_list/Other_DB/Sarina_BGCI_species_country_df.csv")

# figure out whats missing from the tree
missing_countries <- list1_bgci %>% filter(country_name %in% c("Bolivia", "Tanzania", "Iran", "North Korea", "South Korea", "Taiwan"))
missing_species <- missing_countries %>% filter(!(species %in% tr3$species.list$species))


library(WorldFlora)

WFO.remember("WFO_Backbone/classification.csv")

# first we need a tree of all the tree species on the list - from USA plant list - and international

TranslationTableWFO <- read.csv("Host_list/Output/host_translation_table_all_new_wfo.csv")

# before synonymizing tip lable info from guide tree,
# substitute/fix error where only sub was used instead of gsub in the translation table search
TranslationTableWFO$spec.name.ORIG  <- gsub("_"," ",TranslationTableWFO$spec.name.ORIG)

#tocheck <- which(is.na(TranslationTableWFO$scientificName))

# if this is not done, errors are produced in joining steps
TranslationTableWFO<- distinct(TranslationTableWFO)
TranslationTableWFO <- TranslationTableWFO %>% mutate(scientificName = replace(scientificName, is.na(scientificName), spec.name.ORIG[is.na(scientificName)]))


# synonymize

missing_synonyms <- data.frame(spec.name.ORIG = missing_species$species) %>%
  left_join(select(TranslationTableWFO, c("scientificName","spec.name.ORIG")), by="spec.name.ORIG")%>%
  filter(is.na(scientificName))

xx<-WFO.match(spec.data=missing_synonyms$spec.name.ORIG, WFO.data=WFO.data)
yy<-WFO.one(xx)

write.csv(yy, "Host_list/Output/TranslationTable_WFO/missing_spp_5-12-2023.csv")

missing_translated <- missing_species %>%
  rename(spec.name.ORIG = species) %>%
  left_join(distinct(select(rbind(TranslationTableWFO, yy), c("scientificName","spec.name.ORIG"))), by="spec.name.ORIG")

### add to tree

sum(missing_translated$scientificName %in% (tips.info.WP$species %>% gsub("_"," ",.)))

missing_translated[which(missing_translated$scientificName %in% (tips.info.WP$species %>% gsub("_"," ",.))),]
