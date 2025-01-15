# make phyloseq objects to calculate phylogenetic distances between regions and places, etc.

setwd("C:/Users/GeoffreyWilliams/OneDrive - USDA/Project_folders_Personal/Projects/Sentinels/Retrospective")

library(ape)
#library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(reshape2)
library(foreach)
library(doParallel)
load("Host_list/Output/global_tree_phylogeny.RData")

# first part taken from phyloseq_unifrac_2.R version 6/13/23 ################################
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



#dim(hawaii.untranslated)
#dim(na.omit(hawaii.translated))

hawaii.translated <- na.omit(hawaii.translated)

#hawaii.translated %>% filter(scientificName %in% trees.countries.translated$scientificName)

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

#######################

tree_translation2<-
  read.csv("Host_list/Output/host_translation_table_all_new_wfo.csv")%>%
  rbind(read.csv("Host_list/Output/missing_spp_5-12-2023.csv")[,-1]) %>%          #### added 5-16-2023
  select(spec.name,spec.name.ORIG,scientificName, family, genus, infraspecificEpithet)%>%
  distinct

trees.countries.translated <- left_join(
  trees.countries,
  tree_translation2 %>% select(spec.name,spec.name.ORIG,scientificName, family, genus, infraspecificEpithet),
  by = join_by("species" == "spec.name.ORIG"))

prune_families <-
trees.countries.translated %>%
  select(scientificName, family) %>%
  filter(scientificName %in% gsub("_"," ",tr3$tree.scenario.3$tip.label)) %>%
  group_by(family) %>%
  dplyr::summarise(
    sp_rep = (function (x) {
      x[1] %>% gsub(" ", "_", .)
    })(scientificName)
    ) %>% arrange(family) %>% as.data.frame %>% na.omit

prune_genera <-
  trees.countries.translated %>%
  select(scientificName, genus) %>%
  filter(scientificName %in% gsub("_"," ",tr3$tree.scenario.3$tip.label)) %>%
  group_by(genus) %>%
  dplyr::summarise(
    sp_rep = (function (x) {
      x[1] %>% gsub(" ", "_", .)
    })(scientificName)
  ) %>% arrange(genus) %>% as.data.frame %>% na.omit

tax_table_GBIF_single <- function(x) {
  m <- name_backbone(x)
  g <- names (m)
  tax <- c(kingdom=NA, phylum=NA, class=NA, order=NA, family=NA, genus=NA)
  if (m$matchType != "NONE") for (y in names(tax)) if (y %in% g) tax[y] <- m[1,y]
  tax
}

tax_table_GBIF <- function(x) {
  (foreach(a=x, .combine='rbind') %do% tax_table_GBIF_single(a)) %>% as.data.frame %>%
    dplyr::rename(pest_kingdom=kingdom, pest_phylum=phylum, pest_class=class, pest_order=order, pest_family=family, pest_genus=genus)
}

taxonomy_data <- tax_table_GBIF(prune_families$sp_rep)

taxonomy_data <- data.frame(pest_order = as.character(unlist(taxonomy_data$pest_order)))

prune_order <- prune_families %>%
  cbind(taxonomy_data) %>%
  group_by(pest_order) %>%
  dplyr::summarise(
    sp_rep = (function (x) {
      x[1] %>% gsub(" ", "_", .)
    })(sp_rep)
  ) %>% arrange(pest_order) %>% as.data.frame %>% na.omit

familytree <- keep.tip( tr3$tree.scenario.3, prune_families$sp_rep)
familytree$tip.label <- left_join(data.frame(sp_rep=familytree$tip.label), prune_families)$family

ordertree <- keep.tip(tr3$tree.scenario.3, prune_order$sp_rep)
ordertree$tip.label <- left_join(data.frame(sp_rep=ordertree$tip.label), prune_order)$pest_order

ordertree2 <- ordertree
ordertree2$edge.length[which(ordertree$edge[,2] %in% 1:length(ordertree$tip.label))] <-0

library(ggtree)
windows(); ggtree(ordertree2) + geom_tiplab() +
  xlim(0,400)+
  theme(
    panel.background = element_rect(fill = "transparent", colour=NA),
    plot.background = element_rect(fill = "transparent", colour=NA),
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank() #remove minor gridlines
    )

ggsave(filename="Host_list/Output/Tree_of_trees_orders_transparent.png",
       bg="transparent")

windows(8,12);par(bg=NA);plot(ordertree, cex=1)

dev.copy(png, 'Host_list/Output/Orders_tree.png')
dev.off()
  
plot(familytree)
