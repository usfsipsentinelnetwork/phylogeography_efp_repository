library(dplyr)
library(tidyr)

load("Host_list/Output/ordinated_unifrac_data.Rdata")
load("Host_list/Output/countries_unifrac_clusters1.RData")

trees.countries <- read.csv('Host_list/Output/trees.countries.all.translated.csv') %>%
  left_join(bioregions_vector)

glonaf_taxa<- read.csv("Host_list/Input/Taxon_x_List_GloNAF_vanKleunenetal2018Ecology_nonunicode.csv")

glonaf_regions <- read.csv("Host_list/Input/Region_GloNAF_vanKleunenetal2018Ecology.csv")

australia <- glonaf_taxa %>%
  filter(region_id %in% filter(glonaf_regions, country=="Australia")$region_id) %>%
  select(standardized_name) %>% unlist %>% unique

europe <- glonaf_taxa %>%
  filter(region_id %in% filter(glonaf_regions, tdwg1_name=="Europe")$region_id) %>%
  select(standardized_name)  %>% unlist %>% unique

exotics <- rbind(cbind(country="Australia", sciname=australia), cbind(country="Europe",sciname=europe)) %>%
  as.data.frame

exotic_trees <- expand.grid(c("Europe","Australia"), names(obj) %>% setdiff(c("Hawaii")))
names(exotic_trees)<-c("Invaded.Range.Cluster","Geographic.Origin.Cluster")
exotic_trees$number_of_exotic_trees <- 0

#rownames(exotichost_occurrence)
#df$scientificName

#unique(rownames(exotichost_occurrence))
#unique(df$scientificName)

levels(trees.countries$biogeographic_region)
levels(exotic_trees$Geographic.Origin.Cluster)

trees.countries <- trees.countries %>% filter(biogeographic_region != "Hawaii")
trees.countries$biogeographic_region <- factor(trees.countries$biogeographic_region)
levels(trees.countries$biogeographic_region)
levels(exotic_trees$Geographic.Origin.Cluster)

#write.csv(exotics %>% filter(sciname %in% trees.countries$scientificName), 'Covariates_input/introduced_spp_AU_EU.csv', row.names=F)

north.america.exotic.trees <- glonaf_taxa %>%
  filter(region_id %in% filter(glonaf_regions, country %in% c("United States", "Canada"))$region_id) %>%
  select(standardized_name) %>% unlist %>% unique %>%
  intersect(trees.countries$scientificName) # glonaf only counts 118 (we count 168)

# we need to add to this script to also have a list of all the non-native tree species in all the regions
# that are native to Europe and Australia

trees.list1 <- list()

for (i in 1:dim(exotic_trees)[1]) {
  
  exotic_trees$number_of_exotic_trees[i] <- sum(
    unique(exotics[exotics$country==exotic_trees$Invaded.Range.Cluster[i], 'sciname']) %in%
      unique(trees.countries[trees.countries$biogeographic_region==exotic_trees$Geographic.Origin.Cluster[i], 'scientificName'])
  )
    
  current.origin <- exotic_trees[i,'Geographic.Origin.Cluster']
  
  if(!(current.origin %in% names(trees.list1))) {
    trees.list1[[length(trees.list1) + 1]] <- character()
    names(trees.list1)[[length(trees.list1)]] <- as.character(current.origin)
  }
  
  #add any new trees to the list
  trees.list1[[exotic_trees[i,'Geographic.Origin.Cluster']]] <-
    
    # this list represents the exotic trees from the current Geographic.Origin.Cluster
    # and present in the current location (Europe or Australia) from the expanded grid
    intersect(
      unique(exotics[exotics$country==exotic_trees$Invaded.Range.Cluster[i], 'sciname']),
      unique(trees.countries[trees.countries$biogeographic_region==exotic_trees$Geographic.Origin.Cluster[i], 'scientificName'])
    ) %>%
    
    # now we only want to retain the ones that are new to the list, and add them to the list
    union(trees.list1[[exotic_trees[i,'Geographic.Origin.Cluster']]]) %>% unique()

    #  and then later we can sum them all up
}

write.csv(exotic_trees, "Covariates_Input/AU_EU_exotic_trees.csv")

save(trees.list1, file="Covariates_input/AU_EU_trees_list.RData")
