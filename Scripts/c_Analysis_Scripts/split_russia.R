# Split Russia - this script was used to experiment with the analysis with and without
# Russia being split based on E vs. W - in the final analysis we didn't use it because
# We split based on similarity to eastern bordering countries (like China) and western
# bordering countries. In the analysis, E and W Russia clustered together. Therefore
# it would have been necessary to split China into tropical and temperate component in
# order to complete this in a sensible way... not going to do that...

library(dplyr)

trees_countries <- read.csv("trees.countries.all.translated.csv")

ordinated <- load("C:/Users/GeoffreyWilliams/OneDrive - USDA/Project_folders_Personal/Projects/Sentinels/Biogeographical_Patterns_Retro2/Host_list/Output/ordinated_unifrac_data.Rdata")

bioregions_vector$country_name %>% setdiff(trees_countries$country_name)
trees_countries$country_name %>% setdiff(bioregions_vector$country_name)

sort(filter(trees_countries,country_name=="China")$scientificName) %>%
  grep("Larix",.,value=T)

west <- trees_countries %>%
  filter(
    country_name %in%
      setdiff(
        filter(
          bioregions_vector,
          biogeographic_region %in% c(
            "Eurasian.Palearctic",
            "N.North.America"
          )
        )$country_name,
        grep("Russia", bioregions_vector$country_name, value=T)
      )) %>%
  select(scientificName) %>%
  unlist %>%
  as.vector %>%
  unique %>%
  sort

east <- trees_countries %>%
  filter(
    country_name %in%
      setdiff(
        filter(
          bioregions_vector,
          biogeographic_region %in% c(
            "SE.Asia",
            "N.Pacific"
          )
        )$country_name,
        grep("Russia", bioregions_vector$country_name, value=T)
      )) %>%
  select(scientificName) %>%
  unlist %>%
  as.vector %>%
  unique %>%
  sort

russia.w <-
  intersect(
    filter(trees_countries,
           grepl("Russia", country_name))$scientificName,
    west
  ) %>% sort

russia.e<-
  intersect(
    filter(trees_countries,
           grepl("Russia", country_name))$scientificName,
    east
  ) %>% sort

russia.neither <-
  setdiff(
    filter(trees_countries,
           grepl("Russia", country_name))$scientificName,
    union(russia.w, russia.e)
  ) %>% sort

russia.neither

trees_countries2 <-
  trees_countries %>%
  filter(!grepl("Russia", country_name)) %>%
  rbind(data.frame(
    country_name = "W Russia",
    scientificName = russia.w
  )) %>%
  rbind(data.frame(
    country_name = "E Russia",
    scientificName = russia.e
  )) %>%
  arrange(
    country_name,
    scientificName
  )

save(russia.e, russia.w, russia.neither, trees_countries2,file= "new.tree.country.data.RData")

write.csv(trees_countries2, "trees.countries.all.translated.split.csv", row.names=F)
