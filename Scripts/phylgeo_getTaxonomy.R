library(dplyr)
library(foreach)
library(reshape2)
source("Scripts/Script1_Sentinels_functions_2024_reduced_memory.R")

#if("pest_taxonomy.csv" %in% dir("Data")) {pest_taxonomy <- read.csv("Data/pest_taxonomy.csv")} else pest_taxonomy <- data.frame(name=NULL)

unique_pest_list_and_full <-
  read.csv("Pathogens/Input/Expanded_data_EU_AU_NA_HI_Updated_Feb2025.csv") #%>%
  
unique_pest_list_origins <- unique_pest_list_and_full %>%
  dplyr::select(Pathogen.Species, Geographic.Origin.Cluster) %>% distinct
  
unique_pest_list_destinations <- unique_pest_list_and_full %>%
  dplyr::select(Pathogen.Species, Invaded.Range.Cluster) %>% distinct %>%
  reshape2::dcast(Pathogen.Species~Invaded.Range.Cluster, length) %>% arrange(Pathogen.Species)

#unique_pest_list<-
#  unique_pest_list_and_full %>%
#  dplyr::select(Pathogen.Species) %>% unlist %>% as.character %>% unique %>% sort
##levels(as.factor(unique_pest_list_and_full$Pathogen.Species))

loqfalta <- #setdiff(unique_pest_list, pest_taxonomy$pest_sp_full_accepted) %>% na.omit

  levels(as.factor(unique_pest_list_and_full$Pathogen.Species))

#setdiff(unique_pest_list, (pest_taxonomy %>% filter(!(is.na(pest_genus) & is.na(pest_kingdom))))$pest_sp_full_accepted)

#loqfalta <- #setdiff(unique_pest_list, (pest_taxonomy %>% filter(!is.na(pest_kingdom)))$pest_sp_full_accepted) %>% na.omit

#for (i in 1:((length(loqfalta)%/%100)+1)) {
  
#for (i in loqfalta) {

#  a<- ((i-1)*100+1)
#  b<- min((i)*100, length(loqfalta))
  pestlist_chunk <- loqfalta#[a:b]
  #pestlist_chunk <- loqfalta#[a:b]
  
  taxchunk <- tax_table_GBIF(pestlist_chunk)
  
  taxchunk$pest_kingdom <- unlist(taxchunk$pest_kingdom)
  taxchunk$pest_phylum  <- unlist(taxchunk$pest_phylum)
  taxchunk$pest_class   <- unlist(taxchunk$pest_class)
  taxchunk$pest_order   <- unlist(taxchunk$pest_order)
  taxchunk$pest_family  <- unlist(taxchunk$pest_family)
  taxchunk$pest_genus   <- unlist(taxchunk$pest_genus)
  
  pest_taxonomy <- #rbind(pest_taxonomy, 
    cbind(taxchunk, name=pestlist_chunk)#)
  
#  print(taxchunk)
#  print(paste(a, "to", b, "of", length(loqfalta)))
#}

#str(pest_taxonomy)
#names(pest_taxonomy) <- c('pest_kingdom', 'pest_phylum' , 'pest_class'  ,'pest_order','pest_family' ,'pest_genus', 'pest_sp_full_accepted')
rownames(pest_taxonomy) <- NULL

pest_taxonomy %>% dplyr::rename(pest_sp_full_accepted = name) %>%
  dplyr::select(pest_sp_full_accepted,
         pest_kingdom, pest_phylum , pest_class  ,pest_order,pest_family ,pest_genus, pest_sp_full_accepted) %>%
  arrange(pest_kingdom, pest_phylum , pest_class  ,pest_order,pest_family ,pest_genus, pest_sp_full_accepted) %>%
  dplyr::rename(Pathogen.Species = pest_sp_full_accepted) %>%
  left_join(unique_pest_list_origins) %>%
  left_join(unique_pest_list_destinations) %>%
  write.csv("Final_output/pest_list_taxonomy_updatedMarch2025.csv", row.names=F)
