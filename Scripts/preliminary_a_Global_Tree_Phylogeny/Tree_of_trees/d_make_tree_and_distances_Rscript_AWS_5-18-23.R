setwd("/home/ec2-user/Retrospective/")

#source("Scripts/Script1_Sentinels_functions_2023.R")

#library(doParallel)
library(dplyr)
library(tidyr)
library(V.PhyloMaker2)
library(phyloseq)
library(WorldFlora)

source("Scripts/phylo.maker.verbose.R")
#library(ggtree)

# first we need a tree of all the tree species on the list - from USA plant list - and international

TranslationTableWFO <- read.csv("Host_list/Output/host_translation_table_all_new_wfo.csv")%>%          #### added 5-16-2023
  rbind(read.csv("Host_list/Output/missing_spp_5-12-2023.csv")[,-1])          #### added 5-16-2023

# before synonymizing tip lable info from guide tree,
# substitute/fix error where only sub was used instead of gsub in the translation table search
TranslationTableWFO$spec.name.ORIG  <- gsub("_"," ",TranslationTableWFO$spec.name.ORIG)

#tocheck <- which(is.na(TranslationTableWFO$scientificName))

# if this is not done, errors are produced in joining steps
TranslationTableWFO <- distinct(TranslationTableWFO)
TranslationTableWFO <-
  TranslationTableWFO %>% mutate(scientificName = replace(scientificName, is.na(scientificName), spec.name.ORIG[is.na(scientificName)]))

# read in world plant list and synonymize

##########
# NOTE ########
# input file updated with 6 missing countries (Iran, Koreas, Taiwan, Tanzania, and Bolivia): 5/12/2023
# script re-ran: <not yet as of 5/12/2023>
# therefore, phylogenetic tree built in vphyloplot needs updating

list1_bgci <- read.csv("Host_list/Other_DB/Sarina_BGCI_species_country_df.csv")
list1_bgci <- list1_bgci %>% rename(spec.name.ORIG = species) %>% left_join(distinct(select(TranslationTableWFO, c("scientificName","spec.name.ORIG"))), by="spec.name.ORIG")

## solved the problem below by adding distinct above

# produced a warning but shouldnt based on:
# TranslationTableWFO$spec.name.ORIG %>% table %>% sort(desc=T) %>% head

# Row 151004 
#sum(TranslationTableWFO$spec.name.ORIG == list1_bgci[151004 ,"spec.name.ORIG"])
#TranslationTableWFO[TranslationTableWFO$spec.name.ORIG == list1_bgci[151004 ,"spec.name.ORIG"],]


# read in north america plant list and synonymize
# modified 5-16-23

#na_plantlist <- read.csv( "Host_list/USDA_DB/host_trees_byregion_notsynonymized.csv") %>% dplyr::rename(host_full=host) %>% dplyr::rename(host=host_name)

na_plantlist <- read.csv( "Host_list/Output/host_trees_bystateprovince_notsynonymized.csv") %>% dplyr::rename(host_full=host) %>% dplyr::rename(host=host_name)

na_plantlist <- na_plantlist %>% rename(spec.name.ORIG = host) %>% left_join(select(TranslationTableWFO, c("scientificName","spec.name.ORIG")), by="spec.name.ORIG")

# string two lists together, check that all names are in V.PhyloMaker2 - they are not
all_host_sp_list_translated <- list1_bgci[,c('spec.name.ORIG','scientificName')] %>% rbind(na_plantlist[,c('spec.name.ORIG','scientificName')]) %>% unique

# remove NA matches
all_host_sp_list_translated <- all_host_sp_list_translated[!is.na(all_host_sp_list_translated$scientificName),]

# bind genus and family names to list - we'll need these later
# added distinct 5-16-23
all_host_sp_list_translated <- all_host_sp_list_translated %>% left_join(distinct(select(TranslationTableWFO, c("scientificName","genus","family"))), by="scientificName")

# read in v phyloplot plant list and synonymize
# added distinct 5-16-23
list2_vphyloplot <-tips.info.WP %>% filter(group=="seedplant") %>% dplyr::mutate(host= sub("_", " ", species))
list2_vphyloplot <- list2_vphyloplot %>% rename(spec.name.ORIG = host) %>% left_join(distinct(select(TranslationTableWFO, c("scientificName","spec.name.ORIG"))), by="spec.name.ORIG")

# correct for erroneous use of sub instead of gsub (only replaced first space with underscore) in vphyloplot names
list2_vphyloplot$spec.name.ORIG <- gsub(" ","_",list2_vphyloplot$spec.name.ORIG)

######## now make the tree

# first read in tip info and filter out seed plants for species, genera, and families included in the reference tree
tip.info <- tips.info.WP %>% filter(group=="seedplant") %>% select(-group)

# next make a species info table that has translated names from synonym table but uses genus and family designations consistent with v.phylomaker2
species_table <- all_host_sp_list_translated %>%
			select(-spec.name.ORIG) %>%
			rename(species= scientificName) %>%
			distinct %>%
			left_join(tip.info, by='species') %>%
			rename(genus = genus.y) %>%
			rename(family = family.y) %>%
			mutate(node_inheritance = "tip.info")

species_table[which(is.na(species_table$genus)), c('genus','family', 'node_inheritance')] <-
	mutate(species_table[which(is.na(species_table$genus)), c('genus.x','family.x')], node_inheritance="WFO")

species_table <- select(species_table, c(-genus.x, -family.x))

# family data is needed to bind to tree - and the families have to be present in the node info

species_table <- filter(species_table, !(family == ""))
species_table %>% filter(species == "Mollia nitida")
species_table <-
	species_table %>% 
	mutate(family=replace(family, family == "Viburnaceae", "Adoxaceae")) %>%
	filter(family %in% unique(tip.info$family))

## CHECK FOR DUPLICATES

#length(species_table$species)
#length(unique(species_table$species))
#dim(species_table)
#dim(distinct(species_table))

#make a preliminary sub tree of global tree for seedplants (which were used for synonymy table)
#and substitute spaces for underscores in tip info

# check whats in tip.info thats not in our list
# and filter it out in the preliminary tree
names_tip.info <-
	data.frame(spec.name.ORIG = gsub("_"," ",tip.info$species)) %>% 
	left_join(
		distinct(select(TranslationTableWFO, c("scientificName","spec.name.ORIG"))),
		by="spec.name.ORIG") %>%
		distinct %>%
		filter(!is.na(scientificName))

# ERROR HERE IS THAT SYNONYMIZATION OF TIP NAMES PRODUCES DUPLICATE TIP LABEL NAMES

duplicates <- table(names_tip.info$scientificName) %>% (function(x) x[which(x > 1)] %>% names)

sum(duplicates %in% species_table$species)

duplicates.to.account <- duplicates[duplicates %in% species_table$species]
duplicates.to.remove <- duplicates[!(duplicates %in% species_table$species)]

tt.names_tip.info_duplicates.to.account <- names_tip.info %>% filter(scientificName %in% duplicates.to.account)

original.names.to.remove <- (tt.names_tip.info_duplicates.to.account %>%
	filter(!(spec.name.ORIG %in% tt.names_tip.info_duplicates.to.account$scientificName)))$spec.name.ORIG

names_tip.info <-
	names_tip.info %>%
		filter(!(scientificName %in% duplicates.to.remove)) %>%
		filter(!(spec.name.ORIG %in% original.names.to.remove))

## NOW CHECK AGAIN

#names_tip.info$scientificName %>% length
#names_tip.info$scientificName %>% unique %>% length
#names_tip.info %>% dim
#names_tip.info %>% distinct %>% dim

tr.all <- phylo.maker.verbose(
		tip.info %>%
			filter(gsub("_"," ",species) %in% names_tip.info$spec.name.ORIG),
		tree=GBOTB.extended.WP,
		output.tree=TRUE
	)#$scenario.3

tr<- tr.all$scenario.3

# replace tip label names with synonyms used in species_table dataset from bgci/usda plants

#tr$tip.label %in% nodes.info.1.WP$taxa %>% sum

tr2 <- tr

tr2$tip.label <- names_tip.info$scientificName
tr2.1 <- tr2
tr2$tip.label <- gsub(" ","_",tr2$tip.label)

nodes_mod <- build.nodes.1.verbose(tree=tr2, tips=data.frame(cbind(species=tr2$tip.label, tr.all$species.list[,c("genus", "family")])))

#species_table$species %>% length
#species_table$species %>% unique %>% length


#(species_table %>% mutate(species= gsub(" ","_",species)))$species %>% length
#(species_table %>% mutate(species= gsub(" ","_",species)))$species %>% unique %>% length

# be sure to run / paste this command together with the next(to save)

tr3 <-
#	phylo.maker(
	phylo.maker.verbose(
		species_table %>% select(-node_inheritance) ,#%>% mutate(species= gsub(" ","_",species)),
		tree= tr2,#.1,
		nodes=(nodes_mod),# %>% mutate(taxa = gsub(" ", "_", taxa))),
		output.tree=TRUE)#
		#scenarios="S1")

# save data - calculate unifrac distances

save(tr3, file="Host_list/Output/global_tree_phylogeny.RData")

load("Host_list/Output/global_tree_phylogeny.RData")

##### NOTE NEVER SAVED FINAL TREE DATA AFTER THIS STEP - could be an issue

#nodes_final <- build.nodes.1(tree=tr3$scenario.3, tips=tr3$species.list[,c("species","genus","family")])

#length(tr3$scenario.3$tip.label)
#dim(tr3$species.list)

#> length(tr3$scenario.3$tip.label)
#[1] 56981
#> dim(tr3$species.list)
#[1] 56881     4

#length(unique(tr3$scenario.3$tip.label)) #56881 - there are duplicates
# no duplicates
#reps <- table(tr3$scenario.3$tip.label)[table(tr3$scenario.3$tip.label) > 1]

#x<-NULL
#for (a in 1:length(reps)) x <- c(x, seq(0,reps[a]-1) %>% as.character)

#x[which(x=="0")]<-""

dups <- table(tr3$scenario.3$tip.label)[table(tr3$scenario.3$tip.label) > 1] %>% names

dups.w <- which(tr3$scenario.3$tip.label %in% dups)

dups.o <- order(tr3$scenario.3$tip.label[dups.w])

to_replace <- paste(tr3$scenario.3$tip.label[dups.w][dups.o], x, sep="")

tr3$scenario.3$tip.label[dups.w][dups.o]<-to_replace

tips_to_prune <- setdiff(tr3$scenario.3$tip.label[dups.w], names(reps))

tr3.2 <- drop.tip(tr3$scenario.3, tips_to_prune)

#str(tr3.2)
#str(tr3$species.list)

#setdiff(tr3.2$tip.label, gsub(" ","_",tr3$species.list$species))
#setdiff(gsub(" ","_",tr3$species.list$species),tr3.2$tip.label)

tips<-data.frame(species=gsub(" ","_",tr3$species.list$species), tr3$species.list[,c("genus","family")])

nodes_final <- build.nodes.1.verbose(tree=tr3.2, tips=tips)

save(nodes_final, file="Host_list/Output/nodes_final.Rdata")
save(tr3.2, file="Host_list/Output/tree_final.Rdata")