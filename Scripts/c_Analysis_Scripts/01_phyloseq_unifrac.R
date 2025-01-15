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

# this is just a helper function that converts between long and wide format

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

# read in the country data for the tree species - translated by WFO
# you can read in the version with Russia split E vs W or the version with Russia as one country
# trees.countries.translated<- read.csv("Host_list/Output/trees.countries.all.translated.split.csv")
trees.countries.translated<- read.csv("Host_list/Output/trees.countries.all.translated.csv")

# now make a 'otu table' of countries
# make the matrix
mat.world <- reshape2::dcast(trees.countries.translated, scientificName ~ country_name, fill=0)#, fun.aggegate=length)

# replace spaces with underscores
mat.world$scientificName <- gsub(" ", "_", mat.world$scientificName)

# take out taxa from tree not in the global tree list
tree <- drop.tip(tr3$tree.scenario.3, setdiff(tr3$tree.scenario.3$tip.label, mat.world$scientificName))

save(tree, file="Host_list/Output/tree_of_trees_just_trees.RData")

# take out taxa from matrix not in tree
mat.world <- mat.world[-which(!(mat.world$scientificName %in% tree$tip.label)),]
row.names(mat.world) <- mat.world$scientificName

#check
setdiff(mat.world$scientificName, tree$tip.label)
setdiff(tree$tip.label, mat.world$scientificName)

# need to tet rid of duplicates...
mat.world$scientificName %>% table %>% sort(T) %>% (function(x) x[which(x>1)]) %>% names
tree$tip.label %>% table %>% (function(x) x[which(x>1)])

# need to prune duplicate tips from tree
dups <- tree$tip.label %>% table %>% (function(x) x[which(x>1)])
which.dups <- which(tree$tip.label %in% names(dups))

to.drop <- NULL
for(i in dups) {
  to.drop <- c(to.drop, seq(from=1, to=i))
}

to.drop.index<-data.frame(which.dups, tree$tip.label[which.dups]) %>%
  arrange(tree.tip.label.which.dups.) %>%
  cbind(to.drop) %>%
  filter(to.drop > 1) %>%
  select(which.dups) %>%
  unlist

# do the pruning
tree <- drop.tip(tree, to.drop.index)

#check again
tree$tip.label %>% table %>% sort(T) %>% (function(x) x[which(x>1)])
setdiff(mat.world$scientificName, tree$tip.label)
setdiff(tree$tip.label, mat.world$scientificName)

######### now we are going to make a matrix of tree species ocurrence for the countries of the world
mat.world.2 <- replace(mat.world[,-1], mat.world[,-1] != 0, 1) %>% as.matrix
storage.mode(mat.world.2)<-"numeric"

# now we can combine them into a phyloseq object
OTU <- otu_table(object=mat.world.2,T)
PT <- phy_tree(tree)
ps.countries <- phyloseq(OTU, PT)

# check the names
taxa_names(OTU)

# noew we can calculate the unifrac distances based on our tree
registerDoParallel(8, cores=4)
uf.countries <- UniFrac(ps.countries, parallel=T)
save(uf.countries, file="Host_list/Output/countries_unifrac_splitRussia.RData")

# you can reload the data here so you don't have to wait for the unifrac
#load("Host_list/Output/countries_unifrac.RData")
mean(uf.countries)

sums<-colSums(otu_table(ps.countries))
sum(sums<100)

# find the countries with less than 35 tree species
#- might remove them because its throwing off our unifrac analyses
additional.names <- names(sums[sums<35])

hist(sums)
hist(sums[which(sums<1000)])
colnames(otu_table(ps.countries))

# you can reload the data here so you don't have to wait for the unifrac
# load("Host_list/Output/recluster/countries_unifrac_2.RData")
# check the data out
str(uf.countries)

# rename the matrix with the correct country names
names(uf.countries) <- attr(uf.countries, "Labels")

# check the data out
as.matrix(uf.countries)[1:5,1:5]

# remove countries with not enough species and these island nations just for good measure
# they don't cluster well in our analysis
countries.to.remove <- c(additional.names, "Saint Pierre and Miquelon","United States Minor Outlying Islands","Disputed Territory","Svalbard and Jan Mayen","Saint BarthÃ©lemy", "Ã…land Islands", "Monaco", "Nauru", "Tuvalu", "San Marino", "Liechtenstein", "Marshall Islands", "Saint Kitts and Nevis", "Maldives", "Malta", "Grenada","Saint Vincent and the Grenadines","Barbados","Antigua and Barbuda", "Seychelles", "Palau", "Andorra", "Saint Lucia", "Micronesia, Federated States of", "Singapore", "Tonga", "Dominica", "Bahrain", "Kiribati", "Sao TomÃ© and Principe", "Comoros", "Mauritius", "Luxembourg","RÃ©union")
rc1 <- which(names(uf.countries) %in% countries.to.remove)

# remove the rows and columns for those countries and resave the unifrac object
as.matrix(uf.countries)[-rc1,-rc1] %>% row.names

# convert to a distance object
uf.countries <- as.matrix(uf.countries)[-rc1,-rc1] %>% as.dist

# double check the country names
attr(uf.countries, "Labels")
#uf.countries

# export to an r data object
# save(uf.countries, file="Host_list/Output/countries_unifrac_3_russia_split.RData")
# save(uf.countries, file="Host_list/Output/countries_unifrac_3.RData")
########################################################################################################################
