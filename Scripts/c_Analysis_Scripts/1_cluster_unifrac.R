library(ape)
library(tidyr)
library(dplyr)
library(ggplot2)

###########################
# First, we want to make a hierarchical cluster of all the state-country tree species lists
# The unifrac distances were already calculated based on the lists - so we read in the unifrac object

load("Host_list/Output/countries_unifrac_3.RData")#_russia_split.RData")

hc.countries <- hclust(uf.countries, method="ward.D")

windows()
plot(hc.countries, cex=1)

################################
# now lets make a figure

#pdf("Figures/clustered_countries_northamerica_stateprovince.pdf",width=35,height=18);plot(hc.countries, cex=1);dev.off()

###############
# the next is to arbitrarily cut off the clusters into 18 subclusters
# so that we have a list with the countries belonging to each cluster
clusters <- cutree(hc.countries, k = 18)
obj <- list()
for(i in 1:max(clusters)) {
  obj[[length(obj)+1]] <- names(which(clusters == i))
}

# Give them nicer names later
names(obj) <- paste(x=rep("Cluster"), num=1:length(obj), sep=".")

obj

# now we give them tentative names
names(obj) <-
  c(
    "MegaCluster2",                 #1!!
    "SE.North.America",             #2
    "N.North.America",              #3
    "C.Europe",                     #4
    "Levant.Iberia.N.Africa",       #5
    "MegaCluster1",                      #6
    "Trop.Africa",                  #7
    "Carribbean",                   #8
    "W.North.America",                 #9
    "N.Europe",             #10
    "C.America",                    #11
    "Sahel",                   #12
    "SE.Asia",                       #13
    "Amazonia" ,                    #14
    "S.Africa",                    #15
    "NC.North.America",                    #16
    "Arabia.Sahara",                 #17
    "NE.North.America")


###############
# need to break up the southern hemisphere-pacific megacluster
# (break off japan, taiwan, and korea)
#

obj
rc <- which(attr(uf.countries,"Labels") %in% obj$MegaCluster1)
sub.matrix1 <- as.matrix(uf.countries)[rc,rc] %>% as.dist

sub.cluster1 <- hclust(sub.matrix1, method="ward.D")
clusters1 <- cutree(sub.cluster1, k = 4)
clusters1
plot(sub.cluster1)

obj$N.Pacific <- names(clusters1[clusters1==4])
obj$S.Cone.Pacific <- names(clusters1[clusters1 %in% c(2,3,1)])

obj<- obj[-which(names(obj)=="MegaCluster1")]

obj

###################
## now we need to clean up some additional clusters by hand

rc2 <- which(attr(uf.countries,"Labels") %in% with(obj, c(N.Europe, Levant.Iberia.N.Africa, C.Europe, MegaCluster2)))
sub.matrix2 <- as.matrix(uf.countries)[rc2,rc2] %>% as.dist
sub.cluster2 <- hclust(sub.matrix2, method="ward.D")
clusters2 <- cutree(sub.cluster2, k = 7)

# N.Europe, C.Europe, N.Eurasia = Eurasian Palearctic
clusters2[which(clusters2==5)]
clusters2[which(clusters2==4)]
clusters2[which(clusters2==2)]

# MiddleEast, Iberia, Central Asia
clusters2[which(clusters2==1)]
clusters2[which(clusters2==3)]
clusters2[which(clusters2==6)]
clusters2[which(clusters2==7)]

plot(sub.cluster2)

obj$MiddleEast.Iberia.N.Africa <- names(clusters2[clusters2 %in% c(1,3,6,7)])
obj$Eurasian.Palearctic <- names(clusters2[clusters2 %in% c(2,4,5)])
obj<- obj[-which(names(obj) %in% c('N.Europe', 'Levant.Iberia.N.Africa', 'C.Europe', 'MegaCluster2'))]

length(obj)
obj

# finally cluster Sahel, Trop.Africa, and S.Africa together

obj$Subsaharan.Africa <- with(obj, c(Sahel, Trop.Africa, S.Africa))
obj<- obj[-which(names(obj) %in% c('Sahel', 'Trop.Africa', 'S.Africa'))]

#obj$Subsaharan.Africa <- with(obj, c(Sahel, Trop.Africa, S.Africa))
#obj<- obj[-which(names(obj) %in% c('Sahel', 'Trop.Africa', 'S.Africa'))]

obj$NE.North.America <- with(obj, c(NE.North.America, NC.North.America))
obj <- obj[-which(names(obj) == "NC.North.America")]

length(obj)
obj

# fix guyana
obj$Amazonia
obj$Eurasian.Palearctic %>% sort
# appears to be a mapping issue

# fix spain-middleeast-nafrica
obj$Iberia.N.Africa <- c("Algeria",'Morocco',"Portugal","Spain",'Tunisia')
obj$Middle.East <- setdiff(obj$MiddleEast.Iberia.N.Africa, obj$Iberia.N.Africa)
obj <- obj[-which(names(obj) == "MiddleEast.Iberia.N.Africa")]

##################
# separate Hawaii!
########obj$Hawaii <- c("Hawaii")
########obj$S.Cone.Pacific <- setdiff(obj$S.Cone.Pacific, "Hawaii")

###########
##############save(obj, file="Host_list/Output/countries_unifrac_clusters1.RData" )

save(obj, file="Host_list/Output/countries_unifrac_clusters1_nohawaii.RData" )

########
# maybe make a dendrogram with these clusters highlighted

#library(ggdendro)
#library(zoo)
library(dendextend)
library(magrittr)
library(ggplot2)

length(obj)

c16 <- #pals::polychrome(29)
  c(
    #"dodgerblue2"
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    #"black",
    "gold1",
    #"skyblue2",
    "#FB9A99",
    ##############"blue1", 
    "palegreen2",
    #"#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    #"maroon", "orchid1", "deeppink1", #"steelblue4",
    "darkturquoise",# "green1", "yellow4", #"yellow3",
    "darkorange4", #"brown"#,
    #"#822E1C",
    "#BDCDFF","#FA0087","#2ED9FF","#E31A1C" # red
  )

clust <- NULL; for(i in 1:length(obj)) {
  zz <- rep(i,length(obj[[i]]))
  names(zz) <- obj[[i]]
  clust <- c(clust, zz)
}
hcdendro <- as.dendrogram(hc.countries)
clust_reordered <- clust[labels(hcdendro)]

rotateorder<-NULL
for (j in unique(clust_reordered))
  {rotateorder<-c(
    rotateorder,
    clust_reordered[which(clust_reordered==j)])}

names(c16) <- sort(names(obj))
c16a <- c16[names(obj)]


hcdendro <- as.dendrogram(hc.countries) %>%
  rotate(order=names(rotateorder)) %>%
  set("labels_cex", .5) %>%
  set("labels_colors", c16a[rotateorder]) %>%
  color_branches(
    col=c16a[rotateorder] %>% unique,
    clusters=rotateorder) %>%
  set("labels","") %>%
  set("branches_lwd", 1.25)

hc.dendro.hclust <- as.hclust(hcdendro)

save(hc.dendro.hclust, file="Final_output/country_dendro_data_asplotted.RData")

gg.dendro.countries<-as.ggdend(hcdendro) %>% ggplot(horiz=F) + scale_y_sqrt()
windows(6,4);gg.dendro.countries
save(gg.dendro.countries, file="Figures_July2024/1b_clusters_countries.RData")


###

#colors_ordered15 <- c15[rotateorder] %>% unique
