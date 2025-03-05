library(ape)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)

###########################
# First, we want to make a hierarchical cluster of all the state-country tree species lists
# The unifrac distances were already calculated based on the lists - so we read in the unifrac object

load("Host_list/Output/countries_unifrac_3.RData")#_russia_split.RData")

# now look at number of species
# must have corrected a bunch of errors in the country names in the last script
trees.countries.translated.trees.per.country <-
  read.csv("Host_list/Output/trees.countries.all.translated.csv") %>%
  group_by(country_name) %>%
  summarise(nspp = n()) %>%
  ungroup() %>%
  
  # add data 
  rbind(
    data.frame(country_name="Bolivia",
    nspp=dim(read.csv("Host_list/Input/BGCI_Missing/Bolivia, Plurinational State of - Tree Species List - 2023-05-12.csv"))[1])) %>%
  rbind(
    data.frame(country_name="Iran",
    nspp=dim(read.csv("Host_list/Input/BGCI_Missing/Iran, Islamic Republic of - Tree Species List - 2023-05-12.csv"))[1])) %>%
  rbind(
    data.frame(country_name="North Korea",
    nspp=dim(read.csv("Host_list/Input/BGCI_Missing/Korea, Democratic People's Republic of - Tree Species List - 2023-05-12.csv"))[1])) %>%
  rbind(
    data.frame(country_name="South Korea",
    nspp=dim(read.csv("Host_list/Input/BGCI_Missing/Korea, Republic of - Tree Species List - 2023-05-12.csv"))[1])) %>%
  rbind(
    data.frame(country_name="Taiwan",
    nspp=dim(read.csv("Host_list/Input/BGCI_Missing/Taiwan, Province of China - Tree Species List - 2023-05-12.csv"))[1])) %>%
  rbind(
    data.frame(country_name="Tanzania",
    nspp=dim(read.csv("Host_list/Input/BGCI_Missing/Tanzania, United Republic of - Tree Species List - 2023-05-12.csv"))[1])) %>%

  mutate(nspp.log = log(nspp))%>%
  mutate(nspp.log.scale = scale(nspp.log)%>%as.numeric)

# where is Bolivia and North and South Korea in the trees.countries.translated data? -
# they are in uf.countries but not in the translated species table version here for some reason -
# so we added them above
#setdiff(attr(uf.countries,"Labels"), trees.countries.translated.trees.per.country$country_name)



# these were the countries that got removed from the unifrac dataset
#setdiff(trees.countries.translated.trees.per.country$country_name,
#        attr(uf.countries,"Labels"))

#attr(uf.countries,"Labels")
#trees.countries.translated.trees.per.country$country_name %>% unique

hist(log(trees.countries.translated.trees.per.country$nspp))
hist(scale(log(trees.countries.translated.trees.per.country$nspp)))
hist(trees.countries.translated.trees.per.country$nspp.log.scale)

hist(uf.countries)
hist(log(1-uf.countries))
hist(scale(log(1-uf.countries)))

nsppdiff.countries <-
  trees.countries.translated.trees.per.country %>%
  with(expand_grid(a=country_name,b=country_name)) %>%
  filter(a!=b)%>%
  left_join(trees.countries.translated.trees.per.country, by = join_by(a==country_name))%>%
  left_join(trees.countries.translated.trees.per.country, by = join_by(b==country_name))%>%
  mutate(nsppdiff = abs(nspp.x-nspp.y),
         nsppdiff.log.scale = abs(nspp.log.scale.x-nspp.log.scale.y))

hist(nsppdiff.countries$nsppdiff)
hist(log(nsppdiff.countries$nsppdiff.log.scale))

uf.nspp.df <-
  uf.countries %>% as.matrix %>% as.data.frame %>%
  mutate(a = attr(uf.countries,"Labels")) %>%
  reshape2::melt("a",variable.name="b",value.name="unifrac") %>%
  filter(a!=b)%>%
    left_join(nsppdiff.countries)

uf.nspp.df$unifrac.m1log <- NA
uf.nspp.df$unifrac.m1log[which(uf.nspp.df$unifrac<1)] <-
  log(1-uf.nspp.df$unifrac[which(uf.nspp.df$unifrac<1)])

quantile(na.omit(uf.nspp.df$unifrac.m1log), probs=(0:100)/100)
uf.nspp.df$unifrac.m1log[which(uf.nspp.df$unifrac.m1log< -12 | uf.nspp.df$unifrac==1)] <- -12

hist(uf.nspp.df$unifrac.m1log)

center.unifrac.m1log <-
  ecdf(uf.nspp.df$unifrac.m1log[which(uf.nspp.df$unifrac.m1log>=-12)])((0:120-120)/10) %>%
  (function(x) x[-1] -x[-length(x)]) %>% which.max %>%
  (function(x) ((0:120-120)/10)[x])

uf.nspp.df$unifrac.m1log.scale <- NA
uf.nspp.df$unifrac.m1log.scale[which(uf.nspp.df$unifrac.m1log>=-12)] <-
  scale(uf.nspp.df$unifrac.m1log[which(uf.nspp.df$unifrac.m1log>=-12)],
        center = center.unifrac.m1log) %>% as.numeric

hist(uf.nspp.df$unifrac.m1log.scale)

quantile(na.omit(uf.nspp.df$unifrac.m1log.scale), probs=(0:100)/100)

# plot number of species difference against unifrac distance

uf.nspp.df %>% filter(unifrac.m1log.scale > -1.9) %>%
  with(plot(nsppdiff.log.scale, unifrac.m1log.scale, pch='·', col = 'grey'))

uf.nspp.df.model <-
  uf.nspp.df %>% filter(unifrac.m1log.scale > -1.9) %>%
    with(lm(unifrac.m1log.scale ~ nsppdiff.log.scale))

points(0:35/10, predict(uf.nspp.df.model, data.frame(nsppdiff.log.scale=0:35/10)), type='l', lwd=2)

quantile(residuals(uf.nspp.df.model), probs=(0:100)/100)

uf.nspp.df.model %>%
  with(plot(model$nsppdiff.log.scale, residuals, pch='·', col = 'grey'))

# use residuals as distance measure
uf.nspp.df$unifrac.transformed.residual <-
  quantile(residuals(uf.nspp.df.model), probs=(0:100)/100)[1] %>% as.numeric

uf.nspp.df$unifrac.transformed.residual[which(uf.nspp.df$unifrac.m1log.scale > -1.9)] <-
  uf.nspp.df$unifrac.m1log.scale[which(uf.nspp.df$unifrac.m1log.scale > -1.9)] -
  predict(uf.nspp.df.model,
          data.frame(
            nsppdiff.log.scale=
              uf.nspp.df$nsppdiff.log.scale[which(uf.nspp.df$unifrac.m1log.scale > -1.9)]))

# should look the same as the last plot but with a line of ~ -2.00 at the bottom
uf.nspp.df %>%
  with(plot(nsppdiff.log.scale, unifrac.transformed.residual, pch='·', col = 'grey'))

uf.nspp.df$unifrac.trans.res.std.rev <- NA
#  mutate(uf.nspp.df,
#         
#         unifrac.trans.res.std.rev =
#           
#           (function (x) 1 - ((x - min(x))/(max(x)-min(x))))(unifrac.transformed.residual)
#         
#         )

uf.nspp.df$unifrac.trans.res.std.rev[!is.na(uf.nspp.df$unifrac.transformed.residual)] <-
  (function (x) 1 - ((x - min(x))/(max(x)-min(x))))(uf.nspp.df$unifrac.transformed.residual[!is.na(uf.nspp.df$unifrac.transformed.residual)])

hist(uf.nspp.df$unifrac.trans.res.std.rev)

#points(0:35/10, -0.310906 + -0.223156*(0:35)/10, type='l', lwd=2, col='purple')

#######################

uf.countries.2 <-
  
  reshape2::acast(uf.nspp.df %>% filter(!(
    # for now do this for testing - need to reconcile country names later
    
    a %in% setdiff(attr(uf.countries,"Labels"),
                   trees.countries.translated.trees.per.country$country_name) |
      b %in% setdiff(attr(uf.countries,"Labels"),
                     trees.countries.translated.trees.per.country$country_name)
    
    
  )),
                  a ~ b,
                  value.var = 'unifrac.trans.res.std.rev') %>% as.dist

save(uf.countries.2, file = "Host_list/Output/countries_unifrac_corrected_transformed.RData")

# still use the untransformed Unifrac distances to cluster...

hc.countries <- hclust(uf.countries, method="ward.D")
#hc.countries <- hclust(as.dist(uf.countries.2), method="ward.D")

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
    "MegaCluster2",                 #1!!    - Middleast + Spain, Portugal
    "SE.North.America",             #2      - Still SE but some things transferred over
    "N.North.America",              #3      - more W + N
    "C.Europe",                     #4      - Europe
    "Levant.Iberia.N.Africa",       #5      - Some islands in the Pacific
    "MegaCluster1",                      #6 - Africa - southern
    "Trop.Africa",                  #7      - Caribbean
    "Carribbean",                   #8      - This is now basically all of South America
    "W.North.America",                 #9   - Central America
    "N.Europe",             #10             - Just Australia and New Zealand + Norfolk Island
    "C.America",                    #11     - SE Asia
    "Sahel",                   #12          - Sahel
    "SE.Asia",                       #13    - More of the pacific
    "Amazonia" ,                    #14     - NE North America
    "S.Africa",                    #15      - Northern and Northeast Africa
    "NC.North.America",                    #16 - Japan, Koreas, and Mongolia (oddly)
    "Arabia.Sahara",                 #17       - Madagascar and Mayotte
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
