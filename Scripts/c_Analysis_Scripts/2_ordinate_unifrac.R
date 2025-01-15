# use phyloseq object and resulting clusters to ordinate countries and find distance between centroids

library(reshape2)
library(vegan)
library(doBy)
library(ggplot2)
library(dplyr)
library(cowplot)

load("Host_list/Output/countries_unifrac_3.RData")
load("Host_list/Output/countries_unifrac_clusters1.RData")

# ordinate using NMDS
ordination1<- vegan::metaMDS(uf.countries, k=2, trymax=2000, tidy=T, previous.best=F)#, engine="isoMDS")

# create a vector of country names with their assignment

bioregions_vector <- data.frame(country_name=attr(ordination1$points,"dimnames")[[1]]) %>%
  left_join(melt(obj), by = join_by(country_name==value)) %>%
  rename(biogeographic_region=L1)
bioregions_vector$biogeographic_region <- as.factor(bioregions_vector$biogeographic_region)

## now plot the clusters

#c16 <- #pals::polychrome(29)
#  c(
#    #"dodgerblue2"
#    "green4",
#    "#6A3D9A", # purple
#    "#FF7F00", # orange
#    #"black",
#    "gold1",
#    #"skyblue2",
#    "#FB9A99",
#    "blue1", # lt pink
#    "palegreen2",
#    #"#CAB2D6", # lt purple
#    "#FDBF6F", # lt orange
#    "gray70", "khaki2",
#    #"maroon", "orchid1", "deeppink1", #"steelblue4",
#    "darkturquoise",# "green1", "yellow4", #"yellow3",
##    "darkorange4", #"brown"#,
#    #"#822E1C",
#    "#BDCDFF","#FA0087","#2ED9FF","#E31A1C" # red
#  )


# match color pallette by naming elements
c16 <-
  c(
    Amazonia = "green4",               # 1 Amazonia
    Arabia.Sahara = "#6A3D9A",         # 2 Arabia.Sahara
    C.America = "#FF7F00",             # C.America
    Caribbean = "gold1",
    Eurasian.Palearctic = "#FB9A99",
    Hawaii = "blue1",                   # Hawaii
    Iberia.N.Africa = "palegreen2",     # Iberia.N.Africa
    Middle.East = "#FDBF6F",
    N.North.America = "gray70",               
    N.Pacific = "khaki2",                  # N.Pacific
    NE.North.America = "darkturquoise",
    S.Cone.Pacific = "darkorange4",
    SE.Asia = "#BDCDFF",
    SE.North.America = "#FA0087",
    Subsaharan.Africa = "#2ED9FF",
    W.North.America = "#E31A1C",
    Australia = "yellow" # Australia ADDED ONE COLOR TEMPORARILY
  )


data.scores <- as.data.frame(vegan::scores(ordination1))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores) # create a column of site names, from the rownames of data.scores
data.scores$grp <- bioregions_vector$biogeographic_region  %>%
  gsub("\\.", " ", .) %>% gsub("([A-Z]{1}[a-z]+)([A-Z]{1}[a-z]+)","\\1 \\2", .)#  add the grp variable created earlier
head(data.scores)  #look at the data

levels(as.factor(data.scores$grp))

names(c16) <- gsub("\\."," ",names(c16)) %>% gsub('Carribbean','Caribbean',.)

str(data.scores)

data.scores$grp <- factor(sub('Carribbean','Caribbean', data.scores$grp),
                          levels = c("SE North America",
                                     "NE North America",
                                     "Eurasian Palearctic",
                                     "Iberia N Africa",
                                     "Middle East",
                                     "W North America",
                                     "N North America",
                                     "SE Asia",
                                     "S Cone Pacific",
                                     "N Pacific",
                                     "Caribbean",
                                     "C America",
                                     "Amazonia",
                                     "Subsaharan Africa",
                                     "Arabia Sahara",
                                     "Hawaii"))

gg1<- ggplot() +
  geom_point(data=data.scores, aes(x=NMDS1, y=NMDS2, colour=grp), size=3) +
  ggplot2::stat_ellipse(data=data.scores, aes(x=NMDS1, y=NMDS2, colour=grp, fill=grp), geom="polygon", alpha=0.2, linetype=2, type="t", level=.68) +
#  ggplot2::stat_ellipse(data=data.scores, aes(x=NMDS1, y=NMDS2, colour=grp), linetype=3, type="t", level=.95) +
  scale_fill_manual("Unifrac Floristic Cluster", values=c16) +
  scale_color_manual("Unifrac Floristic Cluster", values=c16) +
  coord_cartesian()+#xlim=c(-.5,.6), ylim=c(-.4,.5)) +
  theme_bw()

windows(6,6)
gg1

save(gg1, c16, file="Figures_July2024/S1_ordination_of_cluster.RData")

# now do some analysis

dats2<-ordination1$points[,1:2] %>% as.data.frame %>% cbind(biogeographic_region=as.character(bioregions_vector[,"biogeographic_region"]))
dat2.summary <- summaryBy(
  cbind(MDS1, MDS2) ~ biogeographic_region,
  data=dats2,
  FUN=c(mean=mean)) %>%
  rbind(
    cbind(biogeographic_region="Australia",dats2["Australia",c("MDS1","MDS2")]) %>%
      rename(MDS1.mean=MDS1) %>% rename(MDS2.mean=MDS2)
  )

save(dats2, dat2.summary, gg1, bioregions_vector, file = "Host_list/Output/ordinated_unifrac_data.Rdata")
