library(dplyr)
library(tidyr)

setwd("Host_list/Invasive_Trees_States/")

hosts <- NULL


  filenames <- grep("Hawaii",dir(), value=T) #># should just return hawaii filename
  for (j in filenames) {
    csv <- read.csv(j, skip=3)$Scientific.Name
    hosts <- union(hosts, csv)
  }

df <- data.frame(host= hosts)

  hosts <- NULL
  
  for (j in filenames) {
    hosts <- read.csv(j, skip=3)$Scientific.Name
    df <- df %>% mutate(n = host %in% hosts)
    names(df)[dim(df)[2]] <- strsplit(j, "\\.")[[1]][1]
  }


dim(df)

df <- arrange(df, host)

df <- tidyr::extract(df, col= host, remove= F, regex = "^([a-zA-Z\\×]+[ ][a-zA-Z\\-\\×]+) (.*)", into = c("host_name", "host_authority"))
df <- tidyr::extract(df, col = host_name, remove = F, regex="^([a-zA-Z]+)[ ]", into = "host_genus")

setwd("../..")

translation_table <- read.csv("Host_list/Input/host_translation_table_all_new_WFO.csv")

df <- df %>%
  left_join(
    dplyr::select(translation_table, spec.name, scientificName) %>% distinct,
    by = join_by("host_name"=="spec.name"))

# seems like we lost about ~20% of host names
# may be good to rerun  plant name indexer
with(df %>% filter(!is.na(host_name)),
     length(which(is.na(scientificName)))/ length(scientificName))

# compared to 3% in the other list
with(read.csv("Host_list/Output/ExoticTrees_bystate.csv") %>% filter(!is.na(host_name)),
     length(which(is.na(scientificName)))/ length(scientificName))

# re-do synonymization

library(WorldFlora)
WFO.remember("../Retrospective/WFO_Backbone/classification.csv")

#xx<-WFO.match(spec.data=names_new[start:stop], WFO.data=WFO.data)
#yy<-WFO.one(xx)
dim(df %>% filter(is.na(scientificName)))
xx<-with(df %>% filter(is.na(scientificName)),
     WFO.match(host_name, WFO.data=WFO.data))
yy<-WFO.one(xx)
yy<- filter(yy, spec.name.ORIG != "")
#write.csv(df, "Host_list/Output/ExoticTrees_bystate.csv", row.names=F)

write.csv(yy,"Host_list/Output/Hawaii_missing_transtable_WFO.csv", row.names=F)

yy<-read.csv("Host_list/Output/Hawaii_missing_transtable_WFO.csv")

#with(df %>% filter(!is.na(host_name)),
#     length(which(is.na(scientificName))))

#dim(yy)
#tail(yy)

#sum(is.na(df %>% filter(!is.na(host_name)) %>% select(scientificName)))

#df %>% filter(!is.na(host_name)) %>% select(scientificName)

#length(which(is.na(yy$scientificName)))

#setdiff(df$host_name[which(is.na(df$scientificName))], )


###### RESTART HERE

#df <- read.csv("Host_list/Output/ExoticTrees_bystate.csv")


# remove empty host_name and add newly translated names
df2 <- df[-which(is.na(df$host_name)),]

df3 <- df2[-which(is.na(df2$scientificName)),]

df4 <- df3 %>%
  bind_rows(distinct(data.frame(Hawaii=TRUE, host_name=yy$spec.name.ORIG, scientificName=yy$scientificName)))

# df4 <- df3

# remove duplicates -this messed up the DF****

df4[df4$host_name %in% names(table(df4$host_name)[table(df4$host_name)>1]),]

for (n in names(table(df4$host_name)[table(df4$host_name)>1])) {
  if (length(which(df4$host_name==n))==0) next
  df4 <- df4[-which(df4$host_name==n),] %>%
    rbind(
      c(
        df4[which(df4$host_name==n)[1],1:4],
        (df[which(df$host_name==n),-c(1:4,dim(df)[2])] %>% colSums)>0,
        #Hawaii=(df4[which(df4$host_name==n),-c(1:4,dim(df4)[2])] %>% sum)>0,
        scientificName=df4[which(df4$host_name==n)[1],'scientificName']
  ))
}
names(table(df4$host_name)[table(df4$host_name)>1])


df<- df4


# now read in bioregion data
# first aggregate by cluster
load("Host_list/Output/countries_unifrac_clusters1.RData")

obj[['N.North.America']] <- obj[['N.North.America']][-which(obj[['N.North.America']] %in% c("Labrador","Nunavut"))]

obj[['NE.North.America']][which(obj[['NE.North.America']]=="Deleware")]<-"Delaware"

exotichost_occurrence <- matrix(data=1, nrow=dim(df)[1], ncol=1, dimnames=list(df$host_name, c("Hawaii")))

########## CAUGHT BUG ^^ using untranslated name to derive ** fixed in both scripts

#for (regs in colnames(exotichost_occurrence)) {
#  for (specs in rownames(exotichost_occurrence)) {
#    for (states in obj[[regs]]) {
#      if (df[which(df$host_name==specs),states]) exotichost_occurrence[specs,regs] <- exotichost_occurrence[specs,regs]+1
      ## need to multiply by relative land area... but can just use if > 0 for now...
#    }
#  }
#}

#write.csv(exotichost_occurrence, 'Host_list/Output/HI_exotic_host_occurrence_matrix.csv', row.names=TRUE)

exotichost_occurrence <- read.csv('Host_list/Output/HI_exotic_host_occurrence_matrix.csv', row.names=1)

trees.countries <- read.csv('Host_list/Output/trees.countries.all.translated.csv')

exotic_trees <- expand.grid(colnames(exotichost_occurrence), names(obj)%>% setdiff(colnames(exotichost_occurrence)))
names(exotic_trees)<-c("Invaded.Range.in.North.America.Cluster","Geographic.Origin.Cluster")
exotic_trees$number_of_exotic_trees <- 0

rownames(exotichost_occurrence)
df$scientificName

unique(rownames(exotichost_occurrence))
unique(df$scientificName)

trees.list2 <- list()

for (i in 1:dim(exotic_trees)[1]) {
  exotic_trees$number_of_exotic_trees[i] <- sum(
    na.omit(unique(df$scientificName[
      exotichost_occurrence[,exotic_trees$Invaded.Range.in.North.America.Cluster[i]]>0])) %in%
    unique((trees.countries %>%
       filter(country_name %in% obj[[as.character(exotic_trees$Geographic.Origin.Cluster[i])]]))$scientificName))

  #add any new trees to the list
  
  current.origin <- exotic_trees[i,'Geographic.Origin.Cluster']
  
  # there is only one per Hawaii
  #if(!(current.origin %in% names(trees.list1))) {
    trees.list2[[length(trees.list2) + 1]] <- character()
    names(trees.list2)[[length(trees.list2)]] <- as.character(current.origin)
  #}
  
  #add any new trees to the list
    trees.list2[[exotic_trees[i,'Geographic.Origin.Cluster']]] <-
    
    # this list represents the exotic trees from the current Geographic.Origin.Cluster
    # and present in the current location (Europe or Australia) from the expanded grid
    intersect(
      na.omit(unique(df$scientificName[
        exotichost_occurrence[,exotic_trees$Invaded.Range.in.North.America.Cluster[i]]>0])),
      unique((trees.countries %>%
                filter(country_name %in% obj[[as.character(exotic_trees$Geographic.Origin.Cluster[i])]]))$scientificName)
    ) #%>%
    
    # there is only one per Hawaii
    # now we only want to retain the ones that are new to the list, and add them to the list
    #union(trees.list1[[exotic_trees[i,'Geographic.Origin.Cluster']]]) %>% unique()
  
  #  and then later we can sum them all up
}

write.csv(exotic_trees, "Covariates_Input/HI_exotic_trees.csv")
save(trees.list2, file="Covariates_input/HI_trees_list.RData")
