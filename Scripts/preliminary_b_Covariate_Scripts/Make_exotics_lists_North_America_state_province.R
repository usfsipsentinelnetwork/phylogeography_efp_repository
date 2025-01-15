library(dplyr)
library(tidyr)

setwd("Host_list/Invasive_Trees_States/")

hosts <- NULL


  filenames <- dir()
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

df <- extract(df, col= host, remove= F, regex = "^([a-zA-Z\\×]+[ ][a-zA-Z\\-\\×]+) (.*)", into = c("host_name", "host_authority"))
df <- extract(df, col = host_name, remove = F, regex="^([a-zA-Z]+)[ ]", into = "host_genus")

setwd("../..")

translation_table <- read.csv("Host_list/Input/host_translation_table_all_new_WFO.csv")

df <- df %>%
  left_join(
    select(translation_table, spec.name, scientificName) %>% distinct,
    by = join_by("host_name"=="spec.name"))

write.csv(df, "Host_list/Output/ExoticTrees_bystate.csv", row.names=F)

####

# now read in bioregion data
# first aggregate by cluster

#df <- read.csv("Host_list/Output/ExoticTrees_bystate.csv")

load("Host_list/Output/countries_unifrac_clusters1.RData")

df <- df[-which(is.na(df$host_name)),]
# remove duplicates

for (n in names(table(df$host_name)[table(df$host_name)>1])) {
  
  df <- df[-which(df$host_name==n),] %>%
    rbind(
      c(
        df[which(df$host_name==n)[1],1:4],
        (df[which(df$host_name==n),-c(1:4,dim(df)[2])] %>% colSums)>0,
        scientificName=df[which(df$host_name==n)[1],'scientificName']
  ))
}
names(table(df$host_name)[table(df$host_name)>1])

obj[['N.North.America']] <- obj[['N.North.America']][-which(obj[['N.North.America']] %in% c("Labrador","Nunavut"))]

obj[['NE.North.America']][which(obj[['NE.North.America']]=="Deleware")]<-"Delaware"
exotichost_occurrence <- matrix(data=0, nrow=dim(df)[1], ncol=4, dimnames=list(df$host_name, names(obj) %>% grep("North\\.America",., value=T)))

# not time to use scientificName yet

for (regs in colnames(exotichost_occurrence)) {
  for (specs in rownames(exotichost_occurrence)) {
    for (states in obj[[regs]]) {
      if (df[which(df$host_name==specs),states]) exotichost_occurrence[specs,regs] <- exotichost_occurrence[specs,regs]+1
      ## need to multiply by relative land area... but can just use if > 0 for now...
    }
  }
}

# this hasn't changed yet either

#write.csv(exotichost_occurrence, 'Host_list/Output/exotic_host_occurrence_matrix.csv')

trees.countries <- read.csv('Host_list/Output/trees.countries.all.translated.csv')

exotic_trees <- expand.grid(colnames(exotichost_occurrence), names(obj)%>% setdiff(colnames(exotichost_occurrence)))
names(exotic_trees)<-c("Invaded.Range.in.North.America.Cluster","Geographic.Origin.Cluster")
exotic_trees$number_of_exotic_trees <- 0

# here can use df$scientificName instead

# added this now but may be unneccessary first time around
df[is.na(df$scientificName),'scientificName'] <- df[is.na(df$scientificName),'host_name']

rownames(exotichost_occurrence)
df$scientificName

unique(exotic_trees$Geographic.Origin.Cluster)

#names(table(df$scientificName)[table(df$scientificName)>1])

trees.list3 <- list()

for (i in 1:dim(exotic_trees)[1]) {
  exotic_trees$number_of_exotic_trees[i] <- sum(
    na.omit(df$scientificName[
      exotichost_occurrence[,exotic_trees$Invaded.Range.in.North.America.Cluster[i]]>0]) %in%
    unique((trees.countries %>%
       filter(country_name %in% obj[[as.character(exotic_trees$Geographic.Origin.Cluster[i])]]))$scientificName))

  
  current.origin <- exotic_trees[i,'Geographic.Origin.Cluster']
  
  if(!(current.origin %in% names(trees.list3))) {
    trees.list3[[length(trees.list3) + 1]] <- character()
    names(trees.list3)[[length(trees.list3)]] <- as.character(current.origin)
  }
  
  #add any new trees to the list
  trees.list3[[exotic_trees[i,'Geographic.Origin.Cluster']]] <-
    
    # this list represents the exotic trees from the current Geographic.Origin.Cluster
    # and present in the current location (Europe or Australia) from the expanded grid
    intersect(
      na.omit(df$scientificName[
        exotichost_occurrence[,exotic_trees$Invaded.Range.in.North.America.Cluster[i]]>0]),
      unique((trees.countries %>%
                filter(country_name %in% obj[[as.character(exotic_trees$Geographic.Origin.Cluster[i])]]))$scientificName)
    ) %>%
    
    # now we only want to retain the ones that are new to the list, and add them to the list
    union(trees.list3[[exotic_trees[i,'Geographic.Origin.Cluster']]]) %>% unique()
  
  #  and then later we can sum them all up
  
}

write.csv(exotic_trees, "Covariates_Input/exotic_trees.csv")
save(trees.list3, file="Covariates_input/NA_trees_list.RData")
