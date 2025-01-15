setwd("C:/Users/GeoffreyWilliams/OneDrive - USDA/Project_folders_Personal/Projects/Sentinels/Retrospective")
library(dplyr)
library(tidyr)

setwd("Host_list/USDA_DB/")
lists <- dir()

#listnames<-strcapture("(.+)\\.csv", lists, data.frame(chr=character(length=length(lists))))$chr

#### US Host Lists

hosts <- NULL

for (i in lists) {
  setwd(i)
  filenames <- dir()
  for (j in filenames) {
    csv <- read.csv(j, skip=3)$Scientific.Name
    hosts <- union(hosts, csv)
  }
  setwd("..")
}

df <- data.frame(host= hosts)

for (n in 1:length(lists)) {
  setwd(lists[n])
  filenames <- dir()
  hosts <- NULL
  
  for (j in filenames) {
    csv <- read.csv(j, skip=3)$Scientific.Name
    hosts <- union(hosts, csv)
  }
  df <- df %>% mutate(n = host %in% hosts)
  names(df)[n+1] <- lists[n]
  setwd("..")
}

df <- arrange(df, host)

df <- extract(df, col= host, remove= F, regex = "^([a-zA-Z\\×]+[ ][a-zA-Z\\-\\×]+) (.*)", into = c("host_name", "host_authority"))
df <- extract(df, col = host_name, remove = F, regex="^([a-zA-Z]+)[ ]", into = "host_genus")

#write.csv(df, "host_trees_byregion_notsynonymized.csv", row.names=F)

##### Extract names, add other relevant host names,

# tree species names from global plant list (BGCI)

# species names for V.PhyloPlot2

##### make a translation table (to be used for searching)

# add annotations to correct incorrect synonymizations - from notes

##### and add column with accepted names

# translation table for fungal names

# invasive list
# ars

