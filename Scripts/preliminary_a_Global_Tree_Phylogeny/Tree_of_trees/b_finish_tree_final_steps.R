
setwd("C:/Users/GeoffreyWilliams/OneDrive - USDA/Project_folders_Personal/Projects/Sentinels/Retrospective")

library(V.PhyloMaker2)

load("Host_list/Output/global_tree_phylogeny.RData")

##### NOTE NEVER SAVED FINAL TREE DATA AFTER THIS STEP - could be an issue

#nodes_final <- build.nodes.1(tree=tr3$scenario.3, tips=tr3$species.list[,c("species","genus","family")])

length(tr3$scenario.3$tip.label)
dim(tr3$species.list)

#> length(tr3$scenario.3$tip.label)
#[1] 56981
#> dim(tr3$species.list)
#[1] 56881     4

length(unique(tr3$scenario.3$tip.label)) #56881 - there are duplicates

reps <- table(tr3$scenario.3$tip.label)[table(tr3$scenario.3$tip.label) > 1]

x<-NULL
for (a in 1:length(reps)) x <- c(x, seq(0,reps[a]-1) %>% as.character)

x[which(x=="0")]<-""

dups <- table(tr3$scenario.3$tip.label)[table(tr3$scenario.3$tip.label) > 1] %>% names

dups.w <- which(tr3$scenario.3$tip.label %in% dups)

dups.o <- order(tr3$scenario.3$tip.label[dups.w])

to_replace <- paste(tr3$scenario.3$tip.label[dups.w][dups.o], x, sep="")

tr3$scenario.3$tip.label[dups.w][dups.o]<-to_replace

tips_to_prune <- setdiff(tr3$scenario.3$tip.label[dups.w], names(reps))

tr3.2 <- drop.tip(tr3$scenario.3, tips_to_prune)

save(tr3.2, file="Host_list/Output/global_tree_phylogeny_pruned.RData")

#str(tr3.2)
#str(tr3$species.list)

#setdiff(tr3.2$tip.label, gsub(" ","_",tr3$species.list$species))
#setdiff(gsub(" ","_",tr3$species.list$species),tr3.2$tip.label)

tips<-data.frame(species=gsub(" ","_",tr3$species.list$species), tr3$species.list[,c("genus","family")])

nodes_final <- build.nodes.1(tree=tr3.2, tips=tips)



save(nodes_final, file="Host_list/Output/nodes_final.Rdata")