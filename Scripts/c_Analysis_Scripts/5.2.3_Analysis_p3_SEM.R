library(lavaan)
library(dplyr)
library(semTools)
library(semPlot)
library(semptools)
library(tidySEM)

library(Cairo)
library(devEMF)

load("Final_output/Final_analysis_data_pathogens-transformations-Sept-2024.RData")

#logNumber.of.pests ~ Climate.Distance_centered + root4_trade_centered + NMDS.Distance_centered + log_numexotics_centered + NMDS.Distance_centered:log_numexotics_centered + continentHII + continentAus + continentEur
#NMDS.Distance_centered ~ Climate.Distance_centered
#log_numexotics_centered ~ Climate.Distance_centered
#root4_trade_centered ~~ Climate.Distance_centered
#continentHII ~~ continentHII
#continentAus ~~ continentAus
#continentEur ~~ continentEur

to_analyze_scaled <-
  scale(dplyr::select(to_analyze,
                      logNumber.of.pests,
                      NMDS.Distance,
                      root4_trade,
                      log_numexotics,
                      Climate.Distance,
                      NMDS.Distance)) %>%
  cbind(dplyr::select(to_analyze,
                      continentNAm,
                      continentHII,
                      continentAus,
                      continentEur,
                      continent)) %>%
  indProd(c(2,4,5), 3, match=F)

model1 <- ' 
  # latent variable definitions
    potential_richness =~ log_numexotics + NMDS.Distance + Climate.Distance

  # regressions
    logNumber.of.pests ~ potential_richness + root4_trade

  # residual correlations
  
  NMDS.Distance ~~ Climate.Distance
  NMDS.Distance ~~ root4_trade
  #log_numexotics ~~ NMDS.Distance
'

model1.1 <- ' 
  # latent variable definitions
    establishment_barrier =~ Climate.Distance + NMDS.Distance + log_numexotics

  # regressions
    logNumber.of.pests ~ establishment_barrier + root4_trade #+ continentHII + continentAus + continentEur

  # residual correlations

  #NMDS.Distance ~~ Climate.Distance
  NMDS.Distance ~~ root4_trade
  log_numexotics ~~ NMDS.Distance
  
'

fit1 <- sem(model=model1, data = to_analyze_scaled,
            orthogonal=T, meanstructure=T)

fit1.1 <- sem(model=model1.1, data = to_analyze_scaled,
            orthogonal=T, meanstructure=T)

summary(fit1, rsquare=T)
s.not <- summary(sem(model=model1, data = to_analyze_scaled %>% filter(continentNAm == 0),
            orthogonal=T, meanstructure=T), rsquare=T)

s.not$pe

summary(fit1.1, rsquare=T)
s.not2 <- summary(sem(model=model1.1, data = to_analyze_scaled %>% filter(continentNAm == 0),
            orthogonal=T, meanstructure=T), rsquare=T)

s.not2$pe

###############
## PLOT?? #####
###############

make_figure_layout <-
  function(fit, mat=1, grouped=F) {
    m1 <- matrix(
      data =
        c("NMDS.Distance",          "Climate.Distance",       "log_numexotics",
          "root4_trade",            "potential_richness",     NA,
          NA,                       "logNumber.of.pests",     NA
        ),
      3,3,byrow=T
    )
    
    m2 <- #matrix(
    #  data =
    #    c("NMDS.Distance",          "log_numexotics",           "Climate.Distance",       
    #      "root4_trade",            "establishment_barrier",     NA,
    #      NA,                       "logNumber.of.pests",     NA
    #    ),
    #  3,3,byrow=T
    #)
    
      matrix(
        data =
          c("NMDS.Distance",          "Climate.Distance",       "log_numexotics",
            "root4_trade",            "establishment_barrier",     NA,
            NA,                       "logNumber.of.pests",     NA
          ),
        3,3,byrow=T
      )
      
    if (mat == 1) m<- m1
    else if (mat == 2) m <- m2
    else m<- mat
    
    if (grouped)
      semp <- semPaths(
        fit,
#        what="stand",
        whatLabels="est",
        sizeMan = 20,
        sizeMan2=20,
        edge.label.cex = 2,
        style="ram",
        nCharNodes = 0,
        nCharEdges = 0,
        shapeMan="circle",
        residuals=F,
        layout=m,
        intercepts=F)

    else
      semp <- semPaths(
        fit,
        what="stand",
        whatLabels="est",
        sizeMan = 20,
        sizeMan2=20,
        edge.label.cex = 1.5,
        style="ram",
        nCharNodes = 0,
        nCharEdges = 0,
        shapeMan="circle",
        residuals=F,
        layout=m,
        intercepts=F)
        
    semp
  }

make_figure_format <-
  function(semp, fit, mat=1, grouped =F, fit_summary) {
    if(mat==1) semp <- semp %>%
        change_node_label(c(
          root4_trade="Imports",
          Climate.Distance="Climatic\nDistance",
          NMDS.Distance="Phylo-\ngeographic\nDistance",
          logNumber.of.pests="Established\nPathogen\nRichness",
          log_numexotics="Non-native\nTree\nRichness",
          potential_richness="Potential\nPathogen\nRichness"
        ))
    
    else if (mat==2) semp <- semp %>%
        change_node_label(c(
          root4_trade="Imports",
          Climate.Distance="Climatic\nDistance",
          NMDS.Distance="Phylo-\ngeographic\nDistance",
          logNumber.of.pests="Established\nPathogen\nRichness",
          log_numexotics="Non-native\nTree\nRichness",
          establishment_barrier="Establish-\nment\nBarrier"
        ))
    
    if ( grouped ) {
      labs.orig <- semp$graphAttributes$Edges$labels
      
      semp$graphAttributes$Edges$labels <-
        paste(semp$graphAttributes$Edges$labels,
              ifelse((fit_summary$pvalue[1:9] <= 0.05) & !is.na(fit_summary$pvalue[1:9]), "*", ""),sep="")
      
      semp$graphAttributes$Edges$labels <-
        paste(semp$graphAttributes$Edges$labels,
              ifelse((fit_summary$pvalue[1:9] <= 0.01) & !is.na(fit_summary$pvalue[1:9]), "*", ""),sep="")
      
      semp$graphAttributes$Edges$labels <-
        paste(semp$graphAttributes$Edges$labels,
              ifelse((fit_summary$pvalue[1:9] <= 0.001) & !is.na(fit_summary$pvalue[1:9]), "*", ""),sep="")
      
      }
    else {semp <- semp %>% mark_sig(fit)}
    
    semp$Edgelist$from[1:3] <- 1:3
    semp$Edgelist$to[1:3] <-c(6,6,6)
    semp$graphAttributes$Nodes$width[6]<-20
    semp$graphAttributes$Nodes$height[6]<-20
    
    if (mat == 1) semp$graphAttributes$Edges$curve[6:7]<-c(4,-4)
    else if (mat == 2) semp$graphAttributes$Edges$curve[6:7]<-c(-4,2.25)
    
    semp$graphAttributes$Edges$lty[6:9]<-3
    
    semp
  }

make_figure <-
  function(fit, mat=1) {

    semp <- make_figure_layout(fit, mat)
    
    semp <- make_figure_format(semp, fit, mat)
    
    semp
    
  }

semp1 <- make_figure(fit1, 1)
semp2 <- make_figure(fit1.1, 2)

windows(12,5)
par(mfrow = c(1,2), oma=c(1,1,1,1), xpd =T)
plot(semp1); title(main = "A", outer=T, adj=0)
plot(semp2); title(main = "B", outer=T, adj=0.5)

save(semp1, file="Figures_Sept2024/Fig6A_semplot1.RData")
save(semp2, file="Figures_Sept2024/Fig6B_semplot2.RData")

svg("../Figures/Fig5_9-18-24.svg", 12, 5, family="Cambria Math");
par(mfrow = c(1,2), oma=c(1,1,1,1), xpd =T)
plot(semp1); title(main = "A", outer=T, adj=0)
plot(semp2); title(main = "B", outer=T, adj=0.5)
dev.off()

emf("../Figures/Fig5_9-18-24.emf", 12, 5, family="Cambria Math", emfPlus=TRUE);
par(mfrow = c(1,2), oma=c(1,1,1,1), xpd =T)
plot(semp1); title(main = "A", outer=T, adj=0)
plot(semp2); title(main = "B", outer=T, adj=0.5)
dev.off()

summary(fit1, fit.measures = TRUE, rsquare = TRUE, standardized = TRUE)$pe
summary(fit1.1, fit.measures = TRUE, rsquare = TRUE, standardized = TRUE)

############################################
### II grouped analyses (panel figures) ###
############################################

model1_grouped <-
  '
  # latent variable definitions
    potential_richness =~ log_numexotics + NMDS.Distance + Climate.Distance

  # regressions
    logNumber.of.pests ~ potential_richness + root4_trade

  # residual correlations
  
  NMDS.Distance ~~ Climate.Distance
  NMDS.Distance ~~ root4_trade
  #log_numexotics ~~ NMDS.Distance 
'
fit1_grouped <- sem(model1_grouped, data=to_analyze_scaled,
                    group='continent', orthogonal=T, meanstructure=T)

fit1_grouped.2 <- sem(model1_grouped, data=to_analyze_scaled,
                        group='continentNAm', orthogonal=T, meanstructure=T)

model1.1_grouped <-
  '
  # latent variable definitions
    establishment_barrier =~ Climate.Distance + NMDS.Distance + log_numexotics

  # regressions
    logNumber.of.pests ~ establishment_barrier + root4_trade #+ continentHII + continentAus + continentEur

  # residual correlations

  #NMDS.Distance ~~ Climate.Distance
  NMDS.Distance ~~ root4_trade
  log_numexotics ~~ NMDS.Distance
'
fit1.1_grouped <- sem(model1.1_grouped, data=to_analyze_scaled,
                    group='continent', orthogonal=T, meanstructure=T)

fit1.1_grouped.2 <- sem(model1.1_grouped, data=to_analyze_scaled,
                      group='continentNAm', orthogonal=T, meanstructure=T)

###############
## PLOT?? #####
###############

make_figure_grouped <-
  function(fit, mat=1, fit_summary) {
    
    semp <- make_figure_layout(fit, mat, T)
    
    l <- length(semp)
    
    for (i in 1:l) {
      semp[[i]] <- make_figure_format(semp[[i]], fit, mat, T, fit_summary$pe %>% filter(group == i))
      #semp[[i]]$graphAttributes$Edges$color <- rep('black',9)
      #semp[[i]]$graphAttributes$Edges$width <- rep(2.5, 9)
    }
    
    semp
  }

s <- summary(fit1, rsquare=T)
s2 <- summary(fit1.1, rsquare=T)
s.grouped <- summary(fit1_grouped, rsquare=T)
s2.grouped <- summary(fit1.1_grouped, rsquare=T)
s.grouped.2 <- summary(fit1_grouped.2, rsquare=T)
s2.grouped.2 <- summary(fit1.1_grouped.2, rsquare=T)
attr(s.grouped$test, 'info')$group.label
attr(s.grouped.2$test, 'info')$group.label

rbind(
  cbind(group='all',latent='rich',s$pe%>% select('lhs','op','rhs','est','se','pvalue')),
  cbind(group='nam',latent='rich',s.grouped.2$pe %>% filter(block==2) %>% select('lhs','op','rhs','est','se','pvalue')),
  cbind(group='hieuau',latent='rich',s.grouped.2$pe %>% filter(block==1)%>% select('lhs','op','rhs','est','se','pvalue')),
  cbind(group='hi',latent='rich',s.grouped$pe %>% filter(block==3)%>% select('lhs','op','rhs','est','se','pvalue')),
  cbind(group='eu',latent='rich',s.grouped$pe %>% filter(block==2)%>% select('lhs','op','rhs','est','se','pvalue')),
  cbind(group='au',latent='rich',s.grouped$pe %>% filter(block==1)%>% select('lhs','op','rhs','est','se','pvalue')),
  cbind(group='all',latent='barrier',s2$pe%>% select('lhs','op','rhs','est','se','pvalue')),
  cbind(group='nam',latent='barrier',s2.grouped.2$pe %>% filter(block==2)%>% select('lhs','op','rhs','est','se','pvalue')),
  cbind(group='hieuau',latent='barrier',s2.grouped.2$pe %>% filter(block==1)%>% select('lhs','op','rhs','est','se','pvalue')),
  cbind(group='hi',latent='barrier',s2.grouped$pe %>% filter(block==3)%>% select('lhs','op','rhs','est','se','pvalue')),
  cbind(group='eu',latent='barrier',s2.grouped$pe %>% filter(block==2)%>% select('lhs','op','rhs','est','se','pvalue')),
  cbind(group='au',latent='barrier',s2.grouped$pe %>% filter(block==1)%>% select('lhs','op','rhs','est','se','pvalue'))
) %>% write.csv("Final_output/SEM_grouped_results.csv", row.names=F)

