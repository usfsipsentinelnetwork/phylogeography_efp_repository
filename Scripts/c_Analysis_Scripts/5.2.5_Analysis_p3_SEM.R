library(lavaan)
library(dplyr)
library(semTools)
library(semPlot)
library(semptools)
library(tidySEM)

library(Cairo)
library(devEMF)

#load("Final_output/Final_analysis_data_pathogens-transformations-Sept-2024.RData")
load("Final_output/Final_analysis_data_pathogens-transformations-March4-2025-pnas.RData" )

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

summary(fit1, rsquare=T, fit.measures = T)
s.not <- summary(sem(model=model1,
                     data = to_analyze_scaled %>% filter(continentNAm == 0),
   #                 maxiter=10000,
            orthogonal=T, meanstructure=T), rsquare=T)

s.not$pe

summary(fit1.1, rsquare=T, fit.measures = T)
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
        sizeMan = 14,
        sizeMan2=14,
        edge.label.cex = 2.25,
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
        sizeMan = 14,
        sizeMan2=14,
        edge.label.cex = 2.25,
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
save(semp1, fit1, semp2, fit1.1, file="Figures_Sept2024/Fig_semplots_march2025.RData")

windows(12,5)
par(mfrow = c(1,2), oma=c(1,1,1,1), xpd =T)
plot(semp1); title(main = "A", outer=T, adj=0); text(1.25,-.375,
   paste(c('\u03c7\u00b2(mod)','p(\u03c7\u00b2mod)',
           '\u03c7\u00b2(mod|base)','p(\u03c7\u00b2mod|base)',
           'CFI','TLI','RMSEA','SRMR'),"=",
         fitMeasures(fit1)[c('chisq','pvalue',
                             'baseline.chisq','baseline.pvalue',
                             'cfi','tli','rmsea','srmr')]%>%round(2)%>%format(2),
         collapse='\n') %>%
     paste("R\u00b2(patho) = ", summary(fit1,rsquare=T)$pe$est[22]%>%round(2)%>%format(2),'\n',.,sep=''),
   cex=.8, adj=c(1,1), family='mono')
plot(semp2); title(main = "B", outer=T, adj=0.5); text(1.25,-.375,
   paste(c('\u03c7\u00b2(mod)','p(\u03c7\u00b2mod)',
           '\u03c7\u00b2(mod|base)','p(\u03c7\u00b2mod|base)',
           'CFI','TLI','RMSEA','SRMR'),"=",
         fitMeasures(fit1.1)[c('chisq','pvalue',
                             'baseline.chisq','baseline.pvalue',
                             'cfi','tli','rmsea','srmr')]%>%round(2)%>%format(2),
         collapse='\n') %>%
     paste("R\u00b2(patho) = ", summary(fit1.1,rsquare=T)$pe$est[22]%>%round(2)%>%format(2),'\n',.,sep=''),
   cex=.8, adj=c(1,1), family='mono')

svg("Final_figures/Fig5_3-4-25.svg", 12, 5, family="Cambria Math");
par(mfrow = c(1,2), oma=c(1,1,1,1), xpd =T)
plot(semp1); title(main = "A", outer=T, adj=0); text(1.25,-.375,
   paste(c('\u03c7\u00b2(mod)','p(\u03c7\u00b2mod)',
           '\u03c7\u00b2(mod|base)','p(\u03c7\u00b2mod|base)',
           'CFI','TLI','RMSEA','SRMR'),"=",
         fitMeasures(fit1)[c('chisq','pvalue',
                             'baseline.chisq','baseline.pvalue',
                             'cfi','tli','rmsea','srmr')]%>%round(2)%>%format(2),
         collapse='\n') %>%
     paste("R\u00b2(patho) = ", summary(fit1,rsquare=T)$pe$est[22]%>%round(2)%>%format(2),'\n',.,sep=''),
   cex=.8, adj=c(1,1), family='mono')
plot(semp2); title(main = "B", outer=T, adj=0.5); text(1.25,-.375,
   paste(c('\u03c7\u00b2(mod)','p(\u03c7\u00b2mod)',
           '\u03c7\u00b2(mod|base)','p(\u03c7\u00b2mod|base)',
           'CFI','TLI','RMSEA','SRMR'),"=",
         fitMeasures(fit1.1)[c('chisq','pvalue',
                             'baseline.chisq','baseline.pvalue',
                             'cfi','tli','rmsea','srmr')]%>%round(2)%>%format(2),
         collapse='\n') %>%
     paste("R\u00b2(patho) = ", summary(fit1.1,rsquare=T)$pe$est[22]%>%round(2)%>%format(2),'\n',.,sep=''),
   cex=.8, adj=c(1,1), family='mono')
dev.off()

emf("Final_figures/Fig5_3-4-25.emf", 12, 5, family="Cambria Math", emfPlus=TRUE);
par(mfrow = c(1,2), oma=c(1,1,1,1), xpd =T)
plot(semp1); title(main = "A", outer=T, adj=0); text(1.25,-.375,
   paste(c('\u03c7\u00b2(mod)','p(\u03c7\u00b2mod)',
           '\u03c7\u00b2(mod|base)','p(\u03c7\u00b2mod|base)',
           'CFI','TLI','RMSEA','SRMR'),"=",
         fitMeasures(fit1)[c('chisq','pvalue',
                             'baseline.chisq','baseline.pvalue',
                             'cfi','tli','rmsea','srmr')]%>%round(2)%>%format(2),
         collapse='\n') %>%
     paste("R\u00b2(patho) = ", summary(fit1,rsquare=T)$pe$est[22]%>%round(2)%>%format(2),'\n',.,sep=''),
   cex=.8, adj=c(1,1), family='mono')
plot(semp2); title(main = "B", outer=T, adj=0.5); text(1.25,-.375,
   paste(c('\u03c7\u00b2(mod)','p(\u03c7\u00b2mod)',
           '\u03c7\u00b2(mod|base)','p(\u03c7\u00b2mod|base)',
           'CFI','TLI','RMSEA','SRMR'),"=",
         fitMeasures(fit1.1)[c('chisq','pvalue',
                             'baseline.chisq','baseline.pvalue',
                             'cfi','tli','rmsea','srmr')]%>%round(2)%>%format(2),
         collapse='\n') %>%
     paste("R\u00b2(patho) = ", summary(fit1.1,rsquare=T)$pe$est[22]%>%round(2)%>%format(2),'\n',.,sep=''),
   cex=.8, adj=c(1,1), family='mono')
dev.off()

summary(fit1, fit.measures = TRUE, rsquare = TRUE, standardized = TRUE)$pe
summary(fit1.1, fit.measures = TRUE, rsquare = TRUE, standardized = TRUE)

# not doing panels

