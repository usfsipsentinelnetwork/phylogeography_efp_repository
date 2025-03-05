library(ggplot2)
library(dplyr)
library(doBy)
library(reshape2)
library(ggh4x)
library(scales)
library(devEMF)

library(tidyverse)
library(sjPlot)
library(car)
library(lme4)
library(cv)
library(qpcR)
library(MuMIn)
library(glmmTMB)
library(DHARMa)
library(see)
library(plsRglm)

############

#load("Final_output/linearmodels.RData")

#load("Final_output/Final_analysis_data_pathogens-transformations-Sept-2024.RData")
#load("Final_output/Final_analysis_data_pathogens-transformations-March1-2025-pnas.RData" )

load("Final_output/Final_analysis_data_pathogens-transformations-March4-2025-pnas.RData" )

###############################
###############################
##                          ###
##   PARTIAL LEAST SQ GLM   ###
##                          ###
###############################
###############################

#model_boot_and_R2 <- function (fmla, nt, dat, lnk = 'poisson', R=1000, nms) {
model_boot_and_R2 <- function (modpls, R=1000, nms) {
    
  #modpls <- plsRglm(
  #  fmla,
  #  data=dat,
  #  nt=nt,
  #  modele="pls-glm-family",
  #  family=lnk,
  #  pvals.expli=T)
  
  rock.bootYX<- bootplsglm(modpls, typeboot="plsmodel", R=1000, sim="antithetic")
  rownames(rock.bootYX$t0)<-nms
  
  results.bootYX <- summary(rock.bootYX, high.moments=T) %>% data.frame
  results.bootYX <-results.bootYX %>% mutate(logmed=exp(bootMed))
  row.names(results.bootYX) <- nms
  
  p.boot.plsyx<- sapply(1:length(nms), FUN= function (x) sum(rock.bootYX$t[,x]*sign(results.bootYX$bootMed)[x] < 0)/1000)
  names(p.boot.plsyx)<-nms
  
  return(list(model=modpls, boot=rock.bootYX, boot_summary=results.bootYX, boot.p = p.boot.plsyx))
}

pls_out_1 <-
  model_boot_and_R2(
    plsRglm(Number.of.pests ~ Climate.Distance_centered + root4_trade_centered + NMDS.Distance_centered * log_numexotics_centered,
            nt = 2, data = to_analyze, family='poisson', modele="pls-glm-family", pvals.expli=T),
    nms = c("Intercept", "Climate", "Trade", "NMDS.Dist", "Non-native trees", "Interaction"))

####################################

# no longer needed needed - RSS/TSS calculations -  its done elsewhere

###########################################################

cv.plsRglm(
  Number.of.pests ~ Climate.Distance_centered + root4_trade_centered + NMDS.Distance_centered * log_numexotics_centered,
  nt = 2, data = to_analyze %>% filter(continentNAm==1), family='poisson', modele="pls-glm-family",
  K=8,NK=100) %>% summary %>% cvtable %>% (function (x) {plot(x, type="CVQ2Chi2"); plot(x, type="CVPreChi2")})

pls_NAM <-
  model_boot_and_R2(
    plsRglm(Number.of.pests ~ Climate.Distance_centered + root4_trade_centered + NMDS.Distance_centered * log_numexotics_centered,
            nt = 2,
            data = to_analyze %>% filter(continentNAm==1),
            family='poisson', modele="pls-glm-family", pvals.expli=T),
    nms = c("Intercept", "Climate", "Trade", "NMDS.Dist", "Non-native trees", "Interaction"))

pls_NAM$model$RSS

#cv.plsRglm(
#  Number.of.pests ~ Climate.Distance_centered + root4_trade_centered + NMDS.Distance_centered * log_numexotics_centered,
#  nt = 5, data = to_analyze %>% filter(continentHII==1), family='poisson', modele="pls-glm-family",
#  K=8,NK=100) %>% summary %>% cvtable %>% (function (x) {plot(x, type="CVQ2Chi2"); plot(x, type="CVPreChi2")})

#pls_HII <-
#  model_boot_and_R2(
#    plsRglm(Number.of.pests ~ Climate.Distance_centered + root4_trade_centered + NMDS.Distance_centered * log_numexotics_centered,
#            nt = 1,
#            data = to_analyze %>% filter(continentHII==1),
#            family='poisson', modele="pls-glm-family", pvals.expli=T),
#    nms = c("Intercept", "Climate", "Trade", "NMDS.Dist", "Non-native trees", "Interaction"))

#cv.plsRglm(
#  Number.of.pests ~ Climate.Distance_centered + root4_trade_centered + NMDS.Distance_centered * log_numexotics_centered,
#  nt = 5, data = to_analyze %>% filter(continentEur==1), family='poisson', modele="pls-glm-family",
#  K=8,NK=100) %>% summary %>% cvtable %>% (function (x) {plot(x, type="CVQ2Chi2"); plot(x, type="CVPreChi2")})

#pls_EU <-
#  model_boot_and_R2(
#    plsRglm(Number.of.pests ~ Climate.Distance_centered + root4_trade_centered + NMDS.Distance_centered * log_numexotics_centered,
#            nt = 4,
#            # nt = 2,
#            data = to_analyze %>% filter(continentEur==1),
#            family='poisson', modele="pls-glm-family", pvals.expli=T),
#    nms = c("Intercept", "Climate", "Trade", "NMDS.Dist", "Non-native trees", "Interaction"))

#cv.plsRglm(
#  Number.of.pests ~ Climate.Distance_centered + root4_trade_centered + NMDS.Distance_centered * log_numexotics_centered,
#  nt = 5, data = to_analyze %>% filter(continentAus==1), family='poisson', modele="pls-glm-family",
#  K=8,NK=100) %>% summary %>% cvtable %>% (function (x) {plot(x, type="CVQ2Chi2"); plot(x, type="CVPreChi2")})

#pls_Au <-
#  model_boot_and_R2(
#    plsRglm(Number.of.pests ~ Climate.Distance_centered + root4_trade_centered + NMDS.Distance_centered * log_numexotics_centered,
#            nt = 2,
#            data = to_analyze %>% filter(continentAus==1),
#            family='poisson', modele="pls-glm-family", pvals.expli=T),
#    nms = c("Intercept", "Climate", "Trade", "NMDS.Dist", "Non-native trees", "Interaction"))

#make_model_predictions

# save results

save(
  pls_out_1,
  pls_NAM,
#  pls_HII,
#  pls_EU,
#  pls_Au,
  file= "Final_output/PLS-models-3-11-25.RData"
)
