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

#load("Final_output/Final_analysis_data_pathogens-transformations-Sept-2024.RData")
#load("Final_output/Final_analysis_data_pathogens-transformations-March1-2025-pnas.RData")
load("Final_output/Final_analysis_data_pathogens-transformations-March4-2025-pnas.RData")

###############################
###############################
##                          ###
##   zero inflated models   ###
##                          ###
###############################
###############################

#models06<-list(
#  model_0_centered,
#  model_1_centered,
#  model_2_centered,
#  model_3_centered,
#  model_4_centered,
#  model_5_centered,
#  model_6_centered,
#  model_final_zinf_centered)

####### BEST MODEL w/ CENTERED DATA

model_final_zinf_centered <-
  glmmTMB(
    Number.of.pests ~
      root4_trade_centered+
      Climate.Distance_centered +
      log_numexotics_centered +
      NMDS.Distance_centered +
      NMDS.Distance_centered:log_numexotics_centered +
      (1|continent),
    ziformula= ~1,
    data=to_analyze,
    family="poisson"
  )

#summary(model_final_zinf_centered)

##### centered

model_0_centered <-
  glmmTMB(
    Number.of.pests ~ 1,
    data=to_analyze,
    family="poisson"
    #family="nbinom1"
  )
model_1_centered <-
  glmmTMB(
    Number.of.pests ~ 
      (1|continent),
    data=to_analyze,
    family="poisson"
    #family="nbinom1"
  )
model_2_centered <-
  glmmTMB(
    Number.of.pests ~ 
      root4_trade_centered +
      (1|continent),
    data=to_analyze,
    family="poisson"
    #family="nbinom1"
  )
model_3_centered <-
  glmmTMB(
    Number.of.pests ~ 
      Climate.Distance_centered +
      (1|continent),
    data=to_analyze,
    family="poisson"
    #family="nbinom1"
  )
model_4_centered <-
  glmmTMB(
    Number.of.pests ~ 
      Climate.Distance_centered +
      root4_trade_centered +
      (1|continent),
    data=to_analyze,
    family="poisson"
    #family="nbinom1"
  )
model_5_centered <-
  glmmTMB(
    Number.of.pests ~ 
      Climate.Distance_centered +
      root4_trade_centered +
      NMDS.Distance_centered+
      (1|continent),
    data=to_analyze,
    family="poisson"
    #family="nbinom1"
  )
model_6_centered <-
  glmmTMB(
    Number.of.pests ~ 
      Climate.Distance_centered +
      root4_trade_centered +
      NMDS.Distance_centered+
      log_numexotics_centered +
      (1|continent),
    data=to_analyze,
    family="poisson"
    #family="nbinom1"
  )

summary(model_6_centered)

# save model RData
#save(model_final_zinf_centered, file='Final_output/zinf_model_final.RData')

#save(model_final_zinf_centered, file='Final_output/zinf_model_final_march_2025.RData') # need to resave this with previous version of the script
save(model_final_zinf_centered, file='Final_output/zinf_model_final_march4_2025.RData')

################################################
######### COMPARISON OF ALL THE MODELS #########
################################################

to_analyze$continent <- as.factor(to_analyze$continent)

################# ALL ######################

models06<-list(
  model_0_centered,          # intercept
  model_1_centered,          # intercept + continent
  model_2_centered,          # intercept + continent + trade
  model_3_centered,          # intercept + continent + climate
  model_4_centered,          # intercept + continent + trade + climate
  model_5_centered,          # intercept + continent + trade + climate + nmds
  model_6_centered,          # intercept + continent + trade + climate + nmds + exotics
  model_final_zinf_centered) # intercept + continent + trade + climate + nmds * exotics

# order: NA,8,7,6,5,4,2

####### ouput table - GLM and ANCOVA

# FORMULAS

formulas <- (sapply(models06, FUN=(function(x) as.character(formula(x))[3]%>%
                                    gsub("root4_trade","Trade",.)%>%
                                    gsub("Climate.Distance","Climate",.)%>%
                                    gsub("\\(1 \\| continent\\) \\+ log_spatauto","SpatAuto + Continent(Random)",.)%>%
                                    #gsub("area_sink","A(sink)",.)%>%
                                    #gsub("log\\(sink.div_to_area\\)", "D:A(sink)", .)%>%
                                    gsub("NMDS.Distance","Phylogeography",.)%>%
                                    gsub("log_numexotics","Exotic Trees",.)),
                   simplify='vector'))[c(NA,8,7,6,5,4,2)]


# AIC & chisq

aov06 <- 
  rbind(
    NA,
    anova(model_0_centered,model_1_centered)[2,c('Chisq','Df','Pr(>Chisq)')],
    anova(model_1_centered,model_2_centered)[2,c('Chisq','Df','Pr(>Chisq)')],
    anova(model_1_centered,model_3_centered)[2,c('Chisq','Df','Pr(>Chisq)')],
    anova(model_3_centered,model_4_centered)[2,c('Chisq','Df','Pr(>Chisq)')],
    anova(model_4_centered,model_5_centered)[2,c('Chisq','Df','Pr(>Chisq)')],
    anova(model_4_centered,model_6_centered)[2,c('Chisq','Df','Pr(>Chisq)')],
    anova(model_6_centered,model_final_zinf_centered)[2,c('Chisq','Df','Pr(>Chisq)')])

car::Anova(model_final_zinf_centered)

AIC06 <- sapply(models06, FUN=AIC)
W06 <- qpcR::akaike.weights(AIC06)

explicit_zinf_minus_trade <-
  glmmTMB(
    Number.of.pests ~
      #root4_trade_centered+
      Climate.Distance_centered +
      log_numexotics_centered +
      NMDS.Distance_centered +
      NMDS.Distance_centered:log_numexotics_centered +
      (1|continent),
    start = list(
      beta = model_final_zinf_centered$fit$par[1:6][-2],
      betazi = model_final_zinf_centered$fit$par[7],
      theta  = model_final_zinf_centered$fit$par[8],
      b = model_final_zinf_centered$fit$parfull[8:11]
    ),
    control = glmmTMBControl(optCtrl = list(iter.max=0)),
    ziformula= ~1,
    data=to_analyze,
    family="poisson"
  ) #%>% (performance::r2_zeroinflated)

explicit_zinf_minus_climate <-
  glmmTMB(
    Number.of.pests ~
      root4_trade_centered+
      #Climate.Distance_centered +
      log_numexotics_centered +
      NMDS.Distance_centered +
      NMDS.Distance_centered:log_numexotics_centered +
      (1|continent),
    start = list(
      beta = model_final_zinf_centered$fit$par[1:6][-3],
      betazi = model_final_zinf_centered$fit$par[7],
      theta  = model_final_zinf_centered$fit$par[8],
      b = model_final_zinf_centered$fit$parfull[8:11]
    ),
    control = glmmTMBControl(optCtrl = list(iter.max=0)),
    ziformula= ~1,
    data=to_analyze,
    family="poisson"
  ) #%>% (performance::r2_zeroinflated)

explicit_zinf_minus_full_interaction <-
  glmmTMB(
    Number.of.pests ~
      root4_trade_centered+
      Climate.Distance_centered +
      #log_numexotics_centered +
      #NMDS.Distance_centered +
      #NMDS.Distance_centered:log_numexotics_centered +
      (1|continent),
    start = list(
      beta = model_final_zinf_centered$fit$par[c(1:3)],
      betazi = model_final_zinf_centered$fit$par[7],
      theta  = model_final_zinf_centered$fit$par[8],
      b = model_final_zinf_centered$fit$parfull[8:11]
    ),
    control = glmmTMBControl(optCtrl = list(iter.max=0)),
    ziformula= ~1,
    data=to_analyze,
    family="poisson"
  ) #%>% (performance::r2_zeroinflated)

explicit_zinf_minus_exotics <-
  glmmTMB(
    Number.of.pests ~
      root4_trade_centered+
      Climate.Distance_centered +
      #log_numexotics_centered +
      NMDS.Distance_centered +
      NMDS.Distance_centered:log_numexotics_centered +
      (1|continent),
    start = list(
      beta = model_final_zinf_centered$fit$par[1:6][-4],
      betazi = model_final_zinf_centered$fit$par[7],
      theta  = model_final_zinf_centered$fit$par[8],
      b = model_final_zinf_centered$fit$parfull[8:11]
    ),
    control = glmmTMBControl(optCtrl = list(iter.max=0)),
    ziformula= ~1,
    data=to_analyze,
    family="poisson"
  ) #%>% (performance::r2_zeroinflated)

explicit_zinf_minus_nmds <-
  glmmTMB(
    Number.of.pests ~
      root4_trade_centered+
      Climate.Distance_centered +
      log_numexotics_centered +
      #NMDS.Distance_centered +
      NMDS.Distance_centered:log_numexotics_centered +
      (1|continent),
    start = list(
      beta = model_final_zinf_centered$fit$par[1:6][-5],
      betazi = model_final_zinf_centered$fit$par[7],
      theta  = model_final_zinf_centered$fit$par[8],
      b = model_final_zinf_centered$fit$parfull[8:11]
    ),
    control = glmmTMBControl(optCtrl = list(iter.max=0)),
    ziformula= ~1,
    data=to_analyze,
    family="poisson"
  ) #%>% (performance::r2_zeroinflated)

explicit_zinf_minus_int <-
  glmmTMB(
    Number.of.pests ~
      root4_trade_centered+
      Climate.Distance_centered +
      log_numexotics_centered +
      NMDS.Distance_centered +
      #NMDS.Distance_centered:log_numexotics_centered +
      (1|continent),
    start = list(
      beta = model_final_zinf_centered$fit$par[c(1:5)],
      betazi = model_final_zinf_centered$fit$par[7],
      theta  = model_final_zinf_centered$fit$par[8],
      b = model_final_zinf_centered$fit$parfull[8:11]
    ),
    control = glmmTMBControl(optCtrl = list(iter.max=0)),
    ziformula= ~1,
    data=to_analyze,
    family="poisson"
  ) #%>% (performance::r2_zeroinflated)

explicit_zinf_minus_random_effects <-
  glmmTMB(
    Number.of.pests ~
      root4_trade_centered+
      Climate.Distance_centered +
      log_numexotics_centered +
      NMDS.Distance_centered +
      NMDS.Distance_centered:log_numexotics_centered,# +
      #(1|continent),
    start = list(
      beta = model_final_zinf_centered$fit$par[c(1:6)],
      betazi = model_final_zinf_centered$fit$par[7]#,
 #     theta  = model_final_zinf_centered$fit$par[8]#,
 #     b = model_final_zinf_centered$fit$parfull[8:11]
    ),
    control = glmmTMBControl(optCtrl = list(iter.max=0)),
    ziformula= ~1,
    data=to_analyze,
    family="poisson"
  ) #%>% (performance::r2_zeroinflated)

only_random_effects <-
  glmmTMB(
    Number.of.pests ~ 1 +
      #root4_trade_centered+
      #Climate.Distance_centered +
      #log_numexotics_centered +
      #NMDS.Distance_centered +
      #NMDS.Distance_centered:log_numexotics_centered,# +
    (1|continent),
    start = list(
      beta = model_final_zinf_centered$fit$par[1],
      betazi = model_final_zinf_centered$fit$par[7],
      theta  = model_final_zinf_centered$fit$par[8],
      b = model_final_zinf_centered$fit$parfull[8:11]
    ),
    control = glmmTMBControl(optCtrl = list(iter.max=0)),
    ziformula= ~1,
    data=to_analyze,
    family="poisson"
  ) #%>% (performance::r2_zeroinflated)

r2_choices <- 
  sapply(
    list(
      model_final_zinf_centered,
      explicit_zinf_minus_full_interaction,
      explicit_zinf_minus_int,
      explicit_zinf_minus_exotics,
      explicit_zinf_minus_nmds,
      explicit_zinf_minus_trade,
      explicit_zinf_minus_climate
    ),
    function (x) c(
      nakagawa_marginal_loo = ifelse('theta' %in% names(x$fit$par),
        as.numeric(performance::r2_nakagawa(x)$R2_marginal), NA),
      nakagawa_conditional_loo = ifelse('theta' %in% names(x$fit$par),
        as.numeric(performance::r2_nakagawa(x)$R2_conditional), NA),
      efron_loo = performance::r2_efron(x),
      zinf_loo = as.numeric(performance::r2_zeroinflated(x)$R2)
    )) %>% t %>%
  (function (x) {
    data.frame(
      x[-1,],
      full_model_efron = x[1,3],
      full_model_zinf = x[1,4],
      row.names = c(
        "Full Interaction",
        "Phylogeo*Exotics",
        "Exotic Trees",
        "Phylogeography",
        "Trade",
        "Climate")
    )
  })%>%
  mutate(efron_R2_partial = full_model_efron- efron_loo) %>%
  mutate(zinf_R2_partial = full_model_zinf- zinf_loo)

r2_choices_2 <-
  r2_choices %>% dplyr::select('efron_R2_partial', 'zinf_R2_partial')

r2_choices_2['Fixed',] <- colSums(r2_choices_2[c(1,5,6),])
r2_choices_2['Random',] <- r2_choices[1,c('full_model_efron','full_model_zinf')] - r2_choices_2['Fixed',]
r2_choices_2['Total',] <- r2_choices_2['Fixed',] + r2_choices_2['Random',]

results.glm <- data.frame(parameter=
                      c("Full Interaction",
                        "Phylogeo*Exotics",
                        "Exotic Trees",
                        "Phylogeography",
                        "Trade",
                        "Climate",
                        "Continent",
                        "Fixed",
                        "Total"),
                   glm.delta.AIC = W06$deltaAIC[c(NA,8,7,6,5,4,2,NA,NA)],
                   chisq         = aov06$Chisq[c(NA,8,7,6,5,4,2,NA,NA)],
                   glm.pchisq    = aov06$`Pr(>Chisq)`[c(NA,8,7,6,5,4,2,NA,NA)]) %>%
                   cbind(r2=r2_choices_2[c(
                     "Full Interaction",
                     "Phylogeo*Exotics",
                     "Exotic Trees",
                     "Phylogeography",
                     "Trade",
                     "Climate",
                     "Random",
                     "Fixed",
                     "Total"
                   ),1])

#write.csv(results.glm[c(1:2,4,3,5:9),], "Final_output/results.glm.centered.stats.Sept_2024_pois_2.csv", row.names=F)
write.csv(results.glm, "Final_output/results.glm.centered.stats.March4_2025_pois_2.csv", row.names=F)

Anova(model_final_zinf_centered)
### OLS ANCOVA

### type I anova

t1.factor <- lm(
  logNumber.of.pests ~
    continent +
    NMDS.Distance_centered *
    log_numexotics_centered +
    root4_trade_centered +
    Climate.Distance_centered,
  data=to_analyze
)

ofs2 <-
  as.numeric((c(0,t1.factor$coefficients[c('continentEurasian.Palearctic','continentHawaii','continentNorth.America')]))[to_analyze$continent])

t1.ofs2 <- lm(
  logNumber.of.pests ~
    NMDS.Distance_centered *
    log_numexotics_centered +
    root4_trade_centered +
    Climate.Distance_centered+
    offset(ofs2),
  data=to_analyze
)

anova_t1.factor <- anova(t1.factor)
anova_t1.ofs2<-anova(t1.ofs2)

results.ancova <-
  data.frame(parameter=
               c(
                 "Full Interaction",
                 "Phylogeo*Exotics",
                 "Phylogeography",
                 "Exotic Trees",
                 "Trade",
                 "Climate",
                 "Continent",
                 "Total",
                 "Fixed"
                 ),
             ancova.F = anova_t1.factor$`F value`[c(NA,6,2:5,1,NA,NA)],
             ancova.pF = anova_t1.factor$`Pr(>F)`[c(NA,6,2:5,1,NA,NA)],
             ancova.SS = anova_t1.factor$`Sum Sq`[c(NA,6,2:5,1)] %>%
               c(sum(anova_t1.factor$`Sum Sq`), sum(anova_t1.factor$`Sum Sq`[-7])))

#1 - anova_t1.factor$`Sum Sq`[7]/sum(anova_t1.factor$`Sum Sq`)

results.ancova$ancova.r2part <- results.ancova$ancova.SS / sum(anova_t1.factor$`Sum Sq`)

#write.csv(results.ancova, "Final_output/results.ancova.offset.stats.Sept_2024_2.csv", row.names=F)
write.csv(results.ancova, "Final_output/results.ancova.offset.stats.March4_2025.csv", row.names=F)

####### BEST MODEL w/ CENTERED DATA

model_progressive_zinf_centered <-
  glmmTMB(
    Number.of.pests ~
      
      root4_trade_centered * (
        log_numexotics_centered +
          NMDS.Distance_centered *
          Climate.Distance_centered),#+
      
      #(1|continent),
    
    ziformula= ~1,
    data=to_analyze,
    family="poisson"
  )

anova(model_progressive_zinf_centered, model_final_zinf_centered)

# that improved model fit
summary(model_progressive_zinf_centered)

##########################################

### INTERACTIVE MODELS

### DO THIS LAST - TOGETHER WITH FIGURE


######################

save.image("Final_output/linearmodels_march4.2025.RData")
