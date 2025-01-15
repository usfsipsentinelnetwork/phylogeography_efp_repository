library(ggplot2)
library(dplyr)
#library(doBy)
#library(reshape2)
#library(ggh4x)
#library(scales)
library(devEMF)
library(tidyverse)
#library(sjPlot)
library(car)
#library(lme4)
library(cv)
#library(qpcR)
#library(MuMIn)
library(glmmTMB)
#library(DHARMa)
#library(see)
library(plsRglm)
library(cowplot)

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

load("Final_output/linearmodels.RData")
load('Final_output/zinf_model_final.RData')
load('Final_output/PLS-models-9-10-23.RData')




pls_not <-
  plsRglm(Number.of.pests ~ Climate.Distance_centered + root4_trade_centered + NMDS.Distance_centered * log_numexotics_centered,
          nt = 2, data = to_analyze %>% filter(continentNAm != 1),
          family='poisson', modele="pls-glm-family", pvals.expli=T)
pls_NAM <-
  plsRglm(Number.of.pests ~ Climate.Distance_centered + root4_trade_centered + NMDS.Distance_centered * log_numexotics_centered,
          nt = 2, data = to_analyze %>% filter(continentNAm == 1),
          family='poisson', modele="pls-glm-family", pvals.expli=T)

pls_not.boot <- model_boot_and_R2(pls_not,nms = c("Intercept", "Climate", "Trade", "NMDS.Dist", "Non-native trees", "Interaction"))
pls_NAM.boot <- model_boot_and_R2(pls_NAM,nms = c("Intercept", "Climate", "Trade", "NMDS.Dist", "Non-native trees", "Interaction"))

c(
  all=pls_out_1[c(3,4)],
  nam=pls_NAM.boot[c(3,4)],
  hieuau=pls_not.boot[c(3,4)],
  hi=pls_HII[c(3,4)],
  eu=pls_EU[c(3,4)],
  au=pls_Au[c(3,4)]
) %>% write.csv("Final_output/pls_params_pvals.csv", row.names=T)



dataNAm_modelNAm_pred <- predict(pls_NAM, type='response', newdata = pls_out_1$model$dataX %>% filter(to_analyze$continentNAm == 1))
dataNAm_modelnot_pred <- predict(pls_not, type='response', newdata = pls_out_1$model$dataX %>% filter(to_analyze$continentNAm == 1))
dataNAm_modelAll_pred <- predict(pls_out_1$model, type='response', newdata = pls_out_1$model$dataX %>% filter(to_analyze$continentNAm == 1))
datanot_modelNam_pred <- predict(pls_NAM, type='response', newdata = pls_out_1$model$dataX %>% filter(to_analyze$continentNAm == 0))
datanot_modelnot_pred <- predict(pls_not, type='response', newdata = pls_out_1$model$dataX %>% filter(to_analyze$continentNAm == 0))
datanot_modelAll_pred <- predict(pls_out_1$model, type='response', newdata = pls_out_1$model$dataX %>% filter(to_analyze$continentNAm == 0))
dataAll_modelNam_pred <- predict(pls_NAM, type='response', newdata = pls_out_1$model$dataX )
dataAll_modelnot_pred <- predict(pls_not, type='response', newdata = pls_out_1$model$dataX )
dataAll_modelAll_pred <- predict(pls_out_1$model, type='response')

dataNAm_modelNAm_pred_pearson <- cor(filter(to_analyze, continentNAm == 1)$Number.of.pests, dataNAm_modelNAm_pred) %>% round(2)
dataNAm_modelnot_pred_pearson <- cor(filter(to_analyze, continentNAm == 1)$Number.of.pests, dataNAm_modelnot_pred)%>% round(2)
dataNAm_modelAll_pred_pearson <- cor(filter(to_analyze, continentNAm == 1)$Number.of.pests, dataNAm_modelAll_pred)%>% round(2)
datanot_modelNam_pred_pearson <- cor(filter(to_analyze, continentNAm == 0)$Number.of.pests, datanot_modelNam_pred)%>% round(2)
datanot_modelnot_pred_pearson <- cor(filter(to_analyze, continentNAm == 0)$Number.of.pests, datanot_modelnot_pred)%>% round(2)
datanot_modelAll_pred_pearson <- cor(filter(to_analyze, continentNAm == 0)$Number.of.pests, datanot_modelAll_pred)%>% round(2)
dataAll_modelNam_pred_pearson <- cor(to_analyze$Number.of.pests, dataAll_modelNam_pred) %>% round(2)
dataAll_modelnot_pred_pearson <- cor(to_analyze$Number.of.pests, dataAll_modelnot_pred)%>% round(2)
dataAll_modelAll_pred_pearson <- cor(to_analyze$Number.of.pests, dataAll_modelAll_pred)%>% round(2)

dataNAm_modelNAm_pred_R2 <- 
  (1 - sum((dataNAm_modelNAm_pred - filter(to_analyze, continentNAm == 1)$Number.of.pests)^2)/
     sum((mean(filter(to_analyze, continentNAm == 1)$Number.of.pests) - filter(to_analyze, continentNAm == 1)$Number.of.pests)^2) )%>% round(2)
dataNAm_modelnot_pred_R2 <- 
  (1 - sum((dataNAm_modelnot_pred - filter(to_analyze, continentNAm == 1)$Number.of.pests)^2)/
     sum((mean(filter(to_analyze, continentNAm == 1)$Number.of.pests) - filter(to_analyze, continentNAm == 1)$Number.of.pests)^2) )%>% round(2)
dataNAm_modelAll_pred_R2 <- 
  (1 - sum((dataNAm_modelAll_pred - filter(to_analyze, continentNAm == 1)$Number.of.pests)^2)/
     sum((mean(filter(to_analyze, continentNAm == 1)$Number.of.pests) - filter(to_analyze, continentNAm == 1)$Number.of.pests)^2) )%>% round(2)
datanot_modelNam_pred_R2 <- 
  (1 - sum((datanot_modelNam_pred - filter(to_analyze, continentNAm == 0)$Number.of.pests)^2)/
     sum((mean(filter(to_analyze, continentNAm == 0)$Number.of.pests) - filter(to_analyze, continentNAm == 0)$Number.of.pests)^2) )%>% round(2)
datanot_modelnot_pred_R2 <- 
  (1 - sum((datanot_modelnot_pred - filter(to_analyze, continentNAm == 0)$Number.of.pests)^2)/
     sum((mean(filter(to_analyze, continentNAm == 0)$Number.of.pests) - filter(to_analyze, continentNAm == 0)$Number.of.pests)^2) )%>% round(2)
datanot_modelAll_pred_R2 <- 
  (1 - sum((datanot_modelAll_pred - filter(to_analyze, continentNAm == 0)$Number.of.pests)^2)/
     sum((mean(filter(to_analyze, continentNAm == 0)$Number.of.pests) - filter(to_analyze, continentNAm == 0)$Number.of.pests)^2) )%>% round(2)
dataAll_modelNam_pred_R2 <- 
  (1 - sum((dataAll_modelNam_pred - to_analyze$Number.of.pests)^2)/
     sum((mean(to_analyze$Number.of.pests) - to_analyze$Number.of.pests)^2) )%>% round(2)
dataAll_modelnot_pred_R2 <- 
  (1 - sum((dataAll_modelnot_pred - to_analyze$Number.of.pests)^2)/
     sum((mean(to_analyze$Number.of.pests) - to_analyze$Number.of.pests)^2) )%>% round(2)
dataAll_modelAll_pred_R2 <- 
  (1 - sum((dataAll_modelAll_pred - to_analyze$Number.of.pests)^2)/
     sum((mean(to_analyze$Number.of.pests) - to_analyze$Number.of.pests)^2) )%>% round(2)

dataNAm_modelNAm_pred_MSE <-
  mse(
    y = to_analyze$Number.of.pests[to_analyze$continentNAm==1],
    yhat = predict(pls_NAM, pls_out_1$model$dataX %>% filter(to_analyze$continentNAm == 1), type='response')
  )%>% round(2)
dataNAm_modelnot_pred_MSE <-
  mse(
    y = to_analyze$Number.of.pests[to_analyze$continentNAm==1],
    yhat = predict(pls_not, pls_out_1$model$dataX %>% filter(to_analyze$continentNAm == 1), type='response')
  )%>% round(2)
dataNAm_modelAll_pred_MSE <-
  mse(
    y = to_analyze$Number.of.pests[to_analyze$continentNAm==1],
    yhat = predict(pls_out_1$model, pls_out_1$model$dataX %>% filter(to_analyze$continentNAm == 1), type='response')
  )%>% round(2)
datanot_modelNAm_pred_MSE <-
  mse(
    y = to_analyze$Number.of.pests[to_analyze$continentNAm==0],
    yhat = predict(pls_NAM, pls_out_1$model$dataX %>% filter(to_analyze$continentNAm == 0), type='response')
  )%>% round(2)
datanot_modelnot_pred_MSE <-
  mse(
    y = to_analyze$Number.of.pests[to_analyze$continentNAm==0],
    yhat = predict(pls_not, pls_out_1$model$dataX %>% filter(to_analyze$continentNAm == 0), type='response')
  )%>% round(2)
datanot_modelAll_pred_MSE <-
  mse(
    y = to_analyze$Number.of.pests[to_analyze$continentNAm==0],
    yhat = predict(pls_out_1$model, pls_out_1$model$dataX %>% filter(to_analyze$continentNAm == 0), type='response')
  )%>% round(2)
dataAll_modelNAm_pred_MSE <-
  mse(
    y = to_analyze$Number.of.pests,
    yhat = predict(pls_NAM, pls_out_1$model$dataX, type='response')
  )%>% round(2)
dataAll_modelnot_pred_MSE <-
  mse(
    y = to_analyze$Number.of.pests,
    yhat = predict(pls_not, pls_out_1$model$dataX, type='response')
  )%>% round(2)
dataAll_modelAll_pred_MSE <-
  mse(
    y = to_analyze$Number.of.pests,
    yhat = predict(pls_out_1$model, pls_out_1$model$dataX , type='response')
  )%>% round(2)

validation.points <-
  data.frame(Predicted=dataNAm_modelNAm_pred, Observed = filter(to_analyze, continentNAm == 1)$Number.of.pests, Train.Data = "North America", Test.Data  = "North America") %>%
    rbind(data.frame(Predicted=  dataNAm_modelnot_pred, Observed = filter(to_analyze, continentNAm == 1)$Number.of.pests, Train.Data = "EU, AU & HI", Test.Data  = "North America")) %>%
    rbind(data.frame(Predicted=  dataNAm_modelAll_pred, Observed = filter(to_analyze, continentNAm == 1)$Number.of.pests, Train.Data = "All", Test.Data  = "North America")) %>%
    rbind(data.frame(Predicted=  datanot_modelNam_pred, Observed = filter(to_analyze, continentNAm == 0)$Number.of.pests, Train.Data = "North America", Test.Data  = "EU, AU & HI")) %>%
    rbind(data.frame(Predicted=  datanot_modelnot_pred, Observed = filter(to_analyze, continentNAm == 0)$Number.of.pests, Train.Data = "EU, AU & HI", Test.Data  = "EU, AU & HI")) %>%
    rbind(data.frame(Predicted=  datanot_modelAll_pred, Observed = filter(to_analyze, continentNAm == 0)$Number.of.pests, Train.Data = "All", Test.Data  = "EU, AU & HI")) %>%
    rbind(data.frame(Predicted=  dataAll_modelNam_pred, Observed = to_analyze$Number.of.pests, Train.Data = "North America", Test.Data  = "All")) %>%
    rbind(data.frame(Predicted=  dataAll_modelnot_pred, Observed = to_analyze$Number.of.pests, Train.Data = "EU, AU & HI", Test.Data  = "All")) %>%
    rbind(data.frame(Predicted=  dataAll_modelAll_pred, Observed = to_analyze$Number.of.pests, Train.Data = "All", Test.Data  = "All"))
  
validation.stats <-
  data.frame(
    r=c(dataNAm_modelNAm_pred_pearson,
      dataNAm_modelnot_pred_pearson,
      dataNAm_modelAll_pred_pearson,
      datanot_modelNam_pred_pearson,
      datanot_modelnot_pred_pearson,
      datanot_modelAll_pred_pearson,
      dataAll_modelNam_pred_pearson,
      dataAll_modelnot_pred_pearson,
      dataAll_modelAll_pred_pearson),
    R2=c(dataNAm_modelNAm_pred_R2,
      dataNAm_modelnot_pred_R2,
      dataNAm_modelAll_pred_R2,
      datanot_modelNam_pred_R2,
      datanot_modelnot_pred_R2,
      datanot_modelAll_pred_R2,
      dataAll_modelNam_pred_R2,
      dataAll_modelnot_pred_R2,
      dataAll_modelAll_pred_R2),
    MSE=c(dataNAm_modelNAm_pred_MSE,
      dataNAm_modelnot_pred_MSE,
      dataNAm_modelAll_pred_MSE,
      datanot_modelNAm_pred_MSE,
      datanot_modelnot_pred_MSE,
      datanot_modelAll_pred_MSE,
      dataAll_modelNAm_pred_MSE,
      dataAll_modelnot_pred_MSE,
      dataAll_modelAll_pred_MSE)) %>%
    cbind(distinct(select(validation.points, -c('Observed','Predicted'))))

windows(8,8); ggplot() +
  
  geom_point(data = validation.points, aes(x = Observed, y = Predicted)) +
  
  geom_text(data = validation.stats, aes(x = 1, y = 32, label = paste("r =", r, "\nR² =", R2, "\nMSE =", MSE), hjust=0, vjust=1)) +
  
  geom_abline() +
  
  xlab("\nObserved Pathogen Richness") +
  ylab("Predicted Pathogen Richness\n") +
  
 # scale_x_continuous(transform = 'log', breaks = c(.125,.25,.5,0,1,2,4,8,16,32)) +
 # scale_y_continuous(transform = 'log', breaks = c(0,1,2,4,8,16,32)) +
  
  facet_grid(Test.Data ~ Train.Data) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Training Set\n", breaks=NULL, labels = NULL))+
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Testing Set\n", breaks=NULL, labels = NULL))+
  theme_cowplot()+
  theme(strip.background=element_blank(),
        strip.placement="outside",
        panel.border=element_rect(fill=NA, colour="black")) 

#####################################
