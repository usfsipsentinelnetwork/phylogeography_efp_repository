library(ggeffects)
library(cowplot)
library(sjPlot)
library(patchwork)
library(glmmTMB)
library(Cairo)
library(devEMF)
library(dplyr)
library(ggplot2)
library(ggnewscale)
library(see)

load('Final_output/Final_analysis_data_pathogens-transformations-March4-2025-pnas.RData')
load('Final_output/zinf_model_final_march4_2025.RData')
load('Final_output/PLS-models-3-11-25.RData')

explicit.model.specification <- function (params, d) {

    glmmTMB(
      Number.of.pests ~
        root4_trade_centered+
        Climate.Distance_centered +
        log_numexotics_centered +
        NMDS.Distance_centered +
        NMDS.Distance_centered:log_numexotics_centered,
      start = list(
        beta = params
      ),
      control = glmmTMBControl(optCtrl = list(iter.max=0)),
      data=d,
      family="poisson"
    )
  
}

#######################
#####             #####
#####   EFFECTS   #####
#####    PLOTS    #####
#####             #####
#######################

newnames <- rownames(pls_out_1$model$Coeffs)[1:5] %>% c("log_numexotics_centered:NMDS.Distance_centered")

boot.params <- pls_out_1$boot$t %>% as.data.frame
colnames(boot.params) <- newnames

head(boot.params %>% tail)

#### EFFECTS PLOT  ####

sjp1<-
  plot_model(model_final_zinf_centered, vline.color="black",
             order.terms=c(5,1,2,4,3),
             show.values=T,
             show.p=F,
             type="std",
             value.size=3,
             #transform=NULL,
             value.offset=-.2,
             dot.size = 4,
             line.size = 2,
             pred.type="fe",
             #axis.lim = c(0.5,10)
             ) +
  sjPlot::theme_sjplot2()+
  geom_violinhalf(inherit.aes=F,
                  data=boot.params %>%
                    dplyr::select(-Intercept) %>%
                    melt %>%
                    left_join(
                      data.frame(
                        variable=newnames %>% setdiff("Intercept"),
                        color.value=as.factor(c("neg","pos","neg","pos","pos"))
                      )) %>%
                    mutate(value = exp(value)),
                  mapping=aes(y=value, x=variable, color=color.value),
                  alpha=0) +
  geom_text(data=pls_out_1$boot_summary[-1,] %>%
              mutate(variables = c(
                "Climate.Distance_centered",
                "root4_trade_centered",
                "NMDS.Distance_centered",
                "log_numexotics_centered",
                "log_numexotics_centered:NMDS.Distance_centered"
              )) %>%
              mutate(color.value=as.factor(c("neg","pos","neg","pos","pos"))),
            inherit.aes=F,
            mapping=aes(y= exp(bootMed), x=variables, label = round(exp(bootMed),2)%>%format(2), col=color.value),
            position = position_nudge(x=.55),
            cex=3)+
  scale_x_discrete(
    labels=c(
      "root4_trade_centered"="Imports (∜)",
      "Climate.Distance_centered"="Climate (Dist.)",
      "NMDS.Distance_centered"="Phylogeography (Dist.)",
      "log_numexotics_centered"="Non-native Trees (log)",
      "log_numexotics_centered:NMDS.Distance_centered"="Phylogeo. \U00D7 Non-native Trees"
    ))+
  scale_y_continuous(limits = c(0.5,10), transform='log', breaks = c(0.5,1,2,5,10)) +
  #ggtitle("Richness of Invasive Pathogens")+
  ylab("\nStandardized effect estimate (Fold \U0394 \U22C5 σ \U207B\U00B9)")+
  theme_cowplot()+
  theme(plot.title = element_blank()) 

#windows(6,4);sjp1

save(sjp1,file="Figures_Sept2024/Fig5_effects_march4.2025.RData")

#svg("../Figures/Fig4_9-11-24.svg", 6, 4, family="Cambria Math"); sjp1; dev.off()
#emf("../Figures/Fig4_9-11-24.emf", 6,4, family="Cambria Math", emfPlus=TRUE); sjp1; dev.off()
