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

load('Final_output/Final_analysis_data_pathogens-transformations-Sept-2024.RData')
load('Final_output/zinf_model_final.RData')
load('Final_output/PLS-models-9-10-23.RData')

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
            mapping=aes(y= exp(bootMed), x=variables, label = round(exp(bootMed),2), col=color.value),
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

save(sjp1,file="Figures_Sept2024/Fig5_effects.RData")

#svg("../Figures/Fig4_9-11-24.svg", 6, 4, family="Cambria Math"); sjp1; dev.off()
#emf("../Figures/Fig4_9-11-24.emf", 6,4, family="Cambria Math", emfPlus=TRUE); sjp1; dev.off()


#############################################################
###############                   ###########################
############### interaction plots ###########################
###############                   ###########################
#############################################################

full_model_fit <- explicit.model.specification(pls_out_1$boot_summary[c(1,3,2,5,4,6),'bootMed'], d = to_analyze)

# multi continent plots

nam_model_fit <- explicit.model.specification(pls_NAM$boot_summary[c(1,3,2,5,4,6),'bootMed'], d = to_analyze)
aus_model_fit <- explicit.model.specification(pls_Au$boot_summary[c(1,3,2,5,4,6),'bootMed'], d = to_analyze)
hii_model_fit <- explicit.model.specification(pls_HII$boot_summary[c(1,3,2,5,4,6),'bootMed'], d = to_analyze)
eur_model_fit <- explicit.model.specification(pls_EU$boot_summary[c(1,3,2,5,4,6),'bootMed'], d = to_analyze)

scale_continents <-
  c(
    Australia = "red",
    Eurasian.Palearctic = "darkgreen",
    Hawaii = "cyan",
    North.America = 'purple'
  )

scale_continents[] <- c('red','black','blue','tan')

# trade layers

points_layer_trade <-
  ggplot() +
    geom_point(
      aes(
        x = root4_trade_centered,
        y = Number.of.pests,
        color = continent,
        pch = continent
      ), data = to_analyze
    ) +
#  scale_shape_manual("Observed", values=c(3,4,8,19), guide= 'legend')+
#  scale_color_manual("Observed", values=scale_continents, guide= 'legend')+
  scale_y_continuous(transform = 'log', breaks = c(1,2,4,8,16,32), limits=c(1,32))+
  scale_x_continuous(
    breaks = (c(0,.1,1.3,6.5,20.6))^.25-attr(to_analyze$root4_trade_centered,"scaled:center"),
    labels = c(0,.1,1.3,6.5,20.6),
    limits = (c(0,20.6))^.25-attr(to_analyze$root4_trade_centered,"scaled:center")
  )+
  theme_cowplot()+theme(strip.background=element_blank(), strip.placement="outside",
                        panel.border=element_rect(fill=NA, colour="black"))
#+theme_cowplot()


range_nam <-
  to_analyze %>% filter(continent == 'North.America') %>%
    with(range(root4_trade_centered)*100) %>% round %>%
    (function(x) seq(x[1], x[2])/100)

trade_predict_nam <- ggpredict(nam_model_fit, terms = list(root4_trade_centered = range_nam))
trade_predict_nam$continent <- as.factor('North.America')
trade_predict_nam$std.error <- pls_NAM$boot_summary['Trade','bootSE'] #sum(pls_NAM$boot_summary$bootSE)
trade_predict_nam$conf.low  <- (trade_predict_nam$predicted - 1.96 * trade_predict_nam$std.error) %>% (function (x) {x[x<1] <- 1;x})
trade_predict_nam$conf.high <- (trade_predict_nam$predicted + 1.96 * trade_predict_nam$std.error) %>% (function (x) {x[x>32] <- 32;x})

range_aus <-
  to_analyze %>% filter(continent == 'Australia') %>%
  with(range(root4_trade_centered)*100) %>% round %>%
  (function(x) seq(x[1], x[2])/100)

trade_predict_aus <- ggpredict(aus_model_fit, terms = list(root4_trade_centered = range_aus))
trade_predict_aus$continent <- as.factor('Australia')
trade_predict_aus$std.error <- pls_Au$boot_summary['Trade','bootSE']#sum(pls_Au$boot_summary$bootSE)
trade_predict_aus$conf.low  <- (trade_predict_aus$predicted - 1.96 * trade_predict_aus$std.error) %>% (function (x) {x[x<1] <- 1;x})
trade_predict_aus$conf.high <- (trade_predict_aus$predicted + 1.96 * trade_predict_aus$std.error) %>% (function (x) {x[x>32] <- 32;x})

range_eu <-
  to_analyze %>% filter(continent == 'Eurasian.Palearctic') %>%
  with(range(root4_trade_centered)*100) %>% round %>%
  (function(x) seq(x[1], x[2])/100)

trade_predict_eur <- ggpredict(eur_model_fit, terms = list(root4_trade_centered = range_eu))
trade_predict_eur$continent <- as.factor('Eurasian.Palearctic')
trade_predict_eur$std.error <- pls_EU$boot_summary['Trade','bootSE']#sum(pls_EU$boot_summary$bootSE)
trade_predict_eur$conf.low  <- (trade_predict_eur$predicted - 1.96 * trade_predict_eur$std.error) %>% (function (x) {x[x<1] <- 1;x})
trade_predict_eur$conf.high <- (trade_predict_eur$predicted + 1.96 * trade_predict_eur$std.error) %>% (function (x) {x[x>32] <- 32;x})

range_hi <-
  to_analyze %>% filter(continent == 'Hawaii') %>%
  with(range(root4_trade_centered)*100) %>% round %>%
  (function(x) seq(x[1], x[2])/100)

trade_predict_hii <- ggpredict(hii_model_fit, terms = list(root4_trade_centered = range_hi))
trade_predict_hii$continent <- as.factor('Hawaii')
trade_predict_hii$std.error <- pls_HII$boot_summary['Trade','bootSE']#sum(pls_HII$boot_summary$bootSE)
trade_predict_hii$conf.low  <- (trade_predict_hii$predicted - 1.96 * trade_predict_hii$std.error) %>% (function (x) {x[x<1] <- 1;x})
trade_predict_hii$conf.high <- (trade_predict_hii$predicted + 1.96 * trade_predict_hii$std.error) %>% (function (x) {x[x>32] <- 32;x})

united_trade_data <-
  rbind(trade_predict_nam,trade_predict_aus, trade_predict_eur, trade_predict_hii)

united_trade_data$group
united_trade_data$continent <- factor(united_trade_data$continent,
                                      levels = levels(to_analyze$continent))
gg_trade <-

  points_layer_trade +
  
  scale_shape_manual("Observed", values=c(3,4,8,19), guide= 'legend')+
  scale_color_manual("Observed", values=scale_continents, guide= 'legend')+
  
  new_scale_color() +
  
    geom_line(aes(
      x = x,
      y = predicted,
      color = continent,
      linetype = continent
    ), data = united_trade_data,
    linewidth=1
    ) +
    geom_ribbon(aes(
      x = x,
      ymin = conf.low,
      ymax = conf.high,
      fill = continent
    ), data = united_trade_data,
    alpha = .1) +
  
#    scale_shape_manual("Observed", values=c(3,4,8,19), guide= 'legend')+
#    scale_color_manual("Observed", values=scale_continents, guide= 'legend')+

    scale_color_manual('Predicted',values=scale_continents, guide= 'legend') +
    scale_fill_manual('Predicted',values=scale_continents, guide= 'legend') +
    scale_linetype_manual('Predicted',values=c(1,2,3,4), guide='legend') +
  
#    scale_shape_manual('Predicted', values=NULL, guide='legend') +

    labs(
      x = '\nImports (Trillion USD) [∜ scale]\n',
      y = '\nEstablished pathogen richness [log scale]\n') #+
 

 
 # guides(
#    fill = 'legend_pred',
#    colour='legend_pred',
#    linetype='legend_pred',
#    shape = 'legend_observed'
#  )
  
  #,
      #fill = "Continent", color = "Observed", linetype= "Continent", pch = "Observed")

# climate layers


points_layer_climate <-
  ggplot() +
  geom_point(
    aes(
      x = Climate.Distance_centered,
      y = Number.of.pests,
      color = continent,
      pch = continent
    ), data = to_analyze
  ) +
#  scale_shape_manual("Continent", values=c(3,4,8,19), guide= 'legend')+
  scale_y_continuous(transform = 'log', breaks = c(1,2,4,8,16,32), limits=c(1,32))+
  scale_x_continuous(
    limits = c(0.1,.65)-attr(to_analyze$Climate.Distance_centered,"scaled:center"),
    breaks = c(.2,.4,.6)-attr(to_analyze$Climate.Distance_centered,"scaled:center"),
    labels = c(.2,.4,.6))+
  theme_cowplot()+theme(strip.background=element_blank(), strip.placement="outside",
                        panel.border=element_rect(fill=NA, colour="black"))

range_nam <-
  to_analyze %>% filter(continent == 'North.America') %>%
  with(range(Climate.Distance_centered)*100) %>% round %>%
  (function(x) seq(x[1], x[2])/100)

climate_predict_nam <- ggpredict(nam_model_fit, terms = list(Climate.Distance_centered = range_nam))
climate_predict_nam$continent <- as.factor('North.America')
climate_predict_nam$std.error <- pls_NAM$boot_summary['Climate','bootSE'] #sum(pls_NAM$boot_summary$bootSE)
climate_predict_nam$conf.low  <- (climate_predict_nam$predicted - 1.96 * climate_predict_nam$std.error) %>% (function (x) {x[x<1] <- 1;x})
climate_predict_nam$conf.high <- (climate_predict_nam$predicted + 1.96 * climate_predict_nam$std.error) %>% (function (x) {x[x>32] <- 32;x})

range_aus <-
  to_analyze %>% filter(continent == 'Australia') %>%
  with(range(Climate.Distance_centered)*100) %>% round %>%
  (function(x) seq(x[1], x[2])/100)

climate_predict_aus <- ggpredict(aus_model_fit, terms = list(Climate.Distance_centered = range_aus))
climate_predict_aus$continent <- as.factor('Australia')
climate_predict_aus$std.error <- pls_Au$boot_summary['Climate','bootSE']#sum(pls_Au$boot_summary$bootSE)
climate_predict_aus$conf.low  <- (climate_predict_aus$predicted - 1.96 * climate_predict_aus$std.error) %>% (function (x) {x[x<1] <- 1;x})
climate_predict_aus$conf.high <- (climate_predict_aus$predicted + 1.96 * climate_predict_aus$std.error) %>% (function (x) {x[x>32] <- 32;x})

range_eu <-
  to_analyze %>% filter(continent == 'Eurasian.Palearctic') %>%
  with(range(Climate.Distance_centered)*100) %>% round %>%
  (function(x) seq(x[1], x[2])/100)

climate_predict_eur <- ggpredict(eur_model_fit, terms = list(Climate.Distance_centered = range_eu))
climate_predict_eur$continent <- as.factor('Eurasian.Palearctic')
climate_predict_eur$std.error <- pls_EU$boot_summary['Climate','bootSE']#sum(pls_EU$boot_summary$bootSE)
climate_predict_eur$conf.low  <- (climate_predict_eur$predicted - 1.96 * climate_predict_eur$std.error) %>% (function (x) {x[x<1] <- 1;x})
climate_predict_eur$conf.high <- (climate_predict_eur$predicted + 1.96 * climate_predict_eur$std.error) %>% (function (x) {x[x>32] <- 32;x})

range_hi <-
  to_analyze %>% filter(continent == 'Hawaii') %>%
  with(range(Climate.Distance_centered)*100) %>% round %>%
  (function(x) seq(x[1], x[2])/100)

climate_predict_hii <- ggpredict(hii_model_fit, terms = list(Climate.Distance_centered = range_hi))
climate_predict_hii$continent <- as.factor('Hawaii')
climate_predict_hii$std.error <- pls_HII$boot_summary['Climate','bootSE']#sum(pls_HII$boot_summary$bootSE)
climate_predict_hii$conf.low  <- (climate_predict_hii$predicted - 1.96 * climate_predict_hii$std.error) %>% (function (x) {x[x<1] <- 1;x})
climate_predict_hii$conf.high <- (climate_predict_hii$predicted + 1.96 * climate_predict_hii$std.error) %>% (function (x) {x[x>32] <- 32;x})

united_climate_data <-
  rbind(climate_predict_nam,climate_predict_aus, climate_predict_eur, climate_predict_hii)

united_climate_data$continent <- factor(united_climate_data$continent,
                                      levels = levels(to_analyze$continent))
gg_climate <-
  
  points_layer_climate +
  scale_shape_manual("Observed", values=c(3,4,8,19), guide= 'legend')+
  scale_color_manual("Observed", values=scale_continents, guide= 'legend')+
  
  new_scale_color() +
  geom_line(aes(
    x = x,
    y = predicted,
    color = continent,
    linetype = continent
  ), data = united_climate_data,
  linewidth=1
  ) +
  geom_ribbon(aes(
    x = x,
    ymin = conf.low,
    ymax = conf.high,
    fill = continent
  ), data = united_climate_data,
  alpha = .1) +
  
  scale_color_manual('Predicted',values=scale_continents, guide= 'legend') +
  scale_fill_manual('Predicted',values=scale_continents, guide= 'legend') +
  scale_linetype_manual('Predicted',values=c(1,2,3,4), guide='legend') +
  labs(
    x = '\nClimatic Distance\n',
    y = '\nEstablished pathogen richness [log scale]\n')#,
  #  fill = "Continent", color = "Continent", linetype= "Continent", pch = "Continent")

# NMDS layers

points_layer_NMDS <-
  ggplot() +
  geom_point(
    aes(
      x = NMDS.Distance_centered,
      y = Number.of.pests,
      color = continent,
      pch = continent
    ), data = to_analyze
  ) +
 # scale_shape_manual("Continent", values=c(3,4,8,19), guide= 'legend')+
  scale_y_continuous(transform = 'log', breaks = c(1,2,4,8,16,32), limits=c(1,32))+
  scale_x_continuous(
    limits = c(0.1,.65)-attr(to_analyze$NMDS.Distance_centered,"scaled:center"),
    breaks = c(.2,.4,.6)-attr(to_analyze$NMDS.Distance_centered,"scaled:center"),
    labels = c(.2,.4,.6))+
  theme_cowplot()+theme(strip.background=element_blank(), strip.placement="outside",
                        panel.border=element_rect(fill=NA, colour="black"))

range_nam <-
  to_analyze %>% filter(continent == 'North.America') %>%
  with(range(NMDS.Distance_centered)*100) %>% round %>%
  (function(x) seq(x[1], x[2])/100)

NMDS_predict_nam <- ggpredict(nam_model_fit, terms = list(NMDS.Distance_centered = range_nam))
NMDS_predict_nam$continent <- as.factor('North.America')
NMDS_predict_nam$std.error <- pls_NAM$boot_summary['NMDS.Dist','bootSE'] #sum(pls_NAM$boot_summary$bootSE)
NMDS_predict_nam$conf.low  <- (NMDS_predict_nam$predicted - 1.96 * NMDS_predict_nam$std.error) %>% (function (x) {x[x<1] <- 1;x})
NMDS_predict_nam$conf.high <- (NMDS_predict_nam$predicted + 1.96 * NMDS_predict_nam$std.error) %>% (function (x) {x[x>32] <- 32;x})

range_aus <-
  to_analyze %>% filter(continent == 'Australia') %>%
  with(range(NMDS.Distance_centered)*100) %>% round %>%
  (function(x) seq(x[1], x[2])/100)

NMDS_predict_aus <- ggpredict(aus_model_fit, terms = list(NMDS.Distance_centered = range_aus))
NMDS_predict_aus$continent <- as.factor('Australia')
NMDS_predict_aus$std.error <- pls_Au$boot_summary['NMDS.Dist','bootSE']#sum(pls_Au$boot_summary$bootSE)
NMDS_predict_aus$conf.low  <- (NMDS_predict_aus$predicted - 1.96 * NMDS_predict_aus$std.error) %>% (function (x) {x[x<1] <- 1;x})
NMDS_predict_aus$conf.high <- (NMDS_predict_aus$predicted + 1.96 * NMDS_predict_aus$std.error) %>% (function (x) {x[x>32] <- 32;x})

range_eu <-
  to_analyze %>% filter(continent == 'Eurasian.Palearctic') %>%
  with(range(NMDS.Distance_centered)*100) %>% round %>%
  (function(x) seq(x[1], x[2])/100)

NMDS_predict_eur <- ggpredict(eur_model_fit, terms = list(NMDS.Distance_centered = range_eu))
NMDS_predict_eur$continent <- as.factor('Eurasian.Palearctic')
NMDS_predict_eur$std.error <- pls_EU$boot_summary['NMDS.Dist','bootSE']#sum(pls_EU$boot_summary$bootSE)
NMDS_predict_eur$conf.low  <- (NMDS_predict_eur$predicted - 1.96 * NMDS_predict_eur$std.error) %>% (function (x) {x[x<1] <- 1;x})
NMDS_predict_eur$conf.high <- (NMDS_predict_eur$predicted + 1.96 * NMDS_predict_eur$std.error) %>% (function (x) {x[x>32] <- 32;x})

range_hi <-
  to_analyze %>% filter(continent == 'Hawaii') %>%
  with(range(NMDS.Distance_centered)*100) %>% round %>%
  (function(x) seq(x[1], x[2])/100)

NMDS_predict_hii <- ggpredict(hii_model_fit, terms = list(NMDS.Distance_centered = range_hi))
NMDS_predict_hii$continent <- as.factor('Hawaii')
NMDS_predict_hii$std.error <- pls_HII$boot_summary['NMDS.Dist','bootSE']#sum(pls_HII$boot_summary$bootSE)
NMDS_predict_hii$conf.low  <- (NMDS_predict_hii$predicted - 1.96 * NMDS_predict_hii$std.error) %>% (function (x) {x[x<1] <- 1;x})
NMDS_predict_hii$conf.high <- (NMDS_predict_hii$predicted + 1.96 * NMDS_predict_hii$std.error) %>% (function (x) {x[x>32] <- 32;x})

united_NMDS_data <-
  rbind(NMDS_predict_nam,NMDS_predict_aus, NMDS_predict_eur, NMDS_predict_hii)

united_NMDS_data$continent <- factor(united_NMDS_data$continent,
                                        levels = levels(to_analyze$continent))
gg_NMDS <-
  
  points_layer_NMDS +
  scale_shape_manual("Observed", values=c(3,4,8,19), guide= 'legend')+
  scale_color_manual("Observed", values=scale_continents, guide= 'legend')+
  
  new_scale_color() +
  
  geom_line(aes(
    x = x,
    y = predicted,
    color = continent,
    linetype = continent
  ), data = united_NMDS_data,
  linewidth=1
  ) +
  geom_ribbon(aes(
    x = x,
    ymin = conf.low,
    ymax = conf.high,
    fill = continent
  ), data = united_NMDS_data,
  alpha = .1) +
  scale_color_manual('Predicted',values=scale_continents, guide= 'legend') +
  scale_fill_manual('Predicted',values=scale_continents, guide= 'legend') +
  scale_linetype_manual('Predicted',values=c(1,2,3,4), guide='legend') +  
  #scale_fill_manual(values=scale_continents, guide= 'legend') +
  #scale_color_manual(values=scale_continents, guide= 'legend') +
  labs(
    x = '\nPhylogeographic Distance\n',
    y = '\nEstablished pathogen richness [log scale]\n')#,
  #  fill = "Continent", color = "Continent", linetype= "Continent", pch = "Continent")

# Trees layers

summary(to_analyze$number_of_exotic_trees)
summary(to_analyze$log_numexotics)
summary(exp(to_analyze$log_numexotics)-1)
summary(exp(to_analyze$log_numexotics_centered +
              attr(to_analyze$log_numexotics_centered, "scaled:center")
              )-1)

points_layer_Trees <-
  ggplot() +
  geom_point(
    aes(
      x = log_numexotics_centered,
      y = Number.of.pests,
      color = continent,
      pch = continent
    ), data = to_analyze) +
#  scale_shape_manual("Continent", values=c(3,4,8,19), guide= 'legend')+
  
  scale_y_continuous(transform = 'log', breaks = c(1,2,4,8,16,32), limits=c(1,32))+
  
  scale_x_continuous(
#    limits = c(log(1), log(256)) - attr(to_analyze$log_numexotics_centered, "scaled:center"),
    breaks = log(c(1,2,4,8,16,32,64,128,256)) - attr(to_analyze$log_numexotics_centered, "scaled:center"),
    labels = function(x) round(exp(x+attr(to_analyze$log_numexotics_centered, "scaled:center")), 0))+
  
  theme_cowplot()+theme(strip.background=element_blank(), strip.placement="outside",
                        panel.border=element_rect(fill=NA, colour="black"))
  
#  scale_x_continuous(
# #   transform = scales::new_transform(
##      'reverse_center_log',
##      transform = function (x) exp(x + 2.663239)-1,
##      inverse = function (x) log1p(x)- 2.663239),
##    breaks = (function (x) {log1p(x)- 2.663239})(c(1,2,4,8,16,32,64)),
#    labels = c(1,2,4,8,16,32,64),
#    limits = (function (x) {log1p(x)- 2.663239})(c(1,64)))+
#  theme_cowplot()

range_nam <-
  to_analyze %>% filter(continent == 'North.America') %>%
  with(range(log_numexotics_centered)*100) %>% round %>%
  (function(x) seq(x[1], x[2])/100)

Trees_predict_nam <- ggpredict(nam_model_fit, terms = list(log_numexotics_centered = range_nam))
Trees_predict_nam$continent <- as.factor('North.America')
Trees_predict_nam$std.error <- pls_NAM$boot_summary['Non-native trees','bootSE'] #sum(pls_NAM$boot_summary$bootSE)
Trees_predict_nam$conf.low  <- (Trees_predict_nam$predicted - 1.96 * Trees_predict_nam$std.error) %>% (function (x) {x[x<1] <- 1;x})
Trees_predict_nam$conf.high <- (Trees_predict_nam$predicted + 1.96 * Trees_predict_nam$std.error) %>% (function (x) {x[x>32] <- 32;x})

range_aus <-
  to_analyze %>% filter(continent == 'Australia') %>%
  with(range(log_numexotics_centered)*100) %>% round %>%
  (function(x) seq(x[1], x[2])/100)

Trees_predict_aus <- ggpredict(aus_model_fit, terms = list(log_numexotics_centered = range_aus))
Trees_predict_aus$continent <- as.factor('Australia')
Trees_predict_aus$std.error <- pls_Au$boot_summary['Non-native trees','bootSE']#sum(pls_Au$boot_summary$bootSE)
Trees_predict_aus$conf.low  <- (Trees_predict_aus$predicted - 1.96 * Trees_predict_aus$std.error) %>% (function (x) {x[x<1] <- 1;x})
Trees_predict_aus$conf.high <- (Trees_predict_aus$predicted + 1.96 * Trees_predict_aus$std.error) %>% (function (x) {x[x>32] <- 32;x})

range_eu <-
  to_analyze %>% filter(continent == 'Eurasian.Palearctic') %>%
  with(range(log_numexotics_centered)*100) %>% round %>%
  (function(x) seq(x[1], x[2])/100)

Trees_predict_eur <- ggpredict(eur_model_fit, terms = list(log_numexotics_centered = range_eu))
Trees_predict_eur$continent <- as.factor('Eurasian.Palearctic')
Trees_predict_eur$std.error <- pls_EU$boot_summary['Non-native trees','bootSE']#sum(pls_EU$boot_summary$bootSE)
Trees_predict_eur$conf.low  <- (Trees_predict_eur$predicted - 1.96 * Trees_predict_eur$std.error) %>% (function (x) {x[x<1] <- 1;x})
Trees_predict_eur$conf.high <- (Trees_predict_eur$predicted + 1.96 * Trees_predict_eur$std.error) %>% (function (x) {x[x>32] <- 32;x})

range_hi <-
  to_analyze %>% filter(continent == 'Hawaii') %>%
  with(range(log_numexotics_centered)*100) %>% round %>%
  (function(x) seq(x[1], x[2])/100)

Trees_predict_hii <- ggpredict(hii_model_fit, terms = list(log_numexotics_centered = range_hi))
Trees_predict_hii$continent <- as.factor('Hawaii')
Trees_predict_hii$std.error <- pls_HII$boot_summary['Non-native trees','bootSE']#sum(pls_HII$boot_summary$bootSE)
Trees_predict_hii$conf.low  <- (Trees_predict_hii$predicted - 1.96 * Trees_predict_hii$std.error) %>% (function (x) {x[x<1] <- 1;x})
Trees_predict_hii$conf.high <- (Trees_predict_hii$predicted + 1.96 * Trees_predict_hii$std.error) %>% (function (x) {x[x>32] <- 32;x})

united_Trees_data <-
  rbind(Trees_predict_nam,Trees_predict_aus, Trees_predict_eur, Trees_predict_hii)

united_Trees_data$continent <- factor(united_Trees_data$continent,
                                     levels = levels(to_analyze$continent))
gg_Trees <-
  
  points_layer_Trees +
  
  scale_shape_manual("Observed", values=c(3,4,8,19), guide= 'legend')+
  scale_color_manual("Observed", values=scale_continents, guide= 'legend')+
  
  new_scale_color() +

  geom_line(aes(
    x = x,
    y = predicted,
    color = continent,
    linetype = continent
  ), data = united_Trees_data,
  linewidth=1) +
  geom_ribbon(aes(
    x = x,
    ymin = conf.low,
    ymax = conf.high,
    fill = continent
  ), data = united_Trees_data,
  alpha = .1) +
  scale_color_manual('Predicted',values=scale_continents, guide= 'legend') +
  scale_fill_manual('Predicted',values=scale_continents, guide= 'legend') +
  scale_linetype_manual('Predicted',values=c(1,2,3,4), guide='legend') +
#  scale_fill_manual(values=scale_continents, guide= 'legend') +
#  scale_color_manual(values=scale_continents, guide= 'legend') +
  labs(
    x = '\nNon-native Tree Richness [log scale]\n',
    y = '\nEstablished pathogen richness [log scale]\n')#,
#    fill = "Continent", color = "Continent", linetype= "Continent", pch = "Continent")

############################################

windows();
gg_trade
gg_climate
gg_NMDS
gg_Trees

composite_layout <-
  c(
    area(t = 0, b = 99, l = 0, r = 99),
    area(t = 0, b = 99, l = 100, r = 199),
    area(t = 100, b = 199, l = 0, r = 99),
    area(t = 100, b = 199, l = 100, r = 199)
  )

model_patch <-
  (gg_trade +
     gg_climate +
     gg_NMDS +
     gg_Trees) +
    plot_layout(design = composite_layout,
                guides='collect',
                axes= 'collect_y'#,
         #       axis_titles = 'collect_y'
                ) +
    plot_annotation(tag_levels = 'A')  &
  
#  scale_linetype_manual(values = c("solid",
#                                   "longdash",
#                                   "dashed",
#                                   "dotted"))&
  
  ylab(NULL) &

      theme(axis.title = element_text(size = 12),
          legend.position = 'bottom',
          plot.margin = margin(10,10,10,10),
          plot.tag.location = 'panel',
          plot.tag.position = c(.025,.95))

model_patch_labelled <-
  wrap_elements(model_patch) +
    labs(tag = "                      Established pathogen richness [log scale]") +
    theme(
      plot.tag = element_text(size = 12, angle = 90),
      plot.tag.position = "left",
      legend.key.size = unit(10,'cm'))

windows(7.3,6.9);model_patch_labelled

save(model_patch_labelled,file="Figures_Sept2024/FigS3.RData")

svg("../Figures/FigS3_9-11-24.svg", 7.3,6.9, family="Cambria Math"); model_patch_labelled; dev.off()
emf("../Figures/FigS3_9-11-24.emf", 7.3,6.9, family="Cambria Math", emfPlus=TRUE); model_patch_labelled; dev.off()

