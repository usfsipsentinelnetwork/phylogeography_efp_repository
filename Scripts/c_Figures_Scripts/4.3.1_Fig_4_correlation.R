library(ggplot2)
library(dplyr)
library(doBy)
library(reshape2)
library(ggh4x)
library(scales)
library(devEMF)
library(cowplot)

load("Final_output/Final_analysis_data_pathogens-Sept-2024.RData")

# make data frame
to_analyze <- to_analyze %>%
  mutate(log_spatauto= log1p(spatauto)) %>%
  mutate(log_trade = log1p(Cumulative.Import.Value.T)) %>%
  mutate(root4_trade = sqrt(sqrt(Cumulative.Import.Value.T)))%>%
#  mutate(log_numexotics = ifelse(number_of_exotic_trees>0,log(number_of_exotic_trees),log(number_of_exotic_trees+0.001))) %>%
#  mutate(log_numexotics = ifelse(number_of_exotic_trees>0,log(number_of_exotic_trees),-1)) %>%
  mutate(log_numexotics = log1p(number_of_exotic_trees)) %>%
  mutate(continent = gsub("([A-Z]+\\.)(North\\.America)","\\2",Invaded.Range.Cluster, perl=T))

to_analyze$continent <- as.factor(to_analyze$continent)
to_analyze$continentNAm <- as.numeric(to_analyze$continent=="North.America")
to_analyze$continentHII <- as.numeric(to_analyze$continent=="Hawaii")
to_analyze$continentAus <- as.numeric(to_analyze$continent=="Australia")
to_analyze$continentEur <- as.numeric(to_analyze$continent=="Eurasian.Palearctic")

#to_analyze$logNumber.of.pests <- ifelse(to_analyze$Number.of.pests>0,log(to_analyze$Number.of.pests),log(to_analyze$Number.of.pests+0.001))
#to_analyze$logNumber.of.pests <- ifelse(to_analyze$Number.of.pests>0,log(to_analyze$Number.of.pests),-1)

to_analyze$logNumber.of.pests <- log1p(to_analyze$Number.of.pests)

to_analyze$logNumber.of.pests_centered <- scale(to_analyze$logNumber.of.pests, T, F)
to_analyze$log_numexotics_centered <- scale(to_analyze$log_numexotics, T, F)
to_analyze$NMDS.Distance_centered <- scale(to_analyze$NMDS.Distance, T, F)
to_analyze$Climate.Distance_centered <- scale(to_analyze$Climate.Distance, T, F)
to_analyze$root4_trade_centered <- scale(to_analyze$root4_trade, T, F)

to_analyze_zeroadjusted <- to_analyze %>% filter(Number.of.pests > 0)

#save(to_analyze,to_analyze_zeroadjusted, file="Final_output/Final_analysis_data_pathogens-transformations-Sept-2024.RData")

####### FIRST MAKE A PLOT ###########################


###############################
###############################
##                          ###
##   correlations facet     ###
##                          ###
###############################
###############################

melted <-
  to_analyze %>% 
  dplyr::select(Number.of.pests, Climate.Distance, root4_trade,
                NMDS.Distance, log_numexotics, Invaded.Range.Cluster,
                Geographic.Origin.Cluster) %>%
  melt(id.vars= c("Invaded.Range.Cluster", "Geographic.Origin.Cluster", "Number.of.pests"))

c15 <- #pals::polychrome(29)
  c(
    #"dodgerblue2"
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    #"black",
    "gold1",
    #"skyblue2",
    "#FB9A99",
    # take out hawaii for this plot "blue1", # lt pink
    "palegreen2",
    #"#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    #"maroon", "orchid1", "deeppink1", #"steelblue4",
    "darkturquoise",# "green1", "yellow4", #"yellow3",
    "darkorange4", #"brown"#,
    #"#822E1C",
    "#BDCDFF","#FA0087","#2ED9FF","#E31A1C" # red
  )

reorder_levels <-
  c(9,11,12,          # Asia and Pacific,
    2,5,6,7,14,       # Europe, Middle East, Africa
    1,3,4,8,10,13,15) # Americas

reorder_levels2<-c(3,1,2,4:7)

melted$Geographic.Origin.Cluster<-
  factor(melted$Geographic.Origin.Cluster)
levels(melted$Geographic.Origin.Cluster)<-
  gsub("North\\.America", "CONUS/Canada", levels(melted$Geographic.Origin.Cluster))%>%
  gsub("\\.", " ", .)

melted$Invaded.Range.Cluster<-
  factor(melted$Invaded.Range.Cluster)
levels(melted$Invaded.Range.Cluster)<-
  gsub("North\\.America", "CONUS/Canada", levels(melted$Invaded.Range.Cluster))%>%
  gsub("\\.", " ", .)

to_analyze$Geographic.Origin.Cluster<-
  factor(to_analyze$Geographic.Origin.Cluster)
levels(to_analyze$Geographic.Origin.Cluster)<-
  gsub("North\\.America", "CONUS/Canada", levels(to_analyze$Geographic.Origin.Cluster))%>%
  gsub("\\.", " ", .)

to_analyze$Invaded.Range.Cluster<-
  factor(to_analyze$Invaded.Range.Cluster)
levels(to_analyze$Invaded.Range.Cluster)<-
  gsub("North\\.America", "CONUS/Canada", levels(to_analyze$Invaded.Range.Cluster))%>%
  gsub("\\.", " ", .)

labelmaker <- function (x) {
  paste(
    "coef =",
    round(x$coefficients[2,1],1),
    c("***",
      "**",
      "*",
      ".",
      "ns"
    )[findInterval(x$coefficients[2,4],
                   c(0,
                     0.001,
                     0.01,
                     0.05,
                     0.1,
                     1
                   ))],
    sep = ' '
  )
}


labelmaker2 <- function (x) {
  
#  bquote(
#    italic("r") ~
#    " = " ~
#    .(round(x$estimate,2)) ~
#    .(c("***",
#      "**",
#      "*",
#      ".",
#      "ns"
#    )[findInterval(x$p.value,
#                   c(0,
#                     0.001,
#                     0.01,
#                     0.05,
#                     0.1,
#                     1
#                   ))]))
  
  stars<-paste(
    round(x$estimate,2),
    c("***",
      "**",
      "*",
      ".",
      "ns"
    )[findInterval(x$p.value,
                   c(0,
                     0.001,
                     0.01,
                     0.05,
                     0.1,
                     1
                   ))])
#  paste(
#    expression(
      paste(
#  #    "Pearson r = ",
      "= ",
#      
#      italic("r"), " = ",
#      
#      
#  #    "\uD835\uDC5F = ",
      stars,
      sep = ''
    )
#  stars, sep='')
}

glm_coefs<-
  c(Climate.Distance=glm(Number.of.pests ~ Climate.Distance, data=to_analyze, family='quasipoisson') %>% summary %>% labelmaker,
    root4_trade=glm(Number.of.pests ~ root4_trade, data=to_analyze, family='quasipoisson') %>% summary%>% labelmaker,
    NMDS.Distance=glm(Number.of.pests ~ NMDS.Distance, data=to_analyze, family='quasipoisson') %>% summary%>% labelmaker,
    log_numexotics=glm(Number.of.pests ~ log_numexotics, data=to_analyze, family='quasipoisson') %>% summary%>% labelmaker)

glm_coefs_data <-
  data.frame(variable = as.factor(names(glm_coefs)), text=glm_coefs, y = 30)

glm_coefs_data$x <-
  c(0,0,0,0
    #with(to_analyze, min(Climate.Distance)),
    #with(to_analyze, min(root4_trade)),
    #with(to_analyze, min(NMDS.Distance)),
    #with(to_analyze, min(log_numexotics))
  )

pearson_cors<-
  c(with(to_analyze, cor.test(root4_trade, Climate.Distance)) %>% labelmaker2,
  with(to_analyze, cor.test(NMDS.Distance, Climate.Distance)) %>% labelmaker2,
  with(to_analyze, cor.test(log_numexotics, Climate.Distance)) %>% labelmaker2,
  with(to_analyze, cor.test(NMDS.Distance, root4_trade)) %>% labelmaker2,
  with(to_analyze, cor.test(log_numexotics, root4_trade)) %>% labelmaker2,
  with(to_analyze, cor.test(log_numexotics, NMDS.Distance)) %>% labelmaker2)

levels(melted$variable)
levels(glm_coefs_data$variable)

fig3_correlations_plot<-
  ggplot()+
  geom_point(data=melted, aes(x=value,y=Number.of.pests, color=Geographic.Origin.Cluster, pch=Invaded.Range.Cluster), size=3) +
  scale_color_manual("Source Region", values=c15[reorder_levels], breaks=levels(as.factor(to_analyze$Geographic.Origin.Cluster))[reorder_levels]) +
  scale_shape_manual("Sink Region", values=c(3,4,8,15,17,18,19)[reorder_levels2], breaks=levels(as.factor(to_analyze$Invaded.Range.Cluster))[reorder_levels2])+
  geom_smooth(data=melted, aes(x=value,y=Number.of.pests), method = "glm", color="grey", fill="lightgrey",
              method.args=list(family="quasipoisson"))+
  scale_y_continuous(trans="log",labels=scales::number_format(accuracy=1), breaks=c(1,2,4,8,16,32))+
  
  facet_wrap(~ variable, ncol=2, nrow=2, scales="free_x",
             labeller=as_labeller(c(root4_trade="Imports (Trillion USD) [∜ scale]",#"Imports (∜ Trillion USD)",
                                    Climate.Distance="Climatic Distance",
                                    NMDS.Distance="Phylogeographic Distance",
                                    log_numexotics="Non-native Tree Richness [log scale]")),
             strip.position="bottom")+
  labs(x= NULL, y = "Established pathogen richness [log scale]\n")+
  theme_cowplot()+theme(strip.background=element_blank(), strip.placement="outside",
                        panel.border=element_rect(fill=NA, colour="black"))

fig3_correlations_plot<-
  fig3_correlations_plot +
  facetted_pos_scales(
    x = list(
      scale_x_continuous(limits = c(0,.65), breaks = c(0,.2,.4,.6)),
      scale_x_continuous(breaks = function(x) seq(x[1],x[2], length.out=5), labels = function(x) round(x^4,1)),
      scale_x_continuous(limits = c(0,.65), breaks = c(0,.2,.4,.6)),
      scale_x_continuous(limits = c(log(1), log(256)), breaks = log(c(1,2,4,8,16,32,64,128,256)), labels = function(x) round(exp(x), 0))
    ))+
  geom_text(data = glm_coefs_data,
            mapping = aes(x = x, y = y, label = text), vjust='middle', hjust='left', fontface="bold")

windows();fig3_correlations_plot

save(fig3_correlations_plot, file="Figures_Sept2024/fig4a_correlations_plot.RData")

svg("Final_figures/Fig4_A.svg", 9.5,6.5,family="Cambria Math")
fig3_correlations_plot
dev.off()

corrplot_panel <- function(ggp, left=NULL, bottom=NULL, cvals=c15) {
  gg_to_return <-
    ggp + geom_point() +
    scale_color_manual("Geographic Origin Cluster", values=cvals) +
    scale_shape_manual("Invaded Range Cluster", values=c(3,4,8,15,17,18,19))+
    geom_smooth(aes(pch=NA, color=NA), method="lm", color="grey", fill=NA)+
    theme_cowplot()+
    theme(legend.position="none",
          strip.background=element_blank(), strip.placement="outside",
          panel.border=element_rect(fill=NA, colour="black"),
          axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
          axis.text.x =element_blank(), axis.text.y =element_blank(),
          axis.title= element_blank())+
    labs(x= NULL, y =NULL)
  if (!is.null(left)) gg_to_return <- gg_to_return  + theme(axis.ticks.y = element_line(), axis.text.y=element_text(size=10), axis.title.y=element_text(size=10))+ ylab(left)
  if (!is.null(bottom)) gg_to_return <- gg_to_return + theme(axis.ticks.x = element_line(), axis.text.x=element_text(size=10), axis.title.x=element_text(size=10))+ xlab(bottom) 
  gg_to_return
}

scale_climate_x <- function() scale_x_continuous(limits = c(0,.65), breaks = c(0,.2,.4,.6))
scale_trade_x <- function() scale_x_continuous(#limits = c(0,1),
                                               breaks = function(x) seq(x[1],x[2], length.out=5), labels = function(x) round(x^4,1))
scale_nmds_x <- function() scale_x_continuous(limits = c(0,.65), breaks = c(0,.2,.4,.6))

scale_trade_y <- function() scale_y_continuous(#limits = c(0,1),
                                               breaks = function(x) seq(x[1],x[2], length.out=5), labels = function(x) round(x^4,1))
scale_nmds_y <- function() scale_y_continuous(limits = c(0,.65), breaks = c(0,.2,.4,.6))
scale_trees_y <- function() scale_y_continuous(limits = c(log(1), log(256)), breaks = log(c(1,2,4,8,16,32,64,128,256)), labels = function(x) round(exp(x), 0))

leg <- ggpubr::get_legend(
  fig3_correlations_plot+
    scale_color_manual("\n\n\n\n\n\n\n\nSource Region", values=c15[reorder_levels], breaks=levels(as.factor(to_analyze$Geographic.Origin.Cluster))[reorder_levels])+
    theme(legend.text = element_text(size=9), legend.title = element_text(size=10, face="bold"), legend.key.size = unit(.4, 'cm'))
)

annotate_custom <- function (gp,
                             pc, 
                             #scale_x, scale_y
                             position = 'bottomleft') {
  scale_x <- ggplot_build(gp)$layout$panel_params[[1]]$x.range
  scale_y <- ggplot_build(gp)$layout$panel_params[[1]]$y.range
  
  #print(scale_x)
  #print(scale_y)
  
  x <-
    scale_x[1] + .05*(scale_x[2]- scale_x[1])
  y <-
    ifelse(position == 'bottomleft',
           scale_y[1] + .125*(scale_y[2]- scale_y[1]),
           scale_y[2] - .125*(scale_y[2]- scale_y[1]))
  
  print(x)
  print(y)
  
  l <- paste("paste(italic(r), \"", pc, "\", sep=\"\")")
  
  gp + annotate("text", x = x, y = y,
                #label = pc,
                
                #label = paste(italic('r'), ' = '),
                
                #label = "paste(italic(r)^{2}, \" = \")" %>% paste(pc), parse = T,
                label = l, parse = T,
                
                hjust = 'left', vjust = 'middle', fontface = 'bold')
}

fig3b_correlations <-
  egg::ggarrange(
    
    (ggplot(to_analyze,aes(x=Climate.Distance, y=root4_trade, color=Geographic.Origin.Cluster, pch=Invaded.Range.Cluster)) %>%
      corrplot_panel(left="Imports\n(Trillion USD)\n[∜ scale]")+
      scale_climate_x() + scale_trade_y() +theme(axis.title.y=element_text(face="plain"))) %>%
      annotate_custom(pearson_cors[1]),
    
    (ggplot(to_analyze,aes(x=Climate.Distance, y=NMDS.Distance, color=Geographic.Origin.Cluster, pch=Invaded.Range.Cluster)) %>% corrplot_panel(left="Phylogeographic\nDistance\n") + scale_climate_x()+ scale_nmds_y())%>%
      annotate_custom(pearson_cors[2]),
    
    (ggplot(to_analyze,aes(x=Climate.Distance, y=log_numexotics, color=Geographic.Origin.Cluster, pch=Invaded.Range.Cluster)) %>% corrplot_panel(left="Non-native\nTree Richness\n[log scale]", bottom="\nClimatic\nDistance") + scale_climate_x() + scale_trees_y())%>%
      annotate_custom(pearson_cors[3]),
    
    ggpubr::as_ggplot(leg$grobs[[2]]),
    #ggplot() + theme_nothing(),
    (ggplot(to_analyze,aes(x=root4_trade, y=NMDS.Distance, color=Geographic.Origin.Cluster, pch=Invaded.Range.Cluster)) %>% corrplot_panel() + scale_trade_x() + scale_nmds_y())%>%
      annotate_custom(pearson_cors[4]),
    
    (ggplot(to_analyze,aes(x=root4_trade, y=log_numexotics, color=Geographic.Origin.Cluster, pch=Invaded.Range.Cluster)) %>% corrplot_panel(bottom="\nImports\n(Trillion USD)\n[∜ scale]") + scale_trade_x() + scale_trees_y())%>%
      annotate_custom(pearson_cors[5]),
    
    ggpubr::as_ggplot(leg$grobs[[1]]), ggplot()+ theme_nothing(),
    #ggplot()+ theme_nothing(), ggplot()+ theme_nothing(),
    (ggplot(to_analyze,aes(x=NMDS.Distance, y=log_numexotics, color=Geographic.Origin.Cluster, pch=Invaded.Range.Cluster)) %>% corrplot_panel(bottom="\nPhylogeographic\nDistance") + scale_nmds_x() + scale_trees_y())%>%
      annotate_custom(pearson_cors[6]),
    ncol=3, nrow=3,
    byrow=F
  )

save(fig3b_correlations, file="Figures_Sept2024/fig4b_colinearity_plot.RData")

svg("Final_figures/Fig4_B.svg", 9.5,6.5,family="Cambria Math")
fig3b_correlations
dev.off()

fig3_composite<-
  
  ggpubr::ggarrange(
    ggplot()+theme_nothing(),
    ggpubr::ggarrange(
      ggplot()+theme_nothing(),
      fig3_correlations_plot + theme(legend.position="none"),
      nrow=1,
      ncol=3,
      widths=c(.01,.98,.01)
    ),
    ggplot()+theme_nothing(),
    ggplot()+theme_nothing(),
    ggplot()+theme_nothing(),
    ggpubr::as_ggplot(fig3b_correlations),
    nrow=3,
    ncol=2,
    heights=c(.49,.02,.49),
    widths=c(.05,.95),
    labels=c("A","",
             "","",
             "B","")
  )

windows(7.5,10); fig3_composite

save(fig3_composite, file="Figures_Sept2024/fig4_composite.RData")

svg("Final_figures/Fig4_AB.svg", 7,10,family="Cambria Math")
fig3_composite
dev.off()

emf(file = "Final_figures/Fig4_AB.emf", 7,10, emfPlus=T)
fig3_composite
dev.off()
