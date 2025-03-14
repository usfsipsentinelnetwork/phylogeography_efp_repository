library(sjPlot)
library(tidyverse)
library(lavaan)
library(dplyr)
library(semTools)
library(semPlot)
library(semptools)
library(tidySEM)
library(Cairo)
library(devEMF)
library(patchwork)
library(glmmTMB)
library(see)
library(gridBase)
library(gridGraphics)
library(grid)
library(ggplotify)
library(gt)

load("Figures_Sept2024/Fig5_effects_march4.2025.RData")
load("Figures_Sept2024/Fig_semplots_march2025.RData")
load("Final_output/PLS-models-3-11-25.RData")
load("Final_output/linearmodels_march4.2025.RData")

summary(model_final_zinf_centered)
car::Anova(model_final_zinf_centered, type='III')

#sjp1 / (plot(semp1) + plot(semp2))

dats <- data.frame(
  ziGLMM = c('0.43 ***', '0.09 **', '0.09 ***', '0.01 ns',  '0.00 ns'), # march 4 2025
  PLSGLM = c('0.12 ***', '0.02 *',  '0.05 **', '0.17 ***', '0.01 .'),  # march 4 2025
#  ziGLMM = c('0.37 ***', '0.09 ***', '0.10 ***', '0.02 ns',  '0.00 ns'), # march 1 2025
#  PLSGLM = c('0.13 ***', '0.03 **',  '0.01 ***', '0.23 **', '0.00 ns'),  # march 1 2025
#  ziGLMM = c('0.46 ***', '0.10 ***', '0.07 .', '0.02 ns',  '0.02 **'),
#  PLSGLM = c('0.13 ***', '0.02 **',  '0.03 .', '0.18 ***', '0.00 *'),
  row.names = c(
    "Cum. Imports (∜)",
    "Climatic Dist.",
    "Phylogeographic Dist.",
    "Estab. Trees & Shrubs",
    "Phylogeo. \U00D7 Estab. Trees & Shrubs"
  )
)

sjp1 <-
  sjp1 +
  scale_x_discrete(
    limits = rev(c(
      "root4_trade_centered",
      "Climate.Distance_centered",
      "NMDS.Distance_centered",
      "log_numexotics_centered",
      "log_numexotics_centered:NMDS.Distance_centered"
    )), labels=c(
      "root4_trade_centered"="Cumulative Imports (∜)",
      "Climate.Distance_centered"="Climatic Dist.",
      "NMDS.Distance_centered"="Phylogeographic Dist.",
      "log_numexotics_centered"="Estab. Trees & Shrubs",
      "log_numexotics_centered:NMDS.Distance_centered"="Phylogeo. \U00D7 Trees & Shrubs"
    ))+coord_flip(clip = 'off')

sjp1$layers[[4]]$aes_params$size <- 2.5
sjp1$layers[[6]]$aes_params$size <- 2.5

sjp_colors <- unique(ggplot_build(sjp1)$data[[2]]$colour)

r2_tabs <-
  dats |>
  gt()|>
  tab_header(title = 
               '\uD835\uDC5F \u00b2'
        #       md(
        #       '*r* \u00B2')
             ) |>
  cols_align("left") |>
  tab_style(
    style = list(
      cell_text(color = sjp_colors[2])
    ),
    locations = cells_body(
      columns = 1:2,
      rows = 2:3
    )) |>
  tab_style(
    style = list(
      cell_text(color = sjp_colors[1])
    ),
    locations = cells_body(
      columns = 1:2,
      rows = c(1,4:5)
    )) |>
  tab_options(
    table.font.size = 12,
    data_row.padding = 18,
    table_body.hlines.style = 'none',
    table.border.top.style = 'none',
    table.border.bottom.style = 'none',
    table_body.border.bottom.style = 'none',
    #table.border.bottom.style = 'none',
    heading.border.bottom.style = 'none',
    heading.border.bottom.color = 'black'
    #heading.border.top.style = 'none'
  )

#windows();sjp1 +
#  r2_tabs +
#  plot_layout(ncol = 2, nrow = 2, heights = c(1,1), widths = c(2,1))

sjp1_table <-
  sjp1 +
#  ylab("Standardized effect (fold \U0394 \U22C5 σ \U207B\U00B9)")+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size =12)) +
    wrap_table(
      (r2_tabs ),
      panel = "body",
      space = "free_x") +
    plot_layout(ncol = 2, nrow = 1, heights = c(NA), widths = c(.6,.4))+
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = 'bold', size =12))

#windows();sjp1_table

semp1$graphAttributes$Nodes$labels$log_numexotics <- "Estab.\nTrees\n& Shrubs"
semp1$graphAttributes$Nodes$labels$root4_trade <- "Cum.\nImports"
semp1$graphAttributes$Nodes$labels$logNumber.of.pests <- "EFP\nRichness"
semp1$graphAttributes$Nodes$labels$potential_richness <- "Potential\nEFP\nRichness"

semp2$graphAttributes$Nodes$labels$log_numexotics <- "Estab.\nTrees\n& Shrubs"
semp2$graphAttributes$Nodes$labels$root4_trade <- "Cum.\nImports"
semp2$graphAttributes$Nodes$labels$logNumber.of.pests <- "EFP\nRichness"

windows(8,8)
par(oma = c(0.5,0.5,0,0.5), xpd=T)
layout(matrix(c(1,1,2,3),byrow=T,nrow=2),
       heights = c(.6,.4))
plot.new();vps <- baseViewports(); pushViewport(vps$figure); vp1 <- plotViewport(c(1,1,0,5))
print(
  #sjp1+coord_flip(clip='off'),
  
  sjp1_table,
  
  vp = vp1)
#title(main = "\nA", outer=T, adj=0)
plot(semp1); text(1.4,-.18,
   paste(c('\u03c7\u00b2(mod)','p(\u03c7\u00b2mod)',
           '\u03c7\u00b2(mod|base)','p(\u03c7\u00b2mod|base)',
           'CFI','TLI','RMSEA','SRMR'),"=",
         lavaan::fitMeasures(fit1)[c('chisq','pvalue',
                             'baseline.chisq','baseline.pvalue',
                             'cfi','tli','rmsea','srmr')]%>%round(2)%>%format(2),
         collapse='\n') %>%
     paste("R\u00b2(patho) = ", lavaan::summary(fit1,rsquare=T)$pe$est[22]%>%round(2)%>%format(2),'\n',.,sep=''),
   cex=.8, adj=c(1,1), family='mono')
title(main = "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n       C", outer=T, adj=0)
plot(semp2); text(1.4,-.18,
   paste(c('\u03c7\u00b2(mod)','p(\u03c7\u00b2mod)',
           '\u03c7\u00b2(mod|base)','p(\u03c7\u00b2mod|base)',
           'CFI','TLI','RMSEA','SRMR'),"=",
         lavaan::fitMeasures(fit1.1)[c('chisq','pvalue',
                             'baseline.chisq','baseline.pvalue',
                             'cfi','tli','rmsea','srmr')]%>%round(2)%>%format(2),
         collapse='\n') %>%
     paste("R\u00b2(patho) = ", lavaan::summary(fit1.1,rsquare=T)$pe$est[22]%>%round(2)%>%format(2),'\n',.,sep=''),
   cex=.8, adj=c(1,1), family='mono')
title(main = "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n                                                                                     D", outer=T, adj=0)
###############

library(systemfonts)
svglite::svglite("Final_figures/Fig_5_ABC_march2025.svg", 8,8, system_fonts = list(sans = "Cambria Math",
                                                                                   mono = "Consolas"))  
par(oma = c(0.5,0.5,0,0.5), xpd=T)
layout(matrix(c(1,1,2,3),byrow=T,nrow=2),
       heights = c(.6,.4))
plot.new();vps <- baseViewports(); pushViewport(vps$figure); vp1 <- plotViewport(c(1,1,0,5))
print(
  #sjp1+coord_flip(clip='off'),
  
  sjp1_table,
  
  vp = vp1)
#title(main = "\nA", outer=T, adj=0)
plot(semp1); text(1.4,-.18,
   paste(c('\u03c7\u00b2(mod)','p(\u03c7\u00b2mod)',
           '\u03c7\u00b2(mod|base)','p(\u03c7\u00b2mod|base)',
           'CFI','TLI','RMSEA','SRMR'),"=",
         lavaan::fitMeasures(fit1)[c('chisq','pvalue',
                             'baseline.chisq','baseline.pvalue',
                             'cfi','tli','rmsea','srmr')]%>%round(2)%>%format(2),
         collapse='\n') %>%
     paste("R\u00b2(patho) = ", lavaan::summary(fit1,rsquare=T)$pe$est[22]%>%round(2)%>%format(2),'\n',.,sep=''),
   cex=.8, adj=c(1,1), family='mono')

title(main = "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n        C", outer=T, adj=0)
plot(semp2); text(1.4,-.18,
   paste(c('\u03c7\u00b2(mod)','p(\u03c7\u00b2mod)',
           '\u03c7\u00b2(mod|base)','p(\u03c7\u00b2mod|base)',
           'CFI','TLI','RMSEA','SRMR'),"=",
         lavaan::fitMeasures(fit1.1)[c('chisq','pvalue',
                             'baseline.chisq','baseline.pvalue',
                             'cfi','tli','rmsea','srmr')]%>%round(2)%>%format(2),
         collapse='\n') %>%
     paste("R\u00b2(patho) = ", lavaan::summary(fit1.1,rsquare=T)$pe$est[22]%>%round(2)%>%format(2),'\n',.,sep=''),
   cex=.8, adj=c(1,1), family='mono')
title(main = "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n                                                                          D", outer=T, adj=0)

dev.off()

