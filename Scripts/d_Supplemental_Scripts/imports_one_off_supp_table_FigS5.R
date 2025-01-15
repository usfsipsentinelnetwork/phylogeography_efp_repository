library(tidyverse)
library(cowplot)
library(ggthemes)
library(patchwork)

trade <- read.csv('Covariates_input/trade_aggregated_PPP_wide.csv') %>%
  
  na.omit %>%
  
  filter(biogeographic_region != Reporter.Name) %>%
  filter(!(grepl('Eur',biogeographic_region) &
             grepl('Eur',Reporter.Name))) %>%
  
  dplyr::select(-cumulativeT) %>%
  
  reshape2::melt(variable.name="Year", value.name='USD') %>%
  
  mutate(Year = as.numeric(gsub("^X","",as.character(Year)))) %>% 
  mutate(USD = USD/1000000000) %>% 
  
  arrange(biogeographic_region, Reporter.Name, Year) %>%
  
  group_by(biogeographic_region, Reporter.Name) %>%
  
  mutate(cum_USD = cumsum(USD)) %>% ungroup()
  
g <-

  ggplot() +
#    geom_line(
  geom_area(
    data = trade,
    mapping = 
        aes(
          x = Year,
          y = cum_USD,
          fill = biogeographic_region,
          color=biogeographic_region
        ), linewidth=1.25) +
  scale_y_sqrt(labels = scales::label_currency(prefix = "$", suffix="T"))+#, accuracy=1)) +
  
  scale_color_discrete('biogeographic_region',
    type=c(
      economist_pal()(9),
      fivethirtyeight_pal()(3)
    )
  )+
  scale_fill_discrete('biogeographic_region',
                       type=c(
                         economist_pal()(9),
                         fivethirtyeight_pal()(3)
                       )
  )+
  #scale_color_economist() +
  facet_wrap('Reporter.Name', ncol=1, scales = 'free_y') +
  theme_economist_white() +
  theme(legend.position = 'right', legend.key.size = unit(.5, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 12),
        panel.grid.major.y = element_line(linewidth=.5, linetype = 2),
        legend.key.spacing.y = unit(.5, 'cm')) +
  ylab('Cumulative Imports\n') +
  xlab('\nYear')
  # facet_grid(biogeographic_region~Reporter.Name)

g2 <-
  
  ggplot() +
  #    geom_line(
  geom_area(
    data = trade %>% group_by(Reporter.Name, Year) %>% mutate(proportion_USD = cum_USD/sum(cum_USD)) %>% ungroup,
    mapping = 
      aes(
        x = Year,
        y = proportion_USD,
        fill = biogeographic_region,
        color=biogeographic_region
      ), linewidth=1.25) +
  #scale_y_sqrt(labels = scales::label_currency(prefix = "$", suffix="T"))+#, accuracy=1)) +
  
  scale_color_discrete('biogeographic_region',
                       type=c(
                         economist_pal()(9),
                         fivethirtyeight_pal()(3)
                       )
  )+
  scale_fill_discrete('biogeographic_region',
                      type=c(
                        economist_pal()(9),
                        fivethirtyeight_pal()(3)
                      )
  )+
  #scale_color_economist() +
  facet_wrap('Reporter.Name', ncol=1, scales = 'free_y') +
  theme_economist_white() +
  theme(legend.position = 'right', legend.key.size = unit(.5, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 12),
        panel.grid.major.y = element_line(linewidth=.5, linetype = 2),
        legend.key.spacing.y = unit(.5, 'cm')) +
  ylab('Share of Cumulative Imports\n') +
  xlab('\nYear')
# facet_grid(biogeographic_region~Reporter.Name)

g + g2 +
plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A') 

windows();g
windows();g2

