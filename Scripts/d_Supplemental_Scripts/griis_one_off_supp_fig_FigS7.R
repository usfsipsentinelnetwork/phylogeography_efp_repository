library(tidyverse)
library(cowplot)
library(ggthemes)

griis <-
  
  read.csv("GRIIS - Country Compendium V1_0.csv") %>%
  
  select(kingdom, species) %>%
  
  distinct %>%
  
  group_by(kingdom) %>%
  
  summarise(total_sp = n()) %>%
  
  mutate(percent_sp = 100 * total_sp/sum(total_sp)) %>%
  
  arrange(-percent_sp) %>%
  
  mutate(kingdom = factor(kingdom, levels = c(
    
    "Plantae",
    "Animalia",
    "Chromista",
    "Fungi",
    "Bacteria",
    "Viruses",
    "Protozoa"
    
  )))

ggplot(data = griis, aes(x = "", y = total_sp, fill = kingdom)) +
  geom_bar(stat = "identity", width=1) +
  coord_polar("y") +
  theme_void() +
  scale_fill_brewer(palette = "Set1")

ggplot(data = griis, aes(x = kingdom, y = percent_sp, fill = kingdom)) +
  geom_bar(stat = "identity", width=1) +
  
  geom_text(aes(x = kingdom, y = percent_sp, label = paste(round(percent_sp,2), "%\n\n\n"))) +
  
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(.05,.1,.5,1,5,10,50,100),
                     limits = c(0,75))+
  theme_cowplot() +
  theme(
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.position = "none"#,
    #panel.grid = element_blank(),
    #panel.background = element_blank()
  )+
  scale_fill_brewer(palette = "Set1")

########


ggplot(data = griis, aes(x = kingdom, y = total_sp, fill = kingdom)) +
  geom_bar(stat = "identity", width=1) +
  
  geom_text(aes(x = kingdom, y = total_sp, label = paste(round(percent_sp,2), "%\n\n\n"))) +
  
  scale_y_continuous(trans = 'log',
                     breaks = 2^(0:14),
                     limits = 2^c(0,16)) +
  
  #scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
  #                   breaks = c(.05,.1,.5,1,5,10,50,100),
  #                   limits = c(0,75))+
  theme_cowplot() +
  theme(
    axis.line.x = element_blank(),
    #axis.text.y = element_blank(),
    axis.ticks.x =element_blank(),
    axis.title.x =element_blank(),
    legend.position = "none"#,
    #panel.grid = element_blank(),
    #panel.background = element_blank()
  )+
  labs(y = 'Number of species [log]')+
  scale_fill_brewer(palette = "Set1")

