
library(dplyr)
library(doBy)
library(reshape2)
library(tidyverse)
library(sjPlot)
library(cowplot)
library(patchwork)

source("Scripts/ggsankey_custom_2.R")

load("Final_output/country_dendro_data_asplotted.RData")
load("Figures_July2024/1b_clusters_countries.RData")
load("Host_list/Output/countries_unifrac_clusters1_nohawaii.RData")
load("Figures_July2024/Map_panels.RData")

load("Host_list/Output/ordinated_unifrac_data.Rdata")
bioregions <- unique(bioregions_vector$biogeographic_region) %>% as.vector
bioregions[which(bioregions=="Carribbean")] <- "Caribbean"
bioregions <- c(bioregions, "Australia")

# match color pallette by naming elements
c16 <-
  c(
    Amazonia = "green4",               # 1 Amazonia
    Arabia.Sahara = "#6A3D9A",         # 2 Arabia.Sahara
    C.America = "#FF7F00",             # C.America
    Caribbean = "gold1",
    Eurasian.Palearctic = "#FB9A99",
    Hawaii = "blue1",                   # Hawaii
    Iberia.N.Africa = "palegreen2",     # Iberia.N.Africa
    Middle.East = "#FDBF6F",
    N.North.America = "gray70",               
    N.Pacific = "khaki2",                  # N.Pacific
    NE.North.America = "darkturquoise",
    S.Cone.Pacific = "darkorange4",
    SE.Asia = "#BDCDFF",
    SE.North.America = "#FA0087",
    Subsaharan.Africa = "#2ED9FF",
    W.North.America = "#E31A1C",
    Australia = "yellow" # Australia ADDED ONE COLOR TEMPORARILY
  )

origin_range_specs_dendro <-
  sapply(obj, length)[c(1,9,12,14,15,4,2,6)] %>%
    c(sapply(obj, length)[c(11)] - 1) %>%
    c(sapply(obj, length)[c(10)]) %>%
    c(Hawaii=1) %>%
    c(sapply(obj, length)[c(3,5,7,13,8)]) %>%
    
    # midpoints of the dendrogram clusters
    
    (function (x) data.frame(lab=names(x), region_l=as.vector(x))) %>%
    dplyr::mutate(lab = gsub('Carribbean','Caribbean',lab))%>%
    dplyr::mutate(region_end=cumsum(region_l)+.5) %>%
    dplyr::mutate(region_start=region_end-region_l) %>%
    dplyr::mutate(region_mid = region_start+region_l/2)

########### Make panels for Fig 2

############################################
########### Panel B - pathogens ############
############################################

# LOAD DATA
load("Final_output/Final_analysis_data_pathogens-Sept-2024.RData")
to_analyze_sankey <- 
  to_analyze %>%
#  filter(Geographic.Origin.Cluster != "Unknown") %>%
  make_long(Geographic.Origin.Cluster, Invaded.Range.Cluster, value=Number.of.pests)

# reorder to match dendrogram
dendro_order_source <- levels(factor(to_analyze_sankey$node))[c(15,12,6,8,9,17,10,14,13,11,5,4,1,16,2,7,3)]
to_analyze_sankey$node <-
  factor(to_analyze_sankey$node, levels= c(dendro_order_source))

#dendro_order_sink <- c("Hawaii", "")

# modify values to be proportional to total pests

pathogens <-
  read.csv("Pathogens/Input/Expanded_data_EU_AU_NA_HI_Sept2024.csv")%>%
  dplyr::select(c("Pathogen.Species", "Geographic.Origin.Cluster", "Invaded.Range.Cluster")) %>%
  filter(Geographic.Origin.Cluster != Invaded.Range.Cluster) 

pathogens$Invaded.Range.Cluster <- replace(pathogens$Invaded.Range.Cluster, pathogens$Invaded.Range.Cluster=="S.Cone.Pacific", "Australia")
pathogens$Geographic.Origin.Cluster<-as.factor(pathogens$Geographic.Origin.Cluster)
pathogens <- pathogens%>%
  arrange(Geographic.Origin.Cluster, Invaded.Range.Cluster)

number_of_pathogens_origin <-
  dplyr::select(pathogens, Pathogen.Species, Geographic.Origin.Cluster) %>% distinct() %>%
  dplyr::group_by(Geographic.Origin.Cluster) %>% 
  dplyr::summarise(number_pathogens_origin = n_distinct(Pathogen.Species)) %>%
  ungroup() %>%
  as.data.frame
  
path.sum <- read.csv("Pathogens/Output/path.sum.expanded.csv") %>%
  melt(value.name='pathogens') %>%
  dplyr::rename(Geographic.Origin.Cluster=variable) %>%
  
  # there were some duplicated origins
  
  filter(Geographic.Origin.Cluster != Invaded.Range.Cluster) %>%
  filter(!(grepl('North\\.America', Geographic.Origin.Cluster) & grepl('North\\.America', Invaded.Range.Cluster)))

#number_of_pathogens_origin

  #path.sum %>%
  #group_by(Geographic.Origin.Cluster) %>%
  #summarise(number_pathogens_origin = sum(pathogens))
 
number_of_pathogens_accumulated <-
#  dplyr::select(pathogens, Pathogen.Species, Invaded.Range.Cluster) %>% #distinct() %>%
#  group_by(Invaded.Range.Cluster) %>% 
#  summarise(number_pathogens_invaded = n_distinct(Pathogen.Species))%>%
#  as.data.frame
  
  path.sum %>%
  dplyr::group_by(Invaded.Range.Cluster) %>%
  dplyr::summarise(number_pathogens_invaded = sum(pathogens))%>%
  ungroup() %>%
  as.data.frame

to_analyze_sankey <-
  to_analyze_sankey %>%
  dplyr::filter(x == "Geographic.Origin.Cluster") %>%
  dplyr::left_join(number_of_pathogens_origin, by = join_by("node" == "Geographic.Origin.Cluster")) %>%
  dplyr::left_join(number_of_pathogens_accumulated, by = join_by("next_node" == "Invaded.Range.Cluster")) %>%
  dplyr::bind_rows(
    to_analyze_sankey %>%
      dplyr::filter(x == "Invaded.Range.Cluster") %>%
      dplyr::select(-value) %>%
      dplyr::left_join(number_of_pathogens_accumulated, by = join_by("node" == "Invaded.Range.Cluster"))
  ) %>% distinct()

# get centerpoints for clusters


##########
########## FROM PLANTS

# scaling factor
scaling_factor_1 <-  .75

#origin_range_specs_trees <-
#  origin_range_specs_trees %>%
#  dplyr::mutate(flow_l = ifelse(is.na(flow_l),0,flow_l* scaling_factor_1)) %>% # scale


origin_range_specs <-
  origin_range_specs_dendro%>%
  
    # midpoints of the total flow for source regions
  
    left_join(number_of_pathogens_origin, by = join_by('lab' == 'Geographic.Origin.Cluster')) %>%
    dplyr::rename(flow_l = number_pathogens_origin) %>%
    dplyr::mutate(flow_l = ifelse(is.na(flow_l),0,flow_l* scaling_factor_1)) %>%
    dplyr::mutate(flow_mid = ifelse(flow_l == 0,NA, region_mid))  %>%
    dplyr::mutate(flow_end= flow_mid+ flow_l/2) %>%
    dplyr::mutate(flow_start= flow_mid- flow_l/2) %>%
    dplyr::arrange(flow_start) %>%
    dplyr::mutate(leading_gap = flow_start - lag(flow_end, 1)) %>%
    dplyr::arrange(region_start) %>%
  
    # add column to make it bindable
    
    cbind(x= 'Geographic.Origin.Cluster', .)

space <- 2

# FROM TREES
#space <- 2

#invaded_range_specs_trees <-
#  invaded_range_specs_trees %>%
#  dplyr::mutate(flow_l = flow_l * scaling_factor_1) %>%
#  dplyr::mutate(flow_end= cumsum(flow_l) +.5 + space*(as.numeric(lab)- 1)) %>%
#  dplyr::mutate(flow_start= flow_end - flow_l) %>%
#  dplyr::mutate(flow_mid = flow_start + flow_l/2) %>%
  
  # midpoints for sink regions are the same
#  dplyr::mutate(region_l = flow_l) %>%

invaded_range_specs <-
  data.frame(
    x = 'Invaded.Range.Cluster',
    lab = c("Hawaii", "W.North.America", "N.North.America", "NE.North.America", "SE.North.America", "Eurasian.Palearctic", "Australia")
  ) %>%
  
  # for accumulated pathogens it should fit together nicely
  # midpoints of the total flow for sink regions
  
  left_join(number_of_pathogens_accumulated, by = join_by('lab' == 'Invaded.Range.Cluster')) %>%
  
  dplyr::mutate(lab = factor(lab, levels=c("Hawaii", "W.North.America",
                                           "N.North.America", "NE.North.America",
                                           "SE.North.America", "Eurasian.Palearctic",
                                           "Australia"))) %>%
  
  dplyr::rename(flow_l = number_pathogens_invaded) %>%
  dplyr::mutate(flow_l = flow_l * scaling_factor_1) %>% # PASTED FROM TREES
  dplyr::mutate(flow_end= cumsum(flow_l) +.5 + space*(as.numeric(lab)- 1)) %>%
  dplyr::mutate(flow_start= flow_end - flow_l) %>%
  dplyr::mutate(flow_mid = flow_start + flow_l/2) %>%

  # midpoints for sink regions are the same
  dplyr::mutate(region_l = flow_l) %>%
  dplyr::mutate(region_end=flow_end) %>%
  dplyr::mutate(region_start=flow_start) %>%
  dplyr::mutate(region_mid = flow_mid) %>%
  dplyr::arrange(flow_start) %>%
  dplyr::mutate(leading_gap = flow_start - lag(flow_end, 1)) %>%
  dplyr::arrange(region_start)

flow_specs <-
  rbind(origin_range_specs, invaded_range_specs) %>%
    relocate(x, lab, region_l, region_mid, region_start, region_end, flow_l, flow_mid, flow_start, flow_end)

total_flow_invaded <- sum(number_of_pathogens_accumulated$number_pathogens_invaded)
total_flow_origin<- sum(number_of_pathogens_origin$number_pathogens_origin)

to_analyze_sankey2 <-
  to_analyze_sankey %>% left_join(
    flow_specs,
    by=join_by('x', 'node' == 'lab')
  )

#################################

to_analyze_sankey2 %>% dplyr::filter(x=='Invaded.Range.Cluster') %>% arrange(region_start) %>% select(node, region_start, region_end) %>% distinct

dendro_order_sink <- to_analyze_sankey2 %>% dplyr::filter(x=='Invaded.Range.Cluster') %>% arrange(region_start) %>% select(node, region_start, region_end) %>% distinct %>%
  select(node) %>% unlist %>% as.vector

to_analyze_sankey$node <- factor(to_analyze_sankey$node, levels=dendro_order_source)
to_analyze_sankey$next_node <- factor(to_analyze_sankey$next_node, levels=dendro_order_sink)
to_analyze_sankey2$node <- factor(to_analyze_sankey2$node, levels=dendro_order_source)
to_analyze_sankey2$next_node <- factor(to_analyze_sankey2$next_node, levels=dendro_order_sink)

#source("Scripts/ggsankey_custom_2.R")

### FROM TREES
#layer_positions_flow_source_trees <-
#  to_analyze_sankey_exotic_plants2 %>%
#  dplyr::filter(x == 'Geographic.Origin.Cluster') %>%
#  dplyr::arrange(as.numeric(node), as.numeric(next_node)) %>%
#  dplyr::group_by(node) %>%
  
#  dplyr::do(arrange_intervals3(interval_lengths = .$value * scaling_factor_1, linked_names = .$next_node, total = unique(.$flow_l) , .$node)) %>% # need to scale again here
  

layer_positions_flow_source <-
  to_analyze_sankey2 %>%
  dplyr::filter(x == 'Geographic.Origin.Cluster') %>%
  dplyr::arrange(as.numeric(node), as.numeric(next_node)) %>%
  dplyr::group_by(node) %>%
  
  #dplyr::do(arrange_intervals2(interval_lengths = .$value, linked_names = .$next_node, total = unique(.$flow_l), .$node)) %>%
  
  # MODIFIED W SCALING FACTOR CODE FROM TREES
  # ORIGINAL:   dplyr::do(arrange_intervals3(interval_lengths = .$value, linked_names = .$next_node, total = unique(.$flow_l), .$node)) %>%
  dplyr::do(arrange_intervals3(interval_lengths = .$value * scaling_factor_1, linked_names = .$next_node, total = unique(.$flow_l), .$node)) %>%
  
  dplyr::mutate(x = 'Geographic.Origin.Cluster') %>%
  
  ### ???
  ungroup() %>%
  tidyr::drop_na(start) %>%
  
  right_join(to_analyze_sankey2 %>%
               dplyr::filter(x == 'Geographic.Origin.Cluster'),
             by = join_by('x', 'node', 'linked_names' == 'next_node'))

###### FROM TREES
#layer_positions_flow_sink_trees <-
#  layer_positions_flow_sink_trees %>%
#  
#  dplyr::group_by(next_node) %>%
#  
#  #dplyr::mutate(end = region_start + cumsum(value* scaling_factor_1)) %>%
#  #dplyr::mutate(start = end - value * scaling_factor_1) %>%
#  
#  dplyr::do(arrange_intervals3(interval_lengths = .$value * scaling_factor_1, linked_names = .$node, total = unique(.$flow_l), .$node)) %>%
#  dplyr::relocate(start, .before = stop) %>%
#  dplyr::mutate(x = 'Geographic.Origin.Cluster') %>%
#
#
### ???
#  ungroup() %>%
#    tidyr::drop_na(start) %>%
# 
#  #layer_positions_flow_sink_trees %>%
#  #  dplyr::mutate(span = stop - start) %>%
#  #  group_by(next_node) %>%
#  #  summarise(total_bandwidth = sum(span))
#  
#  # might need to change where the midpoints go
#  
#  #layer_positions_flow_sink_trees <-
#  #  layer_positions_flow_sink_trees %>%
#  
#  right_join(layer_positions_flow_sink_trees %>%
#               dplyr::filter(x == 'Geographic.Origin.Cluster'),
#             by = join_by('x', 'next_node', 'linked_names' == 'node')) %>%
#  
#  dplyr::rename(node = linked_names)

#dplyr::select(x, next_x, node, next_node,
#              flow_end_ymin, flow_end_ymax)

layer_positions_flow_sink <-
  to_analyze_sankey2 %>%
  dplyr::filter(x == 'Geographic.Origin.Cluster') %>%
  dplyr::select(x, node, next_x, next_node, value) %>%
  
  left_join(flow_specs %>%
              dplyr::filter(x == 'Invaded.Range.Cluster') %>%
              dplyr::mutate(x = 'Geographic.Origin.Cluster'),
            by = join_by('x','next_node' == 'lab')) %>%
  
  mutate(next_node = factor(next_node, levels=dendro_order_sink)) %>%
  
  arrange(as.numeric(next_node), as.numeric(node)) %>%
  
  dplyr::group_by(next_node) %>%
  
  dplyr::mutate(end = region_start + cumsum(value*scaling_factor_1)) %>%
  dplyr::mutate(start = end - value*scaling_factor_1) %>%
  dplyr::ungroup() %>%
  dplyr::relocate(start, .before = end) %>%

  dplyr::rename(flow_end_ymin  = start) %>%
  dplyr::rename(flow_end_ymax = end) %>%

  dplyr::select(x, next_x, node, next_node,
                flow_end_ymin, flow_end_ymax)

layer_positions <-
  layer_positions_flow_source %>%
    dplyr::rename(next_node = linked_names) %>%
    dplyr::mutate(flow_start_ymin = flow_start + flow_l*start) %>%
    dplyr::mutate(flow_start_ymax = flow_start + flow_l*stop) %>%
    dplyr::left_join(layer_positions_flow_sink) %>%
  dplyr::select(
    x,
    node,
    next_x,
    next_node,
    flow_start_ymax,
    flow_start_ymin,
    flow_end_ymax,
    flow_end_ymin
  )

bound_data <-
  left_join(to_analyze_sankey, layer_positions) %>%
    dplyr::mutate(x = as.factor(x)) %>%
    dplyr::mutate(node = factor(node, levels=dendro_order_source)) %>%
    dplyr::mutate(next_node = factor(next_node, levels=dendro_order_sink)) %>%
  
    # add the node information
    # for Geographic.Origin.Cluster we want flow_source region_start and region_end
  
  left_join(
    layer_positions_flow_source %>%
      dplyr::select(x,
                    node,
                    ymin=flow_start,
                    ymax=flow_end,                     # this is if you want to remove source nodes
                    label_y=region_mid) %>% distinct %>%#dplyr::mutate(ymin=NA,ymax=NA) %>%
  
      # for Invaded.Range.Cluster we want to_analyze_sankey2 flow_start and flow_end
    
      bind_rows(
        to_analyze_sankey2 %>%
          dplyr::filter(x == "Invaded.Range.Cluster") %>%
          dplyr::select(x, node, label_y=flow_mid, ymin=flow_start, ymax=flow_end) %>% distinct),
    by = c('x','node'))

# gonna need to center it
# need to adjust center positions for sources and add padding to dendrogram

# adjust source positions
top.width <- length(unlist(obj))
bottom.width <- max(na.omit(bound_data$ymax))

scaling.factor <- bottom.width/top.width
  #1.5

bound_data2 <-
  bound_data %>%
  mutate(flow_start_ymin = flow_start_ymin*scaling.factor,
         flow_start_ymax = flow_start_ymax*scaling.factor,
         ymin = ifelse(x=='Geographic.Origin.Cluster', ifelse(node %in% c('Caribbean','Iberia.N.Africa'), NA, ymin*scaling.factor), ymin),
         ymax = ifelse(x=='Geographic.Origin.Cluster', ifelse(node %in% c('Caribbean','Iberia.N.Africa'), NA, ymax*scaling.factor), ymax))

bound_data2$ymin
bound_data2$ymax


sankeyplot_pathogens <-
  ggplot(data = bound_data2) +
  geom_sankey2(
    mapping = aes(
      node = node,
      next_node = next_node,
      x = x,
      next_x = next_x,
      value= value,
      fill= node,
      label = node,
      flow_start_ymax=flow_start_ymax,
      flow_start_ymin =flow_start_ymin ,
      flow_end_ymax=flow_end_ymax,
      flow_end_ymin =flow_end_ymin,
      ymin=ymin,
      ymax=ymax
    ),
    flow.alpha=0.6,
    node.alpha=0.6) + 
  geom_text(
    data = distinct(select(bound_data2, node, number_pathogens_origin, ymin, ymax)) %>% na.omit(),
    mapping = aes(
      y = ymin + (ymax-ymin)/2,
      x = 'Geographic.Origin.Cluster',
      label= number_pathogens_origin
    ),
    size=2.5,
    nudge_x = .01,
    angle=90
  ) +
  geom_text(
    data = distinct(select(filter(bound_data2, x == 'Invaded.Range.Cluster'), node, number_pathogens_invaded, ymin, ymax)) %>% na.omit(),
    mapping = aes(
      y = ymin + (ymax-ymin)/2,
      x = 'Invaded.Range.Cluster',
      label= number_pathogens_invaded
    ),
    size=2.5,
    nudge_x = .01
  )+
  labs(x=NULL, y=NULL) +
  theme_sankey() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'none') +
  scale_fill_manual(values=c16) +
  scale_x_discrete(limits = c('Invaded.Range.Cluster',
                              'Geographic.Origin.Cluster'),
                   expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_flip()

#windows(); sankeyplot_pathogens

# center dendrogram
# checked and looks good for now...

#windows(); gg.dendro.countries +
#  scale_x_discrete(expand=c(0,0))

# gonna need to add the sink cluster map portions

# windows();WorldMap2

p2 <- gg.dendro.countries +
  scale_x_discrete(expand=c(0,0))

to_loop <- bound_data2 %>% filter(x == 'Invaded.Range.Cluster') %>%
  select(node, ymin, ymax)%>%
  arrange(as.numeric(factor(node, levels =
                              c("Hawaii","W.North.America","N.North.America",
                                "NE.North.America","SE.North.America",
                                "Eurasian.Palearctic","Australia"))))

layout <- c(
  area(t = 1, l = 1, b = 149, r = max(to_loop$ymax)),
  area(t = 150, l = 1, b = 399, r = max(to_loop$ymax))
)

for (i in 1:7) {
  layout <- c(layout,
    with(to_loop[i,], area(t=400, b= 499, l = ymin, r = ymax)))
}

panel3<- (p2 + theme(plot.margin = margin(0,0,0,0))) +
  (sankeyplot_pathogens+ theme(plot.margin = margin(0,0,0,0))) +
  #(sink_map_HI + geom_text(mapping = aes(x=x,y=y,label=label), data=data.frame(x=.5,y=.5,label=16)) + theme(plot.margin = margin(0,0,0,0))) +
  (sink_map_HI + theme(plot.margin = margin(0,0,0,0))) +
  (sink_map_WNA +theme(plot.margin = margin(0,0,0,0))) +
  (sink_map_NNA +theme(plot.margin = margin(0,0,0,0))) +
  (sink_map_NEN +theme(plot.margin = margin(0,0,0,0))) +
  (sink_map_SEN +theme(plot.margin = margin(0,0,0,0))) +
  (sink_map_EU +theme(plot.margin = margin(0,0,0,0))) +
  (sink_map_AU +theme(plot.margin = margin(0,0,0,0))) +
  plot_layout(
    design = layout)

windows(5,5); panel3

#######################################################
########### Panel C - invasive plants #################
#######################################################

load("Covariates_input/AU_EU_trees_list.RData")
load("Covariates_input/HI_trees_list.RData")
load("Covariates_input/NA_trees_list.RData")

AU_EU_exotichost_occurrence <- read.csv('Covariates_input/introduced_spp_AU_EU.csv') %>%
  dplyr::group_by(country) %>%
  dplyr::summarise(numspp = n_distinct(sciname)) %>% ungroup
  
trees.countries <- read.csv('Host_list/Output/trees.countries.all.translated.csv')

HI_exotichost_occurrence <- read.csv('Host_list/Output/HI_exotic_host_occurrence_matrix.csv') %>%
  filter(X %in% trees.countries$scientificName) %>%
  dplyr::filter(Hawaii > 0) %>%
  dplyr::group_by(Hawaii) %>%
  dplyr::summarise(numspp = n_distinct(X)) %>% ungroup %>%
  dplyr::mutate(country = 'Hawaii') %>%
  dplyr::select(-Hawaii)

NA_exotichost_occurrence <- read.csv('Host_list/Output/exotic_host_occurrence_matrix.csv') %>%
  melt %>%
  dplyr::filter(value > 0) %>%
  dplyr::select(-value) %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(numspp = n_distinct(X)) %>% ungroup %>%
  dplyr::rename(country = variable)


NAM.total.exotic.trees <- (read.csv('Host_list/Output/exotic_host_occurrence_matrix.csv', row.names=1) %>%
  (function (x) x[(rowSums(x) > 0),]) %>%
  dim)[1]

read.csv('Host_list/Output/exotic_host_occurrence_matrix.csv', row.names=1) %>%
  #(function (x) x[(rowSums(x) > 0),]) %>%
  dim

number_of_exotics_accumulated <-
  dplyr::bind_rows(
    AU_EU_exotichost_occurrence,
    HI_exotichost_occurrence,
    NA_exotichost_occurrence
  ) %>%
  dplyr::rename(Invaded.Range.Cluster = country) %>%
  dplyr::rename(number_exotics_sink = numspp) %>%
  dplyr::mutate(Invaded.Range.Cluster = sub( 'Europe', 'Eurasian.Palearctic', Invaded.Range.Cluster))

#names(trees.list1)
#names(trees.list2)
#names(trees.list3)

#setdiff(names(trees.list1), names(trees.list2))
#setdiff(names(trees.list2), names(trees.list1))
#setdiff(names(trees.list3), names(trees.list1))
#setdiff(names(trees.list3), names(trees.list2))
#setdiff(names(trees.list1), names(trees.list3))
#setdiff(names(trees.list2), names(trees.list3))

trees.list3 <- trees.list3[-which(names(trees.list3)=='Hawaii')]
trees.list3$NE.North.America <- character()
trees.list3$W.North.America <- character()
trees.list3$N.North.America <- character()
trees.list3$SE.North.America <- character()

trees.list.all <- list()
number_of_exotic_trees_origin <- NULL

for (i in names(trees.list1)) {
  trees.list.all[[length(trees.list.all) + 1]] <- 
    union(
      trees.list1[[i]],
      trees.list2[[i]]) %>%
      union(trees.list3[[i]])
  names(trees.list.all)[length(trees.list.all)] <- i
  number_of_exotic_trees_origin <- rbind(number_of_exotic_trees_origin, data.frame(Geographic.Origin.Cluster = i,
                                                             number_exotics_origin = length(trees.list.all[[i]])))
}

number_of_exotic_trees_origin <-
  number_of_exotic_trees_origin %>%
  dplyr::mutate(Geographic.Origin.Cluster = sub('Carribbean','Caribbean',Geographic.Origin.Cluster))

# just for safekeeping
save(trees.list.all, file='Host_list/Output/exotic_tree_originating_each_region.RData')
write.csv(number_of_exotic_trees_origin, 'Host_list/Output/total_exotic_tree_originating_each_region.csv', row.names=F)

# now on to the next thing
to_analyze_sankey_exotic_plants <- 
  to_analyze %>%
  make_long(Geographic.Origin.Cluster, Invaded.Range.Cluster, value=number_of_exotic_trees)

# reorder to match dendrogram
dendro_order_source_trees <- levels(factor(to_analyze_sankey_exotic_plants$node))[c(15,12,6,8,9,17,10,14,13,11,5,4,1,16,2,7,3)]

to_analyze_sankey_exotic_plants$node <-
  factor(to_analyze_sankey_exotic_plants$node, levels= c(dendro_order_source_trees))

# modify values to be proportional to total pests

to_analyze_sankey_exotic_plants <-
  to_analyze_sankey_exotic_plants %>%
  dplyr::filter(x == "Geographic.Origin.Cluster") %>%
  dplyr::left_join(number_of_exotic_trees_origin, by = join_by("node" == "Geographic.Origin.Cluster")) %>%
  dplyr::left_join(number_of_exotics_accumulated, by = join_by("next_node" == "Invaded.Range.Cluster")) %>%
  dplyr::bind_rows(
    to_analyze_sankey_exotic_plants %>%
      dplyr::filter(x == "Invaded.Range.Cluster") %>%
      dplyr::select(-value) %>%
      dplyr::left_join(number_of_exotics_accumulated, by = join_by("node" == "Invaded.Range.Cluster"))
  ) %>% distinct()

##### now make the layers

# get centerpoints for clusters

origin_range_specs_trees <-
  origin_range_specs_dendro %>%
  
  # midpoints of the total flow for source regions
  
  left_join(number_of_exotic_trees_origin, by = join_by('lab' == 'Geographic.Origin.Cluster')) %>%
  dplyr::rename(flow_l = number_exotics_origin)

#exotics_total_bandwidth<- to_analyze_sankey_exotic_plants %>%
#    filter(x == 'Geographic.Origin.Cluster') %>%
#    group_by(next_node) %>%
#    summarise(total = sum(value))

invaded_range_specs_trees <-
  data.frame(
    x = 'Invaded.Range.Cluster',
    lab = c("Hawaii", "W.North.America", "N.North.America", "NE.North.America", "SE.North.America", "Eurasian.Palearctic", "Australia")
  ) %>%
  
  # for accumulated pathogens it should fit together nicely
  # midpoints of the total flow for sink regions
  
  left_join(number_of_exotics_accumulated, by = join_by('lab' == 'Invaded.Range.Cluster')) %>%
  
  dplyr::mutate(lab = factor(lab, levels=c("Hawaii", "W.North.America",
                                           "N.North.America", "NE.North.America",
                                           "SE.North.America", "Eurasian.Palearctic",
                                           "Australia"))) %>%
  
  dplyr::rename(flow_l = number_exotics_sink)

#invaded_range_specs_trees
#exotics_total_bandwidth

# for Hawaii, its greater
# for W north America, its greater
# for N, its greater
# for NE, its greater
# for SE, its greater
# for EU, its greater
# for AU, its greater - for all its greater

# need to determine scaling based on total flow and width of dendrogram

#origin_range_specs_trees
#invaded_range_specs_trees

dendro.width <- length(unlist(obj))
source.width.max <- origin_range_specs_trees$flow_l %>% na.omit %>% sum
sink.width.max <- invaded_range_specs_trees$flow_l %>% na.omit %>% sum

#dendro.width
#source.width.max
#sink.width.max

# scaling factor
scaling_factor_1 <- #dendro.width/source.width.max
  .075

origin_range_specs_trees <-
  origin_range_specs_trees %>%
  dplyr::mutate(flow_l = ifelse(is.na(flow_l),0,flow_l* scaling_factor_1)) %>% # scale
  dplyr::mutate(flow_mid = ifelse(flow_l == 0,NA, region_mid))  %>%
  dplyr::mutate(flow_end= flow_mid+ flow_l/2) %>%            # scaled
  dplyr::mutate(flow_start= flow_mid- flow_l/2) %>%          # scaled
  dplyr::arrange(flow_start) %>%
  dplyr::relocate(flow_start, .before = flow_end) %>%
  
  dplyr::mutate(leading_gap = flow_start - lag(flow_end, 1)) %>%
  #dplyr::mutate(leading_gap2 = flow_start - lag(flow_start, 1)) %>%
  #dplyr::mutate(leading_gap3 = flow_start - lag(flow_end, 2)) %>%
  #dplyr::mutate(trailing_gap = ) %>%
  #dplyr::mutate(trailing_gap2 = ) %>%
  #dplyr::mutate(trailing_gap3 = ) %>%
  dplyr::arrange(region_start) %>%
  
  # add column to make it bindable
  
  cbind(x= 'Geographic.Origin.Cluster', .)

space <- 2

invaded_range_specs_trees <-
  invaded_range_specs_trees %>%
  dplyr::mutate(flow_l = flow_l * scaling_factor_1) %>%
  dplyr::mutate(flow_end= cumsum(flow_l) +.5 + space*(as.numeric(lab)- 1)) %>%
  dplyr::mutate(flow_start= flow_end - flow_l) %>%
  dplyr::mutate(flow_mid = flow_start + flow_l/2) %>%
  
  # midpoints for sink regions are the same
  dplyr::mutate(region_l = flow_l) %>%
  dplyr::mutate(region_end=flow_end) %>%
  dplyr::mutate(region_start=flow_start) %>%
  dplyr::mutate(region_mid = flow_mid) %>%
  dplyr::arrange(flow_start) %>%
  dplyr::mutate(leading_gap = flow_start - lag(flow_end, 1)) %>%
  dplyr::arrange(region_start)

origin_range_specs_trees
invaded_range_specs_trees

# could center it at the end ## ALSO ABOVE WITH PATHOGENS

flow_specs_trees <-
  rbind(origin_range_specs_trees, invaded_range_specs_trees) %>%
  relocate(x, lab, region_l, region_mid, region_start, region_end, flow_l, flow_mid, flow_start, flow_end)

#total_flow_invaded_trees <- sum(number_of_exotic_trees_origin$number_exotics_origin)
#total_flow_origin_trees  <- sum(number_of_exotics_accumulated$number_exotics_sink)

to_analyze_sankey_exotic_plants2 <-
  to_analyze_sankey_exotic_plants %>% left_join(
    flow_specs_trees,
    by=join_by('x', 'node' == 'lab')
  )

#################################

to_analyze_sankey_exotic_plants2 %>% dplyr::filter(x=='Invaded.Range.Cluster') %>% arrange(region_start) %>% select(node, region_start, region_end) %>% distinct

dendro_order_sink_trees <- to_analyze_sankey_exotic_plants2 %>% dplyr::filter(x=='Invaded.Range.Cluster') %>% arrange(region_start) %>% select(node, region_start, region_end) %>% distinct %>%
  select(node) %>% unlist %>% as.vector

to_analyze_sankey_exotic_plants$node <- factor(to_analyze_sankey_exotic_plants$node, levels=dendro_order_source_trees)
to_analyze_sankey_exotic_plants$next_node <- factor(to_analyze_sankey_exotic_plants$next_node, levels=dendro_order_sink_trees)
to_analyze_sankey_exotic_plants2$node <- factor(to_analyze_sankey_exotic_plants2$node, levels=dendro_order_source_trees)
to_analyze_sankey_exotic_plants2$next_node <- factor(to_analyze_sankey_exotic_plants2$next_node, levels=dendro_order_sink_trees)

#to_analyze_sankey_exotic_plants2 %>% filter(x=='Geographic.Origin.Cluster') %>% arrange(as.numeric(node)) %>% print(n=100)
#to_analyze_sankey_exotic_plants2 %>% filter(x=='Invaded.Range.Cluster') %>% arrange(as.numeric(node)) %>% print(n=100)

# up to now everthing looks ok

layer_positions_flow_source_trees <-
  to_analyze_sankey_exotic_plants2 %>%
  dplyr::filter(x == 'Geographic.Origin.Cluster') %>%
  dplyr::arrange(as.numeric(node), as.numeric(next_node)) %>%
  dplyr::group_by(node) %>%
  
  dplyr::do(arrange_intervals3(interval_lengths = .$value * scaling_factor_1, linked_names = .$next_node, total = unique(.$flow_l) , .$node)) %>% # need to scale again here
  
  dplyr::mutate(x = 'Geographic.Origin.Cluster') %>%
  
  ### ???
  ungroup() %>%
  tidyr::drop_na(start)

layer_positions_flow_source_trees <-
  
  layer_positions_flow_source_trees %>%
  
  right_join(to_analyze_sankey_exotic_plants2 %>%
               dplyr::filter(x == 'Geographic.Origin.Cluster'),
             by = join_by('x', 'node', 'linked_names' == 'next_node'))

#layer_positions_flow_source_trees %>% arrange(as.numeric(node)) %>% print(n=100)
#layer_positions_flow_source_trees %>% arrange(as.numeric(linked_names)) %>% print(n=100)

layer_positions_flow_sink_trees <-
  to_analyze_sankey_exotic_plants2 %>%
  dplyr::filter(x == 'Geographic.Origin.Cluster') %>%
  dplyr::select(x, node, next_x, next_node, value)%>%#, number_exotics_sink) %>%
  
  left_join(flow_specs_trees %>%
              dplyr::filter(x == 'Invaded.Range.Cluster') %>%
              dplyr::mutate(x = 'Geographic.Origin.Cluster'),
            by = join_by('x','next_node' == 'lab')) %>%
  
  dplyr::mutate(next_node = factor(next_node, levels=dendro_order_sink_trees)) %>%
  
  dplyr::arrange(as.numeric(next_node), as.numeric(node))

#layer_positions_flow_sink_trees %>% print(n=100)

layer_positions_flow_sink_trees <-
  layer_positions_flow_sink_trees %>%
  
  dplyr::group_by(next_node) %>%
  
  #dplyr::mutate(end = region_start + cumsum(value* scaling_factor_1)) %>%
  #dplyr::mutate(start = end - value * scaling_factor_1) %>%
  
  dplyr::do(arrange_intervals3(interval_lengths = .$value * scaling_factor_1, linked_names = .$node, total = unique(.$flow_l), .$node)) %>%
  dplyr::relocate(start, .before = stop) %>%
  dplyr::mutate(x = 'Geographic.Origin.Cluster') %>%
  
  ### ???
  ungroup() %>%
  tidyr::drop_na(start) %>%
  
#layer_positions_flow_sink_trees %>%
#  dplyr::mutate(span = stop - start) %>%
#  group_by(next_node) %>%
#  summarise(total_bandwidth = sum(span))

# might need to change where the midpoints go

#layer_positions_flow_sink_trees <-
#  layer_positions_flow_sink_trees %>%
  
  right_join(layer_positions_flow_sink_trees %>%
               dplyr::filter(x == 'Geographic.Origin.Cluster'),
             by = join_by('x', 'next_node', 'linked_names' == 'node')) %>%
  
  dplyr::rename(node = linked_names)
  
  #dplyr::select(x, next_x, node, next_node,
  #              flow_end_ymin, flow_end_ymax)

layer_positions_flow_sink_trees %>% print(n=100)

layer_positions_trees_p1 <-
  layer_positions_flow_source_trees %>%
  dplyr::rename(next_node = linked_names) %>%
  dplyr::mutate(flow_start_ymin = flow_start + flow_l*start) %>%
  dplyr::mutate(flow_start_ymax = flow_start + flow_l*stop)# %>%


layer_positions_trees_p1 %>% print(n=100)

layer_positions_trees_p2 <-
  layer_positions_flow_sink_trees %>%
    dplyr::mutate(flow_end_ymin  = flow_start + flow_l*start) %>%
    dplyr::mutate(flow_end_ymax = flow_start + flow_l*stop) %>%
    dplyr::select(node, next_node, x, next_x, flow_end_ymin, flow_end_ymax)

layer_positions_trees_p2 %>% print(n=100)

layer_positions_trees <-
  dplyr::left_join(layer_positions_trees_p1, layer_positions_trees_p2) %>%
  
  dplyr::select(
    x,
    node,
    next_x,
    next_node,
    flow_start_ymax,
    flow_start_ymin,
    flow_end_ymax,
    flow_end_ymin
  )
#####################################

bound_data_trees <-
  left_join(to_analyze_sankey_exotic_plants, layer_positions_trees) %>%
  dplyr::mutate(x = as.factor(x)) %>%
  dplyr::mutate(node = factor(node, levels=dendro_order_source_trees)) %>%
  dplyr::mutate(next_node = factor(next_node, levels=dendro_order_sink_trees)) %>%
  
  # add the node information
  # for Geographic.Origin.Cluster we want flow_source region_start and region_end
  
  left_join(
    layer_positions_flow_source_trees %>%
      dplyr::select(x,
                    node,
                    ymin=flow_start,
                    ymax=flow_end,                     # this is if you want to remove source nodes
                    label_y=region_mid) %>% distinct %>% #dplyr::mutate(ymin=NA,ymax=NA) %>%
      
      # for Invaded.Range.Cluster we want to_analyze_sankey2 flow_start and flow_end
      
      bind_rows(
        to_analyze_sankey_exotic_plants2 %>%
          dplyr::filter(x == "Invaded.Range.Cluster") %>%
          dplyr::select(x, node, label_y=flow_mid, ymin=flow_start, ymax=flow_end) %>% distinct),
    by = c('x','node'))

#bound_data_trees %>% print(n=100)
#bound_data_trees %>% filter(x=='Geographic.Origin.Cluster') %>%
#  group_by(next_node) %>% summarise(sum(value))

# gonna need to center it
# need to adjust center positions for sources and add padding to dendrogram

# adjust source positions
top.width <- length(unlist(obj))
bottom.width <- max(na.omit(filter(bound_data_trees, x=="Invaded.Range.Cluster")$ymax))

scaling_factor2 <- bottom.width/top.width

bound_data_trees2 <-
  bound_data_trees %>%
  mutate(flow_start_ymin = flow_start_ymin*scaling_factor2,
         flow_start_ymax = flow_start_ymax*scaling_factor2,
         ymin = ifelse(x=='Geographic.Origin.Cluster',
                              ymin*scaling_factor2, ymin),
         ymax = ifelse(x=='Geographic.Origin.Cluster',
                       ymax*scaling_factor2, ymax))


# plot it
sankeyplot_exotic_trees <-
  ggplot(data = bound_data_trees2) +
  geom_sankey2(
    mapping = aes(
      node = node,
      next_node = next_node,
      x = x,
      next_x = next_x,
      value= value,
      fill= node,
      label = node,
      flow_start_ymax=flow_start_ymax,
      flow_start_ymin =flow_start_ymin ,
      flow_end_ymax=flow_end_ymax,
      flow_end_ymin =flow_end_ymin,
      ymin=ymin,
      ymax=ymax
    ),
    flow.alpha=0.6,
    node.alpha=0.6) + 
  geom_text(
    data = distinct(select(bound_data_trees2, node, number_exotics_origin, ymin, ymax)) %>% na.omit(),
    mapping = aes(
      y = ymin + (ymax-ymin)/2,
      x = 'Geographic.Origin.Cluster',
      label= number_exotics_origin
    ),
    size=2.5,
    nudge_y = -1,
    angle=90
  ) +
  geom_text(
    data = distinct(select(filter(bound_data_trees2, x == 'Invaded.Range.Cluster'), node, number_exotics_sink, ymin, ymax)) %>% na.omit(),
    mapping = aes(
      y = ymin + (ymax-ymin)/2,
      x = 'Invaded.Range.Cluster',
      label= number_exotics_sink
    ),
    size=2.5,
    nudge_x = .01  )+
  labs(x=NULL, y=NULL) +
  theme_sankey() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'none') +
  scale_fill_manual(values=c16) +
  scale_x_discrete(limits = c('Invaded.Range.Cluster',
                              'Geographic.Origin.Cluster'),
                   expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_flip()

to_loop_trees <- bound_data_trees2 %>% filter(x == 'Invaded.Range.Cluster') %>%
  select(node, ymin, ymax)%>%
  arrange(as.numeric(factor(node, levels =
                              c("Hawaii","W.North.America","N.North.America",
                                "NE.North.America","SE.North.America",
                                "Eurasian.Palearctic","Australia"))))

layout2 <- c(
  area(t = 1,   l = 1, b = 149, r = max(to_loop_trees$ymax)),
  area(t = 150, l = 1, b = 399, r = max(to_loop_trees$ymax)))

for (i in 1:7) {
  layout2 <- c(layout2,
              with(to_loop_trees[i,], area(t=400,
                                           b= 499,
                                           l = ymin,
                                           r = ymax)))
}

panel2<- (p2 + theme(plot.margin = margin(0,0,0,0))) +
  (sankeyplot_exotic_trees+ theme(plot.margin = margin(0,0,0,0))) +
  (sink_map_HI +theme(plot.margin = margin(0,0,0,0))) +
  (sink_map_WNA +theme(plot.margin = margin(0,0,0,0))) +
  (sink_map_NNA +theme(plot.margin = margin(0,0,0,0))) +
  (sink_map_NEN +theme(plot.margin = margin(0,0,0,0))) +
  (sink_map_SEN +theme(plot.margin = margin(0,0,0,0))) +
  (sink_map_EU +theme(plot.margin = margin(0,0,0,0))) +
  (sink_map_AU +theme(plot.margin = margin(0,0,0,0))) +
  plot_layout(design = layout2)

windows(5,5); panel2

##########################################
########### Panel D - trade ##############
##########################################

to_analyze_sankey_trade <- 
  to_analyze %>%
  make_long(Geographic.Origin.Cluster, Invaded.Range.Cluster, value=Cumulative.Import.Value.T) %>%
  dplyr::mutate(node = factor(node, levels=dendro_order_source))

# get centerpoints for clusters

origin_range_specs_trade <-
  origin_range_specs_dendro %>%
  
  # midpoints of the total flow for source regions
  
  left_join(to_analyze_sankey_trade %>%
              filter(x == 'Geographic.Origin.Cluster') %>%
              dplyr::group_by(node) %>%
              dplyr::summarise(flow_l = sum(value)), by = join_by('lab' == 'node'))

invaded_range_specs_trade <-
  data.frame(
    x = 'Invaded.Range.Cluster',
    lab = c("Hawaii", "W.North.America", "N.North.America", "NE.North.America", "SE.North.America", "Eurasian.Palearctic", "Australia")
  ) %>%
  
  # for accumulated pathogens it should fit together nicely
  # midpoints of the total flow for sink regions
  
  left_join(to_analyze_sankey_trade %>%
              filter(x == 'Geographic.Origin.Cluster') %>%
              dplyr::group_by(next_node) %>%
              dplyr::summarise(flow_l = sum(value)), by = join_by('lab' == 'next_node')) %>%
  
  dplyr::mutate(lab = factor(lab, levels=c("Hawaii", "W.North.America",
                                           "N.North.America", "NE.North.America",
                                           "SE.North.America", "Eurasian.Palearctic",
                                           "Australia")))

# need to determine scaling based on total flow and width of dendrogram

#origin_range_specs_trees
#invaded_range_specs_trees

dendro.width <- length(unlist(obj))
source.width.max <- origin_range_specs_trade$flow_l %>% na.omit %>% sum
sink.width.max <- invaded_range_specs_trade$flow_l %>% na.omit %>% sum

#dendro.width
#source.width.max
#sink.width.max

# scaling factor
scaling_factor_1 <- #dendro.width/source.width.max
  .22

origin_range_specs_trade <-
  origin_range_specs_trade %>%
  dplyr::mutate(flow_l = ifelse(is.na(flow_l),0,flow_l* scaling_factor_1)) %>% # scale
  dplyr::mutate(flow_mid = ifelse(flow_l == 0,NA, region_mid))  %>%
  dplyr::mutate(flow_end= flow_mid+ flow_l/2) %>%            # scaled
  dplyr::mutate(flow_start= flow_mid- flow_l/2) %>%          # scaled
  dplyr::arrange(flow_start) %>%
  dplyr::relocate(flow_start, .before = flow_end) %>%
  
  dplyr::mutate(leading_gap = flow_start - lag(flow_end, 1)) %>%
  #dplyr::mutate(leading_gap2 = flow_start - lag(flow_start, 1)) %>%
  #dplyr::mutate(leading_gap3 = flow_start - lag(flow_end, 2)) %>%
  #dplyr::mutate(trailing_gap = ) %>%
  #dplyr::mutate(trailing_gap2 = ) %>%
  #dplyr::mutate(trailing_gap3 = ) %>%
  dplyr::arrange(region_start) %>%
  
  # add column to make it bindable
  
  cbind(x= 'Geographic.Origin.Cluster', .)

space <- 2

invaded_range_specs_trade <-
  invaded_range_specs_trade %>%
  dplyr::mutate(flow_l = flow_l * scaling_factor_1) %>%
  dplyr::mutate(flow_end= cumsum(flow_l) +.5 + space*(as.numeric(lab)- 1)) %>%
  dplyr::mutate(flow_start= flow_end - flow_l) %>%
  dplyr::mutate(flow_mid = flow_start + flow_l/2) %>%
  
  # midpoints for sink regions are the same
  dplyr::mutate(region_l = flow_l) %>%
  dplyr::mutate(region_end=flow_end) %>%
  dplyr::mutate(region_start=flow_start) %>%
  dplyr::mutate(region_mid = flow_mid) %>%
  dplyr::arrange(flow_start) %>%
  dplyr::mutate(leading_gap = flow_start - lag(flow_end, 1)) %>%
  dplyr::arrange(region_start)

origin_range_specs_trade
invaded_range_specs_trade

# could center it at the end ## ALSO ABOVE WITH PATHOGENS

flow_specs_trade <-
  rbind(origin_range_specs_trade, invaded_range_specs_trade) %>%
  relocate(x, lab, region_l, region_mid, region_start, region_end, flow_l, flow_mid, flow_start, flow_end)

#total_flow_invaded_trees <- sum(number_of_exotic_trees_origin$number_exotics_origin)
#total_flow_origin_trees  <- sum(number_of_exotics_accumulated$number_exotics_sink)

to_analyze_sankey_trade2 <-
  to_analyze_sankey_trade %>% left_join(
    flow_specs_trade,
    by=join_by('x', 'node' == 'lab')
  )

#################################

#to_analyze_sankey_exotic_plants2 %>% dplyr::filter(x=='Invaded.Range.Cluster') %>% arrange(region_start) %>% select(node, region_start, region_end) %>% distinct

dendro_order_sink_trade <- to_analyze_sankey_trade2 %>% dplyr::filter(x=='Invaded.Range.Cluster') %>% arrange(region_start) %>% select(node, region_start, region_end) %>% distinct %>%
  select(node) %>% unlist %>% as.vector

to_analyze_sankey_trade$node <- factor(to_analyze_sankey_trade$node, levels=dendro_order_source)
to_analyze_sankey_trade$next_node <- factor(to_analyze_sankey_trade$next_node, levels=dendro_order_sink)
to_analyze_sankey_trade2$node <- factor(to_analyze_sankey_trade2$node, levels=dendro_order_source)
to_analyze_sankey_trade2$next_node <- factor(to_analyze_sankey_trade2$next_node, levels=dendro_order_sink)

#to_analyze_sankey_exotic_plants2 %>% filter(x=='Geographic.Origin.Cluster') %>% arrange(as.numeric(node)) %>% print(n=100)
#to_analyze_sankey_exotic_plants2 %>% filter(x=='Invaded.Range.Cluster') %>% arrange(as.numeric(node)) %>% print(n=100)

# up to now everthing looks ok

layer_positions_flow_source_trade <-
  to_analyze_sankey_trade2 %>%
  dplyr::filter(x == 'Geographic.Origin.Cluster') %>%
  dplyr::arrange(as.numeric(next_node), as.numeric(node)) %>%
  dplyr::group_by(node) %>%
  
  dplyr::mutate(flow_start_ymax = flow_start + cumsum(value * scaling_factor_1)) %>%
  dplyr::mutate(flow_start_ymin = flow_start_ymax - value * scaling_factor_1) %>%
  
  ungroup()# %>%
  #dplyr::mutate(x = 'Geographic.Origin.Cluster') %>%
  #right_join(to_analyze_sankey_trade2 %>%
  #             dplyr::filter(x == 'Geographic.Origin.Cluster'))

layer_positions_flow_sink_trade <-
  to_analyze_sankey_trade2 %>%
  dplyr::filter(x == 'Geographic.Origin.Cluster') %>%
  dplyr::select(x, node, next_x, next_node, value)%>%#, number_exotics_sink) %>%
  
  left_join(flow_specs_trade %>%
              dplyr::filter(x == 'Invaded.Range.Cluster') %>%
              dplyr::mutate(x = 'Geographic.Origin.Cluster'),
            by = join_by('x','next_node' == 'lab')) %>%
  
  dplyr::mutate(next_node = factor(next_node, levels=dendro_order_sink)) %>%
  dplyr::arrange(as.numeric(next_node), as.numeric(node))%>%

#layer_positions_flow_sink_trade <-
#  layer_positions_flow_sink_trade %>%
  
  dplyr::group_by(next_node) %>%
  dplyr::mutate(flow_end_ymax = flow_start + cumsum(value * scaling_factor_1)) %>%
  dplyr::mutate(flow_end_ymin = flow_end_ymax - value * scaling_factor_1) %>%
  ungroup()

layer_positions_trade_p1 <-
  layer_positions_flow_source_trade

layer_positions_trade_p2 <-
  layer_positions_flow_sink_trade %>%
  dplyr::select(node, next_node, x, next_x, flow_end_ymin, flow_end_ymax)

layer_positions_trade <-
  dplyr::left_join(layer_positions_trade_p1, layer_positions_trade_p2) %>%
  
  dplyr::select(
    x,
    node,
    next_x,
    next_node,
    flow_start_ymax,
    flow_start_ymin,
    flow_end_ymax,
    flow_end_ymin
  )

bound_data_trade <-
  left_join(to_analyze_sankey_trade, layer_positions_trade) %>%
  dplyr::mutate(x = as.factor(x)) %>%
  dplyr::mutate(node = factor(node, levels=dendro_order_source)) %>%
  dplyr::mutate(next_node = factor(next_node, levels=dendro_order_sink)) %>%
  
  left_join(
    layer_positions_flow_source_trade %>%
      dplyr::select(x,
                    node,
                    ymin=flow_start,
                    ymax=flow_end,                     # this is if you want to remove source nodes
                    label_y=region_mid) %>% distinct %>% #dplyr::mutate(ymin=NA,ymax=NA) %>%
      
      bind_rows(
        to_analyze_sankey_trade2 %>%
          dplyr::filter(x == "Invaded.Range.Cluster") %>%
          dplyr::select(x, node, label_y=flow_mid, ymin=flow_start, ymax=flow_end) %>% distinct),
    by = c('x','node'))

top.width <- length(unlist(obj))
bottom.width <- max(na.omit(filter(bound_data_trade, x=="Invaded.Range.Cluster")$ymax))

scaling_factor3 <- bottom.width/top.width

bound_data_trade2 <-
  bound_data_trade %>%
  mutate(flow_start_ymin = flow_start_ymin*scaling_factor3,
         flow_start_ymax = flow_start_ymax*scaling_factor3,
         ymin = ifelse(x=='Geographic.Origin.Cluster',
              ymin*scaling_factor3, ymin),
         ymax = ifelse(x=='Geographic.Origin.Cluster',
              ymax*scaling_factor3, ymax))
# plot it
sankeyplot_trade <-
  ggplot(data = bound_data_trade2) +
  geom_sankey2(mapping = aes(
                 node = node,
                 next_node = next_node,
                 x = x,
                 next_x = next_x,
                 value= value,
                 fill= node,
                 label = node,
                 flow_start_ymax=flow_start_ymax,
                 flow_start_ymin =flow_start_ymin ,
                 flow_end_ymax=flow_end_ymax,
                 flow_end_ymin =flow_end_ymin,
                 ymin=ymin,
                 ymax=ymax
               ),
    flow.alpha=0.6,
    node.alpha=0.6) +
  geom_text(
    data = select(filter(bound_data_trade2, x == 'Geographic.Origin.Cluster'),
                  node, ymin, ymax, value) %>% group_by(node, ymin, ymax) %>%
      summarise(total_trade_flow = round(sum(value))),
    mapping = aes(
      y = ymin + (ymax-ymin)/2,
      x = 'Geographic.Origin.Cluster',
      label= total_trade_flow
    ),
    size=2.5,
    nudge_x = 0.01,
    angle=90
  ) +
  geom_text(
    data = select(filter(bound_data_trade, x == 'Invaded.Range.Cluster'),
                  node, ymin, ymax, value) %>% group_by(node, ymin, ymax) %>%
      summarise(total_trade_flow = round(sum(value))),
    mapping = aes(
      y = ymin + (ymax-ymin)/2,
      x = 'Invaded.Range.Cluster',
      label= total_trade_flow
    ),
    size=2.5,
    nudge_x = .01  )+
  labs(x=NULL, y=NULL) +
  theme_sankey() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'none') +
  scale_fill_manual(values=c16) +
  scale_x_discrete(limits = c('Invaded.Range.Cluster',
                              'Geographic.Origin.Cluster'),
                   expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_flip()

windows();sankeyplot_trade

## put them together

to_loop <- bound_data2 %>% filter(x == 'Invaded.Range.Cluster') %>%
  select(node, ymin, ymax)%>%
  arrange(as.numeric(factor(node, levels =
                              c("Hawaii","W.North.America","N.North.America",
                                "NE.North.America","SE.North.America",
                                "Eurasian.Palearctic","Australia"))))

to_loop_trees <- bound_data_trees2 %>% filter(x == 'Invaded.Range.Cluster') %>%
  select(node, ymin, ymax)%>%
  arrange(as.numeric(factor(node, levels =
                              c("Hawaii","W.North.America","N.North.America",
                                "NE.North.America","SE.North.America",
                                "Eurasian.Palearctic","Australia"))))

to_loop_trade <- bound_data_trade2 %>% filter(x == 'Invaded.Range.Cluster') %>%
  select(node, ymin, ymax)%>% distinct %>%
  arrange(as.numeric(factor(node, levels =
                              c("Hawaii","W.North.America","N.North.America",
                                "NE.North.America","SE.North.America",
                                "Eurasian.Palearctic","Australia"))))

# add numbers to legend
#library(plotscale)

world_map <- WorldMap2 +
  #scale_x_discrete(expand=c(0,0)) +
  theme(plot.margin = margin(0,0,0,0))

#map_width <- layer_scales(world_map)$x$range$range %>% diff
#map_height <- layer_scales(world_map)$y$range$range %>% diff

#map_width/map_height

#(ps <- plotsize(world_map, width=1, height=1))

#(psratio1 <- ps$width/ps$height)

world_map$data$biogeographic_region[world_map$data$biogeographic_region=="Carribbean"]<-"Caribbean"

newdf <-
  data.frame(names = 
               factor(world_map$data$biogeographic_region %>% c("Australia"),
                      levels=dendro_order_source)) %>%
  arrange(as.numeric(names)) %>%
  mutate(extrasp = c(rep(" ", 9), rep("", 8))) %>%
  mutate(numbers = 1:17) %>%
  mutate(annotated = paste(extrasp, numbers, extrasp, "   ", names, sep = ""))

#dendro_order_source

world_map$data$biogeographic_region <- with(newdf[-17,]%>%arrange(as.character(names)),
                                            factor(annotated, levels=c(annotated, "17   Australia")))

c16_b <- c16#[-17]
names(c16_b) <- levels(world_map$data$biogeographic_region)

# read in centroids

centroids <- read.csv('Covariates_input/landarea_and_centroids.csv') %>%
  dplyr::rename(biogeographic_region=X) %>%
  dplyr::mutate(biogeographic_region=sub("Carribbean","Caribbean",biogeographic_region)) %>%
  left_join(newdf, by= join_by('biogeographic_region'=='names'))
  
centroids[which(centroids$biogeographic_region=='Hawaii'),'lon'] <-
  centroids[which(centroids$biogeographic_region=='Hawaii'),'lon'] + 10

centroids[which(centroids$biogeographic_region=='Caribbean'),'lon'] <-
  centroids[which(centroids$biogeographic_region=='Caribbean'),'lon'] + 10

centroids[which(centroids$biogeographic_region=='N.Pacific'),'lon'] <-
  centroids[which(centroids$biogeographic_region=='N.Pacific'),'lon'] + 15

centroids[which(centroids$biogeographic_region=='NE.North.America'),'lon'] <-
  centroids[which(centroids$biogeographic_region=='NE.North.America'),'lon'] - 20
centroids[which(centroids$biogeographic_region=='NE.North.America'),'lat'] <-
  centroids[which(centroids$biogeographic_region=='NE.North.America'),'lat'] - 16

# add centroids to map
world_map2 <- world_map + 
  geom_text(data=centroids, mapping=aes(x=lon, y=lat, label=numbers)) +
  scale_fill_manual(values=c16_b)+
  theme(#legend.position = 'right',
        #legend.key.size = unit(.022, 'npc'), 
        #legend.text = element_text(size=7.5, margin = margin(0,0,0,-9.25)),
        plot.margin = margin(0,0,0,0))

world_map2$data$biogeographic_region <-
  factor(world_map2$data$biogeographic_region,
         levels=newdf$annotated)

windows();world_map2

map_width <- layer_scales(world_map2)$x$range$range %>% diff
map_height <- layer_scales(world_map2)$y$range$range %>% diff

#(ps <- plotsize(world_map, width=1, height=1))
#(psratio2 <- ps$width/ps$height)
#psratio1
#map_width/map_height

#(ps <- plotsize(world_map2, width =map_width, height =map_height))

#(devs <- devsize(world_map2, width = ps$width, height = ps$height))

#windows(); world_map2

#windows(7, 7*devs$height/devs$width); world_map2

#map_width<- devs$width
#map_height<- devs$height

#windows(); world_map2

#windows(7, 3); world_map2

layout_total <- c(
  area(t = 1, l = 1, b = map_height, r = 1+map_width),
  area(t = map_height+1, l = 1, b = map_height+map_height/5, r = map_width/3),
  area(t = map_height+1, l = 1+ map_width/3, b = map_height+map_height/5, r = 2*map_width/3),
  area(t = map_height+1, l = 1+ 2*map_width/3, b = map_height+map_height/5, r = map_width),
  area(t =  map_height+map_height/5+1 , l = 1, b = map_height+4*map_height/5, r = map_width/3),
  area(t =  map_height+map_height/5+1 , l = 1+ map_width/3, b = map_height+4*map_height/5, r = 2*map_width/3),
  area(t =  map_height+map_height/5+1 , l = 1+ 2*map_width/3, b = map_height+4*map_height/5, r = map_width)
)

#layout_p1 <-
#  c(
#    area(t = 1, l = 1, b = map_height, r = 1+map_width)
#  )

##layout_p2 <-
##  c(
#    area(t = 1, l = 1, b = map_height/5, r = map_width/3),
#    area(t = 1, l = 1+ map_width/3, b = map_height/5, r = 2*map_width/3),
#    area(t = +1, l = 1+ 2*map_width/3, b = map_height/5, r = map_width)
#  )

#panel_test <- world_map +
#  (gg.dendro.countries+ theme(plot.margin = margin(0,0,0,0))) +
#  (gg.dendro.countries+ theme(plot.margin = margin(0,0,0,0))) +
#  (gg.dendro.countries+ theme(plot.margin = margin(0,0,0,0))) +
#  plot_layout(
#    design = layout_total)

#windows(map_width/20, 1.2*map_height/20); panel_test

#layout_p3 <- c(
#  area(t =  1, l = 1, b = 4*map_height/5, r = map_width/3)
#)

for (i in 1:7) {
  layout_total <- c(layout_total,
              with(to_loop[i,], area(t=1 + map_height+4*map_height/5,
                                     b= 2*map_height,
                                     #l = (map_width/3) - (map_width/3)*((max(to_loop$ymax) - ymin)/max(to_loop$ymax)),
                                     #r = (map_width/3) - (map_width/3)*((max(to_loop$ymax) - ymax)/max(to_loop$ymax))
                                     #l = (map_width/3)*(ymin/(max(to_loop$ymax)-min(to_loop$ymin))),
                                     #r = (map_width/3)*(ymax/(max(to_loop$ymax)-min(to_loop$ymin)))
                                     
                                     l = 1+(map_width/3-1)*(ymin/(max(to_loop$ymax))),
                                     r = 1+(map_width/3-1)*(ymax/(max(to_loop$ymax)))
                                     
                                     )))
#  layout_p3 <- c(layout_p3,
#                 with(to_loop[i,], area(t=1 + 4*map_height/5,
#                                        b= map_height,
#                                        #l = (map_width/3) - (map_width/3)*((max(to_loop$ymax) - ymin)/max(to_loop$ymax)),
#                                        #r = (map_width/3) - (map_width/3)*((max(to_loop$ymax) - ymax)/max(to_loop$ymax))
#                                        #l = (map_width/3)*(ymin/(max(to_loop$ymax)-min(to_loop$ymin))),
##                                        #r = (map_width/3)*(ymax/(max(to_loop$ymax)-min(to_loop$ymin)))
#                                        
#                                        l = 1+(map_width/3-1)*(ymin/(max(to_loop$ymax))),
#                                        r = 1+(map_width/3-1)*(ymax/(max(to_loop$ymax)))
                                        
#                 )))
}

#layout_p4 <-
#  c(
#    area(t =  map_height+map_height/5+1 , l = 1+ map_width/3, b = map_height+4*map_height/5, r = 2*map_width/3)
#  )

for (i in 1:7) {
  layout_total <- c(layout_total,
                    with(to_loop_trees[i,], area(t=1 + map_height+4*map_height/5,
                                           b= 2*map_height,
                                           l = 1 + map_width/3 + (map_width/3-1)*(ymin/max(to_loop_trees$ymax)),
                                           r = 1 + map_width/3 + (map_width/3-1)*(ymax/max(to_loop_trees$ymax)))))
#  layout_p4 <- c(layout_p4,
#                 with(to_loop_trees[i,], area(t=1 + map_height+4*map_height/5,
#                                              b= 2*map_height,
#                                              l = 1 + map_width/3 + (map_width/3-1)*(ymin/max(to_loop_trees$ymax)),
#                                              r = 1 + map_width/3 + (map_width/3-1)*(ymax/max(to_loop_trees$ymax)))))
}

#layout_p5 <-
#  c(
#    area(t =  map_height+map_height/5+1 , l = 1+ 2*map_width/3, b = map_height+4*map_height/5, r = map_width)
#  )
for (i in 1:7) {
  layout_total <- c(layout_total,
                    with(to_loop_trade[i,],
                         area(
                         #print(c(
                           t=1 + map_height+4*map_height/5,
                                                 b= 2*map_height,
                                                 l = 1 + 2*map_width/3 + (map_width/3-1)*(ymin/max(to_loop_trade$ymax)),
                                                 r = 1 + 2*map_width/3 + (map_width/3-1)*(ymax/max(to_loop_trade$ymax)))))
#  layout_p5 <- c(layout_p5,
#                 with(to_loop_trade[i,],
#                      area(
#                        #print(c(
#                        t=1 + map_height+4*map_height/5,
#                        b= 2*map_height,
#                        l = 1 + 2*map_width/3 + (map_width/3-1)*(ymin/max(to_loop_trade$ymax)),
#                        r = 1 + 2*map_width/3 + (map_width/3-1)*(ymax/max(to_loop_trade$ymax)))))
  
}

gg.dendro.countries2 <-
gg.dendro.countries+ 
  geom_text(data = origin_range_specs_dendro %>%
              left_join(newdf, by = join_by('lab'=='names')) %>%
              filter(lab != 'Hawaii'),
            mapping = aes(
              x = region_mid,
              y = .175,
              label=numbers,
              size=2
            ))+
  scale_x_discrete(expand=c(0,0)) +
  theme(plot.margin = margin(0,0,0,0))

#windows();gg.dendro.countries2

#panel1 <- (world_map2 + 
#             theme(legend.position = 'right',
                   #legend.key.size = unit(.022, 'npc'), 
#                   legend.text = element_text(#size=7.5,
#                     margin = margin(0,0,0,-14)))) #+ plot_layout(layout_p1)
  
#windows();panel1

#panel2 <-   (gg.dendro.countries2+
#  gg.dendro.countries2+
#  gg.dendro.countries2)+
#  plot_layout()

#windows(); panel2

#panel3 <- 
#  (sankeyplot_pathogens + theme(plot.margin = margin(0,0,0,0))) +
#  (sink_map_HI +theme(plot.margin = margin(0,0,0,0))) +
#  (sink_map_WNA +theme(plot.margin = margin(0,0,0,0))) +
#  (sink_map_NNA +theme(plot.margin = margin(0,0,0,0))) +
#  (sink_map_NEN +theme(plot.margin = margin(0,0,0,0))) +
#  (sink_map_SEN +theme(plot.margin = margin(0,0,0,0))) +
#  (sink_map_EU +theme(plot.margin = margin(0,0,0,0))) +
#  (sink_map_AU +theme(plot.margin = margin(0,0,0,0))) +
#  plot_layout(layout_p3)


#(world_map2 + 
#    theme(legend.position = 'right',
#          #legend.key.size = unit(.022, 'npc'), 
#          legend.text = element_text(#size=7.5,
#            margin = margin(0,0,0,-14)))) /
#  (gg.dendro.countries2 + gg.dendro.countries2 + gg.dendro.countries2) /
#  ((sankeyplot_pathogens + theme(plot.margin = margin(0,0,0,0))) +
#  (sankeyplot_exotic_trees+  theme(plot.margin = margin(0,0,0,0))) +
#  (sankeyplot_trade+  theme(plot.margin = margin(0,0,0,0))))

composite_figure <- (world_map2 + 
  theme(legend.position = 'right',
        #legend.key.size = unit(.022, 'npc'), 
        legend.text = element_text(#size=7.5,
                                   margin = margin(0,0,0,-14)))) +
  (gg.dendro.countries2)+
  (gg.dendro.countries2)+
  (gg.dendro.countries2)+
  (sankeyplot_pathogens + theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sankeyplot_exotic_trees+  theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sankeyplot_trade+  theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_HI +theme(plot.margin = margin(0,0,0,0)) + plot_layout(tag_level = 'new')) +
  (sink_map_WNA +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_NNA +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_NEN +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_SEN +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_EU +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_AU +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_HI +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_WNA +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_NNA +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_NEN +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_SEN +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_EU +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_AU +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_HI +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_WNA +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_NNA +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_NEN +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_SEN +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_EU +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_AU +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  plot_layout(
    design = layout_total,
    guides = 'collect') + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size=12, face='bold'))

#sankeyplot_pathogens$data$

#layout_p1<-
  #layout_total
#layout_p1$t <- layout_p1$t[1]
#layout_p1$t <- layout_p1$t[1]
#layout_p1$t <- layout_p1$t[1]
#layout_p1$t <- layout_p1$t[1]

#layout_p2<-
#  layout_total
#layout_p2$t <- layout_p2$t[2:4]
#layout_p2$t <- layout_p2$t[2:4]
#layout_p2$t <- layout_p2$t[2:4]
#layout_p2$t <- layout_p2$t[2:4]


windows(8, 6);composite_figure+theme(legend.position='none')

#######################################################################

###### TRYING TO MOVE LEGEND TO BOTTOM ################

c16_c <- c16_b#[-17]
names(c16_c) <- gsub("\\.", " ", names(c16_c)) %>% paste("   ")

world_map3 <- world_map2+
  scale_fill_manual(values=c16_c)

levels(world_map3$data$biogeographic_region) <-
  gsub("\\.", " ", levels(world_map3$data$biogeographic_region))  %>% paste("   ")

composite_figure2 <- (world_map3 + 
                        theme(legend.position = "bottom",
                              legend.key.size = unit(.03, 'npc'),
                              legend.key.spacing = unit(.0015, 'npc'),
                              legend.text = element_text(
                                size=10,
                                margin = margin(0,0,0,-12)))) +
  (gg.dendro.countries2)+
  (gg.dendro.countries2)+
  (gg.dendro.countries2)+
  (sankeyplot_pathogens + theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sankeyplot_exotic_trees+ theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sankeyplot_trade+  theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_HI +theme(plot.margin = margin(0,0,0,0)) + plot_layout(tag_level = 'new')) +
  (sink_map_WNA +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_NNA +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_NEN +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_SEN +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_EU +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_AU +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_HI +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_WNA +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_NNA +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_NEN +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_SEN +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_EU +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_AU +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_HI +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_WNA +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_NNA +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_NEN +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_SEN +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_EU +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  (sink_map_AU +theme(plot.margin = margin(0,0,0,0))+ plot_layout(tag_level = 'new')) +
  guide_area() +
  plot_layout(
    design = layout_total %>% c(area(t = 279, b = 345, l = 10, r = 350)),
    guides = "collect") + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size=12, face='bold'))

windows(6,7);composite_figure2

save(composite_figure, file='Figures_Sept2024/Fig1_composite_ABCD.RData')

library(devEMF)

emf(file = "../Figures/Fig1_composite_ABCD.emf", 6,7, emfPlus=T)
composite_figure2
dev.off()

svg("../Figures/Fig1_composite_ABCD.svg", 6,7, family="Cambria Math")
composite_figure2
dev.off()


##########################################
