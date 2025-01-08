library(ggsankey)
library(ggplot2)
library(dplyr)
library(gtools)

##### FUNCTIONS

arrange_intervals <- function (interval_lengths, linked_names, total, source_node=NA) {
  
  if (length(interval_lengths) == 0) return (NULL)
  if (total == 0) return (data.frame(segment=1:length(interval_lengths), layer=1, start=NA, stop=NA))
  
  standard_lengths <- interval_lengths/ total # now we can look for things that sum to 1 and then multiply back by total
  lengths_remaining <- 1:length(standard_lengths)
  layers <- list()
  i <- 1
  greatest_span <- 0
  current.combo <- NA
  
  while (length(lengths_remaining) > 0) {
    
    combs.indices <- combinations(length(lengths_remaining), i, lengths_remaining)
    combs <- standard_lengths[combs.indices] %>% matrix(ncol=i)
    
    if(length(lengths_remaining) == i) {
      if (sum(standard_lengths[lengths_remaining]) > greatest_span & sum(standard_lengths[lengths_remaining]) <= 1) {
        greatest_span <- sum(standard_lengths[lengths_remaining])
        current.combo <- lengths_remaining
      } else if (sum(standard_lengths[lengths_remaining]) == 0) {
        layers[[length(layers) + 1]] <- lengths_remaining
        lengths_remaining <- NULL
        break
      }
    } else {
      
      select.1 <- rowSums(combs) <= 1
      select.2 <- max(rowSums(combs)[select.1])
      select.3 <- which(rowSums(combs) == select.2)
      
      if (select.2 == 0) {
        
        layers[[length(layers) + 1]] <- lengths_remaining
        lengths_remaining <- NULL
        break
        
      } else if (select.2 == 1) {
        
        j <- 1
        while(sum(combs.indices[select.3[j],] %in% lengths_remaining)==length(combs.indices[select.3[j],]) & j<=length(select.3)) {
          winners <- combs.indices[select.3[j],]
          layers[[length(layers) + 1]] <- winners
          lengths_remaining <- setdiff(lengths_remaining, winners)
          j <- j + 1
        }
        i <- 1
        greatest_span <- 0
        current.combo <- NA
        next
      } else {
        for (j in select.3) {
          if (rowSums(combs)[j] > greatest_span) {
            greatest_span <- rowSums(combs)[j]
            current.combo <- combs.indices[j,]
          }
        }
      }
    }
    
    if (i + 1 > length(lengths_remaining)) {
      layers[[length(layers) + 1]] <- current.combo
      lengths_remaining <- setdiff(lengths_remaining, current.combo)
      greatest_span <- 0
      current.combo <- NA
      i <- 1
    } else {
      i <- i + 1
    }
  }
  results <- data.frame()
  rightleft <- sample(c(0,1), 1)
  
  layerlengths_table <- table(sapply(layers, length))
  number_of_ones <- layerlengths_table[names(layerlengths_table)==1] -
    sum(standard_lengths==1)
  
  if (1 %in% names(layerlengths_table))
    number_of_ones <- number_of_ones - sum(sapply(which(sapply(layers, length)==1), function (x) standard_lengths[layers[[x]]])==0)
  
  one_single <- as.vector(number_of_ones) == 1
  even_odd  <- as.vector(number_of_ones) %% 2  # 0 for even, 1 for odd
  k <- 1
  
  for (i in 1:length(layers)) {
    if (sum(standard_lengths[layers[[i]]])==0) {
      start <- NA
      stop <- NA
    } else if (length(layers[[i]]) == 1) {
      if (standard_lengths[layers[[i]]] == 1) {
        start <- 0
        stop <- 1
      } else if (one_single == 1) {
        start <- .5 - standard_lengths[layers[[i]]]/2
        stop <- .5 + standard_lengths[layers[[i]]]/2
      } else if (even_odd == 0) {
        if (rightleft == 0) {
          start <- 0
          stop <- standard_lengths[layers[[i]]]
          rightleft <- 1
        } else {
          start <- 1 - standard_lengths[layers[[i]]]
          stop <- 1
          rightleft <- 0
        }
      } else if (even_odd == 1) {
        if (k == number_of_ones) {
          start <- .5 - standard_lengths[layers[[i]]]/2
          stop <- .5 + standard_lengths[layers[[i]]]/2
        } else if (rightleft == 0) {
          start < 0
          stop <- standard_lengths[layers[[i]]]
          rightleft <- 1
          k <- k + 1
        } else if (rightleft == 1) {
          start <- 1 - standard_lengths[layers[[i]]]
          stop <- 1
          rightleft <- 0
          k <- k + 1
        }
      }
    } else {
      space <- (1 - sum(standard_lengths[layers[[i]]]))/(length(layers[[i]])-1)
      stop <- cumsum(standard_lengths[layers[[i]]]) + space*(0:(length(layers[[i]])-1))
      start <- stop - standard_lengths[layers[[i]]]
    }
    results <-
      rbind(results, data.frame(segment = layers[[i]], layer = i, start=start, stop=stop))
  }
  
  results %>% arrange(segment) %>% cbind(linked_names=linked_names)
}

arrange_intervals2 <- function (interval_lengths,        # starting lengths
                                linked_names,            # ordered list
                                total,                   # total distance
                                source_node=NA) {
  
  # start in the order you want
  
  if (length(interval_lengths) == 0) return (NULL)
  if (total == 0) return (data.frame(linked_names = linked_names, start=NA, stop=NA))
  
  standard_lengths <- interval_lengths / total # now we can look for things that sum to 1 and then multiply back by total
  
  starting_midpoints <-
    (1:length(interval_lengths) - .5)/length(interval_lengths)
  
  ending_midpoints <- starting_midpoints
  
  for (i in 1:length(interval_lengths)) {
    left_end <- starting_midpoints[i] - .5 * standard_lengths[i]
    right_end <- starting_midpoints[i] + .5 * standard_lengths[i]
    if (left_end < 0) ending_midpoints[i] <- 0 + standard_lengths[i]/2
    else if (right_end > 1) ending_midpoints[i] <- 1 - standard_lengths[i]/2
  }
  
  return(
    data.frame(
      linked_names = linked_names,
      start = ending_midpoints - .5 * standard_lengths,
      stop = ending_midpoints + .5 * standard_lengths
    )
  )
  
}

arrange_intervals3 <- function (interval_lengths,        # starting lengths
                                linked_names,            # ordered list
                                total,                   # total distance
                                source_node=NA) {
  
  if (length(interval_lengths) == 0) return (NULL)
  if (total == 0) return (data.frame(linked_names = linked_names, start=NA, stop=NA))
  
  ends <- cumsum(interval_lengths)
  starts <- ends - interval_lengths
  
  total_length <- sum(interval_lengths)
  
  original_remaining_space <- total_length - interval_lengths
  new_remaining_space <- total - interval_lengths
  
  space_reduction_factor <- new_remaining_space/original_remaining_space
  
  # want to preserve ratios between starts and ends
  #start_to_remaining_space <- starts/(total_length - interval_lengths)
  #end_to_remaining_space   <- (total_length - ends)/(total_length - interval_lengths)
  
  #print(cbind(starts, interval_lengths, ends, new_remaining_space, start_to_remaining_space, end_to_remaining_space))
  
  # effectively reduced start and end space by ratio of whats left from
  # end-to-end cumulative arrangement and whats left in 1 unit length arrangement
  
  starts_compressed <- starts * space_reduction_factor
  ends_compressed   <- starts_compressed + interval_lengths
  
  # at this point, starts_compressed + interval_lengths + ends_compressed should sum to totals
   print(cbind(starts,
               interval_lengths,
               ends,
               total_length,
               original_remaining_space,
               new_remaining_space,
               space_reduction_factor,
               starts_compressed,
               ends_compressed,
               starts_compressed + interval_lengths + (total_length - ends_compressed),
               total))
  
  return(
    data.frame(
      linked_names = linked_names,
      start = starts_compressed/total,
      stop = ends_compressed/total
    )
  )
  
}

##############

StatSankeyFlow2 <- ggplot2::ggproto("StatSankeyFlow", ggplot2::Stat,
                                    extra_params = c("n_grid", "na.rm", "type", "width", "space", "smooth"),
                                    
                                    setup_data = function(data, params) {
                                      purrr::map_dfr(unique(data$PANEL),
                                                     ~{
                                                       
                                                       #print('1')
                                                       #print(tail(data %>% dplyr::filter(next_node %in% c('Hawaii', 'Australia')), n=20))
                                                       
                                                       df <- 
                                                         data %>%
                                                         #tidyr::drop_na(flow_start_ymin , flow_start_ymax ,flow_end_ymin  ,flow_end_ymax ) %>%
                                                         dplyr::filter(PANEL == .x) %>%
                                                         dplyr::mutate(dplyr::across(c(x, next_x), ~as.numeric(.), .names = ("n_{.col}")))%>%
                                                         #dplyr::mutate(xmin = 1,#n_x - params$width/2,
                                                        #               xmax = 2)#n_x + params$width/2)
                                                         dplyr::mutate(xmin = n_x - .5 - params$width/.2,
                                                                     xmax = n_x  - .5+ params$width/.2)
                                                       
                                                      
                                                       print('2')
                                                       print(tail(df %>% dplyr::filter(next_node %in% c('Hawaii', 'Australia')), n=20))
                                                       
                                                       flows <- df %>%
                                                         dplyr::left_join(df %>%
                                                                            dplyr::select(n_x, node)%>%#, xmin_end = xmin, xmax_end = xmax) %>%
                                                                            dplyr::distinct(),
                                                                          by = c("n_next_x" = "n_x", "next_node" = "node"))
                                                       
                                                       print('3')
                                                       print(tail(flows %>% dplyr::filter(next_node %in% c('Hawaii', 'Australia')), n=20))
                                                       
                                                       flows <- flows%>%
                                                         tidyr::drop_na(n_x, node, next_node, n_next_x)#,
                                                                        #xmax_end, xmin_end)
                                                 
                                                       flows <- flows %>%
                                                         dplyr::select(-n_x, -node,
                                                                       #-xmin,
                                                                       -n_next_x, -next_node)%>%#,-xmax_end) %>%
                                                         dplyr::mutate(group = dplyr::row_number())%>%
                                                         dplyr::mutate(smooth = params$smooth) %>%
                                                         as.data.frame()
                                                       
                                                       print('4')
                                                       print(flows,20)
                                                       
                                                       flows
                                                     })},
                                    
                                    compute_group = function(data, scales) {
                                      out1 <- sigmoid(data$xmax, data$xmin,#_end,
                                                      data$flow_start_ymax, data$flow_end_ymax,
                                                      smooth = data$smooth)
                                      out2 <- sigmoid(data$xmin,#_end,
                                                      data$xmax, data$flow_end_ymin, data$flow_start_ymin,
                                                      smooth = data$smooth)
                                      dplyr::bind_rows(out1, out2)
                                    }
)


# NODE LAYER -------
StatSankeyNode2 <- ggplot2::ggproto("StatSankeyNode", ggplot2::Stat,
                                   extra_params = c("n_grid", "na.rm", "type", "width", "space", "smooth"),
                                   
                                   setup_data = function(data, params) {
                                     
                                     purrr::map_dfr(unique(data$PANEL),
                                                    ~{
                                                      
                                                      data <- 
                                                        data %>% 
                                                        #tidyr::drop_na(ymin, ymax) %>%
                                                        dplyr::filter(PANEL == .x) %>%
                                                        dplyr::mutate(dplyr::across(c(x, next_x), ~as.numeric(.), .names = ("n_{.col}")))%>%
                                                          dplyr::mutate(group = 1) %>%
                                                          dplyr::select(-n_next_x, -next_node, -next_x) %>%
                                                        dplyr::mutate(xmin = n_x - params$width/2,
                                                                      xmax = n_x + params$width/2)
                                                      data <- as.data.frame(data)
                                                      #print(tail(data, n=20))
                                                      return(data)
                                                    })
                                     
                                   },
                                   
                                   compute_group = function(data, scales) {
                                     data
                                   }
)

# TEXT LAYER -------
StatSankeyText2 <- ggplot2::ggproto("StatSankeyText", ggplot2::Stat,
                                   extra_params = c("n_grid", "na.rm", "type", "width", "space"),
                                   
                                   setup_data = function(data, params) {
                                     
                                     purrr::map_dfr(unique(data$PANEL),
                                                    ~{
                                                      data <- data %>% dplyr::filter(PANEL == .x) %>%
                                                        dplyr::mutate(dplyr::across(c(x, next_x), ~as.numeric(.), .names = ("n_{.col}")))%>%
                                                        dplyr::mutate(group = 1) %>%
                                                        dplyr::select(-n_next_x, -next_node, -next_x) %>%
                                                        dplyr::mutate(xmin = n_x - params$width/2,
                                                                      xmax = n_x + params$width/2)
                                                      data<- as.data.frame(data)
                                                      print(data)
                                                      return(data)
                                                    })
                                   },
                                   
                                   compute_group = function(data, scales) {
                                     data
                                   }
)


geom_sankey2 <- function(mapping = NULL,
                         data = NULL,
                         position = "identity",
                         na.rm = FALSE,
                         show.legend = NA,
                         space = NULL,
                         type = "sankey",
                         width = .1,
                         smooth = 8,
                         inherit.aes = TRUE,
                         source.width = NULL,
                         ...
) {
  params_list <- prepare_params(...)

  list(
    flow = ggplot2::layer(
      stat = StatSankeyFlow2,
      data = data,
      mapping = mapping,
      geom = "polygon",
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = purrr::flatten(
        list(
          na.rm = na.rm,
          width = width,
          space = space,
          smooth = smooth,
          type = type,
          params_list[[1]]
        )
      )
    ),
    
    node = ggplot2::layer(
      stat = StatSankeyNode2,
      data = data,
      mapping = mapping,
      geom = ggplot2::GeomRect,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = purrr::flatten(
        list(
          na.rm = na.rm,
          width = width,
          space = space,
          smooth = smooth,
          type = type,
          params_list[[2]]
        )
      )
    )
  )
  
  
}

geom_sankey_text2 <- function(mapping = NULL,
                             data = NULL,
                             position = "identity",
                             na.rm = FALSE,
                             show.legend = NA,
                             space = NULL,
                             type = "sankey",
                             width = .1,
                             inherit.aes = TRUE,
                             ...) {
  # Prepare aesthics for label
  label.aes <- list(...)
  
  list(
    label = ggplot2::layer(
      stat = StatSankeyText2,
      data = data,
      mapping = mapping,
      geom = "text",
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = purrr::flatten(
        list(
          na.rm = na.rm,
          width = width,
          space = space,
          type = type,
          label.aes
        )
      )
    )
  )
  
  
}

prepare_params <- function(...) {
  # Prepare aesthics for flow lines
  flow.aes <- list(...)
  removes <- names(flow.aes) %>%
    stringr::str_extract_all(., "(?<=flow.).*") %>% unlist()
  removes2 <- names(flow.aes) %>%
    stringr::str_subset(., "node") %>% unlist()
  flow.aes[c(removes, removes2)] <- NULL
  names(flow.aes) <- names(flow.aes) %>%
    stringr::str_replace_all("flow.", "")
  
  # Prepare aesthics for node boxes
  node.aes <- list(...)
  removes <- names(node.aes) %>%
    stringr::str_extract_all(., "(?<=node.).*") %>% unlist()
  removes2 <- names(node.aes) %>%
    stringr::str_subset(., "flow") %>% unlist()
  node.aes[c(removes, removes2)] <- NULL
  names(node.aes) <- names(node.aes) %>%
    stringr::str_replace_all(., "node.", "")
  
  #flow.aes <- list(alpha = 0.6)
  return(list(flow.aes, node.aes))
}

sigmoid <- function(x_from, x_to, y_from, y_to, smooth = 5, n = 300) {
  x <- seq(-smooth, smooth, length = n)
  y <- exp(x) / (exp(x) + 1)
  out <- data.frame(x = (x + smooth) / (smooth * 2) * (x_to - x_from) + x_from,
                    y = y * (y_to - y_from) + y_from)
}
