library(plsRglm)
library(dplyr)
library(glmmTMB)
library(ggplot2)
library(cowplot)

load("Final_output/PLS-models-9-10-23.RData")

pls_out_1$model$FinalModel$coefficients['tt.1']
pls_out_1$model$CoeffCFull
pls_out_1$model$R2
pls_out_1$model$RSSresidY

pls_out_1$model$InfCrit

pls_out_1$model$FinalModel$y

RSS<-sum((pls_out_1$model$residYChapeau)^2)
TSS<-mean(log())

performance::r2_efron(pls_out_1$model$FinalModel)

pls_out_1$model$FinalModel

pls_out_1$model$tt[,'Comp_1']

names(pls_out_1$model$dataX)

tss.t1 <- 

    with(pls_out_1$model,
         
         sum((
           
           (mean(pls_out_1$model$tt[,'Comp_1']) -
             
              pls_out_1$model$tt[,'Comp_1'])^2))
         
    )


rss.t1.all<-

  with(pls_out_1$model,
       
       sum((
         
         (dataX$Climate.Distance_centered *
         pp['Climate.Distance_centered', 'Comp_1'] *
         FinalModel$coefficients['tt.1'] +
           
           dataX$root4_trade_centered *
           pp['root4_trade_centered', 'Comp_1'] *
           FinalModel$coefficients['tt.1'] +
           
           dataX$NMDS.Distance_centered *
           pp['NMDS.Distance_centered', 'Comp_1'] *
           FinalModel$coefficients['tt.1'] +
           
           dataX$log_numexotics_centered *
           pp['log_numexotics_centered', 'Comp_1'] *
           FinalModel$coefficients['tt.1'] +
           
           dataX$NMDS.Distance_centered * dataX$log_numexotics_centered *
           pp['NMDS.Distance_centered.log_numexotics_centered', 'Comp_1'] *
           FinalModel$coefficients['tt.1']
           ) -
         
           pls_out_1$model$tt[,'Comp_1'])^2)
       
       )

rss.t1.clim<-

  with(pls_out_1$model,
       
       sum((
         
         (#dataX$Climate.Distance_centered *
          #  pp['Climate.Distance_centered', 'Comp_1'] *
          #  FinalModel$coefficients['tt.1'] +
            
            dataX$root4_trade_centered *
            pp['root4_trade_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1'] +
            
            dataX$NMDS.Distance_centered *
            pp['NMDS.Distance_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1'] +
            
            dataX$log_numexotics_centered *
            pp['log_numexotics_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1'] +
            
            dataX$NMDS.Distance_centered * dataX$log_numexotics_centered *
            pp['NMDS.Distance_centered.log_numexotics_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1']
         ) -
           
           pls_out_1$model$tt[,'Comp_1'])^2)
       
  )

rss.t1.trade<-

  with(pls_out_1$model,
       
       sum((
         
         (dataX$Climate.Distance_centered *
             pp['Climate.Distance_centered', 'Comp_1'] *
             FinalModel$coefficients['tt.1'] +
           
           #dataX$root4_trade_centered *
          #   pp['root4_trade_centered', 'Comp_1'] *
          #   FinalModel$coefficients['tt.1'] +
             
             dataX$NMDS.Distance_centered *
             pp['NMDS.Distance_centered', 'Comp_1'] *
             FinalModel$coefficients['tt.1'] +
             
             dataX$log_numexotics_centered *
             pp['log_numexotics_centered', 'Comp_1'] *
             FinalModel$coefficients['tt.1'] +
             
             dataX$NMDS.Distance_centered * dataX$log_numexotics_centered *
             pp['NMDS.Distance_centered.log_numexotics_centered', 'Comp_1'] *
             FinalModel$coefficients['tt.1']
         ) -
           
           pls_out_1$model$tt[,'Comp_1'])^2)
       
  )

rss.t1.nmds<-
  
  with(pls_out_1$model,
       
       sum((
         
         (dataX$Climate.Distance_centered *
            pp['Climate.Distance_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1'] +
            
            dataX$root4_trade_centered *
               pp['root4_trade_centered', 'Comp_1'] *
               FinalModel$coefficients['tt.1'] +
          #######NO EXPERIMENT
          
          #  dataX$NMDS.Distance_centered *
          #  pp['NMDS.Distance_centered', 'Comp_1'] *
          #  FinalModel$coefficients['tt.1'] #+
            
            dataX$log_numexotics_centered *
            pp['log_numexotics_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1'] +
            
            dataX$NMDS.Distance_centered * dataX$log_numexotics_centered *
            pp['NMDS.Distance_centered.log_numexotics_centered', 'Comp_1'] *
           FinalModel$coefficients['tt.1']
         ) -
           
           pls_out_1$model$tt[,'Comp_1'])^2)
       
  )

rss.t1.trees<-
    
  with(pls_out_1$model,
       
       sum((
         
         (dataX$Climate.Distance_centered *
            pp['Climate.Distance_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1'] +
            
            dataX$root4_trade_centered *
            pp['root4_trade_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1'] +
            
            dataX$NMDS.Distance_centered *
            pp['NMDS.Distance_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1'] +
            
           # dataX$log_numexotics_centered *
          #  pp['log_numexotics_centered', 'Comp_1'] *
          #  FinalModel$coefficients['tt.1'] +
            #######NO EXPERIMENT
            
            dataX$NMDS.Distance_centered * dataX$log_numexotics_centered *
            pp['NMDS.Distance_centered.log_numexotics_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1']
         ) -
           
           pls_out_1$model$tt[,'Comp_1'])^2)
       
  )


rss.t1.full.nmds<-
  
  with(pls_out_1$model,
       
       sum((
         
         (dataX$Climate.Distance_centered *
            pp['Climate.Distance_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1'] +
            
            dataX$root4_trade_centered *
            pp['root4_trade_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1'] +
            
            #  dataX$NMDS.Distance_centered *
            #  pp['NMDS.Distance_centered', 'Comp_1'] *
            #  FinalModel$coefficients['tt.1'] #+
            
            dataX$log_numexotics_centered *
            pp['log_numexotics_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1'] #+
            #######EXPERIMENT
            
            #dataX$NMDS.Distance_centered * dataX$log_numexotics_centered *
            #pp['NMDS.Distance_centered.log_numexotics_centered', 'Comp_1'] *
            #FinalModel$coefficients['tt.1']
         ) -
           
           pls_out_1$model$tt[,'Comp_1'])^2)
       
  )

rss.t1.full.trees<-
  
  with(pls_out_1$model,
       
       sum((
         
         (dataX$Climate.Distance_centered *
            pp['Climate.Distance_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1'] +
            
            dataX$root4_trade_centered *
            pp['root4_trade_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1'] +
            
            dataX$NMDS.Distance_centered *
            pp['NMDS.Distance_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1'] #+
            
            # dataX$log_numexotics_centered *
            #  pp['log_numexotics_centered', 'Comp_1'] *
            #  FinalModel$coefficients['tt.1'] +
            ####### EXPERIMENT
            
            #dataX$NMDS.Distance_centered * dataX$log_numexotics_centered *
            #pp['NMDS.Distance_centered.log_numexotics_centered', 'Comp_1'] *
            #FinalModel$coefficients['tt.1']
         ) -
           
           pls_out_1$model$tt[,'Comp_1'])^2)
       
  )


rss.t1.int<-
  
  with(pls_out_1$model,
       
       sum((
         
         (dataX$Climate.Distance_centered *
            pp['Climate.Distance_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1'] +
            
            dataX$root4_trade_centered *
            pp['root4_trade_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1'] +
            
            dataX$NMDS.Distance_centered *
            pp['NMDS.Distance_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1'] +
            
             dataX$log_numexotics_centered *
              pp['log_numexotics_centered', 'Comp_1'] *
              FinalModel$coefficients['tt.1'] #+
            
            #dataX$NMDS.Distance_centered * dataX$log_numexotics_centered *
            #pp['NMDS.Distance_centered.log_numexotics_centered', 'Comp_1'] *
            #FinalModel$coefficients['tt.1']
         ) -
           
           pls_out_1$model$tt[,'Comp_1'])^2)
       
  )

rss.t1.full.int<-
  
  with(pls_out_1$model,
       
       sum((
         
         (dataX$Climate.Distance_centered *
            pp['Climate.Distance_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1'] +
            
            dataX$root4_trade_centered *
            pp['root4_trade_centered', 'Comp_1'] *
            FinalModel$coefficients['tt.1'] #+
            
          #  dataX$NMDS.Distance_centered *
          #  pp['NMDS.Distance_centered', 'Comp_1'] *
          #  FinalModel$coefficients['tt.1'] +
            
           # dataX$log_numexotics_centered *
          #  pp['log_numexotics_centered', 'Comp_1'] *
          #  FinalModel$coefficients['tt.1'] #+
          
          #dataX$NMDS.Distance_centered * dataX$log_numexotics_centered *
          #pp['NMDS.Distance_centered.log_numexotics_centered', 'Comp_1'] *
          #FinalModel$coefficients['tt.1']
         ) -
           
           pls_out_1$model$tt[,'Comp_1'])^2)
       
  )


tss.t2 <- 
  
  with(pls_out_1$model,
       
       sum((
         
         (mean(pls_out_1$model$tt[,'Comp_2']) -
            
            pls_out_1$model$tt[,'Comp_2'])^2))
       
  )


rss.t2.all<-
  
  with(pls_out_1$model,
       
       sum((
         
         (dataX$Climate.Distance_centered *
            pp['Climate.Distance_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] +
            
            dataX$root4_trade_centered *
            pp['root4_trade_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] +
            
            dataX$NMDS.Distance_centered *
            pp['NMDS.Distance_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] +
            
            dataX$log_numexotics_centered *
            pp['log_numexotics_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] +
            
            dataX$NMDS.Distance_centered * dataX$log_numexotics_centered *
            pp['NMDS.Distance_centered.log_numexotics_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2']
         ) -
           
           pls_out_1$model$tt[,'Comp_2'])^2)
       
  )

rss.t2.clim<-
  
  with(pls_out_1$model,
       
       sum((
         
         (#dataX$Climate.Distance_centered *
           #  pp['Climate.Distance_centered', 'Comp_2'] *
           #  FinalModel$coefficients['tt.2'] +
           
           dataX$root4_trade_centered *
             pp['root4_trade_centered', 'Comp_2'] *
             FinalModel$coefficients['tt.2'] +
             
             dataX$NMDS.Distance_centered *
             pp['NMDS.Distance_centered', 'Comp_2'] *
             FinalModel$coefficients['tt.2'] +
             
             dataX$log_numexotics_centered *
             pp['log_numexotics_centered', 'Comp_2'] *
             FinalModel$coefficients['tt.2'] +
             
             dataX$NMDS.Distance_centered * dataX$log_numexotics_centered *
             pp['NMDS.Distance_centered.log_numexotics_centered', 'Comp_2'] *
             FinalModel$coefficients['tt.2']
         ) -
           
           pls_out_1$model$tt[,'Comp_2'])^2)
       
  )

rss.t2.trade<-
  
  with(pls_out_1$model,
       
       sum((
         
         (dataX$Climate.Distance_centered *
            pp['Climate.Distance_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] +
            
            #dataX$root4_trade_centered *
            #   pp['root4_trade_centered', 'Comp_2'] *
            #   FinalModel$coefficients['tt.2'] +
            
            dataX$NMDS.Distance_centered *
            pp['NMDS.Distance_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] +
            
            dataX$log_numexotics_centered *
            pp['log_numexotics_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] +
            
            dataX$NMDS.Distance_centered * dataX$log_numexotics_centered *
            pp['NMDS.Distance_centered.log_numexotics_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2']
         ) -
           
           pls_out_1$model$tt[,'Comp_2'])^2)
       
  )

rss.t2.nmds<-
  
  with(pls_out_1$model,
       
       sum((
         
         (dataX$Climate.Distance_centered *
            pp['Climate.Distance_centered', 'Comp_2'] *
             FinalModel$coefficients['tt.2'] +
           
             dataX$root4_trade_centered *
                pp['root4_trade_centered', 'Comp_2'] *
                FinalModel$coefficients['tt.2'] +
           
           #dataX$NMDS.Distance_centered *
          #   pp['NMDS.Distance_centered', 'Comp_2'] *
          #   FinalModel$coefficients['tt.2'] +
           
             dataX$log_numexotics_centered *
             pp['log_numexotics_centered', 'Comp_2'] *
             FinalModel$coefficients['tt.2'] +
           
           ####### NO EXPERIMENT
             dataX$NMDS.Distance_centered * dataX$log_numexotics_centered *
             pp['NMDS.Distance_centered.log_numexotics_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2']
         ) -
           
           pls_out_1$model$tt[,'Comp_2'])^2)
       
  )

rss.t2.full.nmds<-
  
  with(pls_out_1$model,
       
       sum((
         
         (dataX$Climate.Distance_centered *
            pp['Climate.Distance_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] +
            
            dataX$root4_trade_centered *
            pp['root4_trade_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] +
            
            #dataX$NMDS.Distance_centered *
            #   pp['NMDS.Distance_centered', 'Comp_2'] *
            #   FinalModel$coefficients['tt.2'] +
            
            dataX$log_numexotics_centered *
            pp['log_numexotics_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] #+
            
            #######  EXPERIMENT
            #dataX$NMDS.Distance_centered * dataX$log_numexotics_centered *
            #pp['NMDS.Distance_centered.log_numexotics_centered', 'Comp_2'] *
            #FinalModel$coefficients['tt.2']
         ) -
           
           pls_out_1$model$tt[,'Comp_2'])^2)
       
  )

rss.t2.trees<-
  
  with(pls_out_1$model,
       
       sum((
         
         (dataX$Climate.Distance_centered *
            pp['Climate.Distance_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] +
            
            dataX$root4_trade_centered *
            pp['root4_trade_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] +
            
            dataX$NMDS.Distance_centered *
            pp['NMDS.Distance_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] +
            
            # dataX$log_numexotics_centered *
            #  pp['log_numexotics_centered', 'Comp_2'] *
            #  FinalModel$coefficients['tt.2'] +
          ####### NO EXPERIMENT
          
            dataX$NMDS.Distance_centered * dataX$log_numexotics_centered *
            pp['NMDS.Distance_centered.log_numexotics_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2']
         ) -
           
           pls_out_1$model$tt[,'Comp_2'])^2)
       
  )

rss.t2.full.trees<-
  
  with(pls_out_1$model,
       
       sum((
         
         (dataX$Climate.Distance_centered *
            pp['Climate.Distance_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] +
            
            dataX$root4_trade_centered *
            pp['root4_trade_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] +
            
            dataX$NMDS.Distance_centered *
            pp['NMDS.Distance_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2']# +
            
            # dataX$log_numexotics_centered *
            #  pp['log_numexotics_centered', 'Comp_2'] *
            #  FinalModel$coefficients['tt.2'] +
            ####### EXPERIMENT
            
            #dataX$NMDS.Distance_centered * dataX$log_numexotics_centered *
            #pp['NMDS.Distance_centered.log_numexotics_centered', 'Comp_2'] *
            #FinalModel$coefficients['tt.2']
         ) -
           
           pls_out_1$model$tt[,'Comp_2'])^2)
       
  )

rss.t2.int<-
  
  with(pls_out_1$model,
       
       sum((
         
         (dataX$Climate.Distance_centered *
            pp['Climate.Distance_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] +
            
            dataX$root4_trade_centered *
            pp['root4_trade_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] +
            
            dataX$NMDS.Distance_centered *
            pp['NMDS.Distance_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] +
            
            dataX$log_numexotics_centered *
            pp['log_numexotics_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] #+
          
          #dataX$NMDS.Distance_centered * dataX$log_numexotics_centered *
          #pp['NMDS.Distance_centered.log_numexotics_centered', 'Comp_2'] *
          #FinalModel$coefficients['tt.2']
         ) -
           
           pls_out_1$model$tt[,'Comp_2'])^2)
       
  )

rss.t2.full.int<-
  
  with(pls_out_1$model,
       
       sum((
         
         (dataX$Climate.Distance_centered *
            pp['Climate.Distance_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] +
            
            dataX$root4_trade_centered *
            pp['root4_trade_centered', 'Comp_2'] *
            FinalModel$coefficients['tt.2'] #+
          
          #  dataX$NMDS.Distance_centered *
          #  pp['NMDS.Distance_centered', 'Comp_2'] *
          #  FinalModel$coefficients['tt.2'] +
          
          # dataX$log_numexotics_centered *
          #  pp['log_numexotics_centered', 'Comp_2'] *
          #  FinalModel$coefficients['tt.2'] #+
          
          #dataX$NMDS.Distance_centered * dataX$log_numexotics_centered *
          #pp['NMDS.Distance_centered.log_numexotics_centered', 'Comp_2'] *
          #FinalModel$coefficients['tt.2']
         ) -
           
           pls_out_1$model$tt[,'Comp_2'])^2)
       
  )


r2.t1.all <- 1 - rss.t1.all/tss.t1

r2.t1.clim <-1 - rss.t1.clim/tss.t1
r2.t1.trade <-1 - rss.t1.trade/tss.t1
r2.t1.nmds <-1 - rss.t1.nmds/tss.t1
r2.t1.trees <-1 - rss.t1.trees/tss.t1
r2.t1.int <-1 - rss.t1.int/tss.t1
r2.t1.full.int <-1 - rss.t1.full.int/tss.t1
r2.t1.full.nmds <-1 - rss.t1.full.nmds/tss.t1
r2.t1.full.trees <-1 - rss.t1.full.trees/tss.t1

r2.t2.all <- 1 - rss.t2.all/tss.t2

r2.t2.clim <-1 - rss.t2.clim/tss.t2
r2.t2.trade <-1 - rss.t2.trade/tss.t2
r2.t2.nmds <-1 - rss.t2.nmds/tss.t2
r2.t2.trees <-1 - rss.t2.trees/tss.t2
r2.t2.int <-1 - rss.t2.int/tss.t2
r2.t2.full.int <-1 - rss.t2.full.int/tss.t2
r2.t2.full.nmds <-1 - rss.t2.full.nmds/tss.t2
r2.t2.full.trees <-1 - rss.t2.full.trees/tss.t2

r2.pls <- data.frame(
  
  t1 = c(
    clim = r2.t1.all - r2.t1.clim,
    trade = r2.t1.all - r2.t1.trade,
    nmds = r2.t1.all - r2.t1.nmds,
    trees = r2.t1.all - r2.t1.trees,
    full.nmds = r2.t1.all - r2.t1.full.nmds,
    full.trees = r2.t1.all - r2.t1.full.trees,
    int = r2.t1.all - r2.t1.int,
    full.int = r2.t1.all - r2.t1.full.int,
    fixed      = r2.t1.all ),
  
  t2 = c(
    clim = r2.t2.all - r2.t2.clim,
    trade = r2.t2.all - r2.t2.trade,
    nmds = r2.t2.all - r2.t2.nmds,
    trees = r2.t2.all - r2.t2.trees,
    full.nmds = r2.t2.all - r2.t2.full.nmds,
    full.trees = r2.t2.all - r2.t2.full.trees,
    int = r2.t2.all - r2.t2.int,
    full.int = r2.t2.all - r2.t2.full.int,
    fixed = r2.t2.all)
  
  )

#t(r2.pls)

#with(t(r2.pls) %>% as.data.frame, full.int-full.trees)

r2.pls <-
  r2.pls%>% 
  mutate( total.r2.t1 = 
            
            r2.pls['fixed','t1']
            
#            sum(
#              r2.t1.all - r2.t1.clim,
#              r2.t1.all - r2.t1.trade,
#              #r2.t1.all - r2.t1.nmds,
#              #r2.t1.all - r2.t1.trees,
#              #r2.t1.all - r2.t1.int,
#              r2.t1.all - r2.t1.full.int)
            ) %>%
  
  mutate( total.r2.t2 = 
            
            r2.pls['fixed','t2']
          
#            sum(
#              r2.t2.all - r2.t2.clim,
#              r2.t2.all - r2.t2.trade,
#              #r2.t2.all - r2.t2.nmds,
#              #r2.t2.all - r2.t2.trees,
#              #r2.t2.all - r2.t2.int,
#              r2.t2.all - r2.t2.full.int)
            )


# explained by first axis

r2.pls <-
  r2.pls%>%
  
  # how much is the total explained ?
  mutate( axis.total.t1 = .5819) %>%
  # explained by regression
  
  mutate( rel.r2.t1 = t1 * axis.total.t1)

# explained by two axes

r2.pls <-
  r2.pls%>%
  
  # how much is the total explained ?
  mutate( axis.total.both = .6392) %>%
  # explained by regression
  
  mutate( rel.r2.t1 = t1 * axis.total.t1)

###################

r2.pls <-
  r2.pls%>%
  
  # how much is the total explained ?
  mutate( axis.total.t1 = .5819) %>%
  #mutate( axis.total.t2 = .6392) %>%
  
  mutate( axis.total.t2 = .6392-.5819) %>%
  
  # explained by regression
  
  mutate( rel.r2.t1 = t1 * axis.total.t1) %>%
  mutate( rel.r2.t2 = t2 * axis.total.t2) %>%
  
  # total explained
  mutate(total = rel.r2.t1 + rel.r2.t2)
  

r2.pls['continent',] <- NA
r2.pls['total',] <- NA

r2.pls['total',c('rel.r2.t1','rel.r2.t2','total')] <-
  c(.5819,
    .6392-.5819,
    .6392)

r2.pls[c('full.int',
         'int',
         'nmds',
         'trees',
         'trade',
         'clim',
         'continent',
         'total',
         'fixed'),] %>% 
  write.csv('Final_output/final_pls_R2.csv', row.names=T)

#############

#load('Final_output/zinf_model_final.RData')
load("Final_output/linearmodels.RData")

#anova(t1.factor)

#(1 - sum((predict(t1.factor) - to_analyze$logNumber.of.pests)^2)/sum((mean(to_analyze$Number.of.pests) - to_analyze$Number.of.pests)^2)
#)-
#(1 - sum((predict(t1.factor, newdata=
##          to_analyze %>%
#          mutate(
#            log_numexotics_centered=scale(0, T, F),
#            NMDS.Distance_centered=scale(0, T, F)
#          )) - to_analyze$logNumber.of.pests)^2)/sum((mean(to_analyze$Number.of.pests) - to_analyze$Number.of.pests)^2)
#)

#t1.factor$coefficients['NMDS.Distance_centered:log_numexotics_centered']

zinf_pred <- predict(model_final_zinf_centered, type='response')
ancova_pred <- predict(t1.factor) %>% exp
pls_pred <- predict(pls_out_1$model, type='response')

zinf_pred_pearson <- cor(to_analyze$Number.of.pests, zinf_pred) %>% round(2)
ancova_pred_pearson <- cor(to_analyze$Number.of.pests, ancova_pred)%>% round(2)
pls_pred_pearson <- cor(to_analyze$Number.of.pests, pls_pred)%>% round(2)

zinf_pred_R2 <- 
  (1 - sum((zinf_pred - to_analyze$Number.of.pests)^2)/
     sum((mean(to_analyze$Number.of.pests) - to_analyze$Number.of.pests)^2) )%>% round(2)
ancova_pred_R2 <-
  (1 - sum((ancova_pred - to_analyze$Number.of.pests)^2)/
     sum((mean(to_analyze$Number.of.pests) - to_analyze$Number.of.pests)^2))%>% round(2)
pls_pred_R2 <-
  (1 - sum((pls_pred - to_analyze$Number.of.pests)^2)/
     sum((mean(to_analyze$Number.of.pests) - to_analyze$Number.of.pests)^2))%>% round(2)

to_plot <-
  data.frame(
    Model = "Zero-inflated Poisson GLM",
    r = zinf_pred_pearson,
    R2 = zinf_pred_R2,
    Observed = to_analyze$Number.of.pests,
    Predicted= zinf_pred
  ) %>% rbind(data.frame(
    Model = "Linear Model on Log1p",
    r=ancova_pred_pearson,
    R2 = ancova_pred_R2,
    Observed = to_analyze$Number.of.pests,
    Predicted= ancova_pred
  )) %>% rbind(data.frame(
    Model = "Partial-Least-Sq Poisson GLM",
    r=pls_pred_pearson,
    R2 = pls_pred_R2,
    Observed = to_analyze$Number.of.pests,
    Predicted= pls_pred
  ))

to_plot$Model <- factor(to_plot$Model,
                        levels = levels(as.factor(to_plot$Model))[c(3,1,2)])

windows(8,4); ggplot(data = to_plot) +
  
  geom_point(aes(x = Observed, y = Predicted)) +
  
  geom_text(aes(x = 1, y = 32, label = paste("r =", r, "\nR² =", R2), hjust=0, vjust=1)) +
  
  geom_abline() +
  
  scale_x_continuous(transform = 'log', breaks = c(0,1,2,4,8,16,32)) +
  scale_y_continuous(transform = 'log', breaks = c(.125,.25,.5,0,1,2,4,8,16,32)) +
  
  facet_grid(~ Model) +   theme_cowplot()+theme(strip.background=element_blank(), strip.placement="outside",
                                                          panel.border=element_rect(fill=NA, colour="black"))
#####################################

windows(8,4);to_plot %>%
  mutate(Residuals = Predicted - Observed) %>%
  ggplot() + geom_point(aes(x = Observed, y = Residuals)) + facet_grid(~ Model) +
  scale_x_continuous(transform = 'log', breaks = c(0,1,2,4,8,16,32))


############
