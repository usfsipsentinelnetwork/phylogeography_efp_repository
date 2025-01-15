library("dplyr")
library("tidyr")
library("reshape2")

country_names_to_standardize_table <-
  data.frame(
    matrix(
      data= c(
        "Burma", "Myanmar",
        "Bahamas", "Bahamas",
        "Congo", "Congo",
        "Egypt", "Egypt",
        "Iran",  "Iran",
        "Korea, North|Korea, Dem\\. Rep\\.", "North Korea",
        "Korea, Rep\\.|Korea, South", "South Korea",
        "Kyrgyz", 'Kyrgyzstan',
        "Serbia|Yugoslavia", "Serbia",
        "Sudan|Sahara", "Sudan",
        "Syria","Syrian Arab Republic",
        "Lao","Loas",
        "Hong Kong","Hong Kong",
        "United States","United States of America",
        "Soviet Union|Russia","Russian Federation",
        "Yemen","Yemen",
        "Ethiopia","Ethiopia",
        "Gambia","Gambia",
        "Taiwan","Taiwan",
        "Czech","Czechia",
        "Slovak Republic","Slovakia"
      ),
      ncol=2,
      byrow=T,
      dimnames=list(
        NULL,
        c("searchterm","replacement")
      )
    )
  )

correct_name_usatrade_census <- function (x, name="Hawaii", country_names_to_standardize = country_names_to_standardize_table) {
  x <- x %>% select(colnames(x)[grepl("Country|[0-9]{4}",colnames(x))]) %>% rename("Partner Name"=Country)
  x <- cbind("Reporter Name"=name, x)
  # hawaii not in thousands
  x[,-c(1,2)] <- as.numeric(as.matrix(x[,-c(1,2)]))/1000
  
  melt_format_wits(x)
}

melt_format_wits <- function (x, country_names_to_standardize = country_names_to_standardize_table) {
  x <- x %>% select(colnames(x)[grepl("Reporter Name|Partner Name|[0-9]{4}",colnames(x))]) %>% replace(is.na(.), 0)
  for (i in 1:dim(country_names_to_standardize)[1])
  {
    p <- country_names_to_standardize$searchterm[i]
    r <- country_names_to_standardize$replacement[i]
    j <- which(grepl(pattern=p, x$`Partner Name`))
    
    print(paste("Searching for",p))
    
    if (length(j) > 0) {
      
   #   print(paste("Replacing", p, "with", r))
      
      if(length(j) ==1) {
       # print("One row")
        x$`Partner Name`[j] <- r
      } else {
    #    print("Two rows")
        #print(country_names_to_standardize$searchterm[i])
        #print(y[,-c(1:2)])
        #print(colSums(y[,-c(1:2)]))
        y <- x[j,]
        x[dim(x)[1]+1,1:2] <- c("Reporter Name"=x[1,1],
                     "Partner Name"=r)
        x[dim(x)[1],3:dim(x)[2]]<-
                     colSums(y[,-c(1:2)])
        x <- x[-j,]
      }
    }
  }
  x <- x %>% melt(value.name="Dollars") %>% rename(Year=variable)
  x$Year <- as.integer(as.character(x$Year))
  x
}

setwd("Covariates_input/New_Trade_Data")

eu_not_eu_countries <- c(
  "Cyprus",
  "Greenland",
  "Kazakhstan",
  "Kyrgzystan",
  "Tajikistan",
  "Turkmenistan",
  "Uzbekistan"
)

eu_trade_raw <- read.csv("Europe_CentralAsia_WITS-Partner-Timeseries.csv", check.names=F)%>% melt_format_wits
na_trade_raw <- read.csv("NorthAmerica-WITS-Partner-Timeseries.csv", check.names=F)%>% melt_format_wits

#melt_format_wits(read.csv("Europe_CentralAsia_WITS-Partner-Timeseries.csv", check.names=F))

eu_not_eu <- NULL
for (f in grep(paste(eu_not_eu_countries,collapse="|"),
               grep("\\.csv",dir(), value=T), value=T)) {
  eu_not_eu <- bind_rows(eu_not_eu, read.csv(f, check.names=F))
}
eu_not_eu <- eu_not_eu %>% melt_format_wits

na_not_na <- read.csv("Bermuda-WITS-Partner-Timeseries.csv", check.names=F)%>%
  melt_format_wits

eu_trade_corrected<-
  eu_not_eu %>% 
    group_by(Year, `Partner Name`) %>%
    summarise(total_to_subtract = sum(Dollars)) %>% as.data.frame() %>%
    left_join(eu_trade_raw,., by = c("Year","Partner Name")) %>% replace(is.na(.),0) %>%
    mutate(region_only=Dollars-total_to_subtract)

na_trade_corrected<-
  left_join(na_trade_raw,
            na_not_na %>%
              select("Partner Name","Year","Dollars") %>%
              rename(Dollars_B=Dollars), by=c("Year","Partner Name")) %>%
  replace(is.na(.),0) %>%mutate(region_only=Dollars-Dollars_B)

au_trade <- read.csv("Australia-WITS-Partner-Timeseries.csv", check.names=F)%>% melt_format_wits %>%
  rename(region_only=Dollars)

conversion<- read.csv("Dollars.2023.csv")

PPP<-
  bind_rows(na_trade_corrected,eu_trade_corrected, au_trade) %>%
    select("Reporter Name","Partner Name",Year,region_only) %>%
    left_join(conversion) %>% mutate(Dollars=region_only*Dollars.2023)%>%
  select(-c(region_only, Dollars.2023))

hi_trade_ppp<- read.csv("Hawaii-89-21-PPP.csv", check.names=F) %>%correct_name_usatrade_census #%>%
  #select("Reporter Name","Partner Name",Year,region_only,Dollars)

setdiff(hi_trade_ppp$`Partner Name`, PPP$`Partner Name`) %>% sort
setdiff(PPP$`Partner Name`, hi_trade_ppp$`Partner Name`) %>% sort

#notfilledin.all.wide <- bind_rows(PPP,hi_trade_ppp) %>%
#  rename(country_name="Partner Name")%>%
#  dcast(country_name + `Reporter Name` ~ Year, value.var = "Dollars")

load('../../Host_list/Output/ordinated_unifrac_data.Rdata')

bioregions_vector <- bioregions_vector%>%
  filter(!(
    biogeographic_region %in%
      c("N.North.America",
        "W.North.America",
        "NE.North.America",
        "SE.North.America"))) %>%
  rbind(data.frame(country_name=c("United States of America","Canada"), biogeographic_region="North America"))

trade_aggregated <- bind_rows(PPP,hi_trade_ppp) %>%
  rename(country_name="Partner Name") %>% 
  right_join(
    bioregions_vector ) %>%
  dplyr::select(-country_name) %>%
  aggregate(Dollars~biogeographic_region+Year+`Reporter Name`, ., sum)

trade_aggregated %>%
  write.csv("../trade_aggregated_PPP_long.csv",row.names=F)

trade_aggregated_wide <-
  trade_aggregated %>%
    dcast(biogeographic_region + `Reporter Name` ~ Year, value.var = "Dollars")#%>%
    #replace(is.na(.),0)

trade_aggregated_wide$`Reporter Name` <-
  replace(trade_aggregated_wide$`Reporter Name`, trade_aggregated_wide$`Reporter Name`=="Europe & Central Asia", "Europe")

trade_aggregated_wide$sum_1989_to_2007T <- trade_aggregated_wide %>%
  select(colnames(trade_aggregated_wide)[-c(1:2)] %>%
           (function(x) x[which(as.integer(x)<2008)])) %>% 
  (function(x) rowSums(x)/1000000000)
    
trade_aggregated_wide$sum_2008_to_2021T <- trade_aggregated_wide %>%
  select(colnames(trade_aggregated_wide)[-c(1:2)] %>%
           (function(x) x[which(as.integer(x)>2007)])) %>% 
  (function(x) rowSums(x)/1000000000)

trade_aggregated_wide$sum_1989_to_2007T <-
  replace(trade_aggregated_wide$sum_1989_to_2007T,
          trade_aggregated_wide$`Reporter Name`=='Hawaii',
          
          with(trade_aggregated_wide,
               sum_2008_to_2021T[which(`Reporter Name`=='Hawaii')]*
                 (sum_1989_to_2007T[which(`Reporter Name`=='North America')])/
                 (sum_2008_to_2021T[which(`Reporter Name`=='North America')]))
          
  )

trade_aggregated_wide <-
  trade_aggregated_wide %>% dplyr::mutate(cumulativeT = sum_1989_to_2007T + sum_2008_to_2021T)

trade_aggregated_wide%>%
  write.csv("../trade_aggregated_PPP_wide.csv",row.names=F)

#library(nlme)
#library(lme4)

#names(trade_aggregated)[3] <- "sink_region"

#model.trade <- nlme(Dollars ~ (A + biogeographic_region + sink_region)*exp(Year),
#                    data=trade_aggregated,
#                    fixed = A + Year ~ 1,
#                    random = biogeographic_region * sink_region ~ 1,
#                    start= c(A=1))


#model.trade <- lm(log(Dollars+.001) ~ biogeographic_region * sink_region + Year, data=trade_aggregated)

#model.gam <- filter(trade_aggregated, biogeographic_region=="SE.Asia" & sink_region=="Europe & Central Asia") %>%
#  gam(Dollars ~ s(Year, bs="cr"), data=., family=Gamma(link="inverse"))

#summary(model.gam)
#plot(model.gam)

#predicted.trade.values <-
#  expand.grid(Year=1800:2021,
#              biogeographic_region=unique(trade_aggregated$biogeographic_region),
#              sink_region=unique(trade_aggregated$sink_region)) %>%
#    (function(x) cbind(x, Dollars=exp(predict(model.trade, x))))

#predicted.trade.values.gam <-
#  expand.grid(Year=1800:2021#,
#              #biogeographic_region=unique(trade_aggregated$biogeographic_region),
#              #sink_region=unique(trade_aggregated$sink_region)
#              ) %>%
#  (function(x) cbind(x, Dollars=predict(model.gam, x)))

#predicted.trade.values_wide <- predicted.trade.values %>%
#  dcast(biogeographic_region + sink_region ~ Year, value.var = "Dollars") %>%
#  replace(is.na(.), 0)

#predicted.trade.values$`2021`
#trade_aggregated_wide$`2021`

#ggplot() +
#  geom_line(aes(x=Year, y=Dollars, color=biogeographic_region, lty=sink_region), data=predicted.trade.values) +
#  geom_point(aes(x=Year, y=Dollars, color=biogeographic_region, pch=sink_region), data=trade_aggregated) +
#  xlim(1980,2020) + ylim(0,1000000000) #+ scale_y_continuous(trans="log")



#summary(model.trade)
