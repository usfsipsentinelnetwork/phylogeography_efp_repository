# Packages
#library(plyr)
library(dplyr)
#library(httr)
library(jsonlite)
#library(tidyr)
library(rgbif)
#library(sp)
#library(maps)
#library(foreach)

#############
# Functions #
#############

get_rgbif_occurrences <- function(pestname, translation, summarize_region_codes=NULL, trees.all, trees.focal.region) {
  
  x<-NA
  while(class(x)=="logical") x<-tryCatch(key_table_GBIF_single(x=pestname), error=function(cond) return(NA), warning=function(cond) return(NA))
  
  if(is.null(x)) {
    #skipped <- c(skipped, pestname)
    print(paste("No match for", pestname))
    return("No match")
  }
  if(x[2] %in% c("HIGHER RANK", "NO MATCH")) {
    #skipped <- c(skipped, pests[i])
    print(paste("No match for", pestname))
    return("No match")
  }
  
  pest.x <- x[1]
  key.x <-  x[2]
  
  if(!(
    occ_count(
    taxonKey = key.x,
    occurrenceStatus="PRESENT",
    hasGeospatialIssue=FALSE) > 0
  )) {
    print(paste("No records for", pestname))
    return("No records")
  }
  
  print(paste("Prepping", pestname))
  
  p1<-occ_download_prep(
    #pred_default(),
    pred('hasGeospatialIssue', "false"),
    pred('occurrenceStatus', "PRESENT"),
    pred_not(pred('datasetKey', "da3542b4-9a73-4054-b9a3-2d762e172199")),
    pred('taxonKey', key.x))
  
  #if (length(p1)==0) next
  
  q1<-NA
  #number_of_records<- NA
  #while(is.na(number_of_records)) number_of_records <- tryCatch(occ_download_meta(q1)$totalRecords, error=function(cond) return(NA), warning=function(cond) return(NA))
  while(class(q1)=="logical" & length(q1)>0) q1<-tryCatch(occ_download_queue(.list=list(p=p1)), error=function(cond) return(NA), warning=function(cond) return(NA))
  
  #print(q1)
  
    z<-NA
    
    print(paste("Downloading", pestname))
    
    while(class(z)=="logical") z<-tryCatch(occ_download_get(q1,path="Output/chronologies-GBIF/occ_downloads"), error=function(cond) return(NA), warning=function(cond) return(NA))
    
    print(paste("Importing", pestname))
    
    zzz <- tryCatch(occ_download_import(z), error=function(cond) return(list(NA,cond)))#, warning=function(cond) return(NA))
    if ("list" %in% class(zzz)) { if (is.na(zzz[[1]])) {
      skipped <- c(skipped, pestname)
      print(paste("Could not download ", pestname))
      print(zzz[[2]])
      return("Download error")
    }}
    
    hosts.df <- gethosts(zzz, translation=translation, return_data_frame=T)%>%
      select(species, associatedTaxa_translated, year,
             countryCode, level0Name, stateProvince, level1Name, locality,
             decimalLatitude, decimalLongitude, level1Gid, level2Gid, level2Name,	level3Gid, level3Name,
             eventDate, verbatimEventDate, associatedTaxa_verbatim, scientificName.x, verbatimScientificName,
             institutionCode, catalogNumber, associatedSequences) %>%
      arrange(year, eventDate, countryCode, level1Name, stateProvince, level2Name)
  
    if(!is.null(summarize_region_codes)) {
      number_years_focal_region_after_elsewhere<- NA
      number_records_pre<-NA
      first_year_focal_region<-NA
      first_year_elsewhere <-NA
      records_focal_region<-NA
      records_elsewhere<-NA
      records_focal_region_trees<-NA
      records_elsewhere_trees<-NA
      records_focal_region_focal_region_trees<-NA
      records_elsewhere_focal_region_trees<-NA
      
      records_focal_region <- dim(hosts.df %>% filter(countryCode %in% summarize_region_codes))[1]
      records_elsewhere <- dim(hosts.df %>% filter(!(countryCode %in% c(NA, ""))) %>% filter(!(countryCode %in% summarize_region_codes)))[1]
      
      print(paste("Getting hosts for", pestname))
      
      hosts1 <- (hosts.df %>% filter(countryCode %in% summarize_region_codes))$associatedTaxa_translated
      hosts2 <- (hosts.df %>% filter(!(countryCode %in% c(NA, ""))) %>% filter(!(countryCode %in% summarize_region_codes)))$associatedTaxa_translated
      
      print(paste("Counting records for", pestname))
      
      records_focal_region_trees<- sum(hosts1 %in% trees.all)
      records_elsewhere_trees<- sum(hosts2 %in% trees.all)
      records_focal_region_focal_region_trees<- sum(hosts1 %in% trees.focal.region)
      records_elsewhere_focal_region_trees<- sum(hosts2 %in% trees.focal.region)
      
      print(paste("Getting years for", pestname))
      
      first_year_focal_region <- min(na.omit((hosts.df %>% filter(countryCode %in% summarize_region_codes))$year))
      first_year_elsewhere <- min(na.omit((hosts.df %>% filter(!(countryCode %in% c(NA, ""))) %>% filter(!(countryCode %in% summarize_region_codes)))$year))
      
      number_years_focal_region_after_elsewhere <- first_year_elsewhere - first_year_focal_region
      number_records_pre <- dim((hosts.df %>% filter(!(countryCode %in% c(NA, ""))) %>% filter(!(countryCode %in% summarize_region_codes))) %>% filter(year < first_year_focal_region))[1]
      
      pest_search <- 
                           data.frame(
                             pest=pest.x,
                             number_years_focal_region_after_elsewhere=number_years_focal_region_after_elsewhere,
                             number_records_pre=number_records_pre,
                             first_year_focal_region=first_year_focal_region,
                             first_year_elsewhere = first_year_elsewhere,
                             records_focal_region=records_focal_region,
                             records_elsewhere=records_elsewhere,
                             records_focal_region_trees=records_focal_region_trees,
                             records_elsewhere_trees=records_elsewhere_trees,
                             records_focal_region_focal_region_trees,
                             records_elsewhere_focal_region_trees)
      print("Returning data")
      return(list(hosts.df, pest_search))
    } else {
      return(hosts.df)
    }
}

gethosts <- function(x, full_names_only = T, translation=NULL, epithet_only=TRUE, return_data_frame=FALSE) {
  
  # need to adapt to no associatedTaxa in names
  if ("associatedTaxa" %in% names(x)) {
    
#    if ("associatedOrganisms" %in% names(x)) {
#      x<-filter(x, !(associatedOrganisms %in% c("", NA, "NA", NULL, "NULL")) | !(associatedTaxa %in% c("", NA, "NA", NULL, "NULL"))) %>%
#        filter(!(references %in% c("", NA)) | !(catalogNumber %in% c("", NA)) | !(occurrenceID %in% c("", NA))) %>%
#        filter(!grepl("^vector", associatedOrganisms))
#    } else {
#      x<-filter(x,  !(associatedTaxa %in% c("", NA, "NA", NULL, "NULL"))) %>%
#        filter(!(references %in% c("", NA)) | !(catalogNumber %in% c("", NA)) | !(occurrenceID %in% c("", NA)))
#    }
    
    x<- x %>%
      mutate(associatedTaxa, associatedTaxa=sub("^(infests\\: )(.*$)", "\\2", associatedTaxa)) %>%
      mutate(associatedTaxa, associatedTaxa=sub("^(host\\: )(.*$)", "\\2", associatedTaxa)) %>%
      mutate(associatedTaxa, associatedTaxa=sub("^(on )(.*$)", "\\2", associatedTaxa)) %>%
      mutate(associatedTaxa, associatedTaxa=sub("^(Host of\\: )(.*$)", "\\2", associatedTaxa)) %>%
      mutate(associatedTaxa, associatedTaxa=sub("^(Host\\: )(.*) \\|.*", "\\2", associatedTaxa)) %>%
      mutate(associatedTaxa, associatedTaxa=sub("^(Host\\: )(.*)", "\\2", associatedTaxa))
    
    if ("associatedOrganisms" %in% names(x)) {
      x<- #filter(x, grepl("^[A-Z]", associatedTaxa) | !(associatedOrganisms %in% c("", NA))) %>%
        x%>%mutate(associatedOrganisms, associatedOrganisms=sub("(^.*) spec\\.$", "\\1 sp.", associatedOrganisms))
    } #else {
      #x<- filter(x, grepl("^[A-Z]", associatedTaxa))
    #}
    
    x<-
      mutate(x, associatedTaxa, associatedTaxa=sub("^(.*)(\\([A-Z][a-z]* [a-z]\\))", "\\2", associatedTaxa)) %>%
      mutate(associatedTaxa, associatedTaxa=sub("((^[A-Z][a-z]* [a-z]*)|(^[A-Z][a-z]* [a-z]* (subsp|cv|var)\\. [a-z]*)) [A-Za-z0-9\\.\\(\\) ]*", "\\1", associatedTaxa))
    
    if ("associatedOrganisms" %in% names(x))
      x<- mutate(x, associatedTaxa = replace(associatedTaxa, which(associatedTaxa %in% c("", NA)), associatedOrganisms[which(associatedTaxa %in% c("", NA))]))
    
    x <- mutate(x, associatedTaxa, associatedTaxa=sub("^([A-Z][a-z]+ [a-z]+)[\\|\\,\\&].*", "\\1", associatedTaxa)) %>%
      mutate(associatedTaxa, associatedTaxa=sub("^([A-Z][a-z]+ [a-z]+)\\.", "\\1", associatedTaxa))
    
    if (full_names_only)
      x <- tidyr::separate(x, associatedTaxa, c("hostgenus", "hostspecies"), sep=" ", remove=F, extra="merge") #%>%
#      filter(!(hostspecies %in% c("sp", "sp."))) %>%
#      filter(!(is.na(hostspecies))) %>%
#      filter(!grepl("^sp\\.", hostspecies)) %>%
#      filter(!grepl("^[A-Z]([a-z]*)\\.", hostspecies))
    
    if (epithet_only)
      x <- tidyr::separate(x, hostspecies, c("sp_epithet", "subtaxa"), sep=" ", remove=F, extra="merge") %>%
      mutate(associatedTaxa = paste(hostgenus, sp_epithet))
    
    if(!is.null(translation)) {
      x<-left_join(
        x,
        select(translation, spec.name, scientificName),
        by=join_by("associatedTaxa"=="spec.name")
      ) %>%
        dplyr::rename(associatedTaxa_translated = scientificName.y) %>%
        dplyr::rename(associatedTaxa_verbatim   = associatedTaxa)
      if(return_data_frame) {return(x)} else {return(x$associatedTaxa_translated)}
    }  else {
      if(return_data_frame) {return(x)} else {return(x$associatedTaxa)}
    }
  } else {
      if(return_data_frame) {return(x)} else {return(NULL)}
  }
}

## output summary tables
## 1 present in country 1
## 2 present in both country 1 and country 2
## 3 present in country 1, but not in country 2
## 4 not present in country 1
## 5 not present in country 2
## 6 not present in country 1 or 2
pest_distr_summary <- function(oi_hosts,
                               distributions,
                               associations,
                               taxonomy,
                               host_tax_tab,
                               country_ref1,
                               country_ref2,
                               confidence_level = "confident",
                               option=c(1,2,3,4,5)) {
  
  ## begin data process
  
  host_full <- filter(host_tax_tab) %>% filter(host_sp %in% oi_hosts$host_sp_full_accepted) %>% select(host_sp_full)
  assoc_oi <- associations %>% filter(host_sp_full_accepted %in% host_full$host_sp_full)
  pests_oi <- unique(assoc_oi$pest_sp_full_accepted)
  
  ## filter out what is wanted
  
  if (option == 1) {
  # 1 present in C1
    oi1 <- (distributions %>% filter(pest_sp_full_accepted %in% pests_oi) %>% filter(country == country_ref1))
    oi2 <- (distributions %>% filter(pest_sp_full_accepted %in% pests_oi) %>% filter(country == country_ref2))
    pests_oi <- oi1$pest_sp_full_accepted[which(oi1[,confidence_level])] # present in Poland

  } else if (option == 2) {
  # 2 present in both C1 and C2
    oi1 <- (distributions %>% filter(pest_sp_full_accepted %in% pests_oi) %>% filter(country == country_ref1))
    oi2 <- (distributions %>% filter(pest_sp_full_accepted %in% pests_oi) %>% filter(country == country_ref2))
    pests_oi <- oi1$pest_sp_full_accepted[which(oi1[,confidence_level])] %>%
    intersect(oi2$pest_sp_full_accepted[which(oi2[,confidence_level])])

  } else if (option == 3) {
  # 3 present in C1 not in C2
  oi1 <- (distributions %>% filter(pest_sp_full_accepted %in% pests_oi) %>% filter(country == country_ref1))
  oi2 <- (distributions %>% filter(pest_sp_full_accepted %in% pests_oi) %>% filter(country == country_ref2))
  pests_oi <- oi1$pest_sp_full_accepted[which(oi1[,confidence_level])] %>%
    setdiff(oi2$pest_sp_full_accepted[which(oi2[,confidence_level])])# present in Poland not in US

  } else if (option == 4) {
  # 4 not present in C1
  oi1 <- (distributions %>% filter(pest_sp_full_accepted %in% pests_oi) %>% filter(country == country_ref1))
  pests_oi <- oi1$pest_sp_full_accepted[-which(oi1[,confidence_level])]
  
  } else if (option == 5) {
  # 5 not present in C2
  oi2 <- (distributions %>% filter(pest_sp_full_accepted %in% pests_oi) %>% filter(country == country_ref2))
  pests_oi <- oi2$pest_sp_full_accepted[-which(oi2[,confidence_level])]
  
  } else if (option == 6) {
  # 6 not present in C1 or C2
  oi1 <- (distributions %>% filter(pest_sp_full_accepted %in% pests_oi) %>% filter(country == country_ref1))
  oi2 <- (distributions %>% filter(pest_sp_full_accepted %in% pests_oi) %>% filter(country == country_ref2))
  pests_oi <- oi1$pest_sp_full_accepted[-which(oi1[,confidence_level])] %>%
    intersect(oi2$pest_sp_full_accepted[-which(oi2[,confidence_level])])# present in Poland not in US
  }
  
  tax <- taxonomy %>% filter(name %in% pests_oi)
  assoc <- assoc_oi %>% filter(pest_sp_full_accepted %in% pests_oi)
  
  return (list(tax=tax, assoc=assoc))
}

join_tax_assoc <- function (pt, ht, a) {
  out <- ht %>%
    dplyr::rename(host_sp_full_accepted = host_sp_full) %>%
    plyr::join(a, by = "host_sp_full_accepted", type="inner", match="all") %>%
    dplyr::rename(name = pest_sp_full_accepted) %>%
    plyr::join(pt, by = "name", type="inner", match="all") %>%
    dplyr::rename(pest_sp_full_accepted= name)
  return(out)
}

summary_pests <- function(joined_list, pest_levels=c("pest_kingdom","pest_phylum","pest_class","pest_order","pest_family","pest_genus","pest_sp_full_accepted"), host_level, host_taxon) {
  out<- joined_list[which(joined_list[,host_level] == host_taxon), c("pest_kingdom","pest_phylum","pest_class","pest_order","pest_family","pest_genus","pest_sp_full_accepted")] %>%
    distinct %>%
    mutate(identity = 1)
  if (dim(out)[1]>0)
    return(out %>%
    summaryBy(as.formula(paste("identity ~", paste(pest_levels, collapse="+"))), ., FUN=sum))
  else return (NULL)
}

summary_pests_long <- function(joined_list, host_level, host_taxon) {
  
  kingdoms <- c("Fungi", "Bacteria", "Viruses")
  phyla    <- c("Basidiomycota", "Ascomycota", "Oomycota", "Nematoda")
  classes  <- c("Insecta", "Agaricomycetes", "Pucciniomycetes")
  orders   <- c("Lepidoptera", "Hemiptera", "Coleoptera",
                "Diaporthales", "Capnodiales", "Botryosphaeriales",
                "Pleosporales", "Helotiales", "Hypocreales",
                "Amphisphaeriales", "Ophiostomatales",
                "Xylariales", "Erysiphales", "Microascales","Rhytismatales",
                "Glomerellales", "Santalales", "Xanthomonadales", "Pseudomonadales")
  families <- c("Curculionidae", "Cerambycidae", "Buprestidae")
  taxa.order<-c("Insecta","Lepidoptera","Hemiptera","Coleoptera","Curculionidae","Cerambycidae","Buprestidae","Fungi","Basidiomycota","Agaricomycetes","Pucciniomycetes","Ascomycota","Rhytismatales","Diaporthales","Capnodiales","Botryosphaeriales","Pleosporales","Helotiales","Hypocreales","Amphisphaeriales","Ophiostomatales","Xylariales","Erysiphales","Microascales","Glomerellales","Oomycota","Nematoda","Santalales","Bacteria","Xanthomonadales","Pseudomonadales","Viruses")
  
  total<-summary_pests(joined_list, "pest_sp_full_accepted", host_level, host_taxon)$identity.sum %>% sum

  if (total>0) {
  
    fams<- summary_pests(joined_list, c("pest_family"), host_level, host_taxon) %>% dplyr::filter(pest_family %in% families) %>% select(pest_family, identity.sum) %>% dplyr::rename(Taxa=pest_family)
    ords<- summary_pests(joined_list, c("pest_order"), host_level, host_taxon)%>% dplyr::filter(pest_order %in% orders) %>% select(pest_order, identity.sum)%>% dplyr::rename(Taxa=pest_order)
    clas<- summary_pests(joined_list, c("pest_class"), host_level, host_taxon) %>% dplyr::filter(pest_class %in% classes) %>% select(pest_class,identity.sum)%>% dplyr::rename(Taxa=pest_class)
    phyls<-summary_pests(joined_list, c("pest_phylum"), host_level, host_taxon) %>% dplyr::filter(pest_phylum %in% phyla) %>% select(pest_phylum, identity.sum)%>% dplyr::rename(Taxa=pest_phylum)
    kings<-summary_pests(joined_list, c("pest_kingdom"), host_level, host_taxon) %>% dplyr::filter(pest_kingdom %in% kingdoms) %>% select(pest_kingdom, identity.sum)%>% dplyr::rename(Taxa=pest_kingdom)
  
    joined_taxa_summary <- rbind(fams,ords,clas,phyls,kings)
    rownames(joined_taxa_summary) <- joined_taxa_summary$Taxa
    names(joined_taxa_summary)[2] <- host_taxon
  
    k<-setdiff(taxa.order, row.names(joined_taxa_summary))
    to.add <- data.frame(Taxa=k, rep(0, length(k)))
    names(to.add)[2] <- host_taxon
    row.names(to.add) <- k
    joined_taxa_summary <- rbind(joined_taxa_summary, to.add)
    tt <- data.frame(Taxa="Total", identity.sum=total)
    names(tt)[2]<- host_taxon
    row.names(tt) <-"Total"
    joined_taxa_summary <- rbind(joined_taxa_summary,tt)
  
    return(joined_taxa_summary[c(taxa.order,"Total"),])
  } else {return(NULL)}
}

## Pests for group of interest
get_pest_host_taxon_EPPO <- function(hostlist, taxon, db, t) {
	sp <- hostlist[hostlist$Grouping %in% taxon, "sp"]
	hosts_names_tables <- eppo_names_tables(sp, db)
	hosts_pests <- eppo_tabletools_pests(hosts_names_tables, t)
	hosts_pests
	# species level only
	genus_fam_levels <- grep("\\(as", hosts_pests$long_table$fullname)
	pests_sp_level <- hosts_pests$long_table[-genus_fam_levels,]
	#pests_sp_level
	out<-pests_sp_level[pests_sp_level$labelclass != "Doubtful host",]
	# fill in subsp
	# spvars <- which.na()
	out
}

## EPPO codes for pests
get_pest_eppocodes <- function(pest_spp, db) {
	res <- eppo_names_tables(pest_spp, db)$pref_names
	hyb <- grep(" x ", res$fullname)
	matches <- res[-hyb,2:3]
}

get_host_eppocodes <- function(host_spp, db) {
	res <- eppo_names_tables(host_spp, db)$all_associated_names[,c("fullname","eppocode")]
	names(res)[1] <- "sp"
	host_spp_df <- as.data.frame(host_spp)
	names(host_spp_df)[1] <- "sp"
	(plyr::join(host_spp_df, res, "sp", type="left", match="first"))$eppocode
}

get_pest_eppocodes <- function(pest_spp, db) {
	res <- eppo_names_tables(pest_spp, db)$all_associated_names[,c("fullname","eppocode")]
	names(res)[1] <- "sp"
	pest_spp_df <- as.data.frame(pest_spp)
	names(pest_spp_df)[1] <- "sp"
	(plyr::join(pest_spp_df, res, "sp", type="left", match="first"))$eppocode
}

## Check against reference country to find pests present in U.S.

distributions_sentinels <- function(code, table, country="United States of America")
	{country %in% table[table$eppocode == code, "country"]}

ex_geo_search <- function (pests_geo, ref_loc="United States of America") {
	ex_outside <- sapply(unique(pests_geo$long_table$eppocode), function(x) distributions_sentinels(x, pests_geo$long_table, ref_loc), simplify="matrix")
	ex_outside
}

get_assoc_GLOBI <- function(host_tax="Pinus taeda", type="killedBy") {
	pests.globi <- rglobi::get_interactions(taxon=host_tax, interaction.type=type)
	if (dim(pests.globi)[1]==0) return (NULL)
	pests.globi[,c("source_taxon_external_id","source_taxon_name","interaction_type","target_taxon_external_id","target_taxon_name")]
#	cbind(
#		name_translation_table_sponly[
#			which(name_translation_table_sponly$in_NA=="FALSE"),
#			c("kingdom","phylum","order","family","species","speciesKey")],
#		data.frame(Host=host_tax, Interaction=paste(type, sep="-")))
}

get_taxonomy_GBIF <- function(pests_tax) {
	name_translation_table <- sapply(X=unique(pests_tax), FUN=function (x) c(fullname=x, rgbif::name_backbone(x)))
	if (is.null(dim(name_translation_table))){
		name_translation_table <- dplyr::bind_rows(name_translation_table)
	}else{
		name_translation_table <- as.data.frame(t(name_translation_table))}
	w <- which(is.na(name_translation_table$speciesKey))
	if (length(w)>1) name_translation_table <- name_translation_table[-w,]
	if (dim(name_translation_table)[1] == 0) return(NULL)
	name_translation_table
}

# occ_data
get_distri_GBIF <- function(p.spkey) {
	m_c <- occ_data(p.spkey)$meta$count
	occ_data(taxonKey=p.spkey, limit=m_c)$data
}

plot.pest <- function(c, p.spkey, occ=occ_data(p.spkey), p.name) {
	plot(c)
	distri <- occ$data
	points(distri$decimalLongitude, distri$decimalLatitude, col="red", pch=3)
	text(-180,-90, p.name, adj=0, col="red")
}

pest_in_country <- function (sp_key, occ=occ_data(sp_key), country_codes, country_polys, country_name=NULL) {
	distr <- occ$data
	if ("country" %in% names(distr))
		if (sum(country_name %in% distr$country) > 0) return (T)
	if(!("decimalLongitude" %in% names(distr))) return (NA)
	if(length(na.omit(distr$decimalLongitude))==0 | length(na.omit(distr$decimalLatitude))==0) return(NA)
	sp_sp <- sp::SpatialPoints(na.omit(cbind(Longitude=distr$decimalLongitude, Latitude=distr$decimalLatitude)), proj4string=sp::CRS("+proj=longlat +datum=WGS84"))
	p_in_p<- sp::over(sp_sp, country_polys)
	total_p_in_p <- sum(dplyr::select(dplyr::mutate( p_in_p , in_c = geoNameId %in% country_codes), in_c))
	total_p_in_p > 0
}

pest_in_US_GBIF <- function(x) 
	occ_data(taxonKey=x, country="US", basisOfRecord="HUMAN_OBSERVATION;MATERIAL_SAMPLE;PRESERVED_SPECIMEN;LIVING_SPECIMEN;LITERATURE;OBSERVATION")$meta$count > 0

pest_in_USCA_GBIF <- function(x) 
	(occ_data(taxonKey=x, country="US", basisOfRecord="HUMAN_OBSERVATION;MATERIAL_SAMPLE;PRESERVED_SPECIMEN;LIVING_SPECIMEN;LITERATURE;OBSERVATION")$meta$count > 0) |
		(occ_search(taxonKey=x, country="CA", basisOfRecord="HUMAN_OBSERVATION;MATERIAL_SAMPLE;PRESERVED_SPECIMEN;LIVING_SPECIMEN;LITERATURE;OBSERVATION")$meta$count > 0)

pest_in_NA_GBIF <- function(x)
	occ_data(taxonKey=x, continent="north_america", basisOfRecord="HUMAN_OBSERVATION;MATERIAL_SAMPLE;PRESERVED_SPECIMEN;LIVING_SPECIMEN;LITERATURE;OBSERVATION")$meta$count > 0

reformat_synonym_single <- function(nm, rk, codes=c(FORM="f.",VARIETY="var.",SUBSPECIES="subsp.")) {
  if (rk %in% names(codes)) {
    z <- strsplit(nm, " ")[[1]]
    nm <- z[1:2] %>% paste(collapse=" ") %>% paste(codes[[rk]], z[3])
  }
  return (nm)
}

synonym_table_GBIF <- function (x) {
  out.table <- NULL
  for (i in 1:length(x)) {
    y<- name_backbone(x[i])
    if(is.null(y)) {
      next
    } else if(y$matchType=="NONE") {
      out.table <- rbind.fill(out.table, data.frame(scientificName=i, canonicalName=i, species=i, verbatim_name=i, matchType="NONE", status="NONE"))
    } else if(sum(c("scientificName","canonicalName","species","verbatim_name", "matchType","status") %in% names(y))==6) {
      out.table <- rbind.fill(out.table, y)
    }
  }
  return(out.table)
}

key_table_GBIF_single <- function(x) {
  m <- name_backbone(x)
  if (m$matchType == "NONE") return (c(x, "NO MATCH"))
  if (m$matchType == "HIGHERRANK") return (c(x, "HIGHER RANK"))
  if (m$rank == "SPECIES") return (c(x, m$speciesKey))
}
key_table_GBIF <- function (x) foreach(a=x, .combine='rbind') %do% key_table_GBIF_single(a)

tax_table_GBIF_single <- function(x) {
  m <- name_backbone(x)
  g <- names (m)
  tax <- c(kingdom=NA, phylum=NA, class=NA, order=NA, family=NA, genus=NA)
  if (m$matchType != "NONE") for (y in names(tax)) if (y %in% g) tax[y] <- m[1,y]
  tax
}

tax_table_GBIF <- function(x) {
  (foreach(a=x, .combine='rbind') %do% tax_table_GBIF_single(a)) %>% as.data.frame %>%
    dplyr::rename(pest_kingdom=kingdom, pest_phylum=phylum, pest_class=class, pest_order=order, pest_family=family, pest_genus=genus)
}


############### added from script4 get geography gbif 3/6/23

# Note that currently the function treats all errors as curl timeouts (10 sec max default?)

try.pest.chunk.gbif <- function (pred.list) {
  pest.chunk <- tryCatch({
    q<-occ_download_queue(.list=pred.list);
    return(list(ok=TRUE, contents=q))
  }, error=function(cond) {
    message(paste(cond, "\n"))
    return(list(ok=FALSE))
  })
  return(pest.chunk)
}

get.pest.chunk.gbif <- function (pred.list, a, b, times=5, secs=30) {
  count <- 0L
  while (count < times) {
    g <- try.pest.chunk.gbif(pred.list)
    count <- count + 1L
    OK <- g$ok
    if(OK) break
    paste(paste("Trying chunk ", a, "-", b," again (", count+1, "/", times, ")", sep=""))
    Sys.sleep(time=secs)
  }
  
  if(OK) return(g$contents)
  else {
    print(paste("Chunk" , a, "-", b, " failed", sep=""))
    return(NULL)
  }
}