library(cancensus)
library(sf)
library(dplyr)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### INTRODUCTION ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Dans ce script, nous allons télécharger un ensemble de données de recensement
#' de 2021 à l'échalle des secteurs de recensement pour la métropole de Montréal.
#' Nous nous intéressons au ratio entre la proportion des individus utilisant 
#' principalement un mode motorisé individuel ou un mode durable (TC et actif)
#' pour leur déplacement de travail.
#' 
#' Nous modéliserons ensuite cette variable avec un ensemble de variables de 
#' recensement : 
#' 1) les variables d'accessibilité proposées par StatCan
#' 2) plusieurs variables socio-économiques
#' 3) l'espace
#' 


setwd('E:/Cours/Cours/Analyse spatiale/TP/Tp4/Data2025')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Téléchargement des données de recensement ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



## Afficher la liste des données pour le recensement de 2021
Tableau.Variables2021 <- list_census_vectors('CA21')
View(Tableau.Variables2021)


## Liste des variables sur les modes de transport
transport_vars <- c('mode_tc' = 'v_CA21_7644',
  'total_commuters' = 'v_CA21_7632',
  'mode_auto' = 'v_CA21_7635',
  'mode_pieton' = 'v_CA21_7647',
  'mode_velo' = 'v_CA21_7650'
  )

## Liste des variables socio-demo
socio_demo_vars <- c(
  'prt_personnes_faibles_revenu' = 'v_CA21_1085',
  'nb_non_minorite_visible_1' = 'v_CA21_4914',
  'total_minorite_visible_1' = 'v_CA21_4872',
  'nb_menages_monoparentaux_2' = 'v_CA21_548',
  'total_menage_2' = 'v_CA21_544',
  'nb_sans_emplois_3' = 'v_CA21_6501',
  'total_actif_3' = 'v_CA21_6495',
  'densite_population_km2' = 'v_CA21_6',
  'revenu_median' = 'v_CA21_983'
)


SR2021.MTL <- get_census(dataset='CA21',
                         regions=list(CMA="24462"),
                         level='CT',
                         vectors=c(transport_vars, socio_demo_vars),
                         quiet = TRUE,
                         geo_format = 'sf', labels = 'short')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Structuration des données d'accessibilité ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

acces_tc_peak <- read.csv('statcan_access/acs_public_transit_peak.csv')
access_velo <- read.csv('statcan_access/acs_cycling.csv')
acces_tc_offpeak <- read.csv('statcan_access/acs_public_transit_offpeak.csv')
access_pieton <- read.csv('statcan_access/acs_walking.csv')
ilots <- st_read('statcan_access/lid_000b21a_f.shp')

acces_tc_peak$DBUID <- as.character(acces_tc_peak$DBUID)
access_velo$DBUID <- as.character(acces_tc_peak$DBUID)
acces_tc_offpeak$DBUID <- as.character(acces_tc_offpeak$DBUID)
access_pieton$DBUID <- as.character(access_pieton$DBUID)

# on mets les deux jeux de données dans le même système de référence
SR2021.MTL <- st_transform(SR2021.MTL, st_crs(ilots))

# on ne garde que les ilots dans notre territoire d'étude
ilots <- subset(
  ilots,
  lengths(st_intersects(st_centroid(ilots), SR2021.MTL)) > 0
)

# on effectue les différentes jointures pour avoir nos données spatiales à l'échelle
# des ilots
ilots_complet <- ilots %>% 
  left_join(acces_tc_peak[c('DBUID', 'acs_idx_emp')], by = c('IDIDU' = 'DBUID')) %>% 
  rename(any_of(c('acs_idx_emp_tc_peak' = 'acs_idx_emp'))) %>% 
  left_join(acces_tc_offpeak[c('DBUID', 'acs_idx_emp')], by = c('IDIDU' = 'DBUID')) %>% 
  rename(any_of(c('acs_idx_emp_tc_offpeak' = 'acs_idx_emp'))) %>% 
  left_join(access_velo[c('DBUID', 'acs_idx_emp')], by = c('IDIDU' = 'DBUID')) %>% 
  rename(any_of(c('acs_idx_emp_velo' = 'acs_idx_emp'))) %>% 
  left_join(access_pieton[c('DBUID', 'acs_idx_emp')], by = c('IDIDU' = 'DBUID')) %>% 
  rename(any_of(c('acs_idx_emp_pieton' = 'acs_idx_emp')))

# on ajoute la population à chaque ilot pour pouvoir ensuite calculer des moyennes pondérées
# au niveau des SR
population_ilots <- read.csv('statcan_access/2021_92-151_X.csv')
ilots_complet$population_2021 <- population_ilots$DBPOP2021_IDPOP2021[match(ilots_complet$IDUGD, population_ilots$DBDGUID_IDIDUGD)]

# on trouve pour chaque ilot dans que SR il se trouve
ilots_complet <- st_point_on_surface(ilots_complet)
ilots_complet <- st_intersection(ilots_complet, SR2021.MTL[c('GeoUID')])

# on calcule finalement nos moyennes pondérées par SR
data_SR_access <- ilots_complet %>% 
  st_drop_geometry() %>% 
  group_by(GeoUID) %>%
  summarise(
    acs_idx_emp_tc_peak = weighted.mean(acs_idx_emp_tc_peak, population_2021),
    acs_idx_emp_tc_offpeak = weighted.mean(acs_idx_emp_tc_offpeak, population_2021),
    acs_idx_emp_velo = weighted.mean(acs_idx_emp_velo, population_2021),
    acs_idx_emp_pieton = weighted.mean(acs_idx_emp_pieton, population_2021),
  )


final_data <- left_join(SR2021.MTL, data_SR_access, by = c('GeoUID'))

final_data$prt_monoparental <- final_data$nb_menages_monoparentaux_2 / final_data$total_menage_2
final_data$prt_minorite_vis <- 1 - (final_data$nb_non_minorite_visible_1 / final_data$total_minorite_visible_1)
final_data$prt_chomage <- final_data$nb_sans_emplois_3 / final_data$total_actif_3


st_write(final_data[c('prt_monoparental',
                      'prt_minorite_vis',
                      'prt_chomage',
                      'prt_personnes_faibles_revenu',
                      'revenu_median',
                      'densite_population_km2', 
                      'mode_tc', 
                      'mode_auto', 
                      'mode_pieton', 
                      'mode_velo',
                      "total_commuters",
                      'acs_idx_emp_velo',
                      'acs_idx_emp_tc_peak',
                      'acs_idx_emp_tc_offpeak',
                      'acs_idx_emp_pieton',
                      'Households',
                      'Dwellings',
                      'Population',
                      'GeoUID'
                      )], 'data_sr_access.gpkg', delete_layer = T)



