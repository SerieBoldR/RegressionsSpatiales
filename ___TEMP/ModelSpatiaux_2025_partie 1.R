# Code pour l'analyse par regression spatiale dans R
# Auteurs : Jérémy Gelb et Philippe Apparicio

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### installation et chargement des packages a utiliser ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(sf)
library(foreign)		## Read Data Stored by Minitab, S, SAS, SPSS, Stata, Systat, dBase.
library(nortest)      	## Tests de normalit? supp.
library(spdep)		## Spatial dependence: weighting schemes, statistics and models
library(spatialreg) ## pour les modeles SAR et CAR
library(car)            ## options r?gression (pour VIF)
library(spgwr)       	## GWR
library(moments)
library(tidyverse)
library(performance) # pour calculer de r2 sur les modeles GLM
library(RColorBrewer)
library(classInt)
library(nlme) #pour les modele GLS
library(ggpubr)
library(tmap)
library(corrplot) # pour réaliser des matrices de corrélation sous forme graphique

rm(list=ls())


#' objectif de l'analyse : Nous tentons de modéliser les prix des appartement 
#' Airbnb sur l'île de Montréal à partir d'un ensemble de caractéristiques 
#' propres à chaque logement. Nous introduirons ensuite un terme spatial pour
#' faire ressortir l'impact de l'espace sur le prix des logements (quartiers en demande)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Definition de fonctions utilitaires ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TestNormality <- function(x,showHelp=T){
  #calcul du test shapiro.test pour evaluer la normalite
  if(length(x) > 5000){
    Test <- shapiro.test(sample(x, size = 5000, replace = FALSE))
  }else{
    Test <- shapiro.test(x)
  }
  
  #affichage de la distribution et de la courbe normale pour reference
  Mean <- mean(x)
  Sd <- sd(x)
  df <- data.frame(x = x)
  Plot <- ggplot(df) + 
    geom_histogram(aes(x=x,y =..density..),bins=30,fill="white",color="black")+
    stat_function(aes(x=x),fun = dnorm, args = list(mean = Mean, sd = Sd),color="red",size=1)
  #calcul des descripteurs de la distribution
  Values <- list("mean"=Mean,
                 "sd"=Sd,
                 "skewness"=skewness(x),
                 "kurtosis"=kurtosis(x),
                 "shapiro.p"=Test$p.value,
                 "plot"=Plot)
  print(Plot)
  if (showHelp){
    print("------  interpreting guide  ------")
    print("the variable is close to normality if")
    print("the skweness value is around 0, this describe the symetric aspect of the distribution")
    print("the kurtosis value is around 3, this describe the spiky aspect of the distribution")
    print("this shapiro.p value is significant (above 0.05)")
  }
  return(Values)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### preparation de l'environnement de travail ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setwd("C:/Users/APPP2302/OneDrive - USherbrooke/Articles Publies ne pas supprimer/_Livres/BolR_RegressionsSpatiales/___TEMP")

# on charge les données à partir du fichier CSV
data_airbnb <- read.csv('Data_2025/airbnb_data.csv')

# nous allons ensuite spatialiser cette données
data_airbnb <- st_as_sf(data_airbnb,
                        coords = c('longitude', 'latitude'),
                        crs = 4326) %>% 
  st_transform(32188)


MTL <- st_read('Data_2025/ILE_MTL.shp')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### PRÉPARATION DES DONNÉES ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_airbnb$host_is_superhost <- ifelse(data_airbnb$host_is_superhost == 't', 'YES', 'NO')


hist(data_airbnb$price)

#' les prix par nuit varient entre presque 10 et 250$
#' Nous allons tenter de modéliser les prix en fonction d'un ensemble 
#' de variables : 
#' - private : binaire, est-ce que le logement est partagé ou privé
#' - Free_street_parking : binaire, présence d'un stationnement gratuit
#' - Garden_or_backyard : binaire, présence d'un jardin
#' - bedrooms : nombre de chambres
#' - host_is_superhost : binaire, est-ce que l'hôte a le statut de super-hote
#' - number_of_reviews : nombre d'évaluation
#' - review_scores_value : la note du logement
#' - has_metro_500m : binaire, présence d'une station de métro à 500m
#' - prt_veg_500m : pourcentage, couverture végétale dans un rayon de 500m

#' nous allons recoder certaines variables catégorielles afin de déterminer 
#' la catégorie de référence pour faciliter l'interprétation

data_airbnb$private <- factor(data_airbnb$private, levels = c('Chambre', 'Entier'))
data_airbnb$Free_street_parking <- factor(data_airbnb$Free_street_parking, levels = c('NO', 'YES'))
data_airbnb$Garden_or_backyard <- factor(data_airbnb$Garden_or_backyard, levels = c('NO', 'YES'))
data_airbnb$host_is_superhost <- factor(data_airbnb$host_is_superhost, levels = c('NO', 'YES'))
data_airbnb$has_metro_500m <- factor(data_airbnb$has_metro_500m, levels = c('NO', 'YES'))


# Enfin, nous allons seulement garder les données pour lesquelles nous avons
# toutes les informations

model_data <- data_airbnb %>%
  select(any_of(c('price','private','Free_street_parking', 'Garden_or_backyard', 
                    'bedrooms', 'host_is_superhost', 
                    'number_of_reviews','review_scores_rating','has_metro_500m',
                    'prt_veg_500m'))) %>% 
  filter(complete.cases(st_drop_geometry(.)))


# on va regarder l'apport de 10 review supplémentaire
data_airbnb$number_of_reviews <- data_airbnb$number_of_reviews / 10



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### VISUALISATION DE LA VARIABLE Y ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

summary(model_data$price)

# il sera difficile de cartographier notre variable Y car de nombreux Airbnb se 
# superposent. Nous allons donc calculer la moyenne des valeurs des logements
# dans une maille régulière

# on calcule la superficie totale de la zone d'étude
S <- as.numeric(st_area(MTL))
# on compte la nombre de logement
N <- nrow(model_data)
# on détermine la superficie idéale d'un quadra
A <- (2 * S) / N # formule dans le cours

quadra_hexa <- st_make_grid(
  MTL, # nous créeons une grille englobant notre objet spatial Quartiers
  cellsize = 2 * sqrt(A/((3*sqrt(3)/2))) * sqrt(3)/2, # pour des quadra hexagonaux 
  what = "polygons",
  square = FALSE
)

quadra_hexa <- st_as_sf(
  df <- data.frame(OID = 1:length(quadra_hexa)),
  geometry = quadra_hexa,
  crs = st_crs(quadra_hexa)
)

# nous allons ensuite déterminer dans quel quadra tombe chaque airbnb
# puis calculer la moyenne des prix par quadra
inter <- st_intersection(model_data, quadra_hexa) %>% 
  st_drop_geometry() %>% 
  group_by(OID) %>% 
  summarise(
    mean_price = mean(price)
  )

quadra_hexa <- left_join(quadra_hexa, inter, by = 'OID')
quadra_hexa <- subset(quadra_hexa, !is.na(quadra_hexa$mean_price))

tm_shape(quadra_hexa) + 
  tm_fill('mean_price',title = 'prix ($)', 
          breaks = c(15, 40, 80, 120, 160, 220, 260)
          ) + 
  tm_shape(model_data) + 
  tm_dots('black', size = 0.01, alpha = 0.1) +
  tm_compass(position = c("right", "bottom"), 
             size = 0.5)+
  tm_legend(position = c("left", "top"),
            frame = FALSE, bg.color = "white")+
  tm_scale_bar(breaks  = c(0, 5, 10),
               position = c("left", "bottom"))+
  tm_layout(main.title = "Prix moyen des Aibnb par nuit",
            legend.outside = TRUE, 
            attr.outside = TRUE,
            inner.margins = 0,
            frame = FALSE, 
            legend.format = list(text.separator = "-",
                                 digits = 0,
                                 text.or.more = 'ou plus'
            )
  ) 


# nous pouvons calculer l'autocorrélation spatiale de la variable en utilisant 
# une matrice des N plus proches voisin
XY <- st_coordinates(model_data)
nearest5 <- knearneigh(XY,k=5)
nearest5_nb <- knn2nb(nearest5, row.names = NULL, sym = FALSE)
nearest5_listw <- nb2listw(nearest5_nb, style = "W")
moran.mc(model_data$price, listw = nearest5_listw, nsim = 999)

# nous pouvons constater que l'autocorrélation spatiale est relativement modérées
# (i = 0.13), mais tout de même significative pour la variable Y


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### PRÉPARATION DES DONNÉES ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Le premier modèle que nous réalisons est un simple modèle GLM utilisant une
#' distribution Gaussienne des résidus.

model_GLM <- glm(price ~ 
                private + Free_street_parking + Garden_or_backyard + 
                bedrooms + host_is_superhost +  
                number_of_reviews + review_scores_rating + has_metro_500m + 
                prt_veg_500m,
              family = gaussian,
              data = model_data
                )

# nous pouvons afficher les coefficients ainsi obtenus : 

summary(model_GLM)

#' Interprétation : 
#' 
#' - Le prix pour un logement privé entier est 40$ plus élevé en moyenne que pour 
#' un logement partagé
#' - Chaque chambre supplémentaire augment le prix moyen de 21$
#' - si le host est un super-host, le prix augmente en moyenne de 6.74 $
#' - le nombre de review ne semble pas avoir d'impact
#' - une augmentation de 1% du couvert végétal à 500m d'un airbnb tend 
#' à réduire son coût de 61 centimes
#' - Le fait de se situer à moins de 500 d'un métro tend à augmenter de 2.13$ le
#' coût d'un airbnb, mais cet effet n'est significatif qu'au seuil 0.05
#' - Une augmentation d'un point de pourcentage pour le score de review est 
#' associé avec un prix moyen plus élevé de 62 centimes


# Verification de la correlation entre les variables

vif(model_GLM) # aucun VIF ne dépasse 5, donc absence de multicolinéarité

# On peut regarder le R2 du modèle pour cerner sa capacité de prédiction
r2(model_GLM)

# nous pouvons regarder la distribution des résidus (pour ce modèle, elle devrait 
# être proche d'une distribution normale

TestNormality(residuals(model_GLM))

# Ces résidus sont asymétrique et pointus. Il serait donc pertinent de modéliser
# nos pas le prix des airbnb, mais le log du prix.
# Ainsi, nous obtenons un modèle non plus additif, mais multiplicatif

model_data$log_price <- log(model_data$price)

model_GLM2 <- glm(log_price ~ 
                   private + Free_street_parking + Garden_or_backyard + 
                   bedrooms + host_is_superhost +  
                   number_of_reviews + review_scores_rating + has_metro_500m + 
                   prt_veg_500m,
                 family = gaussian,
                 data = model_data
)

TestNormality(residuals(model_GLM2))

# Les résidus sont bien mieux distribués, et il serait plus cohérents d'analyser
# les résultats de ce modèle.

r2(model_GLM2)

# Nous avons aussi nettement augmenté le R2

# Verification de l'autocorrélation spatiale des résidus
moran.mc(model_data$log_price, listw = nearest5_listw, nsim = 999)

moran.mc(residuals(model_GLM2), listw = nearest5_listw, nsim = 999)

# Nous avons une autocorrélation spatiale résiduelle non négligeable !

#' Interprétation : 
#' ATTENTION : le fait de modéliser le log du prix du airbnb change l'interprétation des coefficients.
#' En effet, on ne modélise plus le prix dans son échelle linéaire ($), mais dans son échelle logarithmique 
#' (log $). Additionner dans l'échelle logarithmique revient à additionner dans l'échelle linéaire. Donc les 
#' effets des coefficients ne sont plus additifs, mais multiplicatifs dans l'échelle linéaire.
#' 
#' - Le prix pour un logement privé entier est 0.64 log($) plus élevé en moyenne que pour 
#' un logement partagé. Ceci signifie qu'avoir un logement entier privé plutôt que 
#' partagé multiplie le prix moyen d'un logement par exp(0.6419), soit 1.9. Un logement privé 
#' tend donc à être presque deux fois plus cher qu'un logement partagé.
#' - Chaque chambre supplémentaire multiplie le prix moyen par exp(0.20), soit 1.22. Ainsi, une 
#' chambre supplémentaire augmente le prix d'un logement de 22%.
#' -  Chaque review supplémentaire contribue à augmenter le prix moyen d'un logement de 0.000329 log($).
#' Donc 50 review supplémentaires multiplie le prix par  exp(50 * 0.000329), soit 1.016. Donc 50 review 
#' supplémentaires augmentent le prix d'un airbnb de en moyenne 1.6%.
#' - si le host est un super-host, le prix est multiplié en moyenne par 1.07, soit une augmentation de 7%
#' - une augmentation de 10% du couvert végétal à 500m d'un airbnb tend à multiplier le prix par 0.937, soit une 
#' diminution de 6.3% du prix du airbnb.


# Analyse des outliers via la distance de Cook
cookD <- cooks.distance(model_GLM2)
plot(cookD)

# on observe clairement deux airbnb problématiques que l'on devrait retirer des données
View(subset(model_data, cookD > 0.1))
# Il s'agit de deux auberges avec 20 chambres !

model_data2 <- subset(model_data, cookD < 0.1)

# on doit recalculer notre matrice de voisinage
XY <- st_coordinates(model_data2)
nearest5 <- knearneigh(XY,k=5)
nearest5_nb <- knn2nb(nearest5, row.names = NULL, sym = FALSE)
nearest5_listw <- nb2listw(nearest5_nb, style = "W")

model_GLM2 <- glm(log_price ~ 
                    private + Free_street_parking + Garden_or_backyard + 
                    bedrooms + host_is_superhost +  
                    number_of_reviews + review_scores_rating + has_metro_500m + 
                    prt_veg_500m,
                  family = gaussian,
                  data = model_data2
)

# on constate que la nouvelle distribution des distances de Cook est beaucoup plus
# acceptable
cookD <- cooks.distance(model_GLM2)
plot(cookD)

# nous pouvons tenter de visualiser les résidus, mais une fois encore, 
# nous allons avoir quelques difficultés car de nombreux points vont se 
# superposer. Nous allons donc représenter en premier les résidus avec des 
# valeurs très faibles, puis très fortes

model_data2$resid_glm <- residuals(model_GLM2)

tmap_mode('view')

map1 <- tm_shape(model_data2[order(model_data2$resid_glm),]) + 
  tm_dots('resid_glm', n = 7, style = 'fisher', 
          palette = 'RdBu'
          )

map2 <- tm_shape(model_data2[order(-1 * model_data2$resid_glm),]) + 
  tm_dots('resid_glm', n = 7, style = 'fisher', 
          palette = 'RdBu'
  )

tmap_arrange(map1, map2, sync = TRUE)
moran.mc(residuals(model_GLM2), listw = nearest5_listw, nsim = 999)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Modeles SAR LAG et SAR ERROR ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Nous allons à présent appliquer des modèles économétriques à notre 
# jeu de donnée.


#### Diagnostic des test de Lagrange sur le mod?le OLS
ols_log <- 	lm(log_price ~ 
                 private + Free_street_parking + Garden_or_backyard + 
                 bedrooms + host_is_superhost +  
                 number_of_reviews + review_scores_rating + has_metro_500m + 
                 prt_veg_500m,
               data = model_data2)

# on effectue un test de Lagrange pour voir quel type de modèle nous devrions privilégier
test <- lm.LMtests(ols_log, listw = nearest5_listw, test=c("LMerr","RLMerr","LMlag","RLMlag"))

test$RSerr
test$RSlag
test$adjRSerr
test$adjRSlag

# NB : Ajuster un modèle SAR avec 9000 observations peut-être relativement long
# pour que l'on puisse avancer pendant le TD, nous allons tricher un peu et ajuster
# le modèle sur un sous-échantillon des données avec 2500 observations

sub_data <- model_data2[sample(1:nrow(model_data2), size = 2500, replace = FALSE),]
XY <- st_coordinates(sub_data)
nearest5 <- knearneigh(XY,k=5)
nearest5_nb <- knn2nb(nearest5, row.names = NULL, sym = FALSE)
nearest5_listw <- nb2listw(nearest5_nb, style = "W")


#ok il semble que l'on doive se tenter un model SARLAG (plus haute statistique pour le test ajusté et non ajusté)
SarLag_log <- lagsarlm(log_price ~ 
                         private + Free_street_parking + Garden_or_backyard + 
                         bedrooms + host_is_superhost +  
                         number_of_reviews + review_scores_rating + has_metro_500m + 
                         prt_veg_500m,
                      data = sub_data,
                      listw = nearest5_listw,
                      tol.solve=1e-15)

summary(SarLag_log,Nagelkerke=T) #

# Le R2 est un peu plus élevé que tout à l'heure (0.45)
# le Rho est le paramètre spatial. Sa valeur de 0.17 peut être interprêtée ainsi : 
# si la moyenne du log du prix des 5 plus proches voisins augmente de 1, alors 
# le log du prix du logement analysé augmenterait de 0.17. L'effet du voisinage se
# reporte de 17% sur le prix de chaque observation


TestNormality(SarLag_log$residuals) # très belle distribution des résidus

qqPlot(SarLag_log$residuals, distribution = 'norm')

# Absence d'autocorélation spatiale dans les résidus
moran.mc(SarLag_log$residuals,nearest5_listw,nsim = 999)


#verifions l'homeoscedasticite
DF <- st_drop_geometry(sub_data)
DF$Pred <- c(predict(SarLag_log))
DF$Resid <- residuals(SarLag_log)

ggplot(DF)+
  geom_point(aes(x=Pred,y=Resid))+
  geom_smooth(aes(x=Pred,y=Resid))

# pas trop mal, mais on observe quand même deux groupes de points qui se distinguent
# assez nettement. Il serait pertinent d'essayer de comprendre ce qui cause cette distinction.

#pour le fun, on va le comparer a un SARERROR
SarErr_log <- errorsarlm(log_price ~ 
                          private + Free_street_parking + Garden_or_backyard + 
                          bedrooms + host_is_superhost +  
                          number_of_reviews + review_scores_rating + has_metro_500m + 
                          prt_veg_500m,
                        data = sub_data,
                        listw = nearest5_listw,
                        tol.solve=1e-15)

summary(SarErr_log,Nagelkerke=T)
AIC(SarErr_log,SarLag_log)

TestNormality(residuals(SarErr_log))

# On a peu de différence entre les deux modèles


#et finalement un model DURBIN complet
# Ce modèle ajoute en plus une version spatialement laguée des prédicteurs.
# Il s'agit ici d'une version "full", car on permet a chaque prédicteur d'avoir
# une version spatialement décalée. Nous pouvons ensuire refaire un modèle avec 
# uniquement les prédicteurs utiles
SarDurbin_log_full <- lagsarlm(log_price ~ 
                                 private + Free_street_parking + Garden_or_backyard + 
                                 bedrooms + host_is_superhost +  
                                 number_of_reviews + review_scores_rating + has_metro_500m + 
                                 prt_veg_500m,
                               data = sub_data,
                               listw = nearest5_listw,
                             tol.solve=1e-15,
                             Durbin=T)

summary(SarDurbin_log_full)


#on constate que plusieurs WX sont significatifs, on ne va conserver que ceux-ci
SarDurbin_log <- lagsarlm(log_price ~ 
                           private + Free_street_parking + Garden_or_backyard + 
                           bedrooms + host_is_superhost +  
                           number_of_reviews + review_scores_rating + has_metro_500m + 
                           prt_veg_500m,
                         data = sub_data,
                         listw = nearest5_listw,
                         tol.solve=1e-15,
                         Durbin = ~ private + Free_street_parking + bedrooms + number_of_reviews)

summary(SarDurbin_log)
TestNormality(residuals(SarDurbin_log))
AIC(SarDurbin_log_full,SarDurbin_log)
moran.mc(residuals(SarDurbin_log), nearest5_listw, nsim = 999)


## Question intéressante : pourquoi celon-vous observe-t-on que le fait d'avoir
## une chambre en plus augmente le prix du logement, mais si la moyenne du nombre 
## de chambre chez les voisins est plus élevé, le prix du logement diminue ?
## En revanche, que se passe-t-il pour les review des voisins ?

#' guide d'interprétation pour les lag : 
#' Si tous les voisins ont le type logement privé entier, alors on observe en 
#' moyenne une diminution du log du prix du logement de 0.2211, soit multiplication 
#' du prix par 0.801 (exp(-0.2211)) ce qui correspond à une baisse de 20% du prix !
#' 
#' Si la moyenne du nombre de review des voisins est augmentée de 50, alors on a 
#' une augmentation du log du prix de 50 * 0.001386, ce qui correspond à une multiplication 
#' du prix par 1.07, soit une augmentation de 7% du prix du logement.



#observation des impacts directs et indirects du model SARLAG
ImpSarLag <- impacts(SarLag_log,listw = nearest5_listw,R=999)

summary(ImpSarLag,zstats=T,short=T)

#' l'interprétation des effets directs et indirects est un peu difficile, on doit 
#' commencer par les valeurs de p pour déterminer lesquels sont significatifs ou non.
#' 
#' Une fois que c'est fait, on peut se pencher sur l'effet direct qui peut s'interprêter comme un simple coefficient classique
#' 
#' 1) Le fait d'avoir un logement entier privé augmente en moyenne le prix de son logement
#' de 0.62 log($) en incluant les effets de rétro-action propres. Ainsi, si mon logement "devient" privé, 
#' son prix va augmenter, ce qui va contaminer le prix des voisins, puis recontaminer le prix de mon propre
#' logement. Ces effets combinés sur mon logement vont multiplier son prix par exp(0.62), soit 1.85, ce qui correspond à une 
#' augmentation de 85% de son prix !
#' 
#' 2) L'effet indirect est de 0.123, il mesure l'impact cumulatif sur les autres unités voisines. Donc si mon propre 
#' logement devient privé, on devrait observer une augmentation de la somme des prix de l'ensemble des autres logements de 0.123 log($).
#' Notez bien que ce 0.123 est une somme de l'effet de propagation sur tous les autres logements. Il serait difficile d'interprêter
#' directement cette valeur. Mais nous pouvons noter qu'elle est relativement faible ici et que l'effet direct représente 84% de 
#' l'effet total.
#' 
#' 3) L'effet total correspond à la somme de l'effet indirect et direct. Il s'agit donc de la taille de l'effet total
#' que l'on peut attendre sur tous notre système spatial.


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Modele CAR ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Un modèle CAR doit être réalisé avec une matrice symétrique ! Ceci pose problème pour notre 
#' matrice des 5 plus proches voisins. Nous allons donc forcer cette matrice à être symétrique, 
#' ce qui va augmenter le nombre de voisin
#' Notez que si l'on standardise la matrice de voisinage, on risque de se retrouver avec une 
#' matrice non symétrique car nous n'aurons pas le même nombre de voisin "

sym_nearest5_nb <- make.sym.nb(nearest5_nb)
sym_nearest5_listw <- nb2listw(sym_nearest5_nb, style = "B")

CAR_log <- spautolm(log_price ~ 
                     private + Free_street_parking + Garden_or_backyard + 
                     bedrooms + host_is_superhost +  
                     number_of_reviews + review_scores_rating + has_metro_500m + 
                     prt_veg_500m,
                   data = sub_data,
                  listw=sym_nearest5_listw,
                  zero.policy=TRUE,
                  tol.solve=1e-15, 
                  family = "CAR")

summary(CAR_log,Nagelkerke=T)

#' Dans ce modèle, le terme lambda est assez faible (0.084), indiquant une autocorrélation 
#' spatiale positive faible, mais tout de même significativement différente de 0

TestNormality(residuals(CAR_log))

moran.mc(residuals(CAR_log),sym_nearest5_listw,nsim=999)

# Nous pouvons maintenant regarder la contribution de notre "intercept spatial"

sub_data$car_intercept <- CAR_log$fit$signal_stochastic

# visualisation du terme spatial
# Une fois encore, nous devons essayer de représenter les logements en 
# montrant d'abord les valeurs fortes, puis les valeurs faible

tmap_mode('view')

map1 <- tm_shape(sub_data[order(sub_data$car_intercept),]) + 
  tm_dots('car_intercept', n = 7, style = 'fisher', 
          palette = 'RdBu'
  )

map2 <- tm_shape(sub_data[order(-1 * sub_data$car_intercept),]) + 
  tm_dots('car_intercept', n = 7, style = 'fisher', 
          palette = 'RdBu'
  )

tmap_arrange(map1, map2, sync = TRUE)



AIC(CAR_log,SarLag_log,SarDurbin_log)
# En terme de AIC, notre meilleur modèle pour le moment et le model SAR-DURBIN



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Modeles GLS ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

XY <- st_coordinates(sub_data)
sub_data$X <- XY[,1]
sub_data$Y <- XY[,2]

# NB : puisque nous avons des observations qui ont exactement la même localisation
# le modèles GLS ne vont pas pouvoir s'ajuster. Nous pouvons contourner le problème
# en ajouter de façon aléaoire un bruit léger dans les coordonnées X et Y. Si cette 
# distance est assez faible, elle n'aura pas d'impact sur l'ajustement du modèle

# Nous rajoutons ici un bruit de maximum 5 mètres pour chaque coordonées en X et en Y
sub_data$X <- sub_data$X + runif(n = nrow(sub_data), min = -5, max = 5)
sub_data$Y <- sub_data$Y + runif(n = nrow(sub_data), min = -5, max = 5)


#Mod?le Gls (Exponentiel)
modelGlsExp <- gls(log_price ~ 
                     private + Free_street_parking + Garden_or_backyard + 
                     bedrooms + host_is_superhost +  
                     number_of_reviews + review_scores_rating + has_metro_500m + 
                     prt_veg_500m,
                   data = sub_data,
                   correlation=corExp(form=~X+Y, nugget=T))

summary(modelGlsExp)

#Mod?le Gls (gaussian)
modelGlGaussian <- gls(log_price ~ 
                         private + Free_street_parking + Garden_or_backyard + 
                         bedrooms + host_is_superhost +  
                         number_of_reviews + review_scores_rating + has_metro_500m + 
                         prt_veg_500m,
                       data = sub_data,
                       correlation=corGaus(form=~X+Y, nugget=T))
summary(modelGlGaussian)

#Mod?le Gls (spher)
modelGlsSpher <- gls(log_price ~ 
                       private + Free_street_parking + Garden_or_backyard + 
                       bedrooms + host_is_superhost +  
                       number_of_reviews + review_scores_rating + has_metro_500m + 
                       prt_veg_500m,
                     data = sub_data,
                     correlation=corSpher(form=~X+Y, nugget=T))

#Mod?le Gls (Rationel Quadratic)
modelGlsRatio <- gls(log_price ~ 
                       private + Free_street_parking + Garden_or_backyard + 
                       bedrooms + host_is_superhost +  
                       number_of_reviews + review_scores_rating + has_metro_500m + 
                       prt_veg_500m,
                     data = sub_data,
                     correlation=corRatio(form=~X+Y, nugget=T))

#comparaison des 4 formes
AIC(modelGlsExp,modelGlGaussian,modelGlsSpher,modelGlsRatio)

# meilleur ajustement pour le model Exponentiel (AIC le plus bas)

#' L'interprétation d'un modèle GLS est très facile car les coefficients peuvent être
#' lus comme pour un simple GLM
summary(modelGlsExp)

# En revanche, nous devons analyser les résidus du modèle qui tiennent comptent 
# de la structure d'autocorrélation spatiale que nous avons rajouté
TestNormality(residuals(modelGlsExp,type="normalized"))

# et nous pouvons tester l'autocorrélation spatiale des résidus
moran.mc(residuals(modelGlsExp,type="normalized"),sym_nearest5_listw,nsim=999) # autocorrelation resolue

## Nous pouvons à présent afficher la structure de corrélation ainsi calculée
struc <- modelGlsExp$modelStruct$corStruct

# En allant sur la documentation de la fonction corExp, on peut trouver son 
# équation, permettant de calculer la corrélation attendue selon la distance, 
# et les deux paramètres nugget et range
exp_cor <- function(n,r,d){
  val <- (1-n) * exp( (-1*r) / d)
}

# On peut ensuite l'appliquer pour pouvoir le représenter avec un graphique
distances <- seq(0,25000,100)

df_corr <- data.frame(
  distance = distances, 
  correlation = exp_cor(n = 0.8422976, 
                        d = 1706.2655321, 
                        r = distances
                        )
)

ggplot(df_corr) + 
  geom_line(aes(x = distance, y = correlation)) + 
  theme_bw() + 
  labs(title = 'Structure de corrélation du modèle GLS')


#' on voit ainsi que selon le modèle GLS exponentiel, la corrélation entre deux 
#' airbnb est assez faible lorsqu'ils sont proches l'un de l'autre (maximum 0.15)
#' et diminue très fortement pour atteindre 0 à 10 km. 

#verifions l'homeoscedasticite
DF <- sub_data
DF$Pred <- modelGlsExp$fitted
DF$Resid <- residuals(modelGlsExp,type="normalized")

ggplot(DF)+
  geom_point(aes(x=Pred,y=Resid))+
  geom_smooth(aes(x=Pred,y=Resid))

# Une fois encore, on observe une forme de regroupement entre nos observations
# il semble donc clairement manquer une variable catégorielle qui permettrait
# d'expliquer efficacement le prix des airbnb. 
# Aussi la variance dans les résidus semble plus grande pour les airbnb avec des 
# prix plus faibles.


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Comparaison finale de tous les models ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Pour pouvoir comparer tous nos modèles, nous devons nous assurer qu'ils ont tous
# bien été ajusté sur le même jeu de données

GLM_log <- glm(log_price ~ 
                 private + Free_street_parking + Garden_or_backyard + 
                 bedrooms + host_is_superhost +  
                 number_of_reviews + review_scores_rating + has_metro_500m + 
                 prt_veg_500m,
               family = gaussian,
               data = sub_data)



# on enregistre nos modèles et données pour la prochaine séance
save(GLM_log, SarLag_log ,SarDurbin_log, SarErr_log, CAR_log,modelGlsExp,
     sub_data, model_data2,sym_nearest5_listw,
     file="Results/ModelSeance11.rda")

### Comparaison des ajustements

# nous commençons par créer une liste avec tous nos modèles
Names <- c("GLM", "SAR-LAG" ,"SAR-DURBIN","SARERROR", "CAR" ,"GLS-Exp")
Models <- list(GLM_log, SarLag_log ,SarDurbin_log, SarErr_log, CAR_log ,modelGlsExp)

# pour chaque modèle nous calculons sont AIC
AICs <- sapply(Models,FUN = AIC)

# pour chaque modèle, nous calculons le I de Moran sur ses résidus
ResidTypes <- c("pearson","pearson","pearson","pearson","pearson","normalized") # calcul des I de moran
MoranIs <- sapply(1:length(ResidTypes),FUN = function(i){
  Resid <- residuals(Models[[i]],type=ResidTypes[[i]])
  MoranI <- moran(Resid,sym_nearest5_listw,n=nrow(sub_data),Szero(sym_nearest5_listw))
  return(MoranI$I)
})



# nous calculons pour chaque modèle le RMSE, soit l'erreur quadatrique moyenne.
Y <- sub_data$price

pred_glm <- exp(predict(GLM_log))
pred_sarlag <- exp(fitted(SarLag_log))
pred_sardurbin <- exp(fitted(SarDurbin_log))
pred_sarerror <- exp(fitted(SarErr_log))
pred_CAR <- exp(CAR_log$fit$fitted.values)
pred_gls <- exp(predict(modelGlsExp))

RMSEs <- c(
  sqrt(mean((Y-pred_glm)**2)),
  sqrt(mean((Y-pred_sarlag)**2)),
  sqrt(mean((Y-pred_sardurbin)**2)),
  sqrt(mean((Y-pred_sarerror)**2)),
  sqrt(mean((Y-pred_CAR)**2)),
  sqrt(mean((Y-pred_gls)**2))
)



DataPlot1 <- data.frame(AIC=AICs,
                        MoranI=MoranIs,
                        Model=Names,
                        RMSE=RMSEs)

ggplot(DataPlot1)+
  geom_point(aes(x=AIC,y=MoranI,size=RMSEs, color = Model))
#bon les modeles CAR et GLM sont definitivement les moins performants ici

# Il semblerait que le meilleur choix soit clairement le SARLAG-Durbin qui combine à la fois
# le meilleur RMSE, AIC et a une autocorrélation spatiale résiduelle négligeable


### Pour finir notre analyse, nous allons comparer les coefficients obtenus dans 
### les différents modèles. Pour faire cela, nous allons calculer des intervalles
### de confiance basés sur les erreurs standards des coefficients. Nous utiliserons
### une approximation simple selon laquelle l'incertitude autours des coefficients
### suit une distribution de student. Une approche plus robuste serait d'effectuer
### une analyse de type boostraping. Cependant ce type d'approche nécessite plus
### de temps de calcul. Vous pouvez me poser une question à ce sujet si vous voulez
### plus de détails.
### NB : cette approche par loi de student n'est valide que pour un N>300

### Comparaison des coefficients !
CoefNames <- names(coef(GLM_log))

# on choisit des intervalles de confiance à 95%
conf_int <- 0.975

Confs <- lapply(1:length(Models),FUN = function(i){
  # pour chaque modèle, nous allons extraire sa table de coefficients
  x <- Models[[i]]
  Name <- Names[[i]]
  if(Name == "GLM"){
    Table <- summary(x)$coef
  } else if (Name %in% c("SAR-LAG" ,"SAR-DURBIN","SARERROR", "CAR")){
    Table <- summary(x)$Coef
  } else if (Name == "GLS-Exp"){
   Table <-  summary(x)$tTable
  }
  
  # on calcule ensuite le nombre de degré de liberté du modèle
  # NB : il s'agit d'une approximation, mais au delà de 30, la différence n'a 
  # qu'un impact anecdotique
  df <- nrow(sub_data)-length(Table[,1])
  
  # on peut ensuite calculer l'intervalle de confiance
  lower <- Table[,1] - qt(p = conf_int,df = df)*Table[,2]
  upper <- Table[,1] + qt(p = conf_int, df = df)*Table[,2]
  
  # on créé pour finir un beau dataframe
  Intervals <- data.frame(cbind(Table[,1],lower,upper))
  names(Intervals) <- c("estimate","lower","upper")
  Intervals$Param <- names(Table[,1])
  Intervals$Model <- Name
  
  return(subset(Intervals,Intervals$Param %in% CoefNames))
})

AllConfs <- do.call(rbind,Confs)

AllConfs$Param <- as.factor(AllConfs$Param)

ggplot(data=AllConfs)+
  geom_errorbar(aes(x=Model,ymin=lower,ymax=upper),color="black",width=0.2)+
  geom_point(aes(x=Model,y=estimate),color="blue")+
  geom_hline(yintercept = 0,color="red",linetype="dashed")+
  facet_wrap(vars(Param), scales = 'free', ncol = 3)

# Dans la figure ci-dessus, chaquef facette du graphique représente un coefficient.
# La ligne rouge représente l'absence d'effet (b = 0) et les lignes noires représentent
# l'intervalle de confiance de chaque coefficient. Si une ligne noire dépasse la ligne 
# rouge, on peut conclure que le coefficient n'est pas significativement différent de 0.
# Si nos différents modèles arrivent au mêmes conclusions, ceci devrait nous conforter dans 
# l'interprétation de nos résultats.


# interprétation finale ?