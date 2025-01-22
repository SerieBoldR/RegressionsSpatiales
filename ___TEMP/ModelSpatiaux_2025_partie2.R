# Code pour l'analyse par regression spatiale dans R
# Auteurs : J?r?my Gelb et Philippe Apparicio

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### installation et chargement des packages a utiliser ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(foreign)		## Read Data Stored by Minitab, S, SAS, SPSS, Stata, Systat, dBase.
library(nortest)      	## Tests de normalit? supp.
library(spdep)		## Spatial dependence: weighting schemes, statistics and models
library(spatialreg) ## pour les modeles SAR et CAR et SEVM
library(mgcv) #pour les modeles GAM
library(spgwr)
library(car)            ## options r?gression (pour VIF)
library(GWmodel)       	## GWR
library(moments)
library(tidyverse)
library(rcompanion) # pour calculer de r2 sur les modeles GLM
library(RColorBrewer)
library(classInt)
library(nlme) #pour les modele GLS
library(ggpubr)
library(adespatial) # pour visualiser les MEM
library(terra) # pour manipuler des raster
library(tmap)

rm(list=ls())

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

#chargement des modeles et des données de la seance precedente
load("Data_2025/ModelSeance11.rda")

MTL <- st_read('Data_2025/ILE_MTL.shp')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Modeles GAM ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

XY <- st_coordinates(sub_data)
nearest5 <- knearneigh(XY,k=5)
nearest5_nb <- knn2nb(nearest5, row.names = NULL, sym = FALSE)
nearest5_listw <- nb2listw(nearest5_nb, style = "W")
sym_nearest5_nb <- make.sym.nb(nearest5_nb)
sym_nearest5_listw <- nb2listw(sym_nearest5_nb, style = "B")


# nous allons réaliser un premier modèle GAM utilisant les coordonnées géographiques
# de nos airbnb pour estimer son terme spatial.

GAM_log <- gam(log_price ~ 
                private + Free_street_parking + Garden_or_backyard + 
                bedrooms + host_is_superhost +  
                number_of_reviews + review_scores_rating + has_metro_500m + 
                prt_veg_500m +
                s(X,Y),
              data = sub_data,
              family = gaussian)

summary(GAM_log) # sans suprise, la spline spatiale est significative

# On peut valider les résidus
TestNormality(residuals(GAM_log,type="pearson"))
moran.mc(residuals(GAM_log,type="response"),nearest5_listw,nsim = 999, alternative = 'greater')


####
##observons la forme que prends cette spline (resolution 250 m)
####

## application de la methode des effets marginaux
## on va predire des valeurs avec notre model a different
## endroits de l'espace tout en maintenant tous les autres
## parametres avec des valeurs fixes. On va ainsi voir se
## dessiner l'effet spatial, toute chose etant egale par ailleurs

## Etape 1 : creer une grille de prediction (resolution 250m)
Xcoords <- seq(min(sub_data$X-1000),max(sub_data$X+1000),by=250)
Ycoords <- seq(min(sub_data$Y-1000),max(sub_data$Y+1000),by=250)
PredDF <- expand.grid(X = Xcoords, Y = Ycoords)

PredDF$private <- 'Entier'
PredDF$Free_street_parking <- 'YES'
PredDF$Garden_or_backyard <- 'NO'
PredDF$bedrooms <- 2
PredDF$host_is_superhost <- 'YES'
PredDF$number_of_reviews <- 10
PredDF$review_scores_rating <- 100
PredDF$has_metro_500m <- 'NO'
PredDF$prt_veg_500m <- 35

##Etape 3 : effectuer la prediction
PredDF$pred_log_price_gam <- predict(GAM_log,newdata=PredDF)

##Etape 4 :centrer la prediction 
#(pour ne voir que l'effet spatial debarasse de la constante des
#autre parametres)
PredDF$CenterPred <- PredDF$pred_log_price_gam - mean(PredDF$pred_log_price_gam)

###Etape 5 : construire un raster
rasterGAM <- rast(as.matrix(PredDF[, c("X", "Y", "CenterPred")]), type = 'xyz')
crs(rasterGAM) <- crs(sub_data)

tm_shape(mask(rasterGAM,MTL)) + 
  tm_raster(palette = '-Spectral', n = 8, style = 'jenks') + 
  tm_shape(MTL) + 
  tm_borders('black') + 
  tm_layout(legend.outside = T)
  
# Cette carte se présente pour le moment sous la forme d'une prédiction sur l'échelle log.
# Nous sommes dans cette échelle car nous avons transformé la variable prix avec la fonction
# log dans chacun de nos modèles. Nous pouvons prendre cette valeur et la convertir avec 
# la fonction exp. En faisant aisni, nous verrons l'effet multiplicatif de l'espace sur le prix
# directement dans son échelle originale ($)

rasterGAM_dollar <- exp(rasterGAM)

tm_shape(mask(rasterGAM_dollar,MTL)) + 
  tm_raster(palette = '-Spectral', n = 8, style = 'jenks') + 
  tm_shape(MTL) + 
  tm_borders('black') + 
  tm_layout(legend.outside = T)

# Dans les secteurs en rouge, tout chose étant égale par ailleur, on s'attend à voir 
# le prix des airbnb entre 30 et 50% plus cher.


####
##Forcons la spline a avoir plus de noeuds !
## on aura ainsi une spline mieux ajustee localement
####

GAM_log_150 <- gam(log_price ~ 
                   private + Free_street_parking + Garden_or_backyard + 
                   bedrooms + host_is_superhost +  
                   number_of_reviews + review_scores_rating + has_metro_500m + 
                   prt_veg_500m +
                   s(X,Y, k = 150),
                 data = sub_data,
                 family = gaussian)

summary(GAM_log_150)
# on voit que le modèle utilise environs 55 degrés de libertés pour ajuster son terme spatial.

# Les résidus se comportent toujours relativement bien
moran.mc(residuals(GAM_log_150,type="response"),nearest5_listw,nsim = 999)
TestNormality(residuals(GAM_log_150,type="pearson")) 


##Etape 3 : effectuer la prediction
PredDF$pred_log_price_gam_150 <- predict(GAM_log_150,newdata=PredDF)

##Etape 4 :centrer la prediction
PredDF$CenterPred_150 <- PredDF$pred_log_price_gam_150 - mean(PredDF$pred_log_price_gam_150)

###Etape 5 : construire un raster
rasterGAM_150 <- rast(as.matrix(PredDF[, c("X", "Y", "CenterPred_150")]), type = 'xyz')
crs(rasterGAM_150) <- crs(sub_data)
rasterGAM_150 <- exp(rasterGAM_150)

tm_shape(mask(rasterGAM_150,MTL)) + 
  tm_raster(palette = '-Spectral', n = 8, style = 'jenks') + 
  tm_shape(MTL) + 
  tm_borders('black') + 
  tm_layout(legend.outside = T)

#' Le terme spatial ainsi obtenu est un peu plus précis. Nous devons cependant nous
#' assurer que l'augmentation de la complexité de ce terme spatial en valait la peine 
#' en comparant les AIC des deux modèles

AIC(GAM_log_150,GAM_log) # il semblerait que oui car nous obtenons un meilleur AIC avec 
# le second modèle


#verifions l'homeoscedasticit?
DF <- data.frame(Pred=predict(GAM_log_150),
                 Resid=residuals(GAM_log_150))
ggplot(DF)+
  geom_point(aes(x=Pred,y=Resid))+
  geom_smooth(aes(x=Pred,y=Resid))
# pas mal ici !


####
## tentons une spline avec un random markov field
####

# nos données se présentent sous forme de point. Cependant, nous pouvons utiliser une grille
# plus ou moins grossière sur le territoire et ajuster un Markov Random Field dessus.
# Pour cet exemple, nous allons utiliser les secteurs de recensement, mais une tesselation 
# régulière, une mesh ou des polygones de voronois auraient aussi pu être utilisés

SR <- st_read('Data_2025/SR_2021_rmr_mtl.gpkg')

# nous ne devons garder que les SR dans lesquels nous avons au moins 1 airbnb
ok_SR <- subset(SR,
                lengths(st_intersects(SR, sub_data)) > 0
                )

# nous attributons ensuite un ID spécial à ces SR
ok_SR$GEOID <- as.factor(1:nrow(ok_SR))

# et déterminer pour chaque airbnb dans quel SR il tombe
sub_data2 <- st_join(sub_data, ok_SR)

# on a ici un airbnb qui ne tombe pas dans un SR... On va le forcer à tomber dans 
# le SR le plus proche
tmap_mode('view')
error <- subset(sub_data2, is.na(sub_data2$GEOID))
tm_shape(ok_SR) + tm_fill('white') + tm_borders('black') + tm_shape(error) + tm_dots('black')

# en effet ! Ce petit malin est tombé dans le canal... on devrait lui associer le SR GEOID = 250
sub_data2$GEOID <- factor(ifelse(is.na(sub_data2$GEOID), 
                          250, sub_data2$GEOID))


nb <- poly2nb(ok_SR, row.names = ok_SR$GEOID)
names(nb) <- attr(nb, "region.id")

GAMN_log_rmf <- gam(log_price ~ 
                      private + Free_street_parking + Garden_or_backyard + 
                      bedrooms + host_is_superhost +  
                      number_of_reviews + review_scores_rating + has_metro_500m + 
                      prt_veg_500m+
                  s(GEOID, bs = 'mrf', xt = list(nb = nb),k=10)
                ,data = sub_data2,family = gaussian)

summary(GAMN_log_rmf)

# Nous pouvons visualiser l'effet obtenu, pour cela nous devons une fois encore
# effectuer une prédiction "toutes choses étant égales par ailleurs"
df_pred <- data.frame(
  private = 'Entier',
  Free_street_parking ='YES',
  Garden_or_backyard = 'NO',
  bedrooms = 2,
  host_is_superhost = 'YES',
  number_of_reviews = 10,
  review_scores_rating = 100,
  has_metro_500m = 'NO',
  prt_veg_500m = 35,
  GEOID = unique(ok_SR$GEOID)
)

# puis retirer la moyenne de nos prédictions et retransformer dans l'échelle linéaire
df_pred$prediction <- predict(GAMN_log_rmf, df_pred)
df_pred$prediction <- exp(df_pred$prediction - mean(df_pred$prediction))

# on peut ensuite joindre le résultats à notre mesh et visualiser le résultat
df_pred2 <- left_join(ok_SR, df_pred, by = 'GEOID')

tm_shape(df_pred2) + 
  tm_fill('prediction', n = 7, style = 'jenks')

# Pour rappel, dans cette première version, nous avons seulement autorisé 10 degrés 
# de liberté pour le terme spatial


##PB : besoin de trouver K : option : tester toutes les valeurs possible

## on va tester des valeur de K entre 25 et 200 en faisant des bonds de 25
Ks <- seq(25,200,by=25)

AssessmentValues <- lapply(Ks,FUN = function(x){
  
  ## calculons le model
  print(paste("Calculating model for k = ",x,sep=""))
  Model <- gam(log_price ~ 
                private + Free_street_parking + Garden_or_backyard + 
                bedrooms + host_is_superhost +  
                number_of_reviews + review_scores_rating + has_metro_500m + 
                prt_veg_500m+
                s(GEOID, bs = 'mrf', xt = list(nb = nb),k=x)
              ,data = sub_data2,family = gaussian)
  
  df_pred$prediction <- predict(Model, df_pred)
  df_pred$prediction <- exp(df_pred$prediction - mean(df_pred$prediction))
  df_pred2 <- left_join(ok_SR, df_pred, by = 'GEOID')
  
  print(tm_shape(df_pred2) + 
    tm_fill("prediction", palette = "-Spectral", n = 7, style = 'jenks', midpoint = 1) + 
    tm_layout(paste("prix des airbnb tendance spatiale avec k = ",x,sep="")))
  
  
  Values <- c(
    "k"=x,
    "aic" = AIC(Model),
    "moranI" = moran(residuals(Model),nearest5_listw,n=nrow(sub_data2),Szero(nearest5_listw))$I,
    "R2" = summary(Model)$r.sq
  )
  return(Values)
})

#convertissons ces indicateurs en un dataframe
GRMFValues <- data.frame(do.call(rbind,AssessmentValues))

##looking at the AIC
ggplot(GRMFValues)+
  geom_line(aes(x=k,y=aic))+
  geom_point(aes(x=k,y=aic),color="red")

##looking at the Moran I
ggplot(GRMFValues)+
  geom_line(aes(x=k,y=moranI))+
  geom_point(aes(x=k,y=moranI),color="red")

#clairement un resultat interessant pour k = 150

GamGRMF <- gam(log_price ~ 
                 private + Free_street_parking + Garden_or_backyard + 
                 bedrooms + host_is_superhost +  
                 number_of_reviews + review_scores_rating + has_metro_500m + 
                 prt_veg_500m+
                 s(GEOID, bs = 'mrf', xt = list(nb = nb),k=150)
               ,data = sub_data2,family = gaussian)

AIC(GamGRMF,GAM_log, GAM_log_150)



#verifions l'homeoscedasticit?
DF <- data.frame(Pred=predict(GamGRMF),
                 Resid=residuals(GamGRMF))
ggplot(DF)+
  geom_point(aes(x=Pred,y=Resid))+
  geom_smooth(aes(x=Pred,y=Resid))
# tres bon ici !

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Modeles SEVM ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# NB : nous montrons ici un exemple simple basé sur notre matrice spatiale des 5 plus proches voisins.
# pour des données d'évènement ponctuel, je recommande de regarder le package spmoran 
# qui fournit notamment des fonctions d'interpolation (meigen0) pour prédire les vecteurs 
# spatiaux sur de nouvelles observations. Il propose aussi des méthdodes de calculs plus 
# adaptée à de gros jeux de données.

#calcul des vecteurs spatiaux
Vectors <- mem(nearest5_listw)
Vectors

#affichage des vecteurs spatiaux
barplot(attr(Vectors, "values"), 
        main = "Eigenvalues of the spatial weighting matrix", cex.main = 0.7)


#faisons une carte du 1er vecteur
sub_data$MEM1 <- Vectors$MEM1

# ce premier vecteur spatial correspond au patron spatial avec la plus haute valeur de I de Moran
# possible
tm_shape(sub_data[order(sub_data$MEM1),]) + 
  tm_dots("MEM1", palette = "-Spectral", n = 7, style = 'jenks')

moran.test(sub_data$MEM1,nearest5_listw)

#faisons une carte du 25 vecteur
sub_data$MEM25 <- Vectors$MEM25

tm_shape(sub_data[order(sub_data$MEM25),]) + 
  tm_dots("MEM25", palette = "-Spectral", n = 7, style = 'fisher')

moran.test(sub_data$MEM25,nearest5_listw)

# comme promis ces vecteur sont orthogonaux (pas de correlation entre eux)
cor.test(sub_data$MEM1,sub_data$MEM25)


#faisons une carte du MEM1500 vecteur
sub_data$MEM1500 <- Vectors$MEM1500

tm_shape(sub_data[order(sub_data$MEM1500),]) + 
  tm_dots("MEM1500", palette = "-Spectral", n = 7, style = 'fisher')

moran.test(sub_data$MEM1500,nearest5_listw)


# Maintenant que nous avons joué un peu avec ces vecteurs, nous pouvons les intégrer
# dans notre modèle grâce à la fonction de sélection ME

SEVM <- ME(log_price ~ 
             private + Free_street_parking + Garden_or_backyard + 
             bedrooms + host_is_superhost +  
             number_of_reviews + review_scores_rating + has_metro_500m + 
             prt_veg_500m,
           data = sub_data,
           family = "gaussian",
           nsim=999,
           listw=nearest5_listw, alpha=0.05, stdev=TRUE, verbose=TRUE)

SEVM # on a donc retenu 28 eigen vectors

glmSEVM_log <- glm(log_price ~ 
                     private + Free_street_parking + Garden_or_backyard + 
                     bedrooms + host_is_superhost +  
                     number_of_reviews + review_scores_rating + has_metro_500m + 
                     prt_veg_500m + fitted(SEVM),
               family = "gaussian", data=sub_data)

summary(glmSEVM_log)

# on valide comme toujours les résidus et on regarde le R2
TestNormality(residuals(glmSEVM_log)) # jolis residus
moran.mc(residuals(glmSEVM_log,type="response"),nearest5_listw,nsim = 999) # et plus d'autocorrelation spatiale !
nagelkerke(glmSEVM_log)

## faisons une carte du predicteur spatial
#preparons le jeu de donnees pour la prediction (marginal effect)

df_pred <- data.frame(
  private = 'Entier',
  Free_street_parking ='YES',
  Garden_or_backyard = 'NO',
  bedrooms = 2,
  host_is_superhost = 'YES',
  number_of_reviews = 10,
  review_scores_rating = 100,
  has_metro_500m = 'NO',
  prt_veg_500m = rep(35, nrow(sub_data))
)

df_pred <- cbind(df_pred, fitted(SEVM))


#effectuons la prediction et retirons la moyenne
sub_data$SEVMPred <- predict(glmSEVM_log,newdata=df_pred)
sub_data$SEVMPred <- exp(sub_data$SEVMPred - mean(sub_data$SEVMPred))

## Nous pouvons ensuite afficher le patron spatial capturé par le modèle

map1 <- tm_shape(sub_data[order(sub_data$SEVMPred),]) + 
  tm_dots('SEVMPred', n = 7, style = 'fisher', 
          palette = 'RdBu'
  )

map2 <- tm_shape(sub_data[order(-1 * sub_data$SEVMPred),]) + 
  tm_dots('SEVMPred', n = 7, style = 'fisher', 
          palette = 'RdBu'
  )

tmap_arrange(map1, map2, sync = TRUE)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Modeles GWR ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(spgwr)
###
#Etape 1 : trouver une bandwidth optimale
###

# nous allons ici utiliser une bandwidth adaptative, 
# avec une fonction de pondération biquarée. La bandwidth sera optimisée
# par cross-validation 

bwaCV.voisins  <- gwr.sel(log_price ~ 
                            private + Free_street_parking + Garden_or_backyard + 
                            bedrooms + host_is_superhost +  
                            number_of_reviews + review_scores_rating + has_metro_500m + 
                            prt_veg_500m ,
                          data = sub_data,
                          method = "cv",          # Méthode cv ou AIC
                          gweight=gwr.bisquare,   # gwr.gauss ou gwr.bisquare
                          adapt=TRUE,
                          verbose = TRUE,
                          RMSE = TRUE,
                          longlat = FALSE,
                          coords=cbind(sub_data$X,sub_data$Y))


bwaCV.voisins
# on voit ainsi que le modèle recommande d'utiliser les 20% d'observation les 
# plus proches pour faire les ajustement locaux.


###
#Etape 2 : calculer la GWR
###

GWR1 <- gwr(log_price ~ 
              private + Free_street_parking + Garden_or_backyard + 
              bedrooms + host_is_superhost +  
              number_of_reviews + review_scores_rating + has_metro_500m + 
              prt_veg_500m ,
            data = sub_data,
            adapt=bwaCV.voisins,
            gweight=gwr.bisquare,
            hatmatrix=TRUE,
            se.fit=TRUE,
            coords=cbind(sub_data$X,sub_data$Y),
            longlat=F)


###
#Etape 3 : Validation de la pertinence d'utiliser une GWR
###


# F1 : test de la qualite d'ajustement
# : repond a la question : est-ce que la GWR est significativement mieux
#ajustee au donnees que le model OLS de base
#en se basant sur les RSS (residual sum of squares)
# dans notre cas : oui (p<0.01)

test1 <- LMZ.F1GWR.test(GWR1)
print(test1)

# F2 : test de la qualite d'ajustement
# : repond a la question : est-ce que la GWR est significativement mieux
#ajustee au donnees que le model OLS de base
#en se basant sur une analyse de la variance
# dans notre cas oui : (p<0.01)

test2 <- LMZ.F2GWR.test(GWR1)
print(test2)

# F3 : test de la variabilite spatiale des coefficients
# : repond a la question : est-ce que mes predicteurs varient spatialement de facon significative ?
# dans notre cas, ils semblent tous varier spatialement (p<0.05)
# sauf peut etre Pct_Etrang (p<0.025)

test3 <- LMZ.F3GWR.test(GWR1)
print(test3)


# pour regarder les résidus du modèle, nous devons ici les calculer manuellement

rez <- GWR1$SDF$pred - sub_data$log_price

TestNormality(rez) # beau residus au global !

moran.mc(rez,nearest5_listw,nsim = 999) # et plus d'autocorrelation spatiale !


###
#Etape 4 : Parcours des résultats
###

#' Dans cette dernière partie, nous allons parcourir les résultats de notre GWR
#' en cartographiant ses différentes composantes.

# Extraction de l'objet résultat sous forme d'un objet sf
gw_results <- st_as_sf(GWR1$SDF)
st_crs(gw_results) <- st_crs(sub_data)

#_______________________________________
# Analyse des R2

summary(gw_results$localR2)
quantile(gw_results$localR2, probs = c(0.025, 0.975))
# Pour 95% de nos observations, Le R2 se situe entre 0.30 et 0.66


## Nous pouvons ensuite afficher le patron spatial capturé par le modèle. Pour rappel
## le R2 mesure la part de la variance de la variable Y qui a été capturée par le modèle

tm_shape(gw_results[order(gw_results$localR2),]) + 
  tm_dots('localR2', n = 7, style = 'fisher', 
          palette = 'YlOrBr'
)

#_______________________________________
# Analyse des coefficients

#' Nous pouvons ensuite cartographier les coefficients de notre modèle pour 
#' voir lesquels ont des impacts qui changent effectivement dans l'espace.
#' Avec le code ci-dessous, nous réalisons une boucle qui va créer une carte par
#' coefficient et stocker toutes ces cartes dans une liste.

vars <- c("privateEntier", "Free_street_parkingYES", "Garden_or_backyardYES",
          "bedrooms", "host_is_superhostYES", "number_of_reviews", "review_scores_rating",
          "has_metro_500mYES", 'prt_veg_500m')


coeff_maps <- lapply(vars, function(col){
  
  # convertion en exponentiel des coefficients du fait de la fonctions log
  gw_results$x <- exp(gw_results[[col]])
  
  # les valeurs de t sont obtenues en divisant le coefficient par son erreur type
  t_val <- gw_results[[col]] / gw_results[[paste0(col,'_se')]]
  
  # identification des observations avec une valeur significative au seuil 0.05
  pvals <- round(2 * (1 - pt(abs(t_val), GWR1$results$edf)), 3)
  
  sign_gw <- subset(gw_results, pvals < 0.05)
  non_sign_gw <- subset(gw_results, pvals >= 0.05)
  
  if(nrow(non_sign_gw) > 0){
    map <- tm_shape(sign_gw) + 
      tm_dots("x", palette = "-Spectral", style = 'jenks', n = 7, title = col, midpoint = 1) + 
      tm_shape(non_sign_gw) + 
      tm_dots("grey") + 
      tm_layout(legend.outside = TRUE)
  }else{
    map <- tm_shape(sign_gw) + 
      tm_dots("x", palette = "-Spectral", style = 'jenks', n = 7, title = col, midpoint = 1) + 
      tm_layout(legend.outside = TRUE)
  }
  
  
  return(map)
})

names(coeff_maps) <- vars

#' on peut comibner toutes les carte d'un coup, 
#' ou alors les afficher les unes après les autre
tmap_arrange(coeff_maps, ncol = 2)


coeff_maps$privateEntier
coeff_maps$Garden_or_backyardYES

#' Il est aussi pertinent d'extraire la distribution de chacun des coefficients
#' et de les présenter dans un tableau de résumé : 

qt025 <- function(x){quantile(x, probs = c(0.025))}
qt50 <- function(x){quantile(x, probs = c(0.50))}
qt975 <- function(x){quantile(x, probs = c(0.975))}


tableau_coeff <- gw_results %>% 
  st_drop_geometry() %>%  # on retire les géométries
  select(any_of(vars)) %>% # on garde seulement les colonnes des coefficients
  mutate_all(exp)%>% # on les passe tous en exponentiel pour avoir les effets multiplicatifs (modèle log)
  summarise_all(
    .funs = list(qt025, qt50, qt975)
  ) %>% 
  pivot_longer(cols = everything(), names_sep = '_fn', names_to = c('coeff', 'qtl')) %>% 
  pivot_wider(names_from = 'qtl', values_from = 'value')

names(tableau_coeff) <- c('coefficient (exp)', 'lb', 'median', 'up')


#' Finalement, il est aussi souvent intéressant de représenter la variable qui 
#' a localement le plus d'impact sur la prédiction. Pour faire cela, nous devons 
#' calculer pour chaque observation une nouvelle colonne étant le produit entre
#' la valeur de cette variable localement et son coefficient local. On appelle
#' ces valeurs les contributions de chaque terme du modèle.
#' on doit donc se créer une extraction des données de base pour calculer ces
#' contributions

contrib_df <- data.frame(
  privateEntier = ifelse(sub_data$private == 'Entier', 1, 0),
  Free_street_parkingYES = ifelse(sub_data$Free_street_parking == 'YES', 1, 0),
  Garden_or_backyardYES = ifelse(sub_data$Garden_or_backyard == 'YES', 1, 0),
  bedrooms = sub_data$bedrooms,
  host_is_superhostYES = ifelse(sub_data$host_is_superhost == 'YES', 1, 0),
  number_of_reviews = sub_data$number_of_reviews,
  review_scores_rating = sub_data$review_scores_rating,
  has_metro_500mYES = ifelse(sub_data$has_metro_500m == 'YES', 1, 0),
  prt_veg_500m = sub_data$prt_veg_500m
)

# on peut maintenant calculer les contributions
for(var in vars){
  gw_results[paste0('contrib_',var)] <- contrib_df[[var]] * gw_results[[var]]
}

# il ne nous reste plus qu'à trouver la colonne avec le plus gros effet
Vars <- names(gw_results)[grepl("contrib",names(gw_results),fixed=TRUE)]


#' on va ensuite choisir pour chaque ligne de notre jeu de donnée, laquelle de 
#' ces colonnes à la valeur la plus élevée
AbsDF <- abs(st_drop_geometry(gw_results[Vars]))

gw_results$highest_contrib <- as.factor(Vars[apply(AbsDF[Vars],1,which.max)])


tm_shape(gw_results) + 
  tm_dots('highest_contrib',
          title = 'variable avec la plus grande contribution',
          size = 0.05
          ) 


## Export de la GWR pour de la belle carto !
st_write(gw_results,"Results/GWR1.gpkg",layer="GWR1",driver="GPKG")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Comparaison finale de tous les models ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Dans cette dernière section, nous allons comparer de façon assez générale les 
#' différents modèles que nous avons ajustés pendant les deux dernières séances. 
#' L'objectif est de comparer les qualités d'ajustement des modèles et est ce que
#' ces derniers nous amènent vers les mêmes conclusions ou non.

### Comparaison des ajustements
Names <- c("GLM", "SAR-LAG" ,"SAR-DURBIN","SARERROR", "CAR" ,"GLS-Exp","GAM","GAM_rmf","SEVM","GWR")
Models <- list(GLM_NO2, SarLagNO2 ,SarDurbinNO2Comp, SarErrNO2, CARNO2 ,modelGlsExp,GAMNo2,GamGRMF,glmSEVM,GWR1)
AICs <- sapply(1:length(Names),FUN = function(x){ #calcul des AIC
  if (Names[[x]]=="GWR"){
    return(Models[[x]]$GW.diagnostic$AICc)
  }else{
    return(AIC(Models[[x]]))
  }
}) 

ResidTypes <- c("pearson","pearson","pearson","pearson","pearson","normalized","pearson","pearson","pearson","GWR") # calcul des I de moran
MoranIs <- sapply(1:length(ResidTypes),FUN = function(i){
  if (ResidTypes[[i]]=="GWR"){
    Resid <- Models[[i]]$SDF$Stud_residual
  }else {
    Resid <- residuals(Models[[i]],type=ResidTypes[[i]])
  }
  MoranI <- moran(Resid,Queen_W,n=nrow(LyonIris),Szero(Queen_W))
  return(MoranI$I)
  
})

RMSEs <- sapply(1:length(Names),FUN=function(i){ # calcul des Root mean square error (quantite absolue de residus)
  if(Names[[i]]=="GWR"){
    Resid <- Models[[i]]$SDF$residual
  } else{
    Resid <- residuals(Models[[i]],type="response")
  }
  return(sqrt(mean(Resid**2)))
})

DataPlot1 <- data.frame(AIC=AICs,
                        MoranI=MoranIs,
                        Model=Names,
                        RMSE=RMSEs)

ggplot(DataPlot1)+
  geom_point(aes(x=AIC,y=MoranI,size=RMSEs))+
  geom_label(aes(x=AIC,y=MoranI,label=Model),position=position_dodge(0.5))
#bon les modeles CAR, GAM et GLM sont definitivement les moins performants ici


ggplot(DataPlot1)+
  geom_point(aes(x=AIC,y=MoranI,size=RMSEs,color=Model))+
  geom_text(aes(x=AIC,y=MoranI,label=Model),nudge_y = 0.005)+
  ylim(c(-0.1,0.1)) + 
  xlim(-800,-650)


### Comparaison des coefficients !
CoefNames <- names(coef(GLM_NO2))

conf_int <- 0.975

## NB : cette approche n'est valide que pour un N>300
Confs <- lapply(1:length(Models),FUN = function(i){
  x <- Models[[i]]
  Name <- Names[[i]]
  print(paste("doing model : ",Name,sep=""))
  if(Name != "GWR") {
    if(Name %in% c("GLM","SEVM")){
      Table <- summary(x)$coef
      Table <- Table[1:length(CoefNames),]
    } else if (Name %in% c("SAR-LAG" ,"SAR-DURBIN","SARERROR", "CAR")){
      Table <- summary(x)$Coef
      Table <- Table[1:length(CoefNames),]
    } else if(Name %in% c("GAM","GAM_rmf")){
      Summ <- summary(x)
      Table <- cbind(Summ$p.coeff,Summ$se[1:length(Summ$p.coeff)])
    } else {
      Table <-  summary(x)$tTable
    }
    
    lower <- Table[,1] - qt(conf_int,nrow(LyonIris)-length(Table[,1]))*Table[,2]   #calcul de l'interval ? 95%
    upper <- Table[,1] + qt(conf_int,nrow(LyonIris)-length(Table[,1]))*Table[,2]
    estimate <- Table[,1]
    
  } else{
    GWRDF <- x$SDF
    TempCoeffName <- CoefNames
    TempCoeffName[[1]] <- "Intercept"
    lower <- sapply(TempCoeffName,FUN=function(Term){
      return(quantile(GWRDF[[Term]],probs=c(0.025)))
    })
    upper <-  sapply(TempCoeffName,FUN=function(Term){
      return(quantile(GWRDF[[Term]],probs=c(0.975)))
    })
    estimate <-sapply(TempCoeffName,FUN=function(Term){
      return(quantile(GWRDF[[Term]],probs=c(0.5)))
    })
    
  }
  
  Intervals <- data.frame(cbind(estimate,lower,upper))
  names(Intervals) <- c("estimate","lower","upper")
  Intervals$Param <- CoefNames
  Intervals$Model <- Name
  
  return(subset(Intervals,Intervals$Param %in% CoefNames))
})

AllConfs <- do.call(rbind,Confs)

Plots <- list()
for(Name in CoefNames){
  Sub <- subset(AllConfs,AllConfs$Param==Name)
  Plot <- ggplot(data=Sub)+
    geom_errorbar(aes(x=Model,ymin=lower,ymax=upper),color="black",width=0.2)+
    geom_point(aes(x=Model,y=estimate),color="blue")+
    geom_hline(yintercept = 0,color="red",linetype="dashed")+
    ggtitle(Name)
  Plots[[Name]]<-Plot
}


#permet de comparer les degr?s d'incertitude des mod?les sur les diff?rents param?tres
#et de se faire une id?e g?n?rale des directions des diff?rentes variables
ggarrange(plotlist = Plots,ncol=2,nrow=3)

