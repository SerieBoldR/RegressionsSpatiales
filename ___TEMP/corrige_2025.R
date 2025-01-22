# Code pour l'analyse par regression spatiale dans R
# Auteurs : Jeremy Gelb

# correction pour le TP4

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
library(performance) # pour calculer de r2 sur les modeles GLM
library(RColorBrewer)
library(classInt)
library(nlme) #pour les modele GLS
library(ggpubr)
library(adespatial) # pour visualiser les MEM
library(terra) # pour manipuler des raster
library(tmap)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Chargement des donnees ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setwd("C:/Users/APPP2302/OneDrive - USherbrooke/Articles Publies ne pas supprimer/_Livres/BolR_RegressionsSpatiales/___TEMP")

data_sr <- st_read('Data_2025/data_sr_access.gpkg')

# structuration des variables de recensement
data_sr$prt_auto <- data_sr$mode_auto  / data_sr$total_commuters
data_sr$prt_durable <- (data_sr$mode_tc + data_sr$mode_pieton + data_sr$mode_velo)  / data_sr$total_commuters
data_sr$ratio <- data_sr$prt_durable / data_sr$prt_auto
data_sr$men_size <- data_sr$Population / data_sr$Households


TestNormality <- function(x,showHelp=T){
  #calcul du test shapiro.test pour evaluer la normalite
  if(length(x) > 5000){
    x2 <- sample(x, 5000)
    Test <- shapiro.test(x2)
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


# visualisation rapide de la variable dépendante
hist(data_sr$ratio, breaks = 30)
tm_shape(data_sr) + 
  tm_fill('ratio', n = 7, style = 'quantile')


# calcul du I de Moran de la variable Y
data_sr <- subset(data_sr,
                  complete.cases(st_drop_geometry(data_sr[c(
  "ratio", "prt_minorite_vis" , "prt_monoparental" , "prt_chomage" , "revenu_median" , 
  "densite_population_km2" , "acs_idx_emp_tc_peak" , 'men_size',
  "acs_idx_emp_velo" , "acs_idx_emp_pieton")])))


nb <- poly2nb(data_sr, queen = T)
listw <- nb2listw(nb, style = 'W', zero.policy = T)

moran.mc(data_sr$ratio, listw = listw, nsim = 999, zero.policy = T)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### modele non spatial ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vif(glm(ratio ~ 
          prt_minorite_vis + prt_monoparental + prt_chomage + revenu_median + 
          densite_population_km2 + acs_idx_emp_tc_peak + men_size +
          acs_idx_emp_velo + acs_idx_emp_pieton, 
        data = data_sr))


glm_base <- glm(ratio ~ 
                  prt_minorite_vis + prt_monoparental + prt_chomage + revenu_median + 
                  densite_population_km2 + acs_idx_emp_tc_peak + men_size +
                  acs_idx_emp_velo + acs_idx_emp_pieton, 
                data = data_sr, 
                family = gaussian)
  

TestNormality(residuals(glm_base))


moran.mc(residuals(glm_base), listw = listw, nsim = 999)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### modeles econometriques ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#### Diagnostic des test de Lagrange sur le mod?le OLS
ols_base <- 	lm(ratio ~ 
                  prt_minorite_vis + prt_monoparental + prt_chomage + revenu_median + 
                  densite_population_km2 + acs_idx_emp_tc_peak + men_size +
                  acs_idx_emp_velo + acs_idx_emp_pieton, 
                data = data_sr)

# on effectue un test de Lagrange pour voir quel type de modèle nous devrions privilégier
test <-  lm.LMtests(ols_base, listw=listw, test=c("LMerr","RLMerr","LMlag","RLMlag"), zero.policy = TRUE)
test
# il semble que l'on doive faire un modèle SAR-LAG
SarLag <- lagsarlm (ratio ~ 
                      prt_minorite_vis + prt_monoparental + prt_chomage + revenu_median + 
                      densite_population_km2 + acs_idx_emp_tc_peak + men_size +
                      acs_idx_emp_velo + acs_idx_emp_pieton,
                    data = data_sr,
                        listw = listw,zero.policy = TRUE,
                        tol.solve=1e-15)

qqPlot(residuals(SarLag), distribution = 'norm')

TestNormality(residuals(SarLag))

moran.mc(residuals(SarLag,type="response"),listw,nsim = 999, zero.policy = TRUE)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### modeles GLS ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# pour les modèles GLS, on doit récupérer les coordonnées XY de nos observations
coords <- st_coordinates(st_point_on_surface(data_sr))
data_sr$X <- coords[,1]
data_sr$Y <- coords[,2]


#Mod?le Gls (Exponentiel)
modelGlsExp <- gls(ratio ~ 
                     prt_minorite_vis + prt_monoparental + prt_chomage + revenu_median + 
                     densite_population_km2 + acs_idx_emp_tc_peak + men_size +
                     acs_idx_emp_velo + acs_idx_emp_pieton,
                   data = data_sr,
                   correlation=corExp(form=~X+Y, nugget=T))

summary(modelGlsExp)

#Mod?le Gls (gaussian)
modelGlGaussian <- gls(ratio ~ 
                         prt_minorite_vis + prt_monoparental + prt_chomage + revenu_median + 
                         densite_population_km2 + acs_idx_emp_tc_peak + men_size +
                         acs_idx_emp_velo + acs_idx_emp_pieton,
                       data = data_sr,
                       correlation=corGaus(form=~X+Y, nugget=T))
summary(modelGlGaussian)

#Mod?le Gls (spher)
modelGlsSpher <- gls(ratio ~ 
                       prt_minorite_vis + prt_monoparental + prt_chomage + revenu_median + 
                       densite_population_km2 + acs_idx_emp_tc_peak + men_size +
                       acs_idx_emp_velo + acs_idx_emp_pieton,
                     data = data_sr,
                     correlation=corSpher(form=~X+Y, nugget=T))

#Mod?le Gls (Rationel Quadratic)
modelGlsRatio <- gls(ratio ~ 
                       prt_minorite_vis + prt_monoparental + prt_chomage + revenu_median + 
                       densite_population_km2 + acs_idx_emp_tc_peak + men_size +
                       acs_idx_emp_velo + acs_idx_emp_pieton,
                     data = data_sr,
                     correlation=corRatio(form=~X+Y, nugget=T))

save(modelGlsExp, modelGlGaussian, modelGlsSpher, modelGlsRatio, file = 'GLS_models_airbnb.rda')

#comparaison des 4 formes
AIC(modelGlsExp,modelGlGaussian,modelGlsSpher,modelGlsRatio)

## le plus petit AIC est obtenu par le modèle ratio

TestNormality(residuals(modelGlsRatio,type="normalized")) # nice residuals
moran.mc(residuals(modelGlsRatio,type="normalized"),listw,nsim=999) # autocorrelation resolue


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Modele CAR ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Queen_CAR <- nb2listw(nb, style="B", zero.policy = T)

model_CAR <-spautolm(ratio ~ 
                       prt_minorite_vis + prt_monoparental + prt_chomage + revenu_median + 
                       densite_population_km2 + acs_idx_emp_tc_peak + men_size +
                       acs_idx_emp_velo + acs_idx_emp_pieton,
                     data = data_sr,
                  listw=Queen_CAR,zero.policy=TRUE,
                  tol.solve=1e-15, 
                  family = "CAR")

summary(model_CAR,Nagelkerke=T)
TestNormality(residuals(model_CAR)) # un peu meilleure normalite
moran.mc(residuals(model_CAR),Queen_CAR,nsim=999) # ok pour l'autocorrélation spatiale

# visualisation du terme spatial dans le modèle CAR
data_sr$car_spatial_component <- model_CAR$fit$signal_stochastic

tm_shape(data_sr) + 
  tm_fill('car_spatial_component', k = 7, style = 'jenks', title = 'terme spatial',
          palette = '-RdBu'
          ) + 
  tm_compass(position = c("right", "bottom"), 
             size = 0.5)+
  tm_legend(position = c("left", "top"),
            frame = FALSE, bg.color = "white")+
  tm_scale_bar(breaks  = c(0, 10, 20),
               position = c("left", "bottom"))+
  tm_layout(main.title = "Composante spatiale du modèle CAR",
            legend.outside = TRUE, 
            attr.outside = TRUE,
            inner.margins = 0,
            frame = FALSE, 
            legend.format = list(text.separator = "-",
                                 text.or.more = 'ou plus'
            )
  ) 
  

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Modeles GAM ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_sr$GEOID <- as.factor(1:nrow(data_sr))
nb <- poly2nb(data_sr, row.names = data_sr$GEOID)
names(nb) <- attr(nb, "region.id")


model_GAM <- gam(ratio ~ 
                   prt_minorite_vis + prt_monoparental + prt_chomage + revenu_median + 
                   densite_population_km2 + acs_idx_emp_tc_peak + men_size +
                   acs_idx_emp_velo + acs_idx_emp_pieton +
                   s(GEOID, bs = 'mrf', xt = list(nb = nb),k=10),
              data = data_sr,family = gaussian)

##PB : besoin de trouver K : option : tester toutes les valeurs possible

## preparons rapidement nos donnees pour une prediction
PredDF <- data.frame(GEOID = data_sr$GEOID,
                     prt_minorite_vis=mean(data_sr$prt_minorite_vis),
                     prt_monoparental=mean(data_sr$prt_monoparental),
                     prt_chomage=mean(data_sr$prt_chomage),
                     revenu_median=mean(data_sr$revenu_median),
                     densite_population_km2=mean(data_sr$densite_population_km2),
                     acs_idx_emp_tc_peak = mean(data_sr$acs_idx_emp_tc_peak),
                     acs_idx_emp_velo = mean(data_sr$acs_idx_emp_velo),
                     acs_idx_emp_pieton = mean(data_sr$acs_idx_emp_pieton),
                     men_size = mean(data_sr$men_size)
)

Queen_W <- nb2listw(nb, style = 'W', zero.policy = T)

## on va tester des valeur de K entre 10 et 200 en faisant des bonds de 10
Ks <- seq(10,200,by=10)
AssessmentValues <- lapply(Ks,FUN = function(x){
  ## calculons le model
  print(paste("Calculating model for k = ",x,sep=""))
  Model <- gam(ratio ~ 
                 prt_minorite_vis + prt_monoparental + prt_chomage + revenu_median + 
                 densite_population_km2 + acs_idx_emp_tc_peak + men_size +
                 acs_idx_emp_velo + acs_idx_emp_pieton +
                 s(GEOID, bs = 'mrf', xt = list(nb = nb),k=x),
               data = data_sr,family = gaussian)
  
  data_sr$Pred <- predict(Model,newdata=PredDF)
  data_sr$Pred <- data_sr$Pred-mean(data_sr$Pred)
  
  
  print(tm_shape(data_sr) + 
          tm_fill("Pred", palette = "-Spectral", n = 7, style = 'jenks') + 
          tm_layout(paste("ratio transport mode trend with k = ",x,sep="")))
  
  
  Values <- c(
    "k"=x,
    "aic" = AIC(Model),
    "moranI" = moran(residuals(Model),Queen_W,n=nrow(data_sr),Szero(Queen_W))$I,
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



model_GAM <- gam(ratio ~ 
                   prt_minorite_vis + prt_monoparental + prt_chomage + revenu_median + 
                   densite_population_km2 + acs_idx_emp_tc_peak + men_size +
                   acs_idx_emp_velo + acs_idx_emp_pieton +
               s(GEOID, bs = 'mrf', xt = list(nb = nb),k=90),
             data = data_sr,family = scat(link = 'identity'))



TestNormality(residuals(model_GAM))
qqPlot(residuals(model_GAM), distribution = 't', df  = 3)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Modele GWR ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(spgwr)

# nous utilisons une bandwidth adaptative avec optimisation par CV

bwaCV.voisins  <- gwr.sel(ratio ~ 
                            prt_minorite_vis + prt_monoparental + prt_chomage + revenu_median + 
                            densite_population_km2 + acs_idx_emp_tc_peak + men_size +
                            acs_idx_emp_velo + acs_idx_emp_pieton,
                          data = data_sr,
                          method = "cv",          # Méthode cv ou AIC
                          gweight=gwr.bisquare,   # gwr.gauss ou gwr.bisquare
                          adapt=TRUE,
                          verbose = TRUE,
                          RMSE = TRUE,
                          longlat = FALSE,
                          coords=cbind(data_sr$X,data_sr$Y))

bwaCV.voisins

# réalisation de la GWR

Modele.GWR <- gwr(ratio ~ 
                    prt_minorite_vis + prt_monoparental + prt_chomage + revenu_median + 
                    densite_population_km2 + acs_idx_emp_tc_peak + men_size +
                    acs_idx_emp_velo + acs_idx_emp_pieton,
                  data = data_sr,
                  adapt=bwaCV.voisins,
                  gweight=gwr.bisquare,
                  hatmatrix=TRUE,
                  se.fit=TRUE,
                  coords=cbind(data_sr$X,data_sr$Y),
                  longlat=F)

# on vérifie si le modèle GWR se justifie

anova(Modele.GWR)

LMZ.F1GWR.test(Modele.GWR)
LMZ.F2GWR.test(Modele.GWR)
LMZ.F3GWR.test(Modele.GWR)

data_sr$R2 <- Modele.GWR$SDF$localR2

tm_shape(data_sr) + 
  tm_fill('R2', k = 7, style = 'jenks', title = 'R2'
  ) + 
  tm_compass(position = c("right", "bottom"), 
             size = 0.5)+
  tm_legend(position = c("left", "top"),
            frame = FALSE, bg.color = "white")+
  tm_scale_bar(breaks  = c(0, 10, 20),
               position = c("left", "bottom"))+
  tm_layout(main.title = "R2 local du modèle GWR",
            legend.outside = TRUE, 
            attr.outside = TRUE,
            inner.margins = 0,
            frame = FALSE, 
            legend.format = list(text.separator = "-",
                                 text.or.more = 'ou plus'
            )
  ) 

# on vérifie l'absence d'autocorrélation spatiale dans les résidus
data_sr$resid_gwr <- (data_sr$ratio - Modele.GWR$SDF$pred)
TestNormality(data_sr$resid_gwr)

moran.mc(data_sr$resid_gwr, Queen_W,nsim=999) # ok pour l'autocorrélation spatiale
