setwd("C:/Users/appariciop/OneDrive - USherbrooke/Articles Publies ne pas supprimer/_Livres/BolR_RegressionsSpatiales/data")


library(sf)
library(tmap)

rm(list = ls())

load("Bcn.Rdata")
tm_shape(Bcn)+tm_fill(col="pm_10")

load("Lyon.Rdata")
head(Bcn)


load("Density.RData")
head(Barcelona)
tm_shape(Barcelona)+tm_fill(col="Density")


head(LyonIris)


load("ProbitData.RData")
head(BMA)
head(MMA)

tmap_mode("")
tm_shape(BMA)+tm_fill(col="density")
tm_shape(MMA)+tm_fill(col="density")


load("Bcn.Rdata")
tm_shape(BMA)+tm_fill(col="dist")
