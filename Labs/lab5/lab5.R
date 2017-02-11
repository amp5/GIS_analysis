setwd("/Users/alexandraplassaras/Desktop/Columbia_Courses/Spring_2016/QMSS_G4071/QMSS_G4071_Spatial_Analysis/Lab5_GWR")

library(maptools)
library (spdep)
library(spgwr)
library (rgdal)


NY_GWR <- readShapeSpatial("lab5")

plot(NY_GWR)

names(NY_GWR)

spplot(NY_GWR, "Z")

spplot(NY_GWR, "PEXPOSURE")

spplot(NY_GWR, "PCTAGE65P")

spplot(NY_GWR, "PCTOWNHOME")

bwG <- gwr.sel(Z~PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=NY_GWR, gweight=gwr.Gauss, verbose=TRUE)

gwrG <- gwr(Z~PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=NY_GWR, bandwidth=bwG, gweight=gwr.Gauss)

gwrG

names(gwrG)

names (gwrG$SDF)

spplot (gwrG$SDF, "localR2")


writeSpatialShape(gwrG$SDF, "GWR_Results")
