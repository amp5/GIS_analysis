setwd("/Users/alexandraplassaras/Desktop/Columbia_Courses/Spring_2016/QMSS_G4071/QMSS_G4071_Spatial_Analysis/Lab4_PointPatternAnalysis")

########Load Lab PAckages###########

library(maptools)
library(rgdal)
library(shapefiles)
library(spatstat)
library(splancs)
library(sp)

########Read .shp file into R for splancs#############


TRI<- readShapePoints("MS_TRI")
TRI_Coords<-coordinates(TRI)
border <- readShapePoly(paste("ms_dissolve", sep=""))
MSbord<-border@polygons[[1]]@Polygons[[1]]@coords

plot(TRI)
plot(border)



########Read .shp file into R for spatstat#############

# sets boundaries for each function, set TRI points as axes not the polygons

border2<-readShapePoly("ms_dissolve")
# can't get to here. No spatstat installed
boundry<-as(border2,"owin")
TRI2<-readShapePoints("MS_TRI")
TRIpts<-as(TRI2,"ppp")
#ppp is an indicator of points. notates TRI as points
TRI_border<-ppp(TRIpts$x,TRIpts$y,window=boundry)

plot(TRI_border,axes=T)


#########Quadrat Method of Exploring Randomness#########
# run quadrat test. 10 by 10 cols
qt<-quadrat.test(TRI_border, nx=10, ny=10)
qt

plot(TRI_border)
plot(qt, add=T, cex=.5)




##########G Estimate in R###############################

G<-Gest(as(TRI_border,"ppp"))

plot(G)

G1 <- envelope(TRI_border, Gest, nsim = 100, rank = 2)
plot(G1)

# if clustered should see rapid increase
#black - obs
#red - Theo
#green (random probability) 
#envelope should be around blue line
#blue - confidence intervals
#looks like no clustering


#############F Estimate in R#############################
# should be a slow rise
r=seq(0,10,by=0.1)
F<-envelope(TRI_border, Fest, r=r, nsim=10, rank=2)

plot(F)



################K Function in R ########################


L<-envelope(TRI_border, Lest, nsim=10, rank=2, global=T)

plot(L)



