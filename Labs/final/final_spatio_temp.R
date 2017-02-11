# Final Project
setwd("/Users/alexandraplassaras/Desktop/Columbia_Courses/Spring_2016/QMSS_G4071/QMSS_G4071_Spatial_Analysis/Final_Project")
library(sp)
library(xts)
library(rgdal)
library(maps)
library(spacetime)
library(RColorBrewer)
library(maptools)

install.packages("reshape2")
library(reshape2)


#***********import x,y (long/lat) of observed locations***********

death<- read.csv("A4_coords_filt_only.csv", sep=",")
death<-as.matrix(death)
class(death)

death<-SpatialPoints(death[,c(2,1)])
class(death)


proj4string(death)<-CRS("+proj=longlat +datum=WGS84")
plot(death)




##### Questions about arbitarily setting years for this? How can I get csv file w/ lat lon, date and death attribute for final project
#*********Specify Time and set as Date Frame*******************
death.yrs<-1989:2014
head(death.yrs)
death.y<-as.Date(paste(death.yrs,"-01-01",sep=""),"%Y-%m-%d")
head(death.y)

#*************Read in Data set, Set Missing Values, Create STFDF Object******
ddata<- read.csv("A4coords_matrix.csv", sep=",")
class(ddata)
head(ddata)


####### Having a problem here! ######
ddata.st<-STFDF(death,death.y,data.frame(counts=as.vector(as.matrix(ddata))))
#head(bdata.st)
plot(ddata.st)




#*************************check dimenstionality of ST Object*********************
dim(bdata.st)

#space - location
#time - num of yrs
#variables - count of birds

#***********prepare data for aggregation example to state level**************

m <- map('state',region=c('florida','georgia'), fill=TRUE,plot=FALSE)
FLGA<-map2SpatialPolygons(m,c("FL","GA"))
proj4string(FLGA)<-proj4string(bdata.st)
plot(FLGA)

mw <- map('world',region=c('myanmar','china', 'india', 'thailand'), fill=TRUE,plot=FALSE)
# won't let me plot four countries together....
MCIT<-map2SpatialPolygons(m,c("MM","CH", "IN", "TH"))
# must change bdata to ddata
proj4string(MCIT)<-proj4string(bdata.st)
plot(MCIT)

m<-map("state","florida",fill=TRUE,plot=FALSE)
#m
FL<-map2SpatialPolygons(m,"FL")
#FL
proj4string(FL)<-proj4string(bdata.st)
plot(FL)


#**************check dimensionality by aggregation categories****************
dim(bdata.st[FL,])
dim(bdata.st[,"1998::2003"])
# across entire space and tim period for counts
dim(bdata.st[,,"counts"])
dim(bdata.st[FL,"1998::2003","counts"])



#**********add square root of counts (comment says square though....)******************************
bdata.st$sqrtcounts<-sqrt(bdata.st$counts)



#*****************aggregate bird counts to two year periods and identify matching space time geometries*********
bb<-STF(FL,bird.y[c(4,6,8,10,12)])
# years 83, 85, 97 and something else.... aggregate for 2 year increments
# we create 2 yr aggregates.... only florida for bird years 4,6,8,10,12.....
over(bb,bdata.st,fn=sum,na.rm=TRUE)



#***********Assign Attribute table to creat new STFDF object****************
b.counts<-new("STFDF",bb,data=over(bb,bdata.st,fn=sum,na.rm=TRUE))




#************alternative aggregation syntax****************************
aggregate(bdata.st,bb,sum,na.rm=TRUE)




#********************Aggregate data to 5-year periods********************
bird.5y<-aggregate(bdata.st,"5 years",mean,na.rm=TRUE)
# head(bird.5y)

bird.3y<-aggregate(bdata.st,"3 years",mean,na.rm=TRUE)
head(bird.3y)

#*******************Select only cases in the state of Florida**********
bird.fl<-bird.5y[FL, ,"sqrtcounts"]
bird.fl3<-bird.3y[FL, ,"sqrtcounts"]

# head(bird.fl)
x<-as(bird.fl,"xts")
x[is.na(x)]<-0

x3<-as(bird.fl3,"xts")
x3[is.na(x3)]<-0


#*************Plot florida 5 year results in lattice and animated form******************************
stplot(bird.fl,cuts=5,animate=0)
stplot(bird.fl,cuts=5,animate=1)


stplot(bird.fl3,cuts=5,animate=0)



#****************************compute count weighted mean time to measure abundance and early presence************************
o<-order(as.vector(1:4 %*% x)/apply(x,2,sum))
head(o)
# weighted mean time....

o3<-order(as.vector(1:6 %*% x3)/apply(x3,2,sum))

#*****************************create the Hovmoller Diagram of weighted mean time count*****************
pal<-brewer.pal(6,"Blues")
cuts<-c(0,2,4,6,8,10,12)
ck<-list(at=cuts, labels=as.character(cuts^2))
stplot(bird.fl[o,], mode="xt", col.regions=pal, cuts=6, asp=0.5, xlab="Sites", colorkey=ck)


pal<-brewer.pal(6,"Blues")
cuts<-c(0,2,4,6,8,10,12)
ck<-list(at=cuts, labels=as.character(cuts^2))
stplot(bird.fl3[o3,], mode="xt", col.regions=pal, cuts=6, asp=0.5, xlab="Sites", colorkey=ck)



#*************************plot counts over time*****************************************
stplot(bird.fl, mode="ts")

stplot(bird.fl3, mode="ts")


#**************************plot counts over time by location*****************************
stplot(bird.fl, mode="tp")
stplot(bird.fl3, mode="tp")
