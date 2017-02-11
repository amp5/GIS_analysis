# FINAL EXAM!

# Question 1
setwd("/Users/alexandraplassaras/Desktop/Columbia_Courses/Spring_2016/QMSS_G4071/QMSS_G4071_Spatial_Analysis/G4071_FinalExam_")

library(maptools)
library(rgdal)
library(spdep)

# 174  (now 159...) features are counties
#25 feilds are variable
# GA education = ga_edu
ga_edu<-readOGR(dsn="/Users/alexandraplassaras/Desktop/Columbia_Courses/Spring_2016/QMSS_G4071/QMSS_G4071_Spatial_Analysis/G4071_FinalExam_",layer="GeorgiaEduc")

# tells us this is a spatial polygon dataframe
class(ga_edu)

# shows our diff vars
names(ga_edu)

#[1] "FID_1"     "AREA"      "PERIMETER" "G_UTM_"    "G_UTM_ID"  "AREANAME"  "AREAKEY"   "X_COORD"  
#[9] "Y_COORD"   "KEY_VAL"   "FID_2"     "AreaKey_1" "Latitude"  "Longitud"  "TotPop90"  "PctRural" 
#[17] "PctBach"   "PctEld"    "PctFB"     "PctPov"    "PctBlack"  "ID"        "X"         "Y"        
#[25] "Distance" 


######## creates initial contiguity based weight matrix ############
# polygons to neighborhood. creates by default queens matrix
# any border = neighbors
ga_edu_nbq<-poly2nb(ga_edu)



# will produce a rook matrix
# must have a certain level of border to = neighbors

ga_edu_nbr<-poly2nb(ga_edu, queen=FALSE)

coords<-coordinates(ga_edu)

plot(ga_edu)
plot(ga_edu_nbq, coords, add=T)
plot(ga_edu)
plot(ga_edu_nbr, coords, add=T)

############ this is for nearest neighbors - first, second and fourth ##########
IDs<-row.names(as(ga_edu, "data.frame"))
ga_edu_kn1<-knn2nb(knearneigh(coords, k=1), row.names=IDs)
ga_edu_kn2<-knn2nb(knearneigh(coords, k=2), row.names=IDs)
ga_edu_kn4<-knn2nb(knearneigh(coords, k=4), row.names=IDs)

plot(ga_edu)
plot(ga_edu_kn1, coords, add=T)
plot(ga_edu)
plot(ga_edu_kn2, coords, add=T)
plot(ga_edu)
plot(ga_edu_kn4, coords, add=T)

############ distance based  matrices ################ distance space approach
dist<-unlist(nbdists(ga_edu_kn1, coords))
summary(dist)

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#431.3 20010.0 24950.0 23470.0 28110.0 40690.0 





# max distance == 40690 - largest min distance
max_k1<-max(dist)

#creating distance matricies - 3/4, ==, 1/2 times max distance
ga_edu_kd1<-dnearneigh(coords, d1=0, d2=0.75*max_k1, row.names=IDs)
ga_edu_kd2<-dnearneigh(coords, d1=0, d2=1*max_k1, row.names=IDs)
ga_edu_kd3<-dnearneigh(coords, d1=0, d2=1.5*max_k1, row.names=IDs)

# as distance increases, network increases
plot(ga_edu)
# many neighborless counties
plot(ga_edu_kd1, coords, add=T)
plot(ga_edu)
# everyone will have at least one neighbor
plot(ga_edu_kd2, coords, add=T)
plot(ga_edu)
# now network is much denser. 
plot(ga_edu_kd3, coords, add=T)


####### Weight Presentation/ Computation #########
# work with contiguity queens matrix
ga_edu_nbq_w<-nb2listw(ga_edu_nbq)
ga_edu_nbq_w

#Characteristics of weights list object:
#  Neighbour list object:
#  Number of regions: 174 
#Number of nonzero links: 934 
#Percentage nonzero weights: 3.084952 
#Average number of links: 5.367816 
#1 region with no links:
#  129

#Weights style: W 
#Weights constants summary:
#  n    nn  S0       S1       S2
#W 173 29929 173 70.35158 720.5274


# binary??
ga_edu_nbq_wb<-nb2listw(ga_edu_nbq, style="B")
ga_edu_nbq_wb

#Characteristics of weights list object:
#  Neighbour list object:
#  Number of regions: 174 
#Number of nonzero links: 934 
#Percentage nonzero weights: 3.084952 
#Average number of links: 5.367816 
#1 region with no links:
#  129

#Weights style: B 
#Weights constants summary:
#  n    nn  S0   S1    S2
#B 173 29929 934 1868 22608



###### inverse distance weighting #########
dist<-nbdists(ga_edu_nbq, coordinates(ga_edu))
idw<-lapply(dist, function(x) 1/(x/1000))
ga_edu_nbq_idwb<-nb2listw(ga_edu_nbq, glist=idw, style="B")
summary(unlist(ga_edu_nbq_idwb$weights))

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.01461 0.02568 0.03065 0.04017 0.03633 2.31800 



############ Moran's I test ############
# from here we run moran's test
# interested in global clustering in sidr79
# significant clustering in our data
# weak to moderate ???
moran.test(ga_edu$PctBach, listw=ga_edu_nbq_w)

#Moran's I test under randomisation

#data:  ga_edu$PctBach  
#weights: ga_edu_nbq_w  

#Moran I statistic standard deviate = 5.7139, p-value = 5.519e-09
#alternative hypothesis: greater
#sample estimates:
#Moran I statistic       Expectation          Variance 
#      0.261890634      -0.005813953       0.002195023 



# producing a plot of moran's I results
# our moran's I is not likely to have occured by chance -> on the far right
set.seed(1234)
perm<-moran.mc(ga_edu$PctBach,listw=ga_edu_nbq_w,nsim=999)
perm


#Monte-Carlo simulation of Moran's I

#data:  ga_edu$PctBach 
#weights: ga_edu_nbq_w  
#number of simulations + 1: 1000 

#statistic = 0.26189, observed rank = 1000, p-value = 0.001
#alternative hypothesis: greater


mean(perm$res[1:999])
#[1] -0.006290542

var(perm$res[1:999])
#[1] 0.002069721

summary(perm$res[1:999])
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.145200 -0.038180 -0.008524 -0.006291  0.022470  0.160600


hist(perm$res, freq=TRUE, breaks=20, xlab="Simulated Moran's I")



########### Local Indicators of Spatial Association (LISA) ########
# ID code for FIPS 
fips <- order(ga_edu$AREAKEY)
# look at global morans for sid79
nclocI <- localmoran(ga_edu$PctBach , ga_edu_nbq_w)
printCoefmat(data.frame(nclocI[fips,],row.names=ga_edu$AREAKEY[fips]), check.names=FALSE)
# moran's plot will be produced








############### Problem!!!!! ##########################
nci <- moran.plot(ga_edu$PctBach,ga_edu_nbq_w,labels=as.character(ga_edu$AREANAME),xlab="PctBach Rate",ylab="SL PctBach Rate")
# as your rate of ga_edu increases, your neighbor's also increase

# will plot the local clusters LISA clusters
infl <- apply(nci$is.inf, 1, any)
x <- ga_edu$PctBach
lhx <- cut(x, breaks=c(min(x), mean(x), max(x)), labels=c("L","H"), include.lowest=TRUE)
wx <- lag(ga_edu_nbq_w, ga_edu$PctBach)
lhwx <- cut(wx, breaks=c(min(wx), mean(wx), max(wx)),labels=c("L", "H"), include.lowest=TRUE)
lhlh <- interaction(lhx, lhwx, infl, drop=TRUE)
cols <- rep(1, length(lhlh))
cols[lhlh == "H.L.TRUE"] <- 2
cols[lhlh == "L.H.TRUE"] <- 3
cols[lhlh == "H.H.TRUE"] <- 4
plot(ga_edu, col=grey.colors(4, 0.95, 0.55, 2.2)[cols])
title("LISA Clusters")
legend("topright", legend=c("None", "HL", "LH", "HH"),fill=grey.colors(4, 0.95, 0.55, 2.2), bty="n", cex=0.8,y.intersp=0.8)




##### Question 4 #######
library(maptools)
library (spdep)
library(spgwr)
library (rgdal)


GA_GWR <- readShapeSpatial("GeorgiaEduc")
plot(GA_GWR)

names(GA_GWR)

spplot(GA_GWR, "PctBach")
spplot(GA_GWR, "PctRural")
spplot(GA_GWR, "PctPov")
spplot(GA_GWR, "PctBlack")

bwG <- gwr.sel(PctBach~PctRural + PctPov + 
                 PctBlack, data=GA_GWR, gweight=gwr.Gauss, verbose=TRUE)
gwrG <- gwr(PctBach~PctRural + PctPov + 
              PctBlack, data=GA_GWR, bandwidth=bwG, gweight=gwr.Gauss)
gwrG

names(gwrG)
names (gwrG$SDF)

spplot (gwrG$SDF, "localR2")


writeSpatialShape(gwrG$SDF, "GWR_Results")




f_lab=read.csv(file.choose(),header=FALSE,check.names=FALSE, row.names=NULL)
View(f_lab)

f=read.csv(file.choose(),header=TRUE,check.names=FALSE, row.names=NULL)
View(f)
f <- f[2:4]
write.csv(f, file = "f.csv")

c=read.csv(file.choose(),header=TRUE,check.names=FALSE, row.names=NULL)
View(c)
write.csv(c, file = "c.csv")



