setwd("/Users/alexandraplassaras/Desktop/Columbia_Courses/Spring_2016/QMSS_G4071/QMSS_G4071_Spatial_Analysis/Lab6_SWR")

library(maptools)
library(rgdal)
library(spdep)

# 100 features are counties
#18 feilds are variable
# sudden infant death - sids
sids<-readOGR(dsn="/Users/alexandraplassaras/Desktop/Columbia_Courses/Spring_2016/QMSS_G4071/QMSS_G4071_Spatial_Analysis/Lab6_SWR",layer="sids2")

# tells us this is a spatial polygon dataframe
class(sids)

# shows our diff vars
names(sids)

# polygons to neighborhood. creates by default queens matrix
# any border = neighbors
sids_nbq<-poly2nb(sids)

# will produce a rook matrix
# must have a certain level of border to = neighbors

sids_nbr<-poly2nb(sids, queen=FALSE)

coords<-coordinates(sids)

plot(sids)
plot(sids_nbq, coords, add=T)
plot(sids)
plot(sids_nbr, coords, add=T)

# this is for nearest neighbors - first, second and fourth
IDs<-row.names(as(sids, "data.frame"))
sids_kn1<-knn2nb(knearneigh(coords, k=1), row.names=IDs)
sids_kn2<-knn2nb(knearneigh(coords, k=2), row.names=IDs)
sids_kn4<-knn2nb(knearneigh(coords, k=4), row.names=IDs)

plot(sids)
plot(sids_kn1, coords, add=T)
plot(sids)
plot(sids_kn2, coords, add=T)
plot(sids)
plot(sids_kn4, coords, add=T)

# distance space approach
dist<-unlist(nbdists(sids_kn1, coords))
summary(dist)

# max distance == .4226 - largest min distance
max_k1<-max(dist)

#creating distance matricies - 3/4, ==, 1/2 times max distance
sids_kd1<-dnearneigh(coords, d1=0, d2=0.75*max_k1, row.names=IDs)
sids_kd2<-dnearneigh(coords, d1=0, d2=1*max_k1, row.names=IDs)
sids_kd3<-dnearneigh(coords, d1=0, d2=1.5*max_k1, row.names=IDs)

# as distance increases, network increases
plot(sids)
# many neighborless counties
plot(sids_kd1, coords, add=T)
plot(sids)
# everyone will have at least one neighbor
plot(sids_kd2, coords, add=T)
plot(sids)
# now network is much denser. 
plot(sids_kd3, coords, add=T)

# work with contiguity queens matrix
sids_nbq_w<-nb2listw(sids_nbq)
sids_nbq_w

# binary??
sids_nbq_wb<-nb2listw(sids_nbq, style="B")
sids_nbq_wb


dist<-nbdists(sids_nbq, coordinates(sids))
idw<-lapply(dist, function(x) 1/(x/1000))
sids_nbq_idwb<-nb2listw(sids_nbq, glist=idw, style="B")
summary(unlist(sids_nbq_idwb$weights))

# from here we run moran's test
# interested in global clustering in sidr79
# significant clustering in our data
# weak to moderate ???
moran.test(sids$SIDR79, listw=sids_nbq_w)

# producing a plot of moran's I results
# our moran's I is not likely to have occured by chance -> on the far right
set.seed(1234)
perm<-moran.mc(sids$SIDR79,listw=sids_nbq_w,nsim=999)
perm

mean(perm$res[1:999])

var(perm$res[1:999])

summary(perm$res[1:999])

hist(perm$res, freq=TRUE, breaks=20, xlab="Simulated Moran's I")

# ID code for FIPS 
fips <- order(sids$FIPSNO)
# look at global morans for sid79
nclocI <- localmoran(sids$SIDR79, sids_nbq_w)
printCoefmat(data.frame(nclocI[fips,],row.names=sids$FIPSNO[fips]), check.names=FALSE)
# moran's plot will be produced
nci <- moran.plot(sids$SIDR79,sids_nbq_w,labels=as.character(sids$NAME),xlim=c(-1,6.5),ylim=c(-1,4.5),xlab="SIDS Rate",ylab="SL SIDS Rate")
# as your rate of sids increases, your neighbor's also increase

# will plot the local clusters LISA clusters
infl <- apply(nci$is.inf, 1, any)
x <- sids$SIDR79
lhx <- cut(x, breaks=c(min(x), mean(x), max(x)), labels=c("L","H"), include.lowest=TRUE)
wx <- lag(sids_nbq_w, sids$SIDR79)
lhwx <- cut(wx, breaks=c(min(wx), mean(wx), max(wx)),labels=c("L", "H"), include.lowest=TRUE)
lhlh <- interaction(lhx, lhwx, infl, drop=TRUE)
cols <- rep(1, length(lhlh))
cols[lhlh == "H.L.TRUE"] <- 2
cols[lhlh == "L.H.TRUE"] <- 3
cols[lhlh == "H.H.TRUE"] <- 4
plot(sids, col=grey.colors(4, 0.95, 0.55, 2.2)[cols])
legend("topright", legend=c("None", "HL", "LH", "HH"),fill=grey.colors(4, 0.95, 0.55, 2.2), bty="n", cex=0.8,y.intersp=0.8)

