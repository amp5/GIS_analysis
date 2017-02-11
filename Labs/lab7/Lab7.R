setwd("/Users/alexandraplassaras/Desktop/Columbia_Courses/Spring_2016/QMSS_G4071/QMSS_G4071_Spatial_Analysis/Lab7_Spatial_Reg")

library(spdep)
library(maptools)
library(rgdal)


TRI_poly<-readOGR(dsn="/Users/alexandraplassaras/Desktop/Columbia_Courses/Spring_2016/QMSS_G4071/QMSS_G4071_Spatial_Analysis/Lab7_Spatial_Reg",layer="TRI_REG")
names(TRI_poly)

TRI_nbq<-poly2nb(TRI_poly)

coords<-coordinates(TRI_poly)
plot(TRI_poly)
# plot connectivity matrix of queen on coords
plot(TRI_nbq, coords, add=T)

summary(TRI_nbq)

#Neighbour list object:
#  Number of regions: 82  -> 82 counties
#Number of nonzero links: 432 
#Percentage nonzero weights: 6.424747  -> connection = 1 about 6% account for your neighnors when compared to all the counties
#Average number of links: 5.268293  - bout 5.2 neighbors
#Link number distribution:
  
# 11 counties that have 3 neighbors... on avg its about 5 
#  3  4  5  6  7  8  9 
#11 14 20 20 14  2  1 
#11 least connected regions:
#  0 2 3 5 20 64 67 72 74 79 81 with 3 links
#1 most connected region:
#  38 with 9 links

# ti use weight matrix, create a list. which is what we do here from queen neighborhood
TRI_nbq_w<-nb2listw(TRI_nbq)

# look for spatial dependence in exposure var
# sig with small p value
moran.test(TRI_poly$Exposure, listw=TRI_nbq_w)

# now we have spatial dependece need to figure out which regression to use

fips <- order(TRI_poly$COUNTY)
MSlocI <- localmoran(TRI_poly$Exposure, TRI_nbq_w)
printCoefmat(data.frame(MSlocI[fips,],row.names=TRI_poly$COUNTY[fips]), check.names=FALSE)
LocI <- moran.plot(TRI_poly$Exposure,TRI_nbq_w,labels=as.character(TRI_poly$NAME),xlim=c(-1,6.5),ylim=c(-1,4.5),xlab="Exposure",ylab="Spatial Lag Exposure")

infl <- apply(LocI$is.inf, 1, any)
x <- TRI_poly$Exposure
lhx <- cut(x, breaks=c(min(x), mean(x), max(x)), labels=c("L","H"), include.lowest=TRUE)
wx <- lag(TRI_nbq_w, TRI_poly$Exposure)
lhwx <- cut(wx, breaks=c(min(wx), mean(wx), max(wx)),labels=c("L", "H"), include.lowest=TRUE)
lhlh <- interaction(lhx, lhwx, infl, drop=TRUE)
cols <- rep(1, length(lhlh))
cols[lhlh == "H.L.TRUE"] <- 2
cols[lhlh == "L.H.TRUE"] <- 3
cols[lhlh == "H.H.TRUE"] <- 4
plot(TRI_poly, col=grey.colors(4, 0.95, 0.55, 2.2)[cols])
legend("topright", legend=c("None", "HL", "LH", "HH"),fill=grey.colors(4, 0.95, 0.55, 2.2), bty="n", cex=0.8,y.intersp=0.8)

# high high clusters in mid region. 


# ******SPATIAL REGRESSION LAB***********************
# run linerar model, expo is regressed on other ones. 
lm(Exposure ~ HS_DO + UNEMP + POV + PBLK, data=TRI_poly)
# then save results. 
TRI.lm <- lm(Exposure ~ HS_DO + UNEMP + POV + PBLK, data=TRI_poly)
summary(TRI.lm)
# can see the vars, only unemployment is sig. as unemply incr, exposure increases. HS drop is approaching sig but... the vars are all related a bunch so be wary

# run lagrange test - use reuslts from linaer model and include tri neighborhood weights list 
# SARMA - indicator that both are sig
TRI.lagrange <- lm.LMtests(TRI.lm,TRI_nbq_w, test=c("LMerr","RLMerr","LMlag","RLMlag","SARMA"))
print(TRI.lagrange)
# SARMA, LMlag, LMerr are sig. 
# if lag was not sig we'd use error. 
# for this case, we stick with LM lag. 
# robust tests for lag in absence of error (RLM) and for the other. 
# SARMA tests for both - like interaction test


TRI.lag <- lagsarlm(Exposure ~ HS_DO + UNEMP + POV + PBLK, data=TRI_poly, TRI_nbq_w)
summary(TRI.lag)
# Rho - so spatial sig. 
# as avg amount of neighbors exposure incr by one, local unit increses by .4 unit.

TRI.err <- errorsarlm(Exposure ~ HS_DO + UNEMP + POV + PBLK, data=TRI_poly, TRI_nbq_w)
summary(TRI.err)
# Lambda - spatial sig
# unemploy is sig. 
# avg residuals of neighbors is related to residuls of local units by .51

# used queen here (TRI_nbq_w)
TRI.durbin <- lagsarlm(Exposure ~ HS_DO + UNEMP + POV + PBLK, data=TRI_poly, TRI_nbq_w, type="mixed") 
# 'mixed' will ask for spatial lag and something to be included. it'll produce the lag... manually control for correlated error structure. while stille stimating rho parameter
summary(TRI.durbin)
# Rho - so spatial sig.
# controlling for error and lag
# unemploy is pos rel to toxic releases and rho effect is sig. 


# not best model. we aren't reducing rho model in between two models. 

# want to visualize the data and then looks at sptail dependence and then walks through diff models
#  might need to adjust models as needed. 


# in geoda can only do linear models. in R you can use a lot of other things. 