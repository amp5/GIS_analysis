######## Read .shp file into R #############
setwd("F:/Spring 2015/Week4_Interpolation/Lab4")

library(sp)
library(rgdal)
BMcD <- readOGR(".", "BMcD")
BMcD$Fldf <- factor(BMcD$Fldf)
names(BMcD)
spTransform(BMcD, CRS("+proj=laea +lat_0=51.998807 +lon_0=5.692291"))

########################################################
####### Explore Data  ##################
########################################################

####### Plot the coordinates ##############
plot(BMcD)
bubble(BMcD, "Zn")
# want to know how much zinc in each location

######## Explore statistical data #######################
# set flood factor earlier. look at zince by diff flood categories
boxplot(Zn ~ Fldf, BMcD, width=table(BMcD$Fldf), col="grey")


###### Create SptialPixelsDataFrame #################
BMcD_grid <- as(readGDAL("BMcD_fldf.txt"), "SpatialPixelsDataFrame")
names(BMcD_grid) <- "Fldf"
# creating underlying grid estimates
BMcD_grid$Fldf <- as.factor(BMcD_grid$Fldf)



###### Plot Points and Grid Overlay ################
pts = list("sp.points", BMcD, pch = 4, col = "white")
spplot(BMcD_grid, "Fldf", col.regions=1:3, sp.layout=list(pts))

# use this grid feild to estimate value sin between x's in underlying flood feilds to estimate zinc deposits.


########################################################
####### Aspatial Model - linear model prediction  ##################
########################################################



########## Set parameters for interpolated plots #######
bluepal <- colorRampPalette(c("azure1", "steelblue4"))
#returning a function that interpolates a given set of colors to create new color palettes
#list(bluepal)
#breaks - categorical breaks for legend colors....
brks <- c(0,130,155,195,250,330,450,630,890,1270,1850)
cols <- bluepal(length(brks)-1)
#list(cols)
scols <- c("green", "red")
#using for image for grid



#######Examine residuals and significance differences in Zinc (ANOVA)######
library(ipred)
#prediction error based on resampling, est.para: list of additional parameters that control calculation of estimation
res <- errorest(Zn ~ 1, data = as(BMcD, "data.frame"), model=lm, est.para=control.errorest(k=nrow(BMcD), random=FALSE, predictions=TRUE))
#control.errorest (k -fold cross-validation; nboot - bootstrap; random - cross validation random ordering; prediction returneds?)
round(res$error, 2)
#linear model prediction
fres <- lm(Zn ~ Fldf, data=BMcD)
anova(fres) #ANOVA : between Zinc and Fldf
#prediction error for Zn ~ Fldf
eres <- errorest(Zn ~ Fldf, data = as(BMcD, "data.frame"), model=lm, est.para=control.errorest(k=nrow(BMcD), random=FALSE, predictions=TRUE))
round(eres$error, 2)


####estimate and aspatail model###########
# rel between zinc and flood planes
# take what we know about flood felid and predict what zinc would be there
library(maptools)
#predict floodplain grid using lm(Zn~Fldf)
BMcD_grid$lm_pred <- predict(fres, newdata=BMcD_grid)
#displaying color grid, breaks - finite numeric breakpoint of colors, must be increasing; list of colors)
image(BMcD_grid, "lm_pred", breaks=brks, col=cols)
title("Flood frequency model interpolation")
pe <- BMcD$Zn-eres$predictions #Prediction error: Zn (actual) - predictions (errorest(Zn ~ Fldf))
#bg=scols[(pe < 0)+1] indicator if (pe<0)==1, then scols[2]: red; else scols[1]:green
symbols(coordinates(BMcD), circles=sqrt(abs(pe)), fg="black", bg=scols[(pe < 0)+1], inches=FALSE, add=TRUE)
legend("topleft", fill=cols, legend=leglabs(brks), bty="n", cex=0.8)


########################################################
####### Thin Plate Spline Interpolation - using "Tps"  ##################
########################################################

####Thin Plate Spline Interpolation#######
library(fields)
#estimate prediction error for thin plate spline
pe_tps <- numeric(nrow(BMcD))
cBMcD <- coordinates(BMcD)
for (i in seq(along=pe_tps)) {
  #thin plate; CBMcD - coordinate (everything except row i); BMc$Zn - dependent var (everything except row i) 
  tpsi <- Tps(cBMcD[-i,], BMcD$Zn[-i])
  pri <- predict(tpsi, cBMcD[i,,drop=FALSE])
  pe_tps[i] <- BMcD$Zn[i]-pri
}
round(sqrt(mean(pe_tps^2)), 2)
#thin plate spline actual interpolation
tps <- Tps(coordinates(BMcD), BMcD$Zn)


#####Plot Spline Results########
BMcD_grid$spl_pred <- predict(tps, coordinates(BMcD_grid))
image(BMcD_grid, "spl_pred", breaks=brks, col=cols)
title("Thin plate spline model")
symbols(coordinates(BMcD), circles=sqrt(abs(pe_tps)), fg="black", bg=scols[(pe_tps < 0)+1], inches=FALSE, add=TRUE)
legend("topleft", fill=cols, legend=leglabs(brks), bty="n", cex=0.8)



########################################################
#######  Kriging  - using "gstat(model = fit.variogram)"##################
########################################################

####### Step 1. Ordinary kriging: Using Zn~1 ################

#####Estimate Spatial Dependence in Data##########
#using variogram to estimate spatial dependence
# now starting k-riging
library(gstat)
# new object
cvgm <- variogram(Zn~1, data=BMcD, width=100, cutoff=1000)
# fit line to model
efitted <- fit.variogram(cvgm, vgm(psill=1, model="Exp", range=100, nugget=1))
efitted
#nuggest(starting point), range(distance of spatial dependence), sill(leveling point of semivariance)

###Plot Semivariogram#######
plot(cvgm, model=efitted, plot.numbers=TRUE, col="black")


####Interpolate Using Original Kriging############
# just zinc based on vriogram we created - effitted to predict zinc
#geostat objects, id = identifier of new var, model:fit.variogram
OK_fit <- gstat(id="OK_fit", formula = Zn ~ 1, data = BMcD, model=efitted)
#cross-validation for kriging residuals : prediction error
pe <- gstat.cv(OK_fit, debug.level=0, random=FALSE)$residual
#summary(gstat.cv(OK_fit, debug.level=0, random=FALSE))
round(sqrt(mean(pe^2)), 2)
#kriging prediction
z <- predict(OK_fit, newdata=BMcD_grid, debug.level=0)
BMcD_grid$OK_pred <- z$OK_fit.pred
BMcD_grid$OK_se <- sqrt(z$OK_fit.var)


####Plot OK Model##################
image(BMcD_grid, "OK_pred", breaks=brks, col=cols)
title("Fitted exponential OK model")
symbols(coordinates(BMcD), circles=sqrt(abs(pe)), fg="black", bg=scols[(pe < 0)+1], inches=FALSE, add=TRUE)
legend("topleft", fill=cols, legend=leglabs(brks), bty="n", cex=0.8)

#######Step 2. Universal Kriging model: Zn ~ Fldf ################

#########Add Flood Frequency to the semivariance estimation############
# want to improve on this. use cokriging
cvgm <- variogram(Zn~Fldf, data=BMcD, width=100, cutoff=1000)
uefitted <- fit.variogram(cvgm, vgm(psill=1, model="Exp", range=100, nugget=1))
uefitted


#####Plot Universial Kriging Model#########
plot(cvgm, model=uefitted, plot.numbers=TRUE, col="black")


####Interpolate Using the Flood Frequency UK Model########
UK_fit <- gstat(id="UK_fit", formula = Zn ~ Fldf, data = BMcD, model=uefitted)
pe_UK <- gstat.cv(UK_fit, debug.level=0, random=FALSE)$residual
round(sqrt(mean(pe_UK^2)), 2)
z <- predict(UK_fit, newdata=BMcD_grid, debug.level=0)
BMcD_grid$UK_pred <- z$UK_fit.pred
BMcD_grid$UK_se <- sqrt(z$UK_fit.var)


######Plot UK Model#########
image(BMcD_grid, "UK_pred", breaks=brks, col=cols)
title("Flood frequency UK model")
symbols(coordinates(BMcD), circles=sqrt(abs(pe_UK)), fg="black", bg=scols[(pe_UK < 0)+1], inches=FALSE, add=TRUE)
legend("topleft", fill=cols, legend=leglabs(brks), bty="n", cex=0.8)


########################################################
#######  Comparing interpolations ##################
########################################################


#####Plot all Interpolation for Comparison######
pts = list("sp.points", BMcD, pch = 4, col = "black", cex=0.5)
spplot(BMcD_grid, c("lm_pred", "spl_pred", "OK_pred", "UK_pred"), at=brks, col.regions=cols, sp.layout=list(pts))


#####Write raster to TIF#####
writeGDAL(BMcD_grid["UK_pred"], "UK_pred.tif")




