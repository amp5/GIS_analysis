setwd("/Users/alexandraplassaras/Desktop/Columbia Courses/Spring 2016/QMSS G4071/Lab1_R_QGIS_Materials")

install.packages('maps')
install.packages('maptools')
install.packages('rgdal')
install.packages('RColorBrewer')
install.packages('classInt')

***************************************PROJECTIONS

library(maps) 

oldpar<-par()

world <- map("world", res=0)
str(world)
head(world$names)
plot(world)

states <- map("state", res=0)
str(states)
head(states$names)
plot(states)

library(maptools)

spworld <- map2SpatialLines(world, proj4string = CRS("+proj=longlat"))
spstates <- map2SpatialLines(states, proj4string = CRS("+proj=longlat"))

str(spworld, max.level=2)
str(spstates,max.level=2)

plot(spworld)
plot(spstates)

#Used to keep projections of our map
library(rgdal)

world.laea <- spTransform(spworld, CRS("+proj=laea +lat_0=0 +lon_0=0"))
states.laea <- spTransform(spstates, CRS("+proj=laea +lat_0=43.0758 +lon_0=-89.3976"))


******Run Following Code chunk together


par(mfrow = c(2, 2), pty = "s", cex.axis = 0.5)
plot(spworld, axes = T)
title(main = "Longitude and\nLatitude")
plot(world.laea, axes = T)
title(main = "Lambert Azimuthal\nEqual Area")
plot(spstates, axes = T)
title(main = "Longitude and\nLatitude")
plot(states.laea, axes = T)
title(main = "Lambert Azimuthal\nEqual Area")


***********************************************SPATIAL REFERENCING
#resetting parameters to original paramters we created in the beginning
par(oldpar)


map.states <- map("state", plot = T, fill = T, res=0)

list.names.states <- strsplit(map.states$names,":")
tail(list.names.states)

map.IDs <- sapply(list.names.states, function(x) x[1])
tail(map.IDs)

polystates <- map2SpatialPolygons(map.states, IDs = map.IDs,proj4string = CRS("+proj=longlat"))

summary(polystates)

plot(polystates)

#now we project the data
states.laea <- spTransform(polystates, CRS("+proj=laea +lat_0=43.0758 +lon_0=-89.3976"))
plot(states.laea)

#identified spatial IDs
sp.IDs <- sapply(slot(states.laea, "polygons"), function(x) slot(x,"ID"))
tail(sp.IDs)

#object to data frame
sat <- read.csv("sat.csv", stringsAsFactors = F,row.names = 1)
head(sat)

#creates Data Frame
states.sat <- SpatialPolygonsDataFrame(polystates,sat)
summary(states.sat)

states.sat.laea <- spTransform(states.sat, CRS("+proj=laea +lat_0=43.0758 +lon_0=-89.3976"))

plot(states.sat.laea)
summary(states.sat.laea)



*******************************MAPPING
par(mfrow = c(1,1), pty = "s", cex.axis = 0.5)
library(RColorBrewer)
display.brewer.all()

library(classInt)

# verbal scores
plotvar <- states.sat.laea$verbal
nclr <- 5
plotclr <- brewer.pal(nclr, "Greys")
plotclr
# using quantile distribution
class <- classIntervals(plotvar, nclr, style = "quantile")
class
colcode <- findColours(class, plotclr, digits = 3)
colcode
plot(states.sat.laea, col = colcode)


# math scores
plotvar2 <- states.sat.laea$math
plotclr <- brewer.pal(nclr, "Purples")
class <- classIntervals(plotvar2, nclr, style = "quantile")
colcode <- findColours(class, plotclr, digits = 3)
plot(states.sat.laea, col = colcode, border = "grey",axes = T)
title(main = "SAT math scores in 1999")
legend("bottomleft", legend = names(attr(colcode,"table")), fill = attr(colcode, "palette"))

#exports data into shapefile
writeSpatialShape(states.sat.laea, "sat")




