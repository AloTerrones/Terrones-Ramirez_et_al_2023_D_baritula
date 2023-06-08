# ---- ECOLOGICAL NICHE MODELS 


####  Installing packages
install.packages("dismo") #functions: gbif
install.packages("dplyr") #functions: select
install.packages("CoordinateCleaner") #functions: clean_coordinates
install.packages("devtools") #function install_github
devtools::install_github("luismurao/nichetoolbox") #function clean_dup
install.packages("raster") #functions: getData, extract
install.packages("corrplot") #functions: cor, corrplot (correlation)
install.packages("rgbif") #functions: occ_download
install.packages("remotes") #functions buffer_area 
install.packages("spocc") 
install.packages("scrubr") 
install.packages("rgdal") #functions: readOGR
remotes::install_github("marlonecobos/ellipsenm") #functions buffer_area
install.packages("stars")
devtools::install_github("rsbivand/sp@evolution") 



####  Uploading packages
library(dismo)
library(rgbif)
library(dplyr)
library(CoordinateCleaner)
library(raster)
library(nichetoolbox)
library(devtools)
library(factoextra) 
library(corrplot) 
library(spocc)
library(scrubr)
library(rgdal) 
library(rgeos)
library(ellipsenm) 
library(sf) 
library(terra) 
library(stars) 
library(sp)
library(ggplot2)
library(RColorBrewer) 
library(ecospat)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(dplyr)
library(rgeos)
library(rnaturalearthhires)
library(tiff)



# Read csv files

complejo <- read.csv("./complejo.csv")
baritula <- read.csv("./baritula.csv")
parva <- read.csv("./parva.csv")
montana <- read.csv("./montana.csv")



# Plotting points

world <- ne_countries(scale = "large", returnclass = "sf")
mexico <- ne_states(country="mexico", returnclass= "sf")


p_complejo <- ggplot(data=world)+
  geom_sf(color="gray")+
  geom_sf(data = mexico, fill = NA, color="gray")+
  coord_sf(xlim=c(-108,-84), ylim=c(10,25))+
  annotation_scale(location = "bl", width_hint = 0.2, line_width = 0.1) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.7, "cm"), pad_y = unit(0.7, "cm"),
                         style = north_arrow_fancy_orienteering) +
  geom_point(data=complejo,aes(x=longitude, y=latitude),size=1)+
  xlab("Long")+
  ylab("Lat")+
  theme(panel.grid.major = element_line(color = "lightgray", linetype = "3313", size = 0.2), panel.background = element_rect(fill = "aliceblue"))+
  ggtitle(expression(paste(italic("Diglossa baritula"))))+
  theme(plot.title = element_text(hjust = 0.5))

p_complejo



#### Bioclimatic variables

# Downloaded from https://www.worldclim.org/data/worldclim21.html
# Defining a variable for the 19 bioclimatic layers from worldclim
BIOS <- stack(list.files(path="./wc2.1_30s_bio/", pattern = "*.tif$",full.names = T)) 
plot(BIOS[[1]])

BIOS <- crop(BIOS, c(-108,-84,10,25))
plot(BIOS[[1]])
points(complejo$longitude, complejo$latitude, col="red", cex=.4, pch=16)

# Create directory
#dir.create("./BIOS.c")

# Save cropped bioclimatic layers in ascii format
writeRaster(BIOS, filename=names(BIOS), bylayer=TRUE, format="ascii")
# If I would like to read this table again:
# BIOS <- stack(list.files(path="BIOS.c/", pattern = "*.asc",full.names = T))


#### Extract values from bioclimatic layers 

complejo_coordinates <- subset(complejo, select=c(longitude, latitude))
BIOS_VALUES <- raster::extract(BIOS, complejo_coordinates)
BIOS_VALUES <- cbind.data.frame(complejo, BIOS_VALUES)

#write.csv(BIOS_VALUES_ALL, "./bios_values_all.csv", row.names = FALSE)

# Eliminate records with no values 
BIOS_VALUES <-na.omit(BIOS_VALUES)
#write.csv(BIOS_VALUES_ALL, "./bios_values_all.csv", row.names = FALSE)


#### Correlation and PCA analysis

BIOS_VALUES_ALL <- BIOS_VALUES 
BIOS_VALUES <- subset(BIOS_VALUES_ALL, select=-c(specie, longitude, latitude))
complejo.cor <-cor(BIOS_VALUES) 

corrplot(complejo.cor, type="lower", order="hclust", tl.col = "black", tl.srt = 20, 
         tl.cex = 0.7, col=brewer.pal(n=8, name="RdYlBu")) 

# Create directory
dir.create("./PCA&COR")

#Save correlation analysis in csv format
write.csv(complejo.cor, file="./PCA&COR/complejo.cor.csv", row.names = T)

# PCA
complejo.pca <- prcomp(BIOS_VALUES, scale = TRUE)

# install.packages("factoextra") #si
# Visualization of explained variances (percentages)
fviz_eig(complejo.pca)

fviz_eig(complejo.pca, choice = c("variance", "eigenvalue"),
         addlabels = TRUE)

# Visualization of variables (circular plot)
fviz_pca_var(complejo.pca,
             col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE) +
  theme(text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

# Contributions of variables to PC1
fviz_contrib(complejo.pca, choice = "var", axes = 1, top = 10)

# Contributions of variables to PC2
fviz_contrib(complejo.pca, choice = "var", axes = 2, top = 10)

# Contributions of variables to PC3
fviz_contrib(complejo.pca, choice = "var", axes = 3, top = 10)

# Contributions of variables to PC4
fviz_contrib(complejo.pca, choice = "var", axes = 4, top = 10)


# Summary statistics PCA
summary(complejo.pca)
complejo.pca$sdev^2    
complejo.pca$rotation   


# Save PCA results
write.csv(complejo.pca$sdev^2, file="./PCA&COR/complejo.pca.sdev2.csv", row.names = T)
write.csv(complejo.pca$rotation, file="./PCA&COR/complejo.pca.rotation.csv", row.names = T)



#### Creating a shape that represents the area "M" (accesibility area)

# Hypothesis of areas for model calibration following:
# https://github.com/marlonecobos/ENM_manuals/blob/master/M_hypotheses/Construction_of_simple_Ms.R
# https://github.com/marlonecobos/ENM_manuals/blob/master/M_hypotheses/M_from_polygon_intersection.R

# getting ecorregions
download.file("http://assets.worldwildlife.org/publications/15/files/original/official_teow.zip",
              destfile = file.path(getwd(), "wwf_ecoregions.zip")) 
unzip(file.path(getwd(), "wwf_ecoregions.zip"), 
      exdir = file.path(getwd(), "WWF_ecoregions"))
file.remove(file.path(getwd(), "wwf_ecoregions.zip"))
ecor <- readOGR("WWF_ecoregions/official", layer = "wwf_terr_ecos")


# Areas by buffering records (100 km buffer)
M_buffer <- buffer_area(complejo, longitude = "longitude", latitude = "latitude", 
                        buffer_distance = 100)

# Areas using convex hulls (including 75 km buffer)
M_convex <- convex_area(complejo, longitude = "longitude", latitude = "latitude", 
                        buffer_distance = 75)

# Areas using concave hulls (including 75 km buffer)
M_concave <- concave_area(complejo, longitude = "longitude", latitude = "latitude", 
                          buffer_distance = 75)

# Areas by selecting polygons (including 25 km buffer)
M_ecorreg <- polygon_selection(complejo, longitude = "longitude", latitude = "latitude",
                               polygons = ecor, buffer_distance = 25)

# Intersection
M_intersect <- gIntersection(M_buffer, M_convex)
M_intersect <- gIntersection(M_intersect, M_ecorreg)

# Visualization
par(mfrow = c(2, 2), cex = 0.8, mai=c(0.4,0.4,0.4,0.4))

plot(BIOS[[1]]); plot(M_buffer, add=TRUE); 
points(complejo[, 2:3], col="red", cex=.4, pch=20); 
legend("topleft", legend = "Buffer", bty = "n")

plot(BIOS[[1]]); plot(M_convex, add=TRUE); 
points(complejo[, 2:3], col="red", cex=.4, pch=20); 
legend("topleft", legend = "Convex", bty = "n")

#plot(BIOS[[1]]); plot(M_concave, add=TRUE); 
#	points(complejo[, 2:3], col="red", cex=.4, pch=20); 
#	legend("topleft", legend = "Concave", bty = "n")

plot(BIOS[[1]]); plot(M_ecorreg, add=TRUE); 
points(complejo[, 2:3], col="red", cex=.4, pch=20); 
legend("topleft", legend = "Ecorregions", bty = "n")

plot(BIOS[[1]]); plot(M_intersect, add=TRUE); 
points(complejo[, 2:3], col="red", cex=.4, pch=20); 
legend("topleft", legend = "Intersection", bty = "n")

par(mfrow=c(1,1))

# Saving shapefile
dir.create("./complejo.M")
writeOGR(M_buffer, "complejo.M/", "M_buffer", driver = "ESRI Shapefile") 



#### Cropping "BIOS" using M_intersect 

BIOS.m <- mask(BIOS, M_buffer)
plot(BIOS.m[[1]])

# Plotting M, bio1 and points

plot(BIOS[[1]])
plot(M_buffer, add=TRUE) #se estan graficando dos mapas uno encima del otro
points(complejo$longitude, complejo$latitude, col="red", cex=.4, pch=16)


# Create directory
#dir.create("./BIOS.m")

# Save cropped bioclimatic layers in ascii format
writeRaster(BIOS.m, filename=names(BIOS.m), bylayer=TRUE, format="ascii")

# If I would like to read this table again:
# BIOS.m <- stack(list.files(path="BIOS.m/", pattern = "*.asc",full.names = T))



#### PAST BIOCLIMATIC LAYERS
# https://worldclim.org/data/v1.4/paleo1.4.html
# LIG: Last Interglacial
# MH: Mid Holocene
# LGM: Last Glacial Maximum
# miroc
# ccsm

# Loading data
LIG <- stack(list.files(path="./lig_30s_bio/", pattern = "*.bil$",full.names = T)) 
MH_miroc <- stack(list.files(path="./mrmidbi_30s/", pattern = "*.tif$",full.names = T)) 
MH_ccsm <- stack(list.files(path="./ccmidbi_30s/", pattern = "*.tif$",full.names = T)) 
LGM_miroc <- stack(list.files(path="./mrlgmbi_2-5m/", pattern = "*.tif$",full.names = T)) 
LGM_ccsm <- stack(list.files(path="./cclgmbi_2-5m/", pattern = "*.tif$",full.names = T)) 

plot(LIG[[1]])

# Resampling
LIG <- raster::resample(LIG, BIOS)
MH_miroc <- raster::resample(MH_miroc, BIOS)
MH_ccsm <- raster::resample(MH_ccsm, BIOS)
LGM_miroc <- raster::resample(LGM_miroc, BIOS)
LGM_ccsm <- raster::resample(LGM_ccsm, BIOS)

# Croping
LIG <- crop(LIG, c(-108,-84,10,25))
plot(LIG[[1]])
points(complejo$longitude, complejo$latitude, col="red", cex=.4, pch=16)

MH_miroc <- crop(MH_miroc, c(-108,-84,10,25))
MH_ccsm <- crop(MH_ccsm, c(-108,-84,10,25))
LGM_miroc <- crop(LGM_miroc, c(-108,-84,10,25))
LGM_ccsm <- crop(LGM_ccsm, c(-108,-84,10,25))


# Create directories
dir.create("./LIG") 
dir.create ("./MH_ccsm") 
dir.create ("./MH_miroc") 
dir.create ("./LGM_miroc") 
dir.create ("./LGM_ccsm")

# Save cropped past bioclimatic layers in ascii format
writeRaster(LIG, filename=names(LIG), bylayer=TRUE, format="ascii")
writeRaster(MH_miroc, filename=names(MH_miroc), bylayer=TRUE, format="ascii")
writeRaster(MH_ccsm, filename=names(MH_ccsm), bylayer=TRUE, format="ascii")
writeRaster(LGM_miroc, filename=names(LGM_miroc), bylayer=TRUE, format="ascii")
writeRaster(LGM_ccsm, filename=names(LGM_ccsm), bylayer=TRUE, format="ascii")


# Checking values in past layers for occurrences

BIOS_VALUES_LIG <- raster::extract(LIG, complejo_coordinates)
BIOS_VALUES_LIG <- cbind.data.frame(complejo, BIOS_VALUES_LIG)
write.csv(BIOS_VALUES_LIG, "./bios_values_LIG.csv", row.names = FALSE)

BIOS_VALUES_LGM_ccsm <- raster::extract(LGM_ccsm, complejo_coordinates)
BIOS_VALUES_LGM_ccsm <- cbind.data.frame(complejo, BIOS_VALUES_LGM_ccsm)
write.csv(BIOS_VALUES_LGM_ccsm, "./bios_values_LGM_ccsm.csv", row.names = FALSE)

BIOS_VALUES_LGM_miroc <- raster::extract(LGM_miroc, complejo_coordinates)
BIOS_VALUES_LGM_miroc <- cbind.data.frame(complejo, BIOS_VALUES_LGM_miroc)

BIOS_VALUES_MH_ccsm <- raster::extract(MH_ccsm, complejo_coordinates)
BIOS_VALUES_MH_ccsm <- cbind.data.frame(complejo, BIOS_VALUES_MH_ccsm)

BIOS_VALUES_MH_miroc <- raster::extract(MH_miroc, complejo_coordinates)
BIOS_VALUES_MH_miroc <- cbind.data.frame(complejo, BIOS_VALUES_MH_miroc)


# Values from bio4, bio5, bio6, and bio7 are multiplied by 10 in past layers

bio4<-raster("./bios_cut/bio4.asc")
bio4 <- bio4*10
writeRaster(bio4, filename="./bios_cut_100/bio4", format="ascii", overwrite=TRUE)

bio5<-raster("./bios_cut/bio5.asc")
bio5 <- bio5*10
writeRaster(bio5, filename="./bios_cut_100/bio5", format="ascii", overwrite=TRUE)

bio6<-raster("./bios_cut/bio6.asc")
bio6 <- bio6*10
writeRaster(bio6, filename="./bios_cut_100/bio6", format="ascii", overwrite=TRUE)

bio7<-raster("./bios_cut/bio7.asc")
bio7 <- bio7*10
writeRaster(bio7, filename="./bios_cut_100/bio7", format="ascii", overwrite=TRUE)



####  After running analysis in MAXENT 
### Results from Maxent were saved in a directory named "maxent_results"
####  Uploading median result from maxent

complejo_maxent<-raster("./maxent_results/Diglossa_baritula_bios_cut_100_median.asc") 
plot(complejo_maxent)

complejo_maxent_presente<-raster("./maxent_results/Diglossa_baritula_median.asc") 
plot(complejo_maxent_presente)


# Defining spatial extent (extremes) 
xmin <- extent(complejo_maxent)
xmax <- extent(complejo_maxent)
ymin <- extent(complejo_maxent)
ymax <- extent(complejo_maxent)

# Create a matrix to set the thresholds for the areas
# which are considered "suitable" or "insuitable". 
# The first column in the matrix created below represents 
# the ‘from’ values and the second column represents the 
# ‘to’ values and the third column sets what those values 
# should be set to.
# so here we are setting anything 0-0.1885 to 0 and anything 0.1885-1 as 1 
# 0.1885 corresponds to Maximum training sensitivity plus 
# specificity Logistic threshold

m<-c(0, 0.1885, 0, 0.1885, 1, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rclmat

# Created matrix: 
#        [,1]   [,2] [,3]
# [1,] 0.0000 0.1885    0
# [2,] 0.1885 1.0000    1

# Reclasify the values in the output model 
# as either 0 (species absent) or 1 (species present) 
# We can plot out the resulting
# map of predicted presence and absence.

rc<- reclassify(complejo_maxent, rclmat)
# par(bg = 'gray') 
plot(rc)

rc_presente<- reclassify(complejo_maxent_presente, rclmat)
# par(bg = 'gray') 
plot(rc_presente)

points(complejo$longitude, complejo$latitude, col="red", cex=.4, pch=16)


#### Uploading PAST MODELS 

# LGM CCSM
complejo_LGM_CCSM <- raster("./maxent_results/Diglossa_baritula_LGM_CCSM_median.asc")
plot(complejo_LGM_CCSM)
rc_LGM_CCSM<- reclassify(complejo_LGM_CCSM, rclmat)
par(bg = 'gray') 
plot(rc_LGM_CCSM)

# LGM MIROC
complejo_LGM_MIROC <- raster("./maxent_results/Diglossa_baritula_LGM_MIROC_median.asc")
plot(complejo_LGM_MIROC)
rc_LGM_MIROC<- reclassify(complejo_LGM_MIROC, rclmat)
par(bg = 'gray') 
plot(rc_LGM_MIROC)

# MH CCSM
complejo_MH_CCSM <- raster("./maxent_results/Diglossa_baritula_MH_CCSM_median.asc")
plot(complejo_MH_CCSM)
rc_MH_CCSM<- reclassify(complejo_MH_CCSM, rclmat)
par(bg = 'gray') 
plot(rc_MH_CCSM)

# MH MIROC
complejo_MH_MIROC <- raster("./maxent_results/Diglossa_baritula_MH_MIROC_median.asc")
plot(complejo_MH_MIROC)
rc_MH_MIROC<- reclassify(complejo_MH_MIROC, rclmat)
par(bg = 'gray') 
plot(rc_MH_MIROC)

# LIG
complejo_LIG <- raster("./maxent_results/Diglossa_baritula_LIG_median.asc")
plot(complejo_LIG)
rc_LIG<- reclassify(complejo_LIG, rclmat)
par(bg = 'gray') 
plot(rc_LIG)


#### PLOTS

# background
par(bg = 'white')

#using shape for geographic distribution (dba071dpgw.shp)
distribution <- readOGR(dsn="./dba071dpgw.shp")
plot(distribution)


# plot using distribution layer
plot(rc, col=c("lightgray", "lightgray")); plot(distribution, col="#669966", border="#669966", add=TRUE)
plot(rc, col=c("#CCFFCC", "#CCFFCC")); plot(distribution, col="#669966", border="#669966", add=TRUE)
plot(rc, col=c("#CCCC99", "#CCCC99")); plot(distribution, col="#669966", border="#669966", add=TRUE)
plot(rc, col=c("lightgoldenrod", "lightgoldenrod")); plot(distribution, col="#669966", border="#669966", add=TRUE)
plot(rc, col=c("#CCEEBB", "#CCEEBB")); plot(distribution, col="#669966", border="#669966", add=TRUE)


# plot using present (current BIOS) area
plot(rc, col=c("#CCEEBB", "#CCEEBB"),legend=FALSE); plot(rc_presente, col=c("#CCEEBB", "#669966"), add=TRUE, legend=FALSE)

# plot using M area and occurrence points
plot(rc, col=c("#CCEEBB", "#CCEEBB")); plot(M_buffer, border="#669966", lwd=2, add=TRUE); points(complejo[, 2:3], col="#669966", cex=.4, pch=17)

# Models 
plot(rc_LGM_CCSM, col=c("#CCEEBB", "#669966"), legend=FALSE)
plot(rc_LGM_MIROC, col=c("#CCEEBB", "#669966"), legend=FALSE)
plot(rc_MH_CCSM, col=c("#CCEEBB", "#669966"), legend=FALSE)
plot(rc_MH_MIROC, col=c("#CCEEBB", "#669966"), legend=FALSE)
plot(rc_LIG, col=c("#CCEEBB", "#669966"), legend=FALSE)



