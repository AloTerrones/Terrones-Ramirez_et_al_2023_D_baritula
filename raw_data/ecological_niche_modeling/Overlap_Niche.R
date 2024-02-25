# OVERLAP NICHE TEST Diglossa baritula complex

install.packages("biomod2", type = "binary")
install.packages("ecospat")
install.packages("mda")
install.packages("MODISTools")


library(ecospat) # Spatial Ecology Miscellaneous Methods, CRAN v3.1
#library(mda) # Mixture and Flexible Discriminant Analysis, CRAN v0.5
library(raster) # Geographic Data Analysis and Modeling, CRAN v3.1-5
#library(MODISTools) # Interface to the 'MODIS Land Products Subsets' Web Services, CRAN v1.1.1
library(ade4) #dudi.pca
library(factoextra)


# Following: plantarum.ca/notebooks/ecospat

# Coordinates for subspecies
baritula <- read.csv("./baritula.csv")
parva <- read.csv("./parva.csv")
montana <- read.csv("./montana.csv")


# Environmental data for "M" area -accesible area-
# (or study area)

# -- We need to generate a dataframe with bioclimatic variables for each pixel
# -- This dataframe corresponds to study area (M area, accessible area)
# --- Defining a variable for upload bioclimatic layers already cut to "M" size
# --- This layers come from ENM analysis

varclim <- stack(list.files(path="bios_M_todos/", pattern = "*.asc$",full.names = T)) 

# --- Creating a SpatialPoints DataFrame with values for each pixel in BIO 1 [[1]]
climpunto <- rasterToPoints(varclim[[1]], fun=NULL, spatial=TRUE)

# --- Extracting values from 19 bioclimatic layers for background
clim <- extract(varclim, climpunto)

# --- Adding coordinates to previous dataframe
# --- This will generate a dataframe with coordinates (x,y) and bioclimatic values
clim <- data.frame(coordinates(climpunto),clim)
clim <- na.omit(clim) # --- eliminate NA values (0 NAs)

# --- Change order of columns in "clim"
colnames(clim)
clim <- clim[, c(1,2,3,14:21,4:13)]
colnames(clim)

# --- Save dataframe
#write.csv(clim, file="0_overlap/clim_M.csv", row.names = F) 


### Add environmental data to subspecies data: 

# extract values for each database

baritula_coordinates <- subset(baritula, select=c(longitude, latitude))
parva_coordinates <- subset(parva, select=c(longitude, latitude))
montana_coordinates <- subset(montana, select=c(longitude, latitude))


baritula_v <- raster::extract(varclim, baritula_coordinates)
baritula_v <- cbind(baritula_coordinates, baritula_v)
colnames(baritula_v)
baritula_v <- baritula_v[, c(1,2,3,14:21,4:13)]
colnames(baritula_v)
#write.csv(baritula_v, "0_overlap/baritula_v.csv", row.names = FALSE)

parva_v <- raster::extract(varclim, parva_coordinates)
parva_v <- cbind(parva_coordinates, parva_v)
colnames(parva_v)
parva_v <- parva_v[, c(1,2,3,14:21,4:13)]
colnames(parva_v)
#write.csv(parva_v, "0_overlap/parva_v.csv", row.names = FALSE)

montana_v <- raster::extract(varclim, montana_coordinates)
montana_v <- cbind(montana_coordinates, montana_v)
colnames(montana_v)
montana_v <- montana_v[, c(1,2,3,14:21,4:13)]
colnames(montana_v)
#write.csv(montana_v, "0_overlap/montana_v.csv", row.names = FALSE)


#-----
## pca scores 
# 19 WorldClim variables

Xvar <- c("wc2.1_30s_bio_1", "wc2.1_30s_bio_2", 
          "wc2.1_30s_bio_3", "wc2.1_30s_bio_4", 
          "wc2.1_30s_bio_5", "wc2.1_30s_bio_6", 
          "wc2.1_30s_bio_7", "wc2.1_30s_bio_8",
          "wc2.1_30s_bio_9", "wc2.1_30s_bio_10", 
          "wc2.1_30s_bio_11", "wc2.1_30s_bio_12", 
          "wc2.1_30s_bio_13", "wc2.1_30s_bio_14",
          "wc2.1_30s_bio_15", "wc2.1_30s_bio_16",
          "wc2.1_30s_bio_17", "wc2.1_30s_bio_18",
          "wc2.1_30s_bio_19")
nvar <- length(Xvar)

pca.cal <- dudi.pca(clim[, Xvar], center = TRUE,
                    scale = TRUE, scannf = FALSE, nf = 2)
summary(pca.cal)


# --- Plot variable contribution
ecospat.plot.contrib(contrib=pca.cal$co, eigen=pca.cal$eig)

fviz_eig(pca.cal)

fviz_eig(pca.cal, choice = c("variance", "eigenvalue"),
         addlabels = TRUE)

# Visualization of variables (circular plot)
fviz_pca_var(pca.cal,
             col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE) +
    theme(text = element_text(size = 10),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12))

# Contributions of variables to PC1
fviz_contrib(pca.cal, choice = "var", axes = 1, top = 10)

# Contributions of variables to PC2
fviz_contrib(pca.cal, choice = "var", axes = 2, top = 10)



# subspecies scores
baritula.scores <- suprow(pca.cal, baritula_v[,Xvar])$li
montana.scores <- suprow(pca.cal, montana_v[,Xvar])$li
parva.scores <-suprow(pca.cal, parva_v[,Xvar])$li

scores.clim <- pca.cal$li



# Ordination data to get actual grids:

baritula.grid <- ecospat.grid.clim.dyn(scores.clim, scores.clim,
                            baritula.scores, R = 100)

montana.grid <- ecospat.grid.clim.dyn(scores.clim, scores.clim,
                                       montana.scores, R = 100)

parva.grid <- ecospat.grid.clim.dyn(scores.clim, scores.clim,
                                       parva.scores, R = 100)


par(mfrow = c(1, 1))
#### --- Plot of environmental space 
# -- one by one
ecospat.plot.niche (baritula.grid, title="baritula", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.niche (montana.grid, title="montana", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.niche (parva.grid, title="parva", name.axis1="PC1", name.axis2="PC2", cor=F)


# -- overlapped
ecospat.plot.niche.dyn (baritula.grid, montana.grid, 
                        title="baritula vs montana", quant=0.1, name.axis1="PC1", name.axis2="PC2", 
                        col.unf = "darkorange1", col.exp = "deepskyblue3", col.stab = "lightgray",  
                        transparency = 38)


ecospat.plot.niche.dyn (baritula.grid, parva.grid,
                        title="baritula vs parva", quant=0.1, name.axis1="PC1", name.axis2="PC2", 
                        col.unf="darkorange1", col.exp = "gold", col.stab ="lightgray", 
                        transparency = 38)


ecospat.plot.niche.dyn (montana.grid, parva.grid, 
                        title="montana vs parva", quant=0.1, name.axis1="PC1", name.axis2="PC2", 
                        col.unf="deepskyblue3", col.exp = "gold", col.stab ="lightgray", 
                        transparency = 38)


#### --- Calculate Niche overlap

ecospat.niche.overlap (baritula.grid, montana.grid, cor=T) 
#$D 0.2414861
#$I 0.3560643

ecospat.niche.overlap (baritula.grid, parva.grid, cor=T) 
#$D 0.06433742
#$I 0.2352018

ecospat.niche.overlap (montana.grid, parva.grid, cor=T) 
#$D 0.06209386
#$I 0.2383898


# -- Calculate niche overlap cor=F
# if cor=FALSE, the z$uncor objects created by ecospat.grid.clim are used to calculate the overlap. 
# if cor=TRUE, the z$cor objects are used.
ecospat.niche.overlap (baritula.grid, montana.grid, cor=F) 
#$D 0.4962726
#$I 0.7418024

ecospat.niche.overlap (baritula.grid, parva.grid, cor=F) 
#$D 0.2347161
#$I 0.4571215

ecospat.niche.overlap (montana.grid, parva.grid, cor=F) 
#$D 0.4495401
#$I 0.6638163



#### ---- Perform niche equivalence test
eq.test.bm <-ecospat.niche.equivalency.test(baritula.grid , montana.grid, rep=100) 
eq.test.bp <-ecospat.niche.equivalency.test(baritula.grid , parva.grid, rep=100) 
eq.test.mp <-ecospat.niche.equivalency.test(montana.grid , parva.grid, rep=100) 

#### ---- Plot equivalency test
ecospat.plot.overlap.test(eq.test.bm,"D","Equivalency: baritula vs montana")
ecospat.plot.overlap.test(eq.test.bp,"D","Equivalency: baritula vs parva")
ecospat.plot.overlap.test(eq.test.mp,"D","Equivalency: montana vs parva")



#### ---- Perform niche similarity test
sim.test.bm <-ecospat.niche.similarity.test(baritula.grid , montana.grid, rep=100) 
sim.test.bp <-ecospat.niche.similarity.test(baritula.grid , parva.grid, rep=100) 
sim.test.mp <-ecospat.niche.similarity.test(montana.grid , parva.grid, rep=100) 

sim.test.mb <-ecospat.niche.similarity.test(montana.grid , baritula.grid, rep=100) 
sim.test.pb <-ecospat.niche.similarity.test(parva.grid , baritula.grid, rep=100) 
sim.test.pm <-ecospat.niche.similarity.test(parva.grid , montana.grid, rep=100) 


#### --- Plot similarity test
ecospat.plot.overlap.test(sim.test.mb,"D","Similarity: montana vs baritula")
ecospat.plot.overlap.test(sim.test.pb,"D","Similarity: parva vs baritula")
ecospat.plot.overlap.test(sim.test.pm,"D","Similarity: parva vs montana")

