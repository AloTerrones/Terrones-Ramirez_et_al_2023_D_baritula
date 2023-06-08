##########################################
########## COLORATION ANALYSIS ###########

library(pavo)
library(ggplot2)
require(MCMCglmm)
require(scatterplot3d)
require(gridExtra)
require(vegan)
require(biotools)
require(MASS)
require(pbmcapply)
require(RColorBrewer)
require(dplyr)



Subespecie<-getspec("I:/Mi unidad/Proyecto_Dbaritula_Alo/Paper_multidisciplinario/Coloration/measures/recent_figshare/raw_color_data/Subspecies",                    
                    subdir.names=TRUE, subdir=T,lim = c(300, 700))



# Three measures average
Subespecie<-aggspec(Subespecie, by=3, FUN=mean)

#Male Subset
Male<-subset(Subespecie,subset = "male")
# To smooth and correct the mistake
Male<-procspec(Male,opt=c("smooth",span=0.25),fixneg="zero")
plot(Male)

#males by subspecies
baritula<-subset(Male,subset="baritula")
plot(baritula, main="baritula males")
montana<-subset(Male,subset="montana")
plot(montana, main="montana males")
parva<-subset(Male,subset="parva")
plot(parva, main="parva males")

#males, subespecie, parche
#baritula

esa.bari<-subset(baritula, subset = "esa")
plot(esa.bari, main="baritula male upper back")

esb.bari<-subset(baritula, subset = "esb")
plot(esb.bari, main="baritula male lower back")

gar.bari<-subset(baritula, subset = "gar")
plot(gar.bari, main="baritula male throat")

pec.bari<-subset(baritula, subset = "pec")
plot(pec.bari, main="baritula male breast")

via.bari<-subset(baritula, subset = "via")
plot(via.bari, main="baritula male upper belly")

vib.bari<-subset(baritula, subset = "vib")
plot(vib.bari, main="baritula male lower belly")


#montana

esa.mont<-subset(montana, subset = "esa")
plot(esa.mont, main="montana male upper back")

esb.mont<-subset(montana, subset = "esb")
plot(esb.mont, main="montana male lower back")

gar.mont<-subset(montana, subset = "gar")
plot(gar.mont, main="montana male throat")

pec.mont<-subset(montana, subset = "pec")
plot(pec.mont, main="montana male breast")

via.mont<-subset(montana, subset = "via")
plot(via.mont, main="montana male upper belly")

vib.mont<-subset(montana, subset = "vib")
plot(vib.mont, main="montana male lower belly")


#parva

esa.par<-subset(parva, subset = "esa")
plot(esa.par, main="parva male upper back")

esb.par<-subset(parva, subset = "esb")
plot(esb.par, main="parva male lower back")

gar.par<-subset(parva, subset = "gar")
plot(gar.par, main="parva male throat")

pec.par<-subset(parva, subset = "pec")
plot(pec.par, main="parva male breast")

via.par<-subset(parva, subset = "via")
plot(via.par, main="parva male upper belly")

vib.par<-subset(parva, subset = "vib")
plot(vib.par, main="parva male lower belly")





#### REFLECTANCE CURVES #####

# Subset by patch

# vector to specific the colors
color.sub<-c("darkorange1","deepskyblue3","gold")

#upper back
uback<-subset(Male, subset = "esa")
plot(uback)
male.vec<-do.call(rbind, strsplit(names(uback),"\\/"))[,1]
par(oma=c(0,0,0,0)) # all sides have 3 lines of space
par(mar=c(5,5,4,1) + 0.1)
aggplot(uback,by=male.vec, legend = F, lcol=color.sub, shadecol=color.sub, ylim=c(0, 30),  
        main= "Males upper back", lwd=2, cex.main=4, cex.lab=2, cex.axis=1.7, las=1) 


#lower back
lback<-subset(Male, subset = "esb")
plot(lback)
male.vec<-do.call(rbind, strsplit(names(lback),"\\/"))[,1]
par(oma=c(0,0,0,0)) # all sides have 3 lines of space
par(mar=c(5,5,4,1) + 0.1)
aggplot(lback,by=male.vec, legend = F, lcol=color.sub, shadecol=color.sub, ylim=c(0, 30), 
        main= "Males lower back",  lwd=2, cex.main=4, cex.lab=2, cex.axis=1.7, las=1, ylab="")


#throat
throat<-subset(Male, subset = "gar")
plot(throat)
male.vec<-do.call(rbind, strsplit(names(throat),"\\/"))[,1]
par(oma=c(0,0,0,0)) # all sides have 3 lines of space
par(mar=c(5,5,4,1) + 0.1)
aggplot(throat,by=male.vec, legend = F, lcol=color.sub, shadecol=color.sub, ylim=c(0, 30), 
        main= "Males throat", lwd=2, cex.main=4, cex.lab=2, cex.axis=1.7, las=1, ylab="")

#breast
breast<-subset(Male, subset = "pec")
plot(breast)
male.vec<-do.call(rbind, strsplit(names(breast),"\\/"))[,1]
par(oma=c(0,0,0,0)) # all sides have 3 lines of space
par(mar=c(5,5,4,1) + 0.1)
aggplot(breast,by=male.vec, legend = F, lcol=color.sub, shadecol=color.sub, ylim=c(0, 30), 
        main= "Males breast", lwd=2, cex.main=4, cex.lab=2, cex.axis=1.7, las=1)

#upper belly
ubelly<-subset(Male, subset = "via")
plot(ubelly)
male.vec<-do.call(rbind, strsplit(names(ubelly),"\\/"))[,1]
par(oma=c(0,0,0,0)) # all sides have 3 lines of space
par(mar=c(5,5,4,1) + 0.1)
aggplot(ubelly,by=male.vec, legend = F, lcol=color.sub, shadecol=color.sub, ylim=c(0, 30),  
        main= "Males upper belly", lwd=2, cex.main=4, cex.lab=2, cex.axis=1.7, las=1, ylab="")

#lower belly
lbelly<-subset(Male, subset = "vib")
plot(lbelly)
male.vec<-do.call(rbind, strsplit(names(lbelly),"\\/"))[,1]
par(oma=c(0,0,0,0)) # all sides have 3 lines of space
par(mar=c(5,5,4,1) + 0.1)
aggplot(lbelly,by=male.vec, legend = F, lcol=color.sub, shadecol=color.sub, ylim=c(0, 30), 
        main= "Males lower belly", lwd=2, cex.main=4, cex.lab=2, cex.axis=1.7, las=1, ylab="")



################ JND ####################

Subespecie<-getspec("I:/Mi unidad/Proyecto_Dbaritula_Alo/Paper_multidisciplinario/Coloration/measures/recent_figshare/raw_color_data/Subspecies",                    
                    subdir.names=TRUE, subdir=T,lim = c(300, 700))


# Three measures average
Subespecie<-aggspec(Subespecie, by=3, FUN=mean)



##### upper back males ####

Male<-subset(Subespecie,subset = "male")
Male<-procspec(Male,opt=c("smooth",span=0.25),fixneg="zero")

plot(Male)

esa.bari<-subset(Male, subset = "esa")
plot(esa.bari, main="baritula males upper back")


#### the distances
data(esa.bari)
vis.esa.bari <- vismodel(esa.bari)
cd.esa.bari <- coldist( vis.esa.bari)
jnd_bari_male_esa<-jnd2xyz(cd.esa.bari)


specs.esa <- as.rspec(esa.bari, interp = FALSE)

esa_vis <- sensmodel(c(371, 448, 502, 563),
                     beta = FALSE,
                     lambdacut = c(330, 413, 507, 572),
                     oiltype = c("T", "C", "Y", "R"), om = TRUE
) 

names(esa_vis) <- c("wl", "u", "s", "m", "l")
model.esa <- vismodel(specs.esa, visual = "bluetit", relative = FALSE) 
class(model.esa)
space.esa <- colspace(model.esa)
deltaS.esa <- coldist(model.esa, achro = TRUE, 
                      weber = 0.1, noise = "neural")
mat.esa <- dist(coldist2mat(deltaS.esa)[["dS"]])
df.model.esa<-as.data.frame(model.esa)
df.esa <- tibble::rownames_to_column(df.model.esa, "ID")
x.esa <- df.esa[,1]
pat.esa <- "(baritula|montana|parva)"
x.esa <- transform(x.esa,subspecie=stringr::str_extract(x.esa,pat.esa))
subspecie.esa<-as.vector(x.esa$subspecie)
gr.esa <- subspecie.esa
class(gr.esa)
bdisp.esa <- betadisper(mat.esa, gr.esa, type="centroid")
anova(bdisp.esa)
 
# Difference between dispersions
TukeyHSD(bdisp.esa)

table(gr.esa)
pmanova.esa <- adonis(mat.esa~gr.esa)
pmanova.esa


bootds.esa<-bootcoldist(model.esa, by = gr.esa, n = c(1, 2, 2, 4), weber = 0.1, 
                        weber.achro = 0.1)
bootds.esa
plot.esa<-plot(bootds.esa[,1], ylim=c(0, 8), pch=21, bg=1, cex=2, xaxt="n", xlab="", 
               ylab="", lwd=2, cex.main=4, cex.axis=1.7, las=1)
title(ylab="Chromatic contrast (JND)", line=2, cex.lab=2)
mtext("baritula-montana    baritula-parva    montana-parva", side = 1, font = 3, cex = 2, line = 1)
abline(h=1, lty=3, lwd=2)
segments(1, bootds.esa[1,2], 1, bootds.esa[1,3], lwd=2)
segments(2, bootds.esa[2,2], 2, bootds.esa[2,3], lwd=2)
segments(3, bootds.esa[3,2], 3, bootds.esa[3,3], lwd=2)




#### lower back males #####

Male<-subset(Subespecie,subset = "male")
Male<-procspec(Male,opt=c("smooth",span=0.25),fixneg="zero")
plot(Male)

esb.bari<-subset(Male, subset = "esb")
plot(esb.bari, main="baritula males lower back")

plot(esb.bari, col = spec2rgb(esb.bari))

#### the distances
data(esb.bari)
vis.esb.bari <- vismodel(esb.bari)
cd.esb.bari <- coldist( vis.esb.bari)
jnd_bari_male_esb<-jnd2xyz(cd.esb.bari)


specs.esb <- as.rspec(esb.bari, interp = FALSE)

esb_vis <- sensmodel(c(371, 448, 502, 563),
                     beta = FALSE,
                     lambdacut = c(330, 413, 507, 572),
                     oiltype = c("T", "C", "Y", "R"), om = TRUE
) 

names(esb_vis) <- c("wl", "u", "s", "m", "l")
model.esb <- vismodel(specs.esb, visual = "bluetit", relative = FALSE) 
class(model.esb)
space.esb <- colspace(model.esb)
deltaS.esb <- coldist(model.esb, achro = TRUE, 
                      weber = 0.1, noise = "neural")
mat.esb <- dist(coldist2mat(deltaS.esb)[["dS"]])
df.model.esb<-as.data.frame(model.esb)
df.esb <- tibble::rownames_to_column(df.model.esb, "ID")
x.esb <- df.esb[,1]
pat.esb <- "(baritula|montana|parva)"
x.esb <- transform(x.esb,subspecie=stringr::str_extract(x.esb,pat.esb))
subspecie.esb<-as.vector(x.esb$subspecie)
gr.esb <- subspecie.esb
class(gr.esb)
bdisp.esb <- betadisper(mat.esb, gr.esb, type="centroid")
anova(bdisp.esb)
 


# Difference between dispersions
TukeyHSD(bdisp.esb)

table(gr.esb)
pmanova.esb <- adonis(mat.esb~gr.esb)
pmanova.esb


bootds.esb<-bootcoldist(model.esb, by = gr.esb, n = c(1, 2, 2, 4), weber = 0.1, 
                        weber.achro = 0.1)
bootds.esb
plot.esb<-plot(bootds.esb[,1], ylim=c(0, 8), pch=21, bg=1, cex=2, xaxt="n", xlab="", 
                lwd=2, cex.main=4, cex.lab=2, cex.axis=1.7, las=1, ylab="")
mtext("baritula-montana    baritula-parva    montana-parva", side = 1, font = 3, cex = 2, line = 1)
abline(h=1, lty=3, lwd=2)
segments(1, bootds.esb[1,2], 1, bootds.esb[1,3], lwd=2)
segments(2, bootds.esb[2,2], 2, bootds.esb[2,3], lwd=2)
segments(3, bootds.esb[3,2], 3, bootds.esb[3,3], lwd=2)



#### Males throat ####

Male<-subset(Subespecie,subset = "male")
Male<-procspec(Male,opt=c("smooth",span=0.25),fixneg="zero")
plot(Male)

gar.bari<-subset(Male, subset = "gar")
plot(gar.bari, main="baritula male throat")

#### the distances
data(gar.bari)
vis.gar.bari <- vismodel(gar.bari)
cd.gar.bari <- coldist( vis.gar.bari)
jnd_bari_male_gar<-jnd2xyz(cd.gar.bari)

specs.gar <- as.rspec(gar.bari, interp = FALSE)

gar_vis <- sensmodel(c(371, 448, 502, 563),
                     beta = FALSE,
                     lambdacut = c(330, 413, 507, 572),
                     oiltype = c("T", "C", "Y", "R"), om = TRUE
) 

names(gar_vis) <- c("wl", "u", "s", "m", "l")
model.gar <- vismodel(specs.gar, visual = "bluetit", relative = FALSE) 
class(model.gar)
space.gar <- colspace(model.gar)
deltaS.gar <- coldist(model.gar, achro = TRUE, 
                      weber = 0.1, noise = "neural")
mat.gar <- dist(coldist2mat(deltaS.gar)[["dS"]])
df.model.gar<-as.data.frame(model.gar)
df.gar <- tibble::rownames_to_column(df.model.gar, "ID")
x.gar <- df.gar[,1]
pat.gar <- "(baritula|montana|parva)"
x.gar <- transform(x.gar,subspecie=stringr::str_extract(x.gar,pat.gar))
subspecie.gar<-as.vector(x.gar$subspecie)
gr.gar <- subspecie.gar
class(gr.gar)
bdisp.gar <- betadisper(mat.gar, gr.gar, type="centroid")
anova(bdisp.gar)


# Difference between dispersions
TukeyHSD(bdisp.gar)
table(gr.gar)
pmanova.gar <- adonis(mat.gar~gr.gar)
pmanova.gar


bootds.gar<-bootcoldist(model.gar, by = gr.gar, n = c(1, 2, 2, 4), weber = 0.1, 
                        weber.achro = 0.1)

bootds.gar
plot.gar<-plot(bootds.gar[,1], ylim=c(0, 8), pch=21, bg=1, cex=2, xaxt="n", xlab="",
               lwd=2, cex.main=4, cex.lab=2, cex.axis=1.7, las=1, ylab="")
mtext("baritula-montana    baritula-parva    montana-parva", side = 1, font = 3, cex = 2, line = 1)
abline(h=1, lty=3, lwd=2)
segments(1, bootds.gar[1,2], 1, bootds.gar[1,3], lwd=2)
segments(2, bootds.gar[2,2], 2, bootds.gar[2,3], lwd=2)
segments(3, bootds.gar[3,2], 3, bootds.gar[3,3], lwd=2)



#### breast males ####

Male<-subset(Subespecie,subset = "male")
Male<-procspec(Male,opt=c("smooth",span=0.25),fixneg="zero")
plot(Male)

pec.bari<-subset(Male, subset = "pec")
plot(pec.bari, main="baritula males breast")

#### the distances
data(pec.bari)
vis.pec.bari <- vismodel(pec.bari)
cd.pec.bari <- coldist( vis.pec.bari)
jnd_bari_male_pec<-jnd2xyz(cd.pec.bari)


specs.pec <- as.rspec(pec.bari, interp = FALSE)

pec_vis <- sensmodel(c(371, 448, 502, 563),
                     beta = FALSE,
                     lambdacut = c(330, 413, 507, 572),
                     oiltype = c("T", "C", "Y", "R"), om = TRUE
) 

names(pec_vis) <- c("wl", "u", "s", "m", "l")
model.pec <- vismodel(specs.pec, visual = "bluetit", relative = FALSE) 
class(model.pec)
space.pec <- colspace(model.pec)
deltaS.pec <- coldist(model.pec, achro = TRUE, 
                      weber = 0.1, noise = "neural")
mat.pec <- dist(coldist2mat(deltaS.pec)[["dS"]])
df.model.pec<-as.data.frame(model.pec)
df.pec <- tibble::rownames_to_column(df.model.pec, "ID")
x.pec <- df.pec[,1]
pat.pec <- "(baritula|montana|parva)"
x.pec <- transform(x.pec,subspecie=stringr::str_extract(x.pec,pat.pec))
subspecie.pec<-as.vector(x.pec$subspecie)
gr.pec <- subspecie.pec
class(gr.pec)
bdisp.pec <- betadisper(mat.pec, gr.pec, type="centroid")
anova(bdisp.pec)
 
# Difference between dispersions
TukeyHSD(bdisp.pec)

table(gr.pec)
pmanova.pec <- adonis(mat.pec~gr.pec)
pmanova.pec


bootds.pec<-bootcoldist(model.pec, by = gr.pec, n = c(1, 2, 2, 4), weber = 0.1, 
                        weber.achro = 0.1)
bootds.pec
plot.pec<-plot(bootds.pec[,1], ylim=c(0, 8), pch=21, bg=1, cex=2, xaxt="n", xlab="", 
               ylab="", lwd=2, cex.main=4, cex.axis=1.7, las=1)
title(ylab="Chromatic contrast (JND)", line=2, cex.lab=2)
mtext("baritula-montana    baritula-parva    montana-parva", side = 1, font = 3, cex = 2, line = 1)
abline(h=1, lty=3, lwd=2)
segments(1, bootds.pec[1,2], 1, bootds.pec[1,3], lwd=2)
segments(2, bootds.pec[2,2], 2, bootds.pec[2,3], lwd=2)
segments(3, bootds.pec[3,2], 3, bootds.pec[3,3], lwd=2)



#### upper belly males #####

Male<-subset(Subespecie,subset = "male")
Male<-procspec(Male,opt=c("smooth",span=0.25),fixneg="zero")
plot(Male)

via.bari<-subset(Male, subset = "via")
plot(via.bari, main="baritula males upper belly")


#### the distances
data(via.bari)
vis.via.bari <- vismodel(via.bari)
cd.via.bari <- coldist( vis.via.bari)
jnd_bari_male_via<-jnd2xyz(cd.via.bari)


specs.via <- as.rspec(via.bari, interp = FALSE)

via_vis <- sensmodel(c(371, 448, 502, 563),
                     beta = FALSE,
                     lambdacut = c(330, 413, 507, 572),
                     oiltype = c("T", "C", "Y", "R"), om = TRUE
) 

names(via_vis) <- c("wl", "u", "s", "m", "l")
model.via <- vismodel(specs.via, visual = "bluetit", relative = FALSE) 
class(model.via)
space.via <- colspace(model.via)
deltaS.via <- coldist(model.via, achro = TRUE, 
                      weber = 0.1, noise = "neural")
mat.via <- dist(coldist2mat(deltaS.via)[["dS"]])
df.model.via<-as.data.frame(model.via)
df.via <- tibble::rownames_to_column(df.model.via, "ID")
x.via <- df.via[,1]
pat.via <- "(baritula|montana|parva)"
x.via <- transform(x.via,subspecie=stringr::str_extract(x.via,pat.via))
subspecie.via<-as.vector(x.via$subspecie)
gr.via <- subspecie.via
class(gr.via)
bdisp.via <- betadisper(mat.via, gr.via, type="centroid")
anova(bdisp.via)
 


# Difference between dispersions
TukeyHSD(bdisp.via)

table(gr.via)
pmanova.via <- adonis(mat.via~gr.via)
pmanova.via


bootds.via<-bootcoldist(model.via, by = gr.via, n = c(1, 2, 2, 4), weber = 0.1, 
                        weber.achro = 0.1)
bootds.via
plot.via<-plot(bootds.via[,1], ylim=c(0, 8), pch=21, bg=1, cex=2, xaxt="n", xlab="", 
               lwd=2, cex.main=4, cex.lab=2, cex.axis=1.7, las=1, ylab="")
mtext("baritula-montana    baritula-parva    montana-parva", side = 1, font = 3, cex = 2, line = 1)
abline(h=1, lty=3, lwd=2)
segments(1, bootds.via[1,2], 1, bootds.via[1,3], lwd=2)
segments(2, bootds.via[2,2], 2, bootds.via[2,3], lwd=2)
segments(3, bootds.via[3,2], 3, bootds.via[3,3], lwd=2)



#### lower belly males ####

Male<-subset(Subespecie,subset = "male")
Male<-procspec(Male,opt=c("smooth",span=0.25),fixneg="zero")
plot(Male)

vib.bari<-subset(Male, subset = "vib")
plot(vib.bari, main="baritula males lower belly")

#### the distances
data(vib.bari)
vis.vib.bari <- vismodel(vib.bari)
cd.vib.bari <- coldist( vis.vib.bari)
jnd_bari_male_vib<-jnd2xyz(cd.vib.bari)


specs.vib <- as.rspec(vib.bari, interp = FALSE)

vib_vis <- sensmodel(c(371, 448, 502, 563),
                     beta = FALSE,
                     lambdacut = c(330, 413, 507, 572),
                     oiltype = c("T", "C", "Y", "R"), om = TRUE
) 

names(vib_vis) <- c("wl", "u", "s", "m", "l")
model.vib <- vismodel(specs.vib, visual = "bluetit", relative = FALSE) 
class(model.vib)
space.vib <- colspace(model.vib)
deltaS.vib <- coldist(model.vib, achro = TRUE, 
                      weber = 0.1, noise = "neural")
mat.vib <- dist(coldist2mat(deltaS.vib)[["dS"]])
df.model.vib<-as.data.frame(model.vib)
df.vib <- tibble::rownames_to_column(df.model.vib, "ID")
x.vib <- df.vib[,1]
pat.vib <- "(baritula|montana|parva)"
x.vib <- transform(x.vib,subspecie=stringr::str_extract(x.vib,pat.vib))
subspecie.vib<-as.vector(x.vib$subspecie)
gr.vib <- subspecie.vib
class(gr.vib)
bdisp.vib <- betadisper(mat.vib, gr.vib, type="centroid")
anova(bdisp.vib)
 


# Difference between dispersions
TukeyHSD(bdisp.vib)

table(gr.vib)
pmanova.vib <- adonis(mat.vib~gr.vib)
pmanova.vib


bootds.vib<-bootcoldist(model.vib, by = gr.vib, n = c(1, 2, 2, 4), weber = 0.1, 
                        weber.achro = 0.1)
bootds.vib
plot.vib<-plot(bootds.vib[,1], ylim=c(0, 8), pch=21, bg=1, cex=2, xaxt="n", xlab="", 
               lwd=2, cex.main=4, cex.lab=2, cex.axis=1.7, las=1, ylab="")
mtext("baritula-montana    baritula-parva    montana-parva", side = 1, font = 3, cex = 2, line = 1)
abline(h=1, lty=3, lwd=2)
segments(1, bootds.vib[1,2], 1, bootds.vib[1,3], lwd=2)
segments(2, bootds.vib[2,2], 2, bootds.vib[2,3], lwd=2)
segments(3, bootds.vib[3,2], 3, bootds.vib[3,3], lwd=2)







