##########################################
####### MORPHOLOGY ANALYSIS ##############
##########################################

library(ggfortify)
library(plotly)
library(cluster)
library(ggbiplot)
library(multcompView)
library(lsmeans)


setwd("C:/Users/alote/OneDrive/Documentos/upload_FigShare/morphology/data/mean_three")

setwd("H:/Mi unidad/correcciones_Dbaritula_PeerJ/Data")

setwd("I:/Mi unidad/correcciones_Dbaritula_PeerJ/Data")

morpho<-read.csv(file="morphology_mean.csv", header=TRUE)

summary(morpho)
str(morpho)
class(morpho)

### Values without log transformation
#subsets
Male<-subset(morpho, morpho$Gender=="M")
Female<-subset(morpho, morpho$Gender=="F")  

# Females and Males Wing Chord data

IDs<-morpho[,1:9]
WC<-morpho[,14]
morpho_WC<-cbind(IDs,WC)
Male_WC<-subset(morpho_WC, morpho_WC$Gender=="M")
Female_WC<-subset(morpho_WC, morpho_WC$Gender=="F")

####### log transformed data

morpho$log_BW <- log(morpho$BW)   
morpho$log_BL <- log(morpho$BL) 
morpho$log_BHL <- log(morpho$BHL) 
morpho$log_TL <- log(morpho$TL) 
Male_WC$log_WC <- log(Male_WC$WC)
Female_WC$log_WC <- log(Female_WC$WC) 

#### Data by subspecies

baritula <- subset(morpho, Subspecies == "baritula")
montana <- subset(morpho, Subspecies == "montana")
parva <- subset(morpho, Subspecies == "parva")

baritula_females <- subset(baritula, Gender == "F")
baritula_males <- subset(baritula, Gender == "M")

montana_females <- subset(montana, Gender == "F")
montana_males <- subset(montana, Gender == "M")

parva_females <- subset(parva, Gender == "F")
parva_males <- subset(parva, Gender == "M")



### Check normality of males and females

## Females
# plot
hist(Female$BW)
hist(Female$BL)
hist(Female$BHL)
hist(Female$TL)
hist(Female$WC)

## Shapiro test for females variables

shapiro.test(Female$BW)
shapiro.test(Female$BL)
shapiro.test(Female$BHL)
shapiro.test(Female$TL)
shapiro.test(Female$WC)

## Males
# plot
hist(Male$BW)
hist(Male$BL)
hist(Male$BHL)
hist(Male$TL)
hist(Male$WC)


# Shapiro test for males variables

shapiro.test(Male$BW)
shapiro.test(Male$BL)
shapiro.test(Male$BHL)
shapiro.test(Male$TL)
shapiro.test(Male$WC)


# Test the homocedasticity of variances:

var.test(Female$BW,Male$BW)
var.test(Female$BL,Male$BL)
var.test(Female$BHL,Male$BHL)
var.test(Female$TL,Male$TL)
var.test(Female$WC,Male$WC)

## t test to compare the sexes

t.test(morpho$BW~morpho$Gender)
t.test(morpho$BL~morpho$Gender)
t.test(morpho$BHL~morpho$Gender)
t.test(morpho$TL~morpho$Gender)
t.test(morpho$WC~morpho$Gender)

IDs<-morpho[,1:9]
WC<-morpho[,14]
morpho_WC<-cbind(IDs,WC)
Male_WC<-subset(morpho_WC, morpho_WC$Gender=="M")
Female_WC<-subset(morpho_WC, morpho_WC$Gender=="F")

# Females and males together

shapiro.test(morpho$BW)
shapiro.test(morpho$BL)
shapiro.test(morpho$BHL)
shapiro.test(morpho$TL)

hist(morpho$BW)
hist(morpho$BL)
hist(morpho$BHL)
hist(morpho$TL)

#Test the homocedasticity

bartlett.test(morpho$BW~morpho$Subspecies)
bartlett.test(morpho$BL~morpho$Subspecies)
bartlett.test(morpho$BHL~morpho$Subspecies)
bartlett.test(morpho$TL~morpho$Subspecies)


### Pearson test to explore collinearity between the measurements

morpho1 <- na.omit(morpho)
dat <- morpho1[,10:14]
res <- cor(dat)
round(res, 2)

cor(morpho1$BL, morpho1$BHL)

### ANOVA with the three Subspeciess

# ANOVA Bill weight

anovaBW<-aov(BW~Subspecies,data=morpho)
summary(anovaBW)

# Check residuals distribution
hist(anovaBW$residuals,col="gray",
     main="Bill weight ANOVA residuals distribution",
     xlab="Residuals",ylab="Frequency") 

# Test the normality of residuals
shapiro.test(anovaBW$residuals)

# Tukey test to examine which Subspeciess is different
TukeyHSD(anovaBW)

#Plot Bill weight ANOVA results
plot(TukeyHSD(anovaBW))


lsm<-lsmeans(anovaBW,"morpho$Subspeciess",adjust="tukey")
lsm
cld(lsm, alpha=.05,Letters=letters)


multcompBoxplot(BW~Subspecies,morpho,horizontal = TRUE,
                compFn ="TukeyHSD",sortFn ="mean", 
                decreasing = TRUE, 
                plotList =list(boxplot = list(fig =c(0, 0.75, 0, 1)), 
                               multcompTs = list(fig = c(0.7, 0.85, 0, 1)),
                               multcompLetters = list(fig = c(0.87, 0.97, 0.03, 0.98), 
                                                      fontsize = 20,fontface = "bold")))


###ANOVA BL

anovaBL<-aov(BL~Subspecies,data=morpho)
summary(anovaBL)

#Check residuals distribution
hist(anovaBL$residuals,col="gray",
     main="Bill lenght ANOVA residuals distribution",
     xlab="Residuals",ylab="") 

#Test the normality of residuals
shapiro.test(anovaBL$residuals)

# Tukey test to examine which Subspeciess is different
TukeyHSD(anovaBL)

#Plot Bill lenght ANOVA results
plot(TukeyHSD(anovaBL))


multcompBoxplot(BL~Subspecies,morpho,horizontal = TRUE,
                compFn ="TukeyHSD",sortFn ="mean", 
                decreasing = TRUE, 
                plotList =list(boxplot = list(fig =c(0, 0.75, 0, 1)), 
                               multcompTs = list(fig = c(0.7, 0.85, 0, 1)),
                               multcompLetters = list(fig = c(0.87, 0.97, 0.03, 0.98), 
                                                      fontsize = 20,fontface = "bold")))


###ANOVA BHL
anovaBHL<-aov(BHL~Subspecies,data=morpho)
summary(anovaBHL)

###ANOVA TL
anovaTL<-aov(TL~Subspecies,data=morpho)
summary(anovaTL)

#Check residuals distribution
hist(anovaTL$residuals,col="gray",
     main="Tarsus lenght ANOVA residuals distribution",
     xlab="Residuals",ylab="Frequency") 


#Test the normality of residuals
shapiro.test(anovaTL$residuals)

# Tarsus lenght ANOVA results
TukeyHSD(anovaTL)

#Plot Tarsus lenght ANOVA results
plot(TukeyHSD(anovaTL))


multcompBoxplot(TL~Subspecies,morpho,horizontal = TRUE,
                compFn ="TukeyHSD",sortFn ="mean", 
                decreasing = TRUE, 
                plotList =list(boxplot = list(fig =c(0, 0.75, 0, 1)), 
                               multcompTs = list(fig = c(0.7, 0.85, 0, 1)),
                               multcompLetters = list(fig = c(0.87, 0.97, 0.03, 0.98), 
                                                      fontsize = 20,fontface = "bold")))

### Females ANOVA WC 

anova_Female_WC<-aov(WC~Subspecies,data=Female_WC)
summary(anova_Female_WC)

### Males ANOVA WC 

anova_Male_WC<-aov(WC~Subspecies,data=Male_WC)
summary(anova_Male_WC)

####### ANOVA log transformed data

morpho$log_BW <- log(morpho$BW)   
morpho$log_BL <- log(morpho$BL) 
morpho$log_BHL <- log(morpho$BHL) 
morpho$log_TL <- log(morpho$TL) 
Male_WC$log_WC <- log(Male_WC$WC)
Female_WC$log_WC <- log(Female_WC$WC) 

###ANOVA log BW
anovalog_BW<-aov(morpho$log_BW~morpho$Subspecies)
summary(anovalog_BW)

#Check residuals distribution
hist(anovalog_BW$residuals,col="gray",
     main="Bill weight ANOVA residuals distribution",
     xlab="Residuals",ylab="Frequency") 

#Test the normality of residuals
shapiro.test(anovalog_BW$residuals)

# Tukey test to examine which Subspecies is different
TukeyHSD(anovalog_BW)

#Plot ANOVA results
plot(TukeyHSD(anovalog_BW))

###ANOVA log BL
anovalog_BL<-aov(morpho$log_BL~morpho$Subspecies)
summary(anovalog_BL)

#Check residuals distribution
hist(anovalog_BL$residuals,col="gray",
     main="Bill lenght ANOVA residuals distribution",
     xlab="Residuals",ylab="Frequency") 

#Test the normality of residuals
shapiro.test(anovalog_BL$residuals)

TukeyHSD(anovalog_BL)

plot(TukeyHSD(anovalog_BL))


###ANOVA log BHL
anovalog_BHL<-aov(morpho$log_BHL~morpho$Subspecies)
summary(anovalog_BHL)

###ANOVA log TL
anovalog_TL<-aov(morpho$log_TL~morpho$Subspecies)
summary(anovalog_TL)

#Check residuals distribution
hist(anovalog_TL$residuals,col="gray",
     main="Tarsus lenght ANOVA residuals distribution",
     xlab="Residuals",ylab="Frequency") 

#Test the normality of residuals
shapiro.test(anovalog_TL$residuals)

# Tukey test to examine which Subspecies is different
TukeyHSD(anovalog_TL)

plot(TukeyHSD(anovalog_TL))

###Females ANOVA WC
anovalog_FemaleWC<-aov(Female_WC$log_WC~Female_WC$Subspecies)
summary(anovalog_FemaleWC)


#Check residuals distribution
hist(anovalog_FemaleWC$residuals,col="gray",
     main="Females Wing Chord ANOVA residuals distribution",
     xlab="Residuals",ylab="Frequency") 

#Test the normality of residuals
shapiro.test(anovalog_FemaleWC$residuals)

# Tukey test to examine which Subspecies is different
TukeyHSD(anovalog_FemaleWC)

#Plot Bill weight ANOVA results
plot(TukeyHSD(anovalog_FemaleWC))


###Males ANOVA WC
anovalog_MaleWC<-aov(Male_WC$log_WC~Male_WC$Subspecies)
summary(anovalog_MaleWC)



################ Boxplots by Subspecies ################


# vector to specific the colors
colors.box<-c("darkorange1","deepskyblue3","gold")

windows(5,10)

# Bill width
dev.new(width=5, height=10, unit="in")
boxplot(morpho$BW ~ morpho$Subspecies, las=1, col=colors.box, notch=T, lwd = 3, 
        ylab = "mm", xlab = "", names = labels, cex.lab = 2, xaxt = "n",
        cex.axis = 1.7, boxwex = 0.5, boxlwd = 4)
title("Bill width", cex.main=3.5)
mtext("D. b. baritula                  D. b. montana                 D. b. parva", 
      side = 1, font = 3, cex = 2.7, line = 1)


# Bill length 
dev.new(width=5, height=10, unit="in")
boxplot(morpho$BL ~ morpho$Subspecies, las=1, col=colors.box, notch=T, lwd = 3, 
        ylab = "mm", xlab = "", names = labels, cex.lab = 2, xaxt = "n", 
        cex.axis = 1.7, boxwex = 0.5, boxlwd = 4)
title("Bill length", cex.main=3.5)
mtext("D. b. baritula                  D. b. montana                 D. b. parva", 
      side = 1, font = 3, cex = 2.7, line = 1)

# Bill hook length
dev.new(width=5, height=10, unit="in")
boxplot(morpho$BHL ~ morpho$Subspecies, las=1, col=colors.box, notch=T, lwd = 3, 
        ylab = "mm", xlab = "", names = labels, cex.lab = 2, xaxt = "n", 
        cex.axis = 1.7, boxwex = 0.5, boxlwd = 4)
title("Bill hook length", cex.main=3.5)
mtext("D. b. baritula                  D. b. montana                 D. b. parva", 
      side = 1, font = 3, cex = 2.7, line = 1)

# Tarsus length 
dev.new(width=5, height=10, unit="in")
boxplot(morpho$TL ~ morpho$Subspecies, las=1, col=colors.box, notch=T, lwd = 3, 
        ylab = "mm", xlab = "", names = labels, cex.lab = 2, xaxt = "n", 
        cex.axis = 1.7, boxwex = 0.5, boxlwd = 4)
title("Tarsus length", cex.main=3.5)
mtext("D. b. baritula                  D. b. montana                 D. b. parva", 
      side = 1, font = 3, cex = 2.7, line = 1)

# Females wing chord length
dev.new(width=5, height=10, unit="in")
boxplot(Female$WC ~ Female$Subspecies, las=1, col=colors.box, notch=T, lwd = 3, 
        ylab = "mm", xlab = "", names = labels, cex.lab = 2, xaxt = "n", 
        cex.axis = 1.7, boxwex = 0.5, boxlwd = 4)
title("Females wing chord length", cex.main=3.5)
mtext("D. b. baritula                  D. b. montana                 D. b. parva", 
      side = 1, font = 3, cex = 2.7, line = 1)


# Males wing chord length
dev.new(width=5, height=10, unit="in")
boxplot(Male$WC ~ Male$Subspecies, las=1, col=colors.box, notch=T, lwd = 3, 
        ylab = "mm", xlab = "", names = labels, cex.lab = 2, xaxt = "n", 
        cex.axis = 1.7, boxwex = 0.5, boxlwd = 4)
title("Males wing chord length", cex.main=3.5)
mtext("D. b. baritula                  D. b. montana                 D. b. parva", 
      side = 1, font = 3, cex = 2.7, line = 1)






############### PCA for all variables except WC
### PCA with bbbiplot


biplot = ggbiplot(morpho1.pca,
                  choice=c(1,2),
                  obs.scale = 1, 
                  var.scale = 1,
                  groups=morpho1$Subspecies, 
                  ellipse=TRUE, 
                  colour=colors.pca, 
                  varname.size = 5,
                  points_size = 4,
                  ellipse_size = 0.8,
                  arrows_size = 0.5,
                  labels_textsize = 3,
                  base_textsize = 10)

biplot + scale_fill_manual(values = colors.pca) + 
  scale_color_manual(values = colors.pca) + 
  theme(panel.background = element_blank()) + 
  theme_bw() 




########## Mean and standard deviation ##################

class(morpho$Subspecies)

baritula <- subset(morpho, Subspecies == "baritula")
montana <- subset(morpho, Subspecies == "montana")
parva <- subset(morpho, Subspecies == "parva")

baritula_females <- subset(baritula, Gender == "F")
baritula_males <- subset(baritula, Gender == "M")

montana_females <- subset(montana, Gender == "F")
montana_males <- subset(montana, Gender == "M")

parva_females <- subset(parva, Gender == "F")
parva_males <- subset(parva, Gender == "M")


#### Mean baritula

mean(baritula$BW,na.rm=TRUE)
mean(baritula$BL,na.rm=TRUE)
mean(baritula$BHL,na.rm=TRUE)
mean(baritula$TL,na.rm=TRUE)
mean(baritula_females$WC,na.rm=TRUE)
mean(baritula_males$WC,na.rm=TRUE)


#### Mean montana

mean(montana$BW,na.rm=TRUE)
mean(montana$BL,na.rm=TRUE)
mean(montana$BHL,na.rm=TRUE)
mean(montana$TL)
mean(montana_females$WC)
mean(montana_males$WC)

#### Mean parva

mean(parva$BW,na.rm=TRUE)
mean(parva$BL,na.rm=TRUE)
mean(parva$BHL,na.rm=TRUE)
mean(parva$TL)
mean(parva_females$WC)
mean(parva_males$WC)


#### Standard deviation baritula

sd(baritula$BW,na.rm=TRUE)
sd(baritula$BL,na.rm=TRUE)
sd(baritula$BHL,na.rm=TRUE)
sd(baritula$TL)
sd(baritula_females$WC)
sd(baritula_males$WC)

#### Standard deviation montana

sd(montana$BW,na.rm=TRUE)
sd(montana$BL,na.rm=TRUE)
sd(montana$BHL,na.rm=TRUE)
sd(montana$TL)
sd(montana_females$WC)
sd(montana_males$WC)

#### Standard deviation  parva

sd(parva$BW,na.rm=TRUE)
sd(parva$BL,na.rm=TRUE)
sd(parva$BHL,na.rm=TRUE)
sd(parva$TL)
sd(parva_females$WC)
sd(parva_males$WC)