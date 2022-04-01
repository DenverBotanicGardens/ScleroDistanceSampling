library(devtools)
library(terra)
library(raster)
library(unmarked)
library(dplyr)
library(ggplot2)
library(plyr)

#read in the data
setwd("~/.ssh/ScleroDistanceSampling")
DSdata <- read.csv("DistSamp_ScGl_JuneJuly.csv")

#get only data from certain dates
DSdata$Date<-as.character(DSdata$Date)
DSdata$Site<- as.character(DSdata$Site)
#data cleaning, two spellings for picnic
DSdata$Site[grep("icnic",DSdata$Site)] <- "Picnic"
DSdata$Site[grep("unction", DSdata$Site)] <- "T-Junction"
dates<- c("8/4/21","8/5/21","8/6/21", "8/1/21")
DSdata<- DSdata %>%
  filter(DSdata$Date %in% dates)

#T-Junction n=25
nrow(DSdata[DSdata$Site=="T-Junction",])
#Picnic n=159
nrow(DSdata[DSdata$Site=="Picnic",])

#add unique identifier
DSdata$trans <- paste(DSdata$Site,DSdata$Line, sep = ".")
# We're adding more same plot and same transect number so paste the date on too for a unique Transect identifier
DSdata$trans <- paste(DSdata$trans, DSdata$Date, sep = ".")


################################################# COMBINED DATA SET ##################################################
DS<-DSdata
hist(DS$Dist[DS$Site == "Picnic"], breaks = 10, col = "blue")
hist(DS$Dist[DS$Site == "T-Junction"], breaks = 10, add = TRUE)
#VISUALIZE 
ggplot(DS, aes(Dist, colour = Site))+
  geom_density()+
  theme_bw()

DS$trans <- paste(DS$Site,DS$Line, sep = ".")
#UID
DS$trans <- paste(DS$trans, DS$Date, sep = ".")

DS.ordered <- DS[with(DS, order(DS$Site, DS$Date, DS$Dist)),]
str(DS.ordered) # $trans is currently a character, could be anything. Setting as a factor restricts it to only the current values (character strings)
DS.ordered$trans <- as.factor(DS.ordered$trans)

## TRUNCATION

#COMBINED
DS.trunc <- DS[DS$Dist < 6,]
#this is how it looks with 0.25 cm binning
par(1,2)
ggplot(DS.trunc, aes(Dist, fill =  Site))+
  geom_histogram(bins=24)+
  facet_wrap(~Site)+
  theme_bw()


#format the data for umf, set correct truncation and binning
yDat <- formatDistData(DS.trunc, distCol = "Dist", transectNameCol = "trans", dist.breaks = seq(0,6,0.25))


#COVARIATES aka cover
table(DS$trans, DS$Length)
identical(row.names(table(DS$trans, DS$Length)), row.names(yDat))
# get the column names (the lengths) in the order of the rows, when >0
x <- table(DS$trans, DS$Length)[1,]

siteCovs <- data.frame(Site = unlist(lapply(row.names(yDat), function(x) strsplit(x,"[.]")[[1]][1])),
                       length = 
                         c(apply(table(DS$trans, DS$Length), 1, function(x) as.numeric(as.character(names(x[x>0]))))))
siteCovs$Cover <- ifelse(siteCovs$Site %in% "Picnic", 0.175633,0.480958)

#format into umf format
umf <- unmarkedFrameDS(y = as.matrix(yDat), 
                       siteCovs = siteCovs,
                       survey = "line",
                       dist.breaks = seq(0,6,0.25), 
                       tlength= siteCovs$length, 
                       unitsIn = "m")
summary(umf)
#VISUALIZE
ggplot(DS, aes(Dist, fill =  Site))+
  geom_histogram()+
  facet_wrap(~Site)+
  theme_bw()

###################################################### MODELS #########################################################

m.half <- distsamp(~1 ~1, umf, keyfun="halfnorm", output="density", unitsOut="kmsq")

m.haz <- distsamp(~1 ~1, umf, keyfun="hazard", output="density", unitsOut="kmsq")

m.uni <- distsamp(~1 ~1, umf, keyfun="uniform", output="density", unitsOut="kmsq")
#M.HALF AIC =
m.half
#m.haz AIC = 
m.haz
#m.uni AIC = 
m.uni
library(MuMIn)
model.sel(m.half,m.haz,m.uni)

#USING COVER AS A COVARIATE
m.haz.1.cover <- distsamp(~1 ~Cover, umf, keyfun="hazard", output="density", unitsOut="kmsq") 
m.haz.1.cover2 <- distsamp(~Cover ~1, umf, keyfun="hazard", output="density", unitsOut="kmsq") 
m.haz.1.cover3 <- distsamp(~Cover ~Cover, umf, keyfun="hazard", output="density", unitsOut="kmsq") 
#used 0 here as a start value for cover, 0 behaves the same as -1 and 0.5, they all give the same answer.
m.haz.1.cover #improves here
here<- data.frame(model.sel(m.haz,m.haz.1.cover,m.haz.1.cover2,m.haz.1.cover3, m.uni,m.half))
here$delta

#BEST
m.haz.1.cover
exp(-4.418716)
exp(-1.655486)


##################################################  RESULTS  ##########################################################

#DENSITY BACK TRANSFORM
#when cover is 0.17
lcP <- linearComb(m.haz.1.cover, c(1, 0.176), type="state") 
lcP<-backTransform(lcP)
lcP
#when cover is 0.481
lcJ <- linearComb(m.haz.1.cover, c(1, 0.481), type="state") 
lcJ<-backTransform(lcJ)
lcJ

#THIS IS THE SCALE
backTransform(m.haz.1.cover, type="det")
#THIS IS THE SHAPE
backTransform(m.haz.1.cover, type="scale")

#CI's
#scale
confint(m.haz.1.cover, type="det")
exp(-0.02975884)
exp(0.3342364)
#shape
confint(m.haz.1.cover, type="scale")
exp(0.7466846)
exp(1.152163)

#######  AVERAGE DETECTION PROBABILITY - Formula from Chandler Richard (works with Royle and Kery) #######
halfwidth <- 6 ## Transect half width
integrate(function(x) gxhaz(x, shape=3.17, scale=1.40), lower=0, upper=halfwidth)$value / halfwidth


#not sure what this number means for abundance, total abundance estimated below
m.haz.1.site.abund <- distsamp(~1 ~Cover, umf, keyfun="hazard", output="abund", unitsOut="kmsq",starts=c(0,0,0,0)) 
m.haz.1.site.abund
lcJabund <- linearComb(m.haz.1.site.abund, c(1, 0.481), type="state") 
backTransform(lcJabund)

#PREDICTING and FIGURES
#https://cran.r-project.org/web/packages/unmarked/vignettes/distsamp.pdf
#create sequence of cover levels
head(habConstant <- data.frame(Cover = seq(0, 1, length=25),habitat=factor("A", levels=c("A", "B"))))
#predict density with changing cover
Elambda <- predict(m.haz.1.cover, type="state", newdata=habConstant,appendData=TRUE)


#DENSITY DECREASING WITH COVER
pdf("Plant Cover.pdf", width=6, height=5)
with(Elambda, {x <- barplot(Predicted/1000000, names= Cover*100, xlab="",
                            ylab="", ylim=c(0, .38), cex.names=1,cex.lab=1.5, cex.axis=1)
box()
legend('topright', c("Density at Picnic", "Density at T-junction"),col=c("blue", "dark green"), lty=2, cex=1, lwd=2)
abline(h= .158853,col="blue",lwd=3,lty=2)
abline(h= .063303, col= "dark green", lwd=3,lty=2)
abline(h= 0.0418, col= "dark green", lwd=2,lty=3)
abline(h= 0.0959, col= "dark green", lwd=2,lty=3)
abline(h= 0.130,col="blue",lwd=2,lty=3)
abline(h= 0.197,col="blue",lwd=2,lty=3)
title(ylab="Density (plants / m^2)", line=2.4, cex.lab=1.2)
title(xlab="Plant Cover (%)", line=2.4, cex.lab=1.2)
})
dev.off()


#BINNING OF DATA FIGURE
#this is how it looks with 0.25 cm binning
jpeg("Binning.jpeg", quality=75, res=150)
ggplot(DS.trunc, aes(Dist, fill =  Site))+
  geom_histogram(bins=24)+
  labs(x= "Distance", y="Count")+
  theme_bw()+
  scale_fill_manual(values=c("light green", "dark green"))
dev.off()

jpeg("Binning.jpeg", quality=75, res=150)
ggplot(DS.trunc, aes(Dist, fill =  Site))+
  geom_histogram(bins=24)+
  facet_wrap(~Site)+
  theme_bw()+
  labs(x= "Distance", y="Count")+
  theme(legend.position = "none")+
  scale_fill_manual(values=c("blue", "dark blue"))
dev.off()

par(mfrow=c(1,1))
ggplot(DS.trunc, aes(Dist))+
  geom_histogram(bins=24)+
  theme_bw()+
  labs(x= "Distance", y="Count")+
  theme(legend.position = "none")+
  scale_fill_manual(values = "blue")

#confidence of densities
confint(lcJ, level=.95)
confint(lcP, level=.95)

#abundance conf ints
tarea<-6455.991
tareakm<- tarea/1000000

parea<- 7802.578
pareakm<- parea/1000000

#DESNITIES
#t-junction density is 63303 per sqkm 
ttotalpop<- 63303*tareakm #408.69 individuals
#95%
tupper<-95891.52*tareakm #619.07 individuals
tlower<-41790.22*tareakm #269.80

#picnic density is 159853 per sqkm 
ptotalpop<- 159853*pareakm #1247.27 individuals
#95%
pupper<-196702.6*pareakm #1534.79
plower<-129907*pareakm #1013.61

#ABUNDANCE check
Jcover<- data.frame(Cover=c(0.481))
Jsite.level.density <- predict(m.haz.1.cover, type="state", newdata=Jcover)$Predicted
JplotArea.inkm <- 0.006455991
Jsite.level.abundance <- Jsite.level.density * JplotArea.inkm
JN.hat <- sum(Jsite.level.abundance)

Pcover<- data.frame(Cover=c(0.176))
site.level.density <- predict(m.haz.1.cover, type="state", newdata=Pcover)$Predicted
plotArea.inkm <- 0.0078
site.level.abundance <- site.level.density * plotArea.inkm
N.hat <- sum(site.level.abundance)


#DETECTION PROB
jpeg("Detection Probability.jpeg", quality=100, res=150)
plot(function(x) gxhaz(x, shape=2.58, scale=1.16), 0, 6, xlab="Distance (m)",ylab="Detection Probability", cex.lab=0.7,
     cex.axis=0.7, las=1,lwd=2)
dev.off()

ggsave("Detection Probability.jpg",plot=dp, width=8,height=6,units=c('cm'),dpi= 300)


apply(DSdata[4:8],2,range)










######################################## SPLIT DATA INTO REPRODUCTIVE AND NON ##########################################

#could not do this for selected data because there were no reproductive individuals detected#

#reproductive n = 25
DSR<-DSdata[DSdata$FlFr>0,]
DSR.ordered <- DSR[with(DSR, order(DSR$Site, DSR$Date, DSR$Dist)),]

#nonrepro n = 224
DSN<-DSdata[DSdata$FlFr==0,]
DSN.ordered <- DSN[with(DSN, order(DSN$Site, DSN$Date, DSN$Dist)),]

DSR.ordered$trans <- as.factor(DSR.ordered$trans)
DSN.ordered$trans <- as.factor(DSN.ordered$trans)

#small sample size for reproductive plants

## Graph them ##

ggplot(DSR.ordered, aes(Dist, colour = Site))+
  geom_density()+
  theme_bw()

ggplot(DSN.ordered, aes(Dist, colour = Site))+
  geom_density()+
  theme_bw()

##TRUNCATION to 6m
#REPRO
DSR.trunc<- DSR[DSR$Dist<6,] #one outlier at T-junction with 6 meters away
ggplot(DSR.trunc, aes(Dist, fill =  Site))+
  geom_histogram()+
  facet_wrap(~Site)+
  theme_bw()

#NONREPRO
DSN.trunc<-DSN[DSN$Dist<6,]
ggplot(DSN.trunc, aes(Dist, fill =  Site))+
  geom_histogram()+
  facet_wrap(~Site)+
  theme_bw()

#### FORMAT THE DATA FOR UNMARKED #### #this is where i would change truncation or binning

yDatR <- formatDistData(DSR.trunc, distCol = "Dist", transectNameCol = "trans", dist.breaks =seq(0,6,0.25))
yDatN <- formatDistData(DSN.trunc, distCol = "Dist", transectNameCol = "trans", dist.breaks = seq(0,6,0.25))

nrow(yDatR)
#making sure everything is lining up #all good transect length
identical(row.names(table(DSR$trans, DSR$Length)), row.names(yDatR))
identical(row.names(table(DSN$trans, DSN$Length)), row.names(yDatN))

xR <- table(DSR$trans, DSR$Length)[1,]
xN <- table(DSN$trans, DSN$Length)[1,]

## site covariates: site, cover ##
#repro
siteCovsR<- data.frame(Site = unlist(lapply(row.names(yDatR), function(xR) strsplit(xR,"[.]")[[1]][1])), Cover = NA,length = 
                         c(apply(table(DSR$trans, DSR$Length), 1, function(xR) as.numeric(as.character(names(xR[xR>0]))))))
siteCovsR$Cover <- ifelse(siteCovsR$Site %in% "Picnic", 0.175633,0.480958)
rownames(siteCovsR)<-rownames(yDatR)

#nonrepro
siteCovsN<- data.frame(Site = unlist(lapply(row.names(yDatN), function(xN) strsplit(xN,"[.]")[[1]][1])), Cover = NA,length = 
                         c(apply(table(DSN$trans, DSN$Length), 1, function(xN) as.numeric(as.character(names(xN[xN>0]))))))
siteCovsN$Cover <- ifelse(siteCovsN$Site %in% "Picnic", 0.175633,0.480958)
rownames(siteCovsN)<-rownames(yDatN)

#Very difficult to see a curve here for Picnic
ggplot(DSR, aes(Dist, fill =  Site))+
  geom_histogram(bins=24)+
  facet_wrap(~Site)+
  theme_bw() 

ggplot(DSN, aes(Dist, fill =  Site))+
  geom_histogram(bins=24)+
  facet_wrap(~Site)+
  theme_bw() 

### UNMARKED DATA FRAMES ###
umfR <- unmarkedFrameDS(y = as.matrix(yDatR), 
                        siteCovs = siteCovsR,
                        survey = "line",
                        dist.breaks =seq(0,6,0.25),
                        tlength= siteCovsR$length, 
                        unitsIn = "m")
umfN <- unmarkedFrameDS(y = as.matrix(yDatN), 
                        siteCovs = siteCovsN,
                        survey = "line",
                        dist.breaks = seq(0,6,0.25),
                        tlength= siteCovsN$length, 
                        unitsIn = "m")

#### COMPETE MODELS ####
#fit to either half-normal, hazard or uniform

#REPRO

m.halfR <- distsamp(~1 ~1, umfR, keyfun="halfnorm", output="density", unitsOut="kmsq")

m.hazR <- distsamp(~1 ~1, umfR, keyfun="hazard", output="density", unitsOut="kmsq")

m.uniR <- distsamp(~1 ~1, umfR, keyfun="uniform", output="density", unitsOut="kmsq")

library(MuMIn)
model.sel(m.halfR,m.hazR,m.uniR) #HAZARD
#now lets try covariates

m.hazRsite <- distsamp(~Cover ~1, umfR, keyfun="hazard", output="density", unitsOut="kmsq")
model.sel(m.hazR,m.hazRsite)
#site does not improve model delta AIC 5.19 #on slope or intercept or both


#NONREPRO
m.halfN <- distsamp(~1 ~1, umfN, keyfun="halfnorm", output="density", unitsOut="kmsq")

m.hazN <- distsamp(~1 ~1, umfN, keyfun="hazard", output="density", unitsOut="kmsq")

m.uniN <- distsamp(~1 ~1, umfN, keyfun="uniform", output="density", unitsOut="kmsq")

library(MuMIn)
model.sel(m.halfN,m.hazN,m.uniN) #hazard

m.haz.1.siteN <- distsamp(~1 ~Cover, umfN, keyfun="hazard", output="density", unitsOut="kmsq")
m.haz.1.siteN2<-distsamp(~Cover ~Cover, umfN, keyfun="hazard", output="density", unitsOut="kmsq")
m.haz.1.siteN3<-distsamp(~Cover ~1, umfN, keyfun="hazard", output="density", unitsOut="kmsq")

model.sel(m.hazN,m.haz.1.siteN,m.haz.1.siteN2,m.haz.1.siteN3) #adding cover to just the intercept (density) is best

