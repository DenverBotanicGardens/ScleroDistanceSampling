library(devtools)
#Distance
library(terra)
library(raster)
library(unmarked)
library(dplyr)

#Plotting 
library(ggplot2)

setwd("~/Documents/Regis/R /Externship")
DSdata <- read.csv("DistSamp_ScGl_JuneJuly.csv")

sum(DSdata$Length[,DSdata$Site="T-Junction"])

this<-DSdata[DSdata$Site=="T-Junction",]
sum(this$Length)

summary(DSdata)

#data cleaning, two spellings for picnic
DSdata$Site[grep("icnic",DSdata$Site)] <- "Picnic"
DSdata$Site[grep("unction", DSdata$Site)] <- "T-Junction"

T<- DSdata[DSdata$Site=="Picnic",]

#add unique identifier
DSdata$trans <- paste(DSdata$Site,DSdata$Line, sep = ".")
# We're adding more same plot and same transect number so paste the date on too for a unique Transect identifier
DSdata$trans <- paste(DSdata$trans, DSdata$Date, sep = ".")

#### SPLIT DATA INTO REPRODUCTIVE AND NON ####

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


######### COMBINED DATA SET #########

DS <- read.csv("DistSamp_ScGl_JuneJuly.csv")
summary(DS)

DS$Site[grep("icnic",DS$Site)] <- "Picnic"
DS$Site[grep("unction", DS$Site)] <- "T-Junction"

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

##TRUNCATION

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

#MODELS

m.half <- distsamp(~1 ~1, umf, keyfun="halfnorm", output="density", unitsOut="kmsq")

m.haz <- distsamp(~1 ~1, umf, keyfun="hazard", output="density", unitsOut="kmsq")

m.uni <- distsamp(~1 ~1, umf, keyfun="uniform", output="density", unitsOut="kmsq")
#M.HALF AIC = 1003.13
m.half
#m.haz AIC = 942.5986
m.haz
#m.uni AIC = 1323.66
m.uni

model.sel(m.half,m.haz,m.uni)

m.haz.1.site <- distsamp(~1 ~Cover, umf, keyfun="hazard", output="density", unitsOut="kmsq",starts=c(0,0,0,0)) 
m.haz.1.site2 <- distsamp(~Cover ~1, umf, keyfun="hazard", output="density", unitsOut="kmsq",starts=c(0,0,0,0)) 
m.haz.1.site3 <- distsamp(~Cover ~Cover, umf, keyfun="hazard", output="density", unitsOut="kmsq",starts=c(0,0,0,0,0)) 
#used 0 here as a start value for cover, 0 behaves the same as -1 and 0.5, they all give the same answer.
m.haz.1.site #improves here
here<-model.sel(m.haz,m.haz.1.site, m.haz.1.site2,m.haz.1.site3)

m.haz.1.site <- distsamp(~1 ~1+Cover, umf, keyfun="hazard", output="density", unitsOut="kmsq",starts=c(0,0,0,0))

#BEST
m.hazR
m.haz.1.siteN
m.haz.1.site
backTransform(m.hazR, type= "state")

#TABLE?
library(knitr)
library(kableExtra)
library(magick)
table<- data.frame("Data"= c("All Individuals","    Picnic", "    T-Junction","Reproductive", "Vegetative","    Picnic", 
                                "T-Junction"), "Estimated Density (individuals/m^2)" = c(" " , 1.53, 0.393,0.287," ",1.45,
                                  0.349), "Confidence Interval" = c(" ","1.25-1.87","0.289-0.538", "0.139-0.593"," ","1.18-1.78",
                                  "0.239-0.507"))
colnames(table)[2]<- "Estimated Density (m^2)"
colnames(table)[3]<- "Confidence Interval"
kable(table)

#####RESULTS####
#DENSITY BACK TRANSFORM
#when cover is 0.17
lcP <- linearComb(m.haz.1.site, c(1, 0.176), type="state") 
lcP<-backTransform(lcP)
#when cover is 0.481
lcJ <- linearComb(m.haz.1.site, c(1, 0.481), type="state") 
lcJ<-backTransform(lcJ)

#THIS IS THE SCALE
backTransform(m.haz.1.site, type="det")
#THIS IS THE SHAPE UGH
backTransform(m.haz.1.site, type="scale")
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
Elambda <- predict(m.haz.1.site, type="state", newdata=habConstant,appendData=TRUE)


#DENSITY DECREASING WITH COVER
jpeg("Plant Cover.jpeg", quality=75, res=150)
with(Elambda, {x <- barplot(Predicted/1000000, names= Cover*100, xlab="Plant Cover (%)",
                            ylab=expression("Density (plants / m^2)"), ylim=c(0, .38), cex.names=.8,cex.lab=.9, cex.axis=0.8)
box()
legend('topright', c("Density at Picnic", "Density at T-junction"),col=c("pale green", "dark green"), lty=2, cex=1, lwd=2)
abline(h= .153003,col="light green",lwd=3,lty=2)
abline(h= .039393, col= "dark green", lwd=3,lty=2)
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
  scale_fill_manual(values=c("pale green", "dark green"))
dev.off()

#confidence of densities
confint(lcJ, level=.95)
confint(lcP, level=.95)

#abundance conf ints
tarea<-6455.991
tareakm<- tarea/1000000

parea<- 7802.578
pareakm<- parea/1000000

#total populations
#t-junction density is 39393 per sqkm 
ttotalpop<- 39393*tareakm #254.32 individuals
#95%
tupper<-53777.04*tareakm #347.18 individuals
tlower<-28857*tareakm #186.30

#picnic density is 153003 per sqkm 
ptotalpop<- 153003*pareakm #1193.82 individuals
pa<- 
#95%
pupper<-186587.7*pareakm #1455.86
plower<-125463.8*pareakm #978.94

#ABUNDANCE check
Jcover<- data.frame(Cover=c(0.481))
Jsite.level.density <- predict(m.haz.1.site, type="state", newdata=Jcover)$Predicted
JplotArea.inkm <- 0.006455991
Jsite.level.abundance <- Jsite.level.density * JplotArea.inkm
JN.hat <- sum(Jsite.level.abundance)

Pcover<- data.frame(Cover=c(0.176))
site.level.density <- predict(m.haz.1.site, type="state", newdata=Pcover)$Predicted
plotArea.inkm <- 0.0078
site.level.abundance <- site.level.density * plotArea.inkm
N.hat <- sum(site.level.abundance)

JgetN.hat <- function(fit) {
  d <- predict(m.haz.1.site, type="state", newdata=Jcover)$Predicted
  a <- d * 0.006455991
  N.hat <- c(N.hat = sum(a))
  return(N.hat)
}

#DETECTION PROB
jpeg("Detection Probability.jpeg", quality=100, res=150)
plot(function(x) gxhaz(x, shape=2.43, scale=0.981), 0, 6, xlab="Distance (m)",ylab="Detection Probability", cex.lab=0.7,
     cex.axis=0.7, las=1)
dev.off()
ggsave("Detection Probability.jpg",plot=dp, width=8,height=6,units=c('cm'),dpi= 300)


apply(DSdata[4:8],2,sd)

