library(devtools)
library(terra)
library(raster)
library(unmarked)
library(dplyr)
library(ggplot2)
library(plyr)

#read in the data
setwd("~/.ssh/ScleroDistanceSampling")
setwd("C:/Users/deprengm/ScleroDistanceSampling")
DSdata <- read.csv("DistSamp_ScGl_JuneJuly.csv")

## How big are our T-Junction and Picnic sites?
library(dssd)
picnic.region <- make.region(region.name = "Picnic", shape = "C:/Users/deprengm/OneDrive - Denver Botanic Gardens/P drive/hackathon/ScleroDistanceSampling/Picnic/Tjunction.shp",
                            units = "m")
plot(picnic.region)
## 6459 m^2??

distsamp <- read.csv("C:/Users/deprengm/OneDrive - Denver Botanic Gardens/P drive/hackathon/ScleroDistanceSampling/Data/distanceSampling_1.csv")
head(distsamp)

#get only data from certain dates
DSdata$Date<-as.character(DSdata$Date)
DSdata$Site<- as.character(DSdata$Site)
#data cleaning, two spellings for picnic
DSdata$Site[grep("icnic",DSdata$Site)] <- "Picnic"
DSdata$Site[grep("unction", DSdata$Site)] <- "T-Junction"
dates<- c("8/4/21","8/5/21","8/6/21", "8/1/21")
DSdata<- DSdata %>%
  filter(DSdata$Date %in% dates)

# See only larger individuals at farther distances
ggplot(DSdata, aes(Dist, Width, color = Site ))+
  geom_point()


#T-Junction n=25
nrow(DSdata[DSdata$Site=="T-Junction",])
#Picnic n=159
nrow(DSdata[DSdata$Site=="Picnic",])

## Total length of 180 meters by w = 6 for Picnic, total length of 120 half width is 6 for T-junction
DSdata %>%
  distinct(Site, Line, .keep_all = TRUE) %>%
  group_by(Site) %>%
  dplyr::summarise(TotalLength = sum(Length)*6*2)
  

#add unique identifier
DSdata$trans <- paste(DSdata$Site,DSdata$Line, sep = ".")
# We're adding more same plot and same transect number so paste the date on too for a unique Transect identifier
DSdata$trans <- paste(DSdata$trans, DSdata$Date, sep = ".")


t.test(Height~Site, DSdata)
t.test(Width~Site, DSdata)


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

## TRUNCATION - half-width is 6 meters

#COMBINED
DS.trunc <- DS[DS$Dist < 6,]
#this is how it looks with 0.25 cm binning
par(1,2)
ggplot(DS.trunc, aes(Dist, fill =  Site))+
  geom_histogram(bins=24)+
  facet_wrap(~Site)+
  theme_bw()

dev.off()

scatter.smooth(DS.trunc$Height[DS.trunc$Site == "T-Junction"], family = "gaussian", pch= 20, cex=.9, lpars=list(lwd=3),
               xlab="Plant Height",ylab="Distance (m)")

scatter.smooth(DS.trunc$Height[DS.trunc$Site == "Picnic"], family = "gaussian", pch= 20, cex=.9, lpars=list(lwd=3),
               xlab="Plant Height",ylab="Distance (m)")




#format the data for umf, set correct truncation and binning
yDat <- formatDistData(DS.trunc, distCol = "Dist", transectNameCol = "trans", dist.breaks = seq(0,6,0.25))
yDat.Pic <- formatDistData(DS.trunc[DS.trunc$Site == "Picnic",], distCol = "Dist", transectNameCol = "trans", dist.breaks = seq(0,6,0.25))

#COVARIATES aka cover; percent plant cover from LPI
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


umfPic <- unmarkedFrameDS(y = as.matrix(yDat.Pic), 
                       siteCovs = siteCovs[siteCovs$Site == "Picnic",],
                       survey = "line",
                       dist.breaks = seq(0,6,0.25),
                       tlength= siteCovs$length[siteCovs$Site == "Picnic"], 
                       unitsIn = "m")

###################################################### MODELS #########################################################

m.half <- distsamp(~1 ~1, umf, keyfun="halfnorm", output="density", unitsOut="kmsq")
predict(m.half, type = "state")[1,]/1000000
m.halfPic <- distsamp(~1 ~1, umfPic, keyfun="halfnorm", output="density", unitsOut="kmsq")
predict(m.halfPic, type = "state")[1,]

m.haz <- distsamp(~1 ~1, umf, keyfun="hazard", output="density", unitsOut="kmsq")
predict(m.haz, type = "state")[1,]/1000000
m.hazPic <- distsamp(~1 ~1, umfPic, keyfun="hazard", output="density", unitsOut="kmsq")
predict(m.hazPic, type = "state")[1,]/1000000

m.uni <- distsamp(~1 ~1, umf, keyfun="uniform", output="density", unitsOut="kmsq")
#M.HALF AIC =
m.half # AIC: 630.954 
#m.haz AIC = 
m.haz  # AIC: 593.3394 
#m.uni AIC = 
m.uni  # AIC: 842.1926
# library(MuMIn)
# model.sel(m.half,m.haz,m.uni)
AIC(m.uni, m.half, m.haz)
Distance::summarize_ds_models(m.uni, m.half, m.haz)

#USING COVER AS A COVARIATE; detection and then abundance
m.haz.1.cover <- distsamp(~1 ~Cover, umf, keyfun="hazard", output="density", unitsOut="kmsq")      # AIC: 572.3241
m.haz.1.cover2 <- distsamp(~Cover ~1, umf, keyfun="hazard", output="density", unitsOut="kmsq")     # AIC: 583.327 
m.haz.1.cover3 <- distsamp(~Cover ~Cover, umf, keyfun="hazard", output="density", unitsOut="kmsq") # AIC: 574.307 
#used 0 here as a start value for cover, 0 behaves the same as -1 and 0.5, they all give the same answer.




# here<- data.frame(model.sel(m.haz,m.haz.1.cover,m.haz.1.cover2,m.haz.1.cover3, m.uni,m.half))
# here$delta

#BEST
m.haz.1.cover
exp(-4.418716)
exp(-1.655486)

#USING SITE AS A COVARIATE
m.haz.1.site <- distsamp(~1 ~Site, umf, keyfun="hazard", output="density", unitsOut="kmsq") 
m.haz.1.site2 <- distsamp(~Site ~1, umf, keyfun="hazard", output="density", unitsOut="kmsq") 
m.haz.1.site3 <- distsamp(~Site ~Site, umf, keyfun="hazard", output="density", unitsOut="kmsq") 
m.haz.1.site #improves here
here<- data.frame(model.sel(m.haz,m.haz.1.cover,m.haz.1.cover2,m.haz.1.cover3, m.uni,m.half))
here$delta

library(MuMIn)
model.sel(m.haz,m.haz.1.site,m.haz.1.site2,m.haz.1.site3)
model.sel(m.haz,m.uni,m.half)

exp(((confint(m.half, type="state"))))


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

##SITE##
exp(11.983) #PICNIC = 160011.3
exp(11.983+(-0.928))
#THIS IS THE SCALE
backTransform(m.haz.1.site, type="det")
#THIS IS THE SHAPE
backTransform(m.haz.1.site, type="scale")

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
pdf("Plant Cover.pdf", width=6, height=4)
with(Elambda, {x <- barplot(Predicted/1000000, names= Cover*100, xlab="",
                            ylab="", ylim=c(0, .38), cex.names=1,cex.lab=1.5, cex.axis=1)
box()
legend('topright', c("Density at Picnic", "Density at T-Junction"),col=c("blue", "dark green"), lty=2, cex=1, lwd=2)
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
pdf("detection.pdf", width=4, height=5)
plot(function(x) gxhaz(x, shape=2.58, scale=1.16), 0, 6, xlab="Distance (m)",ylab="Detection Probability", cex.lab=1,
     cex.axis=1.2, las=1,lwd=2)
dev.off()

ggsave("Detection Probability.jpg",plot=dp, width=8,height=6,units=c('cm'),dpi= 300)


apply(DSdata[4:8],2,range)

pdf("counts.pdf", width=5, height=5)
ggplot(DS.trunc, aes(Dist, fill =  Site))+
  geom_histogram(bins=24)+
  facet_wrap(~Site)+
  theme_bw()+
  scale_fill_manual(values = c("blue","dark green"))+
  xlab("Distance")+
  ylab("Number of Detected Individuals")+
  theme(legend.position="none")
dev.off()


ptotalpop<- 159853*pareakm #1247.27 individuals
##ABUNDANCES OF OTHER MODEL
#m.haz.1.site2
exp(11.8)
133252.4*pareakm+133252.4*tareakm
exp(confint(m.haz.1.site2, type="state"))
((111346.8*pareakm)+(111346.8*tareakm))
((168640.7*pareakm)+(168640.7*tareakm))
#m.haz.1.site3
m.haz.1.site3
(exp(11.978))*pareakm
exp(11.978-0.0249)*tareakm

confint(m.haz.1.site3, type="state")
exp(11.760246)*pareakm
exp(12.1967063)*pareakm
exp(11.760246-1.454224 )*tareakm
exp(12.1967063-0.3532933)*tareakm

#m.haz
m.haz
(exp(11.8)*pareakm)+(exp(11.8)*tareakm)
exp(confint(m.haz, type="state"))
((108446.8*pareakm)+(108446.8*tareakm))
((161595*pareakm)+(161595*tareakm))

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
# https://cran.r-project.org/web/packages/unmarked/vignettes/distsamp.html 
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


## Indivdiual covariates height and width have little impact. 
library(Distance)
conversion.factor <- convert_units("Metre", NULL, "square metre")

DS.Picnic <- DS.trunc %>%
  filter(Site == "Picnic") %>%
  dplyr::rename("distance" = "Dist")

scglPic.unif <- ds(DS.Picnic, transect = "line", key = "unif", 
                 adjustment = NULL, convert_units = conversion.factor)

scglPic.hn <- ds(DS.Picnic, transect = "line", key = "hn", 
                   adjustment = NULL, convert_units = conversion.factor)

summary(scglPic.hn)

scglPic.hr <- ds(DS.Picnic, transect = "line", key = "hr", 
                 adjustment = NULL, convert_units = conversion.factor)

scglPic.height <- ds(DS.Picnic, transect = "line", key = "hr", formula = ~ Height,
                     adjustment = NULL, convert_units = conversion.factor)

scglPic.width <- ds(DS.Picnic, transect = "line", key = "hr", formula = ~ Width,
                    adjustment = NULL, convert_units = conversion.factor)


scglPic.widthHeight <- ds(DS.Picnic, transect = "line", key = "hr", formula = ~ Width + Height,
                    adjustment = NULL, convert_units = conversion.factor)

AIC(scglPic.unif,scglPic.hr, scglPic.hn,scglPic.height, scglPic.width, scglPic.widthHeight)

plot(scglPic.height, pdf = TRUE, main = "Hazard rate with height of plant", showpoints=FALSE)

plot(scglPic.width, pdf = TRUE, main = "Hazard rate with width of plant", showpoints=FALSE)
plot(scglPic.unif, pdf = TRUE, main = "Uniform", showpoints=FALSE)

plot(scglPic.hr, pdf = TRUE, main = "Hazard rate", showpoints=FALSE)

Pic.distsamp.hn <- distsamp(~1 ~1, DS.Picnic, keyfun = "halfnorm",
                            output = "abund",
                            unitsOut = "kmsq")

##################################################################################################

### TEST BAYES see if same results  
library(rjags)
library(R2jags)

### When want to use data augmentation 
# Data augmentation: add a bunch of "pseudo-individuals"
# nz <- 500 # Augment by 500
# nind <- nrow(data)
# y <- c(data[,2], rep(0, nz)) # Augmented detection indicator y
# site <- c(data[,1], rep(NA, nz)) # Augmented site indicator,
# # unknown (i.e., NA) for augmented inds.
# d <- c(data[,5], rep(NA,nz)) # Augmented distance data (with NAs)


## Binned distance sampling to compare to ds  
## Evaluate the detection function at the midpoint of each interval 
# log(p[g]) <- -midpt[g]*midpt[g]/(2*sigma*sigma) ## half-normal function
## The probability mass for each distance interval, given the proportion of the bin
# pi[g] <- delta/B

x <- DS.Picnic$distance                    # distance data
nind <- nrow(DS.Picnic)
nz <- 200                                  # Augment observer data with zeros
y <- c(rep(1, nind), rep(0, nz))           # non-detections added
x <- c(x, rep(NA, nz))                     # distances are missing for augmented

B <- 6                                     # Strip half-width; truncation distance.
delta <- 0.25                              # width of the distance bins.
xg <- seq(0,B, delta)                      # interval cut points
dclass <- x %/% delta + 1                  # vector of distance classes; convert distances to cat. distances
midpt <- seq(delta/2, B, delta)            # make mid-points and chop up data
# midpt <- xg[-1] - delta/2        # SAME
nD <- length(xg) - 1                       # number of distance bins.

jags.data <- list(nind = nind, dclass = dclass, midpt = midpt, delta = delta, B = B,
                  nz = nz, y = y,
                  nD = nD)

modelPicnic <- 
  paste("
        model {
          
          # Priors
          psi ~ dunif(0,1)
          
          ## Hazard Rate
          # sigma ~ dunif(0,100)
          # tau ~ dunif(0,10)
          
          ## Half-normal
          sigma ~ dunif(0,1000)
          
          # Conditional detection and Pr(x) for each bin
          for(g in 1:nD){             # just each mid point
            log(p[g]) <- -midpt[g]/(2*sigma*sigma)  # half-normal
            # log(p[g]) <- 1 - exp(-pow(midpt[g]/sigma, -tau) )  # hazard rate??
            pi[g] <- delta/B                        # probability of x in each interval
          }
          
          for(i in 1:(nind+nz)){
            z[i] ~ dbern(psi)          # model for individual covariates
            dclass[i] ~ dcat(pi[])     # population distribution of distance class
            mu[i] <- z[i]*p[dclass[i]] # p depends on distance class
            y[i] ~ dbern(mu[i])
          }
          
          # Derived Population size and density
          N <- sum(z[]) 
          D <- N/1560  ### 2160      ## 180 meters, half width is 6 meters
          Ntot <- D * 6459           # the square meter density by the total m^2 of Picnic
  }") 

writeLines(modelPicnic, "Picnic_hn.jags")

# Inits function
zst <- y # DA variables start at observed value of y
inits <- function(){ list (psi=runif(1), z=zst, sigma=runif(1,40,200)) }
# Parameters to save
params <- c("N", "Ntot", "sigma", "D", "p")

ni <- 11000
nt <- 2
nb <- 1000
nc <- 3
## Call JAGS from R
out <- jags(jags.data, inits, params, "Picnic_hn.jags", 
             n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, 
             working.directory = getwd())

print(out)   
summary(scglPic.hn)

