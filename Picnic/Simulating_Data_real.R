library(spam)
library(DSsim)
library(shapefiles)
library(raster)
library(terra)
library(parallelly)
library(devtools)


setwd("~/.ssh/ScleroDistanceSampling/Picnic")

# Create a polgon
tjunc <- read.shapefile("Tjunction")
picnic <- read.shapefile("Picnic")

library(terra)
x <- vect('Tjunction.shp')
tarea<-expanse(x, unit="m") #6455 m^2
tareakm<- tarea/1000000

x2 <- vect('Picnic.shp')
parea<-expanse(x2) #7802 m^2
pareakm<- parea/1000000

#total populations
#t-junction density is 39393 per sqkm 
ttotalpop<- 39393*tareakm #254.32 individuals

#picnic density is 153003 per sqkm 
ptotalpop<- 153003*pareakm #1193.82 individuals

# Create the survey region
tregion <-make.region(region.name = "Survey Region", units = "m",shapefile = tjunc)
plot(tregion, plot.units = "km")
str(tregion)

pregion <-make.region(region.name = "Survey Region", units = "m",shapefile = picnic)
plot(pregion, plot.units = "km")


# Create the density surface
tdensity <- make.density(region.obj = tregion,
                         x.space = 1, 
                         y.space = 1, 
                         constant = 1)
pdensity <- make.density(region.obj = pregion, 
                         x.space = 1, 
                         y.space = 1, 
                         constant = 1)

# Plot the density surface (can be in 'm' or 'km')
# Note that if converting units the values of the density surface are also converted
# In this example they are plotted in animals per square km instead of per square m.
plot(tdensity, style = "blocks", plot.units = "km")
plot(tregion, add = TRUE)

plot(pdensity, style = "blocks", plot.units = "km")
plot(pregion, add = TRUE)

#DETECTABILITY FUNCTION
tpop.desc <- make.population.description(region.obj = tregion, 
                                        density.obj = tdensity,
                                        N= 254.32)

detect.hz <- make.detectability(key.function = "hr", scale.param =0.981, shape.param = 2.43, truncation = 6)

system.time(summaryoutt<-do.call(rbind,mclapply(6,function(x) {
  design <- make.design(transect.type = "line",    
                        design.details = c("parallel","systematic"),
                        region.obj = tregion,
                        spacing = x)
  ttransects <- generate.transects(design, region = tregion)
  sum(ttransects@sampler.info$length)
  plot(tregion, plot.units = "km")
  plot(ttransects, col = 4, lwd = 2)
  
  #analysis #NEXT TO LOOK AT
  ddf.analyses <- make.ddf.analysis.list(dsmodel = list(~cds(key = "hn", formula = ~1),
                                                        ~cds(key = "hr", formula = ~1)), 
                                         method = "ds",
                                         criteria = "AIC",
                                         truncation = 6)
  
  tsim <- make.simulation(reps = 2, 
                          region.obj = tregion,
                          design.obj = design,
                          population.description.obj = tpop.desc,
                          detectability.obj = detect.hz,
                          ddf.analyses.list = ddf.analyses)
  check.sim.setup(tsim)
  
  tsim <- run(tsim)
  sum<-summary(tsim, description.summary = TRUE)
  data.frame(spacing=x, "Mean Area Covered"= sum@individuals$summary$mean.Cover.Area, 
             "Mean Effort" = sum@individuals$summary$mean.Effort, "True Abundance" = sum@individuals$N$Truth,
             "Mean Abundance Estimate" = sum@individuals$N$mean.Estimate, "Abundance Bias" = sum@individuals$N$percent.bias,
             "Mean Abundance SE" = sum@individuals$N$mean.se, "Mean Abundance SD" = sum@individuals$N$sd.of.means, "True Density"=
              sum@individuals$D$Truth, "Density Estimate" = sum@individuals$D$mean.Estimate, "Density Bias" =
              sum@individuals$D$percent.bias, "Density SE" = sum@individuals$D$mean.se, "Density SD" =sum@individuals$D$sd.of.means)
})))


#RUN THIS, done
#write.csv(summaryoutt,"~/Documents/Regis/R /Externship/Picnic/Tjuncsummary2.csv")


#PICNIC#
ppop.desc <- make.population.description(region.obj = pregion, 
                                        density.obj = pdensity, 
                                        N =1193.82)

detect.hz <- make.detectability(key.function = "hr", scale.param =0.981, shape.param = 2.43, truncation = 6)

summaryoutp<-do.call(rbind,lapply(seq(6,8, by=2), function(x){
  design <- make.design(transect.type = "line",    
                        design.details = c("parallel","systematic"),
                        region.obj = pregion,
                        spacing = x)
  ptransects <- generate.transects(design, region = pregion)
  plot(pregion, plot.units = "km")
  plot(ptransects, col = 4, lwd = 2)
  
  #analysis #NEXT TO LOOK AT
  ddf.analyses <- make.ddf.analysis.list(dsmodel = list(~cds(key = "hn", formula = ~1),
                                                        ~cds(key = "hr", formula = ~1)), 
                                         method = "ds",
                                         criteria = "AIC",
                                         truncation = 6)
  
  psim <- make.simulation(reps = 10, 
                          region.obj = pregion,
                          design.obj = design,
                          population.description.obj = ppop.desc,
                          detectability.obj = detect.hz,
                          ddf.analyses.list = ddf.analyses)
  check.sim.setup(psim)
  
  psim <- run(psim)
  sum<-summary(psim, description.summary = F)
  data.frame(spacing=x, "Mean Area Covered"= sum@individuals$summary$mean.Cover.Area, 
             "Mean Effort" = sum@individuals$summary$mean.Effort, "True Abundance" = sum@individuals$N$Truth,
             "Mean Abundance Estimate" = sum@individuals$N$mean.Estimate, "Abundance Bias" = sum@individuals$N$percent.bias,
             "Mean Abundance SE" = sum@individuals$N$mean.se, "Mean Abundance SD" = sum@individuals$N$sd.of.means, "True Density"=
               sum@individuals$D$Truth, "Density Estimate" = sum@individuals$D$mean.Estimate, "Density Bias" =
               sum@individuals$D$percent.bias, "Density SE" = sum@individuals$D$mean.se, "Density SD" =sum@individuals$D$sd.of.means)
}))

#write.csv(summaryoutp,"~/Documents/Regis/R /Externship/Picnic/Picnicsummary2.csv")
#ran this simulation on different computer, ignore above


krent<- data.frame("habitat area"=c(40542,17450,9069.5,10429.7,51420.3,14298.1,24078.8,14951.9,3038.4,17827.9,17840.2,
                          15122.9,11942.2,4287.9,8408.7,15144.7), "macroplot area"=c(990,2000,560,392,720,960,1400,320,720,
                                                                                     500,420,480,800,216,900,600), "quadrat area"
                          =c(55,100,28,28,40,48,100,32,60,50,42,40,100,18,75,60), "N"= c(18,20,20,14,18,20,14,10,12,10,10,12,
                                                                                         8,12,12,10), "n" = c(8,10,9,6,8,11,7,
               
                                                                                                                                                                                                             5,6,7,5,6,5,8,7,5))
propofmacro<-mean(krent$macroplot.area/krent$habitat.area)
propofquadarea<- krent$quadrat.area/krent$habitat.area

sum(krent$macroplot.area)/ sum(krent$habitat.area)

sum(krent$habitat.area)



#DENSITY OF PICNIC = 0.160 individuals/m2 with a 95% confidence interval of 0.130 - 0.197
16996891*0.160
#lower
16996891*0.130
#upper
16996891*0.197

#his calculation
16996891*.0061

krent2<-data.frame(" "=c("Krening", "Picnic", "Lower Confidence Interval for Picnic", "Upper Confidence Interval for Picnic"), 
                   "Estimated Density"=c(0.0061,0.0633,0.0418,0.0959), "Estimated Abundance"
                   = c(103086,2719503,2209596,3348388))

####################EFFORT##########################

Tjunc<- read.csv("Tjuncsummary2.csv")
Picnic<- read.csv("psummary2.csv")
Tjunc$Site<- "T-Junction"
Picnic$Site<- "Picnic"
both<-rbind(Tjunc,Picnic)

###FIX
library(car)
plot(x=c(Tjunc$Mean.Area.Covered,Picnic$Mean.Area.Covered), y=c(Tjunc$Mean.Abundance.SD,Picnic$Mean.Abundance.SD),
     pch=19,ylim=c(30,200),ylab="SD of Abundance Estimate", 
     xlab=expression("Mean Area Covered (m^2)"), col=rep(c("dark green","pale green"),each=13))
legend('topright', c("Picnic", "T-Junction"),col=c("pale green", "dark green"), lty=1, cex=1, lwd=2)
text(y=Tjunc$Mean.Abundance.SD, x=Tjunc$Mean.Area.Covered, Tjunc$spacing, pos=3, col="dark green")
text(y=Picnic$Mean.Abundance.SD, x=Picnic$Mean.Area.Covered, Picnic$spacing, pos=3, col="green")
abline(lm(Tjunc$Mean.Abundance.SD~Tjunc$Mean.Area.Covered))


#GG
ggplot(both, aes(y=(Mean.Abundance.SD/Mean.Abundance.Estimate), x=Mean.Area.Covered, color=Site)) +
  geom_point()+
  scale_color_manual(values=c("blue","dark green"))+
  ylim(0, 0.5)+
  xlim(1800,17000)+
  geom_text(aes(label=spacing),hjust=-.5, vjust=-.5, show.legend = F)+
  ylab("CV of Abundance Estimate")+
  xlab(expression("Area Covered by Sampling Effort " (~m^2)))+
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE)


