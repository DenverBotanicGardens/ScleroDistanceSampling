install.packages("DSsim")
library

# Create a polgon
poly1 <- data.frame(x = c(0,0,20000,20000,0), y = c(0,5000,5000,0,0))

# Create an empty list
coords <- list()
# Store the polygon inside a list in the first element of the coords list referring to strata 1.
coords[[1]] <- list(poly1)

# Create the survey region
region <- make.region(region.name = "study area", 
                      units = "m",
                      coords = coords)
# The plot function allows plotting in km or m.
plot(region, plot.units = "km")