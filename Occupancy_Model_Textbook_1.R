library(jagsUI)
library(R2WinBUGS)
library(unmarked)

# Sample Simulated Occupancy Model - using unmarked

set.seed(24)
M = 100 # Number of Sites
J = 2   # Number of Pres/Abs measurements
y = matrix(NA, nrow = M, ncol = J)  # To log Observation Data

psi = 0.8 # Probability of Presence
p = 0.5   # Probability of Detection

# Generate True Presence / Absence
z = rbinom(n = M, size = 1, prob = psi) # Random Pres/Ab (1 or 0) with probability of 1 determined by psi

# Generate Detection / NON Detection
for (j in 1:J) {
  y[,j] = rbinom(n = M, size = 1, prob = z*p) # For each count J, what is the detection - Which equals the actual Presence/Absence multiplied by the probability of seeing it if it is Present - the 2 columns represent the results of the 2 counts
}


sum(z) # True number of Occupied Sites (out of 100) - True Probability = 0.80, but realised = 0.86

sum(apply(X = y, MARGIN = 1, FUN = max)) # The sum of the maximum number (ie 1 or 0) per row (Margin =1) of the matrix y - returns how many sites returned at least 1 sighting = 61
# The measurement error is (86-61)/86 = -29%

# We expect a Combined Detection Probability over J number of Surveys 1 - (1-p)^J
# Where (1-p)^J = Probability of Failure (1 - probability of success) to the power of Number of Counts
# Then 1 - above = 1 - Probability of Failure = Probability of Success
# Which in this case is 1 - (1-0.5)^2 = 0.25%

# So our realised Failure was 29% and our expected was 25% - Difference due to stochastic sampling error

#######
# Now Analyse Data Using UNMARKED Package
library(unmarked)
umf = unmarkedFrameOccu(y = y) # Creates an UNMARKED df
summary(umf)                   # Summarises data frame
fm1 = occu(formula = ~1 ~1, data = umf) # Fits Model = Double Right Hand Side Formula for Covariates of Detection and then Occupancy

backTransform(obj = fm1, 'state') # Estimate of True Occurrence (74% despite actual observation of 61%)
backTransform(obj = fm1, 'det')   # Estimate of Detection Rate (58% which is in between the realised of 61% and programmed truth of 50%)

################################
# More Complex Model Occupancy with Covariates
# Occupancy affected by Vegetation Height
# Detection affected by Wind Speed
rm(list=ls())
set.seed(1)
M = 100    # Sites
J = 3      # Observations
y = matrix(NA, nrow = M, ncol = J)


vegHt = sort(runif(n = M, min = -1, max = 1)) # M qty random veg heights, sorted for graphical conveniences

# Parameter Values for Occupancy Model
beta0 = 0 # Intercept Value
beta1 = 3 # Slope of Impact of vegHeight on Occupancy
psi = plogis(beta0 + beta1*vegHt) # Occupancy Probability (True) - Logistic Distribution (Binary outcome and probability between 0 & 1)
plot(vegHt, psi, ylim = c(0,1), type = "l", lwd = 3)

# Create 'True' State of Presence / Absence
z = rbinom(n = M, size = 1, prob = psi) # M number of Random 1 or 0 based on probability psi 

table(z) # 40 x 0, 51 x 1

# Plot the True State
par(mfrow = c(1,3), mar = c(5,5,2,2), cex.axis=1.5, cex.lab=1.5)
plot(vegHt, z, xlab='Veg Height', ylab='True Presence', frame=F, cex=1.5)
plot(function(x)plogis(beta0+beta1*x), -1, 1, add=T, lwd=3, col='red')


# Simulate impact of Wind on Probability of Observation
# Intercept = -2 and slope -3

wind = array(runif(M*J, -1, 1), dim=c(M,J)) # Creates random matrix of 300 numbers between -1 and 1, split into M rows & J columns
# Simulates windspeed at each of the survey sessions

# Choose Observation Parameter Values
alpha0 = -2 # Observation Inercept
alpha1 = -3 # Observation Slope
p = plogis(alpha0 + alpha1*wind) # Detection probability per observation
plot(p~wind, ylim = c(0,1)) # Plots how wind impacts probability of detection

# For each of J survey sessions, get M number of random numbers to represent whether a sighting was made or not.
# Random number based on Probability that there was something to see(z)
# and the (p) probility it was seen if it was there
for (j in 1:J) {
  y[,j] = rbinom(n = M, size = z, prob = p[,j])
}

sum(apply(X = y, MARGIN = 1, FUN = max)) # Apply across every row(margin=1) of matrix y, the maximum number ie. is there a 1...and sum the number of rows with a 1 - this is the actual number of plots where sightings were made = 32

plot(wind, y, xlab='Wind', ylab = 'Detected', frame=F, cex = 1.5)
tab = cbind(psiProbObs=round(psi,2),zTrueProb=z, yObserved1=y[,1],y2=y[,2],y3=y[,3])
tab

#######
# Fit Site Occupancy Models in UNMARKED

# Create 2 variables that are not involved in the probability of occupancy or observation to add to the models

time = matrix(rep(as.character(1:J), M), ncol = J, byrow = T)
hab = c(rep("A",33), rep("B",33), rep("C", 34))

library(unmarked)
# Set up data in correct format
umf = unmarkedFrameOccu( 
  y = y,                 # Observation Data
  siteCovs = data.frame(vegHt = vegHt, hab = hab), # Site specific covariates
  obsCovs = list(wind = wind, time = time)) # observation specific covariates

summary(umf)

# Fit Model & Extract Estimates
fm1.occ = occu(~wind ~vegHt, data = umf) # Detection covariates then occupancy
summary(fm1.occ)


##
# Predict Occupancy and Detection as function of Covariates (with 95% CI)
newdat = data.frame(vegHt = seq(-1, 1, 0.01)) # Create sequence of numbers from -1 to 1 in intervals of 0.01 as sample covariate values to make predictions from
pred.occ = predict(fm1.occ, type = 'state', newdata = newdat, append=T) # predict the "State" values (True occupancy) of the fm1.occ models using input from the newdat data frame

newdat = data.frame(wind = seq(-1, 1, 0.01))
pred.det = predict(fm1.occ, type = 'det', newdata = newdat, append=T) # predict the "Det" values (Detections) of the fm1.occ models using input from the newdat data frame



# Predict using Specific Values of vegHt and Wind
newdat = data.frame(vegHt = c(0.2, 2.1))
predict(fm1.occ, type = "state", newdata = newdat, append=T)


newdat = data.frame(wind = seq(-1, 1, 0.5))
predict(fm1.occ, type = "det", newdata = newdat, append = T)


# Fit GLM observe occurrence that does not take imperfect detection into account
par(mfrow = c(1,1))
(fm.glm = glm(apply(X = y, MARGIN = 1, FUN = max)~vegHt, family = binomial))

# Plot dots of actual Presence/Absence
plot(vegHt, apply(X = y, MARGIN = 1, FUN = max), xlab = 'Veg Height', ylab = 'Observed Occurrence? (3 Observations)', frame =F, cex=1.5)
# Plot Line Linking True Presence with Veg Height
plot(function(x) plogis(beta0+beta1*x), -1, 1, add=T, lwd=3, col = 'red')

# Plot Line linking Predicted Presence with Veg Height
lines(vegHt, predict(fm1.occ, type = 'state')[,1], col='blue', lwd=3)

# Plot Line simple logistic regression with no imperfect detection (p probability)
lines(vegHt, predict(fm.glm,,'response'), type = "l", lwd=3)

ranef(fm1.occ)
# The random effects provides the best guess as to whether a site is occupied of not.
# These are CONDITIONAL OCCUPANCY PROBABILITIES
# And are conditional on the Observed Data

# Check for Site 1 where species was never detected
(psi1 = predict(fm1.occ, type = 'state')[1,1])
# The Probability of Occupancy (Prob of True STATE) = 0.076

(p1 = predict(fm1.occ, type = 'det')[c(1:3),1])
# The 3 probabilities of DETECTION - 1 for each visit
# 0.433, 0.059, 0.065



