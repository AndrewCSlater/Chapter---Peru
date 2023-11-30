# Occupancy Modelling using JAGS
# Tutorial from https://bcss.org.my/tut/

library(jagsUI)
library(wiqid)  # for getGammaPar; also loads mcmcOutput


####################################################
# ----- A SIMPLE MODEL WITH 1 PARAMETER -----------#

# Models need to be written and saved to a text file

write("
model {
  
  # Likelihood:
  lambda <- N * a       # Expected no. of ou in area a
  C ~ dpois(lambda)
  # Prior:
  N ~ dgamma(shape, rate)  # No. of ou in whole area
}
", "orangutan.jags")
######
# JAGS syntax: The code above looks a lot like R code but there are several differences. JAGS does not work through the lines of code in order, they can be in any order; often we prefer to put the likelihood first, then define the necessary priors. There are two kinds of relations in JAGS: the <- symbol defines a deterministic relation, where the left hand side is exactly determined by the values on the right; the ~ symbol is for stochastic relations, and indicates the distribution a parameter is drawn from


# Prepare the data:
# -----------------
priorPars <- getGammaPar(450, 37)

jagsData <- list(shape=priorPars[1], rate=priorPars[2],
                 a = 1/8, C = 45)

# Run the model in JAGS
# ------------------------
jagsout <- jags(jagsData, inits=NULL,
                parameters.to.save="N",
                model.file="orangutan.jags",
                n.chains=3, n.iter=20000, n.burnin=0, DIC=FALSE)
jagsout

diagPlot(jagsout) # Check all went well

# Convert to an mcmcOutput object and plot
ouOut <- mcmcOutput(jagsout, default="N")
dim(ouOut)
# [1] 60000     1
head(ouOut)
#      N
# 1 442.1731
# 2 428.6705
# 3 408.7989
# 4 450.6569
# 5 469.5928
# 6 469.7083
plot(ouOut)

######
# JAGS produces a “mountain of numbers” to describe the posterior distribution of N. With a large enough mountain of draws we can get very close to an accurate representation of the posterior, and can calculate all the summaries we need. A big advantage of this representation is the ease of calculating derived values.






























