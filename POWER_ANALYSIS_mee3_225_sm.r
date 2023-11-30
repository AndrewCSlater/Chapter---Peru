##################################################################
# Power analysis for occupancy studies under imperfect detection
# Examples on how to use the provided functions
# Authors: Gurutzeta Guillera-Arroita & Jose J. Lahoz-Monfort
##################################################################

load("models/1.B-Training_Min11sp_5PCA_50000_25000_50_1.Rdata")
out = out.sfMsPGOcc
rm(out.sfMsPGOcc); gc()
# source('PAfuncts.R') # load the required functions

# Example A: use function 'calcPowerFormula' to calculate the 
# power of a design using the formula (as in Figure 1)
S <- nrow(out$y[1,,])		# number of stations
K <- round(mean(apply(X = !is.na(out$y[1,,]), MARGIN = 1, sum)), digits = 0)   # ncol(out$y[,1,]) # 3		# number of replicates
K=27

R <- 0.5	# proportional change in occupancy We Want to Detect

bet = as.data.frame(out$beta.samples)  # Probabilities of OCCURRENCE
bet = as.data.frame(rowMeans(out$psi.samples, na.rm = T, dims = 2)) # Mean prob of occurrence across all stations per replicate
alp = as.data.frame(out$alpha.samples) # Probabilities of DETECTION

sites.needed = c()
sp.power = c()
for (jj in 1:nrow(out$y[,1,])) { # For each species
  #jj=1
  # psi1 <- plogis(mean(bet[,jj])) # 0.3    # species specific occupancy probability
  psi1 = mean(bet[,jj])
  p <- plogis(mean(alp[,jj]))    # 0.6		# species specific detection probability
  alpha <-0.05	# significance level (FALSE POSITIVE ie detect effect where none exists)
  (powerF <- calcPowerFormula(S1=S,S2=S,K1=K,K2=K,p1=p,p2=p,psi1,R,alpha))
  sp.power = append(sp.power, powerF)
  
  # Example B: use function 'calcSFormula' to calculate the 
  # number of sites needed to achieve a given power (thick lines in Figure 3)
  pow <- 0.7		# target power level i.e. Probability of detecting an effect of magnitude R
  (SS <- calcSFormula(K1=K,K2=K,p1=p,p2=p,psi1,R,alpha,pow))
  sites.needed = append(sites.needed, SS)
}

summary(sites.needed)[]
summary(sp.power)
sd(sp.power)
hist(sp.power)
length(which(sp.power>=0.7)) # 5 species with a 70% probability of detecting a 50% difference in occupancy between groups (with a mean of 3 replicate surveys)
# Example C: use function runPowerSims to assess the actual power 
# using simulations (Wald test on probability scale, red lines in Figure 3) 

nsims<-5000		# number of simulations to run
powerR <- runPowerSims(S1=SS,S2=SS,K1=K,K2=K,p1=p,p2=p,psi1,R,alpha,nsims)
 