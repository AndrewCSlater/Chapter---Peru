##################################################################
# Power analysis for occupancy studies under imperfect detection
# Examples on how to use the provided functions
# Authors: Gurutzeta Guillera-Arroita & Jose J. Lahoz-Monfort
##################################################################

source('PAfuncts.R') # load the required functions

# Example A: use function 'calcPowerFormula' to calculate the 
# power of a design using the formula (as in Figure 1)
S <- nrow(out$y[1,,]) # 50		# number of sites
K <- 3 # ncol(out$y[,1,]) # 3		# number of replicates

sites.needed = c()
for (jj in 1:nrow(out$y[,1,])) {
#jj=1
R <- .5 # (max(plogis(bet[,jj]))-min(plogis(bet[,jj])))/max(plogis(bet[,jj]))  # 0.5	# proportional change in occupancy
psi1 <- plogis(mean(bet[,jj])) # 0.3       # initial occupancy probability
p <- plogis(mean(alp[,jj]))    # 0.6		# detection probability
alpha <-0.05	# significance level
(powerF <- calcPowerFormula(S1=S,S2=S,K1=K,K2=K,p1=p,p2=p,psi1,R,alpha))

# Example B: use function 'calcSFormula' to calculate the 
# number of sites needed to achieve a given power (thick lines in Figure 3)
pow <- 0.8		# target power level
(SS <- calcSFormula(K1=K,K2=K,p1=p,p2=p,psi1,R,alpha,pow))
sites.needed = append(sites.needed, SS)
}

hist(sites.needed, breaks = 20)
summary(sites.needed)

# Example C: use function runPowerSims to assess the actual power 
# using simulations (Wald test on probability scale, red lines in Figure 3) 

nsims<-5000		# number of simulations to run
powerR <- runPowerSims(S1=SS,S2=SS,K1=K,K2=K,p1=p,p2=p,psi1,R,alpha,nsims)
 