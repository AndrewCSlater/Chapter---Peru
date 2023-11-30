##################################################################
# Power analysis for occupancy studies under imperfect detection
# Examples on how to use the provided functions
# Authors: Gurutzeta Guillera-Arroita & Jose J. Lahoz-Monfort
##################################################################

# source('PAfuncts.R') # load the required functions

# Example A: use function 'calcPowerFormula' to calculate the 
# power of a design using the formula (as in Figure 1)

### POWER = Probability of Detecting a Change on Occupancy (The proportion of change required is given by value R below)

alp = as.data.frame(out$alpha.samples)  # Probability of Detection matrix
bet = as.data.frame(out$beta.samples)   # Probability of Occurrence matrix

S <- nrow(out$y[1,,])                                                        # 50		# number of sites
K <- round((length(which(!is.na(out$y[1,,]))))/nrow(out$y[1,,]),digits = 0)  # 3		# number of replicates
R <- 0.5	# proportional change in occupancy
## Get species prevalence out of interest)
prev=rowSums(out$y, dims = 2, na.rm = T)
prev[prev>0]=1
prev = rowSums(prev)
####################
df = matrix(data = NA, nrow = nrow(out$y[,1,]), ncol = 6,  dimnames = list(rownames(out$y[,1,]), c("prob_occ","prob_det","power","sites","prevalence","change%")))
for (jj in 1:nrow(out$y[,1,])) {
#  jj=30

psi1 <- plogis(mean(bet[,jj]))    # 0.3  # initial occupancy probability
p <- plogis(mean(alp[,jj]))      # 0.6  # detection probability
alpha <-0.05	# significance level i.e Probability of a False Positive detection of change
(powerF <- calcPowerFormula(S1=S,S2=S,K1=K,K2=K,p1=p,p2=p,psi1,R,alpha))

# Example B: use function 'calcSFormula' to calculate the 
# number of sites needed to achieve a given power (thick lines in Figure 3)

### POWER = Probability of Detecting a Change on Occupancy (The proportion of change required is given by value R above)
 
pow <- 0.5		# target power level
(SS <- calcSFormula(K1=K,K2=K,p1=p,p2=p,psi1,R,alpha,pow))

df[jj,1] = psi1
df[jj,2] = p
df[jj,3] = powerF
df[jj,4] = SS
df[jj,5] = prev[jj]
df[jj,6] = R

}

length(which(df[,5]<11))

xx=lm(df[,3]~df[,5])
summary(xx)

plot(df[,3]~df[,5])
abline(xx)
ord = order(df[,1], decreasing = T)
df = df[ord,]


# Example C: use function runPowerSims to assess the actual power 
# using simulations (Wald test on probability scale, red lines in Figure 3) 

nsims<-500		# number of simulations to run
powerR <- runPowerSims(S1=SS,S2=SS,K1=K,K2=K,p1=p,p2=p,psi1,R,alpha,nsims)

df=df[df[,5]>32, ]
hist(df[,3])
boxplot(df[,3])
