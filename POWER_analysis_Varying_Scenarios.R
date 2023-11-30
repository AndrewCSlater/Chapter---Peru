##################################################################
# Power analysis for occupancy studies under imperfect detection
# Examples on how to use the provided functions
# Authors: Gurutzeta Guillera-Arroita & Jose J. Lahoz-Monfort
##################################################################

load("models/1.B-Training_Min11sp_5PCA_50000_25000_50_1.Rdata")
out = out.sfMsPGOcc
rm(out.sfMsPGOcc); gc()

### Sumarise Species by Site
y = out$y # 135 Species in 1101 Stations
y2 = as.data.frame(rowSums(y, dims = 2, na.rm = T)) ## Number of replicates each species seen in per station
y2 = y2[colSums(y2)>0]                              ## Remove stations with no sightings
y3 = y2; y3[y3>0] = 1                               ## Species * Station Presence/Absence
summary(colSums(y2))
# zz = which(rowSums(y2)>110)
# zzz = which(rowSums(y3)>54)
# min.sp = unique(c(zz,zzz))

quantiles = c("1st_quantile", "Median", "3rd_Quantile")

bet = as.data.frame(out$beta.samples)  # Probabilities of OCCURRENCE
bet = as.data.frame(rowMeans(out$psi.samples, na.rm = T, dims = 2)) # Mean prob of occurrence across all stations per replicate
alp = out$alpha.samples # Probabilities of DETECTION

rr = c(0.25, 0.5, 0.75) # Range of Proportional Changes in occupancy to be examined

par(mfrow = c(1,1))
par(mar = c(4,4,1,1))

for (iii in 1:length(rr)) {  # Number of proportional changes 
# iii=1
R <- rr[iii]  # 0.5	# proportional change in occupancy

ss = c(1000, 3000, 5000, 7000, 9000) # Range of number of Sites
kk = c(3, 9, 15, 21, 27) # Range of number of Replicates

mat = matrix(data = NA, nrow = length(ss), ncol = length(kk), dimnames = list(ss,kk))

quants = c(2,3,5)   # select 1st quant, median & 3rd Quant from Summary statistic
for (q in 1:3) {    # For Each Quantile
# q=1
qq = quants[q]  
for (i in 1:length(ss) ) { # Range of number of Sites
# i = 1
for (ii in 1:length(kk)) { # Range of number of Replicates
# ii = 1
S = ss[i]
K = kk[ii]

sp.power = c()
# for (jj in 1:nrow(out$y[,1,])) { # For each species
for (jj in 1:nrow(y2)) { # For each species 
# jj=1
  # psi1 <- plogis(mean(bet[,jj]))    # species specific occupancy probability (Intercept)
  psi1 = mean(bet[,jj])
  p <- plogis(mean(alp[,jj]))    		# species specific detection probability (Intercept)
 alpha <-0.05	# significance level (FALSE POSITIVE ie detect effect where none exists)

 # Example A: use function 'calcPowerFormula' to calculate the
# power of a design using the formula (as in Figure 1)
  (powerF <- calcPowerFormula(S1=S,S2=S,K1=K,K2=K,p1=p,p2=p,psi1,R,alpha))
  sp.power = append(sp.power, powerF)
  }  # End Species loop
sp.power = summary(sp.power)[qq]  ## [2]=1stQuant, [3]=Median, [5]=3rdQuant
mat[i,ii] = sp.power
  } # End Replicates loop
  
# Example B: use function 'calcSFormula' to calculate the
# number of sites needed to achieve a given power (thick lines in Figure 3)
#   pow <- 0.7		# target power level i.e. Probability of detecting an effect of magnitude R
#   (SS <- calcSFormula(K1=K,K2=K,p1=p,p2=p,psi1,R,alpha,pow))
#   sites.needed = append(sites.needed, SS)
} # End Sites loop
mat = round(x = mat, digits = 2)

filename = paste0("images/Power_Analysis_R",rr[iii],"_",quantiles[q],  ".png")
png(filename=filename, res=600, width = 11.69, height = 8.27, units = "in")

matplot(as.data.frame(t(mat)), type = "l", ylab = "Probability of Detecting", xaxt="n", xlab = "Replicates", ylim = c(0.05, 1))
title(main = paste0("Proportional Change ",rr[iii]), line = -1)
title(main = paste0(quantiles[q]), line = -2, cex.main=1)
legend(x = 1, y = 1, col = c(1:5), lty = c(1:5),  legend = rownames(mat), title = "Number of Stations", cex = 0.8)
axis(side = 1, at = c(1:5), labels = colnames(mat))

# matplot(as.data.frame(mat), type = "l", ylab = "Power", xaxt="n", xlab = "Sites", ylim = c(0.05, 0.7))
# title(main = paste0("Proportional Change ",rr[iii]))
# legend(x = 1, y = 0.7, col = c(1:5), lty = c(1:5),  legend = colnames(mat), title = "Number of Replicates", cex = 0.8)
# axis(side = 1, at = c(1:5), labels = rownames(mat))
dev.off()
} # End Quantile Loop
} # End Proportional change loop

rr = c(0.25, 0.5, 0.75)

par(mfrow = c(1,3))

for (iii in 1:3) {
# iii=1
R <- rr[iii]  # 0.5	# proportional change in occupancy
# R = 0.5 # Proportional Change in Occupancy

kk = c(6, 12, 24, 48, 96) # Range of number of Replicates
pp = c(0.5, 0.6, 0.7, 0.8, 0.9) # Range of Powers
mat = matrix(data = NA, nrow = 5, ncol = 5, dimnames = list(kk,pp))

for (ii in 1:length(kk)) {
  # ii = 1
  K = kk[ii]

for (iv in 1:5) {
  # iv = 1
pow = pp[iv]
# pow <- 0.7		# target power level i.e. Probability of detecting an effect of magnitude R  

sites.needed = c()
# for (jj in 1:nrow(out$y[,1,])) { # For each species
for (jj in min.sp) { # For each species  
  #jj=1
  psi1 <- plogis(mean(bet[,jj])) # species specific occupancy probability
  p <- plogis(mean(alp[,jj]))    # species specific detection probability
  alpha <-0.05	# significance level (FALSE POSITIVE ie detect effect where none exists)
  
# Example B: use function 'calcSFormula' to calculate the 
# number of sites needed to achieve a given power (thick lines in Figure 3)
(SS <- calcSFormula(K1=K,K2=K,p1=p,p2=p,psi1,R,alpha,pow))
sites.needed = append(sites.needed, SS)
} # End species loop
sites.needed = round(mean(sites.needed), digits = 0)
mat[ii, iv] = sites.needed

testdf = as.data.frame(as.table(mat))
plot(testdf$Freq)
testdf$Var1 = as.character(testdf$Var1)
testdf$Var2 = as.character(testdf$Var2)
plot(testdf$Var1, testdf$Var2, )

} # End Power loop
} # End Replicates loop
  
# plot(mat[5,])
matplot(as.data.frame(t(mat)), type = "l", ylab = "Number of Sites", xaxt="n", xlab = "Power to Observe Change", ylim = c(14000, max(mat)))
title(main = paste0("Proportional Change ",rr[iii]))
legend(x = 1, y = max(mat), col = c(1:5), lty = c(1:5),  legend = rownames(mat), title = "Number of Replicates")
axis(side = 1, at = c(1:5), labels = colnames(mat))

} # End proportional change loop

length(which(sp.power>=0.7)) # 5 species with a 70% probability of detecting a 50% difference in occupancy between groups (with a mean of 3 replicate surveys)
# Example C: use function runPowerSims to assess the actual power 
# using simulations (Wald test on probability scale, red lines in Figure 3) 

nsims<-5000		# number of simulations to run
powerR <- runPowerSims(S1=SS,S2=SS,K1=K,K2=K,p1=p,p2=p,psi1,R,alpha,nsims)
