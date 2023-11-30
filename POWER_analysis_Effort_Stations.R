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

# bet = as.data.frame(out$beta.samples)  # Probabilities of OCCURRENCE
occ = as.data.frame(rowMeans(out$psi.samples, na.rm = T, dims = 2)) # Mean prob of occurrence across all stations per replicate
alp = out$alpha.samples # Probabilities of DETECTION for a SINGLE SURVEY

rr = c(0.1, 0.25, 0.5, 0.75) # Range of Proportional Changes in occupancy to be examined

# par(mfrow = c(1,1))

par(mar = c(2,2,1.5,1))

filename = paste0("images/Power_Analysis_R.png")
png(filename=filename, res=600, width = 7, height = 5, units = "in")

# create layout
par(mar = c(2,2,1.5,1))
layout(matrix(c(1, 2, 3, 4), nrow = 2, 
              ncol = 2, byrow = TRUE))

for (iii in 1:length(rr)) {  # Number of proportional changes 
  # iii=1
  R <- rr[iii]  # 0.5	# proportional change in occupancy
  
  ss = c(100, 500, 1000, 2500, 5000) # Range of number of Sites
  kk = c(1, 3, 5, 10, 15, 20, 25) # Range of number of Replicates
  effort = sort(ss%*%t(kk)) %>% unique()
  # effort = effort[effort<100001]  ## Reduce MAx Effort to 100,000 surveys
  
mat = matrix(data = NA, nrow = length(ss), ncol = length(effort), dimnames = list(ss,effort))
mat.1st = mat
mat.3rd = mat

for (i in 1:length(ss) ) { # Range of number of Sites
  # i = 1
    for (ii in 1:length(effort)) { # Range of Effort (relates back to number of Replicates)
      # ii = 2
        S = ss[i]
        K = round(effort[ii]/S, digits = 0) # kk[ii] Number of Replicates based on Effort & Number of survey stations
        
        sp.power = c()
        for (jj in 1:nrow(y2)) { # For each species 
          # jj=1
          psi1 = mean(occ[,jj])             # species specific occupancy probability
          p <- plogis(mean(alp[,jj]))    		# species specific detection probability (Intercept)
          alpha <-0.05	# significance level (FALSE POSITIVE ie detect effect where none exists)
          
# Example A: use function 'calcPowerFormula' to calculate the
# Power of a design using the formula (as in Figure 1)
(powerF <- calcPowerFormula(S1=S,S2=S,K1=K,K2=K,p1=p,p2=p,psi1,R,alpha))
sp.power = append(sp.power, powerF)
}  # End Species loop

med = summary(sp.power)[3]  ## [2]=1stQuant, [3]=Median, [5]=3rdQuant
frst = summary(sp.power)[2]
thrd = summary(sp.power)[5]
mat[i,ii] = med
mat.1st[i,ii] = frst
mat.3rd[i,ii] = thrd
} # End Replicates loop
      
} # End Sites loop
mat = round(x = mat, digits = 2)
    
matplot(as.data.frame(t(mat)), type = "l", lty = 1,  lwd = 1, xaxt="n", ylab = "", ylim = c(0.05, 1)) # ylab = "Probability of detecting change",  xlab = "Effort - Stations * Replicate surveys"
legend(x = 1, y = 0.9, col = c(1:5), lty = 1,  legend = rownames(mat), title = "Number of Stations", cex = 0.5)
title(main = paste0("Proportional Change ",rr[iii]), line = -1, cex.main=0.7, adj=0.05)
# title(main = "Median", line = -2, cex.main=1)

## !!!!  UNBLANK the following to Add Polygon Between 1st & 3rd Quantiles  !!!!    
# for (gg in 1:nrow(mat)) {
# #gg=5
# xx = c(1:length(mat.3rd[gg,]), length(mat[gg,]):1)
# yy = c(mat.3rd[gg,], rev(mat.1st[gg,]))
# polygon(x = xx, y = yy, col = gg, lty = 1, density = 5+gg, angle = 180)
# polygon(x = xx, y = yy, col = gg, lty = 1, density = 5+gg, angle = 270)
#     
# matplot(as.data.frame(t(mat)), type = "l", lty = 1,  lwd = 2,  ylab = "Probability of detecting change", xaxt="n", xlab = "Effort - Stations * Replicate surveys", ylim = c(0.05, 1), add = T)
# title(main = paste0("Proportional Change ",rr[iii]), line = -1)
# title(main = "Median", line = -2, cex.main=1)

# legend(x = 1, y = 1, col = c(1:5), lty = 1,  legend = rownames(mat), title = "Number of Stations", cex = 0.8)
axis(side = 1, at = seq(1, length(effort), 2), labels = as.character(effort[seq(1, length(effort), 2)]))
    
  } # End Quantile Loop
dev.off()
} # End Proportional change loop

# 
# 
# rr = c(0.25, 0.5, 0.75)
# 
# par(mfrow = c(1,3))
# 
# for (iii in 1:3) {
#   # iii=1
#   R <- rr[iii]  # 0.5	# proportional change in occupancy
#   # R = 0.5 # Proportional Change in Occupancy
#   
#   kk = c(6, 12, 24, 48, 96) # Range of number of Replicates
#   pp = c(0.5, 0.6, 0.7, 0.8, 0.9) # Range of Powers
#   mat = matrix(data = NA, nrow = 5, ncol = 5, dimnames = list(kk,pp))
#   
#   for (ii in 1:length(kk)) {
#     # ii = 1
#     K = kk[ii]
#     
#     for (iv in 1:5) {
#       # iv = 1
#       pow = pp[iv]
#       # pow <- 0.7		# target power level i.e. Probability of detecting an effect of magnitude R  
#       
#       sites.needed = c()
#       # for (jj in 1:nrow(out$y[,1,])) { # For each species
#       for (jj in min.sp) { # For each species  
#         #jj=1
#         psi1 <- plogis(mean(bet[,jj])) # species specific occupancy probability
#         p <- plogis(mean(alp[,jj]))    # species specific detection probability
#         alpha <-0.05	# significance level (FALSE POSITIVE ie detect effect where none exists)
#         
#         # Example B: use function 'calcSFormula' to calculate the 
#         # number of sites needed to achieve a given power (thick lines in Figure 3)
#         (SS <- calcSFormula(K1=K,K2=K,p1=p,p2=p,psi1,R,alpha,pow))
#         sites.needed = append(sites.needed, SS)
#       } # End species loop
#       sites.needed = round(mean(sites.needed), digits = 0)
#       mat[ii, iv] = sites.needed
#       
#       testdf = as.data.frame(as.table(mat))
#       plot(testdf$Freq)
#       testdf$Var1 = as.character(testdf$Var1)
#       testdf$Var2 = as.character(testdf$Var2)
#       plot(testdf$Var1, testdf$Var2, )
#       
#     } # End Power loop
#   } # End Replicates loop
#   
#   # plot(mat[5,])
#   matplot(as.data.frame(t(mat)), type = "l", ylab = "Number of Sites", xaxt="n", xlab = "Power to Observe Change", ylim = c(14000, max(mat)))
#   title(main = paste0("Proportional Change ",rr[iii]))
#   legend(x = 1, y = max(mat), col = c(1:5), lty = c(1:5),  legend = rownames(mat), title = "Number of Replicates")
#   axis(side = 1, at = c(1:5), labels = colnames(mat))
#   
# } # End proportional change loop
# 
# length(which(sp.power>=0.7)) # 5 species with a 70% probability of detecting a 50% difference in occupancy between groups (with a mean of 3 replicate surveys)
# # Example C: use function runPowerSims to assess the actual power 
# # using simulations (Wald test on probability scale, red lines in Figure 3) 
# 
# nsims<-5000		# number of simulations to run
# powerR <- runPowerSims(S1=SS,S2=SS,K1=K,K2=K,p1=p,p2=p,psi1,R,alpha,nsims)
