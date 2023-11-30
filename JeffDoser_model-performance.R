# model-performance.R: script that simulates a spatial factor multi-species
#                      occupancy model and assess the performance of the model
#                      using AUC. Also shows how to extract a confusion matrix
#                      for each iteration of the MCMC, as well as sensitivity
#                      and specificity values.
rm(list = ls())
library(spOccupancy)
# For calculating AUC
library(pROC)

# Simulate Data -----------------------------------------------------------
J.x <- 20
J.y <- 20 
J <- J.x * J.y
n.rep<- sample(2:4, size = J, replace = TRUE)
N <- 5
# Community-level covariate effects
# Occurrence
beta.mean <- c(-1, 1.2, -0.9)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 1.5, 2.3)
# Detection
alpha.mean <- c(0, 1.2, -0.5)
tau.sq.alpha <- c(1, 0.5, 1.3)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list()
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
n.factors <- 1
alpha.true <- alpha
phi <- runif(n.factors, 3 / 0.9, 3 / 0.2)
sigma.sq <- runif(n.factors, 0.5, 2.3)
# phi <- runif(N, 3 / 0.9, 3 / 0.2)
# sigma.sq <- runif(N, 0.5, 2.3)
nu <- rep(2, n.factors)

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
		psi.RE = psi.RE, p.RE = p.RE, phi = phi, sp = TRUE, 
		cov.model = 'exponential', n.factors = n.factors, factor.model = TRUE, 
                sigma.sq = sigma.sq)

y <- dat$y
X <- dat$X
X.p <- dat$X.p
coords <- dat$coords

occ.covs <- cbind(X)
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2')
det.covs <- list(det.cov.1 = X.p[, , 2], 
		 det.cov.2 = X.p[, , 3])
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, 
                  det.covs = det.covs)
# Priors. notice how the half-t priors are specified.
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72), 
		   tau.beta.half.t = list(df = 1, A = 25),
		   # tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
		   tau.alpha.half.t = list(df = 1, A = 25))
		   # tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Starting values
inits.list <- list(alpha.comm = 0, 
		      beta.comm = 0, 
		      beta = 0, 
		      alpha = 0,
		      tau.sq.beta = 1, 
		      tau.sq.alpha = 1, 
		      z = apply(y, c(1, 2), max, na.rm = TRUE)) 
# Tuning
tuning.list <- list(phi = 1)

out <- sfMsPGOcc(occ.formula = ~ occ.cov.1 + occ.cov.2,
                    det.formula = ~ det.cov.1 + det.cov.2,
                    data = data.list,
                    inits = inits.list,
                    n.batch = 100,
		    n.factors = n.factors, 
                    batch.length = 25,
                    accept.rate = 0.43,
                    priors = prior.list,
                    cov.model = "exponential",
                    tuning = tuning.list,
                    n.omp.threads = 1,
                    verbose = TRUE,
                    NNGP = TRUE,
                    n.neighbors = 5,
                    search.type = 'cb',
                    n.report = 10,
                    n.burn = 1500,
		    n.thin = 1,
		    n.chains = 3)

summary(out)

# Posterior predictive check ----------------------------------------------
ppc.out <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out)
# Get the actual Bayesian p-values for individual species and community
# Species-level
species.bpvs <- apply(ppc.out$fit.y.rep > ppc.out$fit.y, 2, mean)
species.bpvs
# Community-level
comm.bpv <- mean(ppc.out$fit.y.rep > ppc.out$fit.y)
comm.bpv

# Calculate AUC using the detection-nondetection values -------------------
# Extract the fitted values
y.rep.samples <- fitted(out)$y.rep.samples
str(y.rep.samples)
# Will calculate an AUC value for each species and each iteration of the 
# MCMC in order to get an AUC estimate with uncertainty.
n.samples <- out$n.post * out$n.chains
N <- nrow(data.list$y)
auc.vals <- matrix(NA, n.samples, N)
for (j in 1:n.samples) {
  print(j)
  for (i in 1:N) {
    auc.vals[j, i] <- auc(response = c(y[i, , ]), predictor = c(y.rep.samples[j, i, , ]))
  } # i (species)
} # j (iteration)

# Calculate quantiles of AUC values for each species
auc.quants <- t(apply(auc.vals, 2, quantile, c(0.025, 0.5, 0.975)))
# The 50% quantile is our median estimate, and the 2.5 and 97.5% quantiles
# form our credible interval for the AUC for each species.
auc.quants

# Calculate species-specific confusion matrices ---------------------------
# Number of non-missing y values for each species
n.obs <- sum(!is.na(c(y[1, , ])))
# True positives
tp.samples <- matrix(NA, n.samples, N)
# True negatives
tn.samples <- matrix(NA, n.samples, N)
# False positives
fp.samples <- matrix(NA, n.samples, N)
# False negatives
fn.samples <- matrix(NA, n.samples, N)
for (j in 1:n.samples) {
  print(j)
  for (i in 1:N) {
    tp.samples[j, i] <- sum((c(y[i, , ]) == 1) & (c(y.rep.samples[j, i, , ]) == 1), na.rm = TRUE)
    tn.samples[j, i] <- sum((c(y[i, , ]) == 0) & (c(y.rep.samples[j, i, , ]) == 0), na.rm = TRUE)
    fp.samples[j, i] <- sum((c(y[i, , ]) == 0) & (c(y.rep.samples[j, i, , ]) == 1), na.rm = TRUE)
    fn.samples[j, i] <- sum((c(y[i, , ]) == 1) & (c(y.rep.samples[j, i, , ]) == 0), na.rm = TRUE)
  }
}
# Calculate Sensitivity
sensitivity.samples <- tp.samples / (tp.samples + fn.samples)
# Calculate Specificity
specificity.samples <- tn.samples / (tn.samples + fp.samples)
# Calculate quantiles of sensitivity and specificity values
sensitivity.quants <- t(apply(sensitivity.samples, 2, quantile, c(0.025, 0.5, 0.975)))
specificity.quants <- t(apply(specificity.samples, 2, quantile, c(0.025, 0.5, 0.975)))
# The 50% quantile is our median estimate, and the 2.5 and 97.5% quantiles
# form our credible interval for the sensitivity and specificity
sensitivity.quants
specificity.quants

