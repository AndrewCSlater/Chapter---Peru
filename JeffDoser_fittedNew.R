fittedNew <- function(object, ...) {
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
      if(! i %in% formal.args)
          warning("'",i, "' is not an argument")
  }
  # Call ----------------------------------------------------------------
  cl <- match.call()
  # Functions -------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

  # Some initial checks -------------------------------------------------
  # Object ----------------------------
  if (missing(object)) {
    stop("error: object must be specified")
  }
  # if (!is(object, c("msPGOcc", "spMsPGOcc", "lfMsPGOcc", "sfMsPGOcc"))) {
  if (!(class(object) %in% c('msPGOcc', 'spMsPGOcc', 'lfMsPGOcc', 'sfMsPGOcc', 
			     'svcMsPGOcc'))) {
    stop("error: object must be of class msPGOcc, spMsPGOcc, lfMsPGOcc, sfMsPGOcc, or svcMsPGOcc\n")
  }
  n.post <- object$n.post * object$n.chains
  X.p <- object$X.p
  y <- object$y
  n.rep <- apply(y[1, , , drop = FALSE], 2, function(a) sum(!is.na(a)))
  K.max <- dim(y)[3]
  J <- dim(y)[2]
  N <- dim(y)[1]
  if (object$pRE) {
    if (nrow(object$X.p.re) == dim(y)[2]) {
      X.p.re <- do.call(rbind, replicate(dim(y)[3], object$X.p.re, simplify = FALSE))
      X.p.re <- X.p.re[!is.na(c(y[1, , ])), , drop = FALSE]
      # Add 1 to get it to R indexing. 
      X.p.re <- X.p.re + 1
    } else {
      # Add 1 to get it to R indexing.
      X.p.re <- object$X.p.re + 1
    }
  }
  if (nrow(X.p) == dim(y)[2]) {
    X.p <- do.call(rbind, replicate(dim(y)[3], X.p, simplify = FALSE))
    X.p <- X.p[!is.na(c(y[1, , ])), , drop = FALSE]
  }
  z.long.indx <- rep(1:J, dim(y)[3])
  z.long.indx <- z.long.indx[!is.na(c(y[1, , ]))]
  psi.samples <- object$psi.samples
  alpha.samples <- object$alpha.samples
  n.obs <- nrow(X.p)
  det.prob.samples <- array(NA, dim = c(n.obs, N, n.post))
  sp.indx <- rep(1:N, ncol(X.p))
  y <- matrix(y, N, J * K.max)
  y <- y[, apply(y, 2, function(a) !sum(is.na(a)) > 0)]
  for (i in 1:N) {
    if (object$pRE) {
      sp.re.indx <- rep(1:N, each = ncol(object$alpha.star.samples) / N)
      tmp.samples <- matrix(0, n.post, n.obs)
      tmp.alpha.star <- object$alpha.star.samples[, sp.re.indx == i]
      for (j in 1:ncol(X.p.re)) {
        tmp.samples <- tmp.samples + tmp.alpha.star[, X.p.re[, j]]
      }
      det.prob.samples[, i, ] <- logit.inv(X.p %*% t(alpha.samples[, sp.indx == i]) + t(tmp.samples))
    } else {
      det.prob.samples[, i, ] <- logit.inv(X.p %*% t(alpha.samples[, sp.indx == i]))
    }
  }

  out <- list()
  # Get detection probability
  # Need to be careful here that all arrays line up. 
  det.prob.samples <- aperm(det.prob.samples, c(3, 2, 1))
  tmp <- array(NA, dim = c(n.post, N, J * K.max))
  names.long <- which(!is.na(c(object$y[1, , ])))
  tmp[, , names.long] <- det.prob.samples
  p.samples <- array(tmp, dim = c(n.post, N, J, K.max))
  out$p.samples <- p.samples
  # Get fitted values
  det.prob.samples <- det.prob.samples * psi.samples[, , z.long.indx]
  y.rep.samples <- array(NA, dim = dim(det.prob.samples))
  for (i in 1:N) {
    y.rep.samples[, i, ] <- apply(det.prob.samples[, i, ], 2, function(a) rbinom(n.post, 1, a))
  }
  tmp <- array(NA, dim = c(n.post, N, J * K.max))
  names.long <- which(!is.na(c(object$y[1, , ])))
  tmp[, , names.long] <- y.rep.samples
  y.rep.samples <- array(tmp, dim = c(n.post, N, J, K.max))
  out$y.rep.samples <- y.rep.samples
  return(out)
}

