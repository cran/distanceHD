############################################
# calculate Centroid/Euclidean distance between two groups in a high dimensional space
# also work for a single member cluster
# also work in low dimension (n > d)
############################################

dist_cen <- function(x, group){
   # x: n by d ( d > n)
   # group is a binary group label with the length of n1, n2
  
   if (!(length(unique(group)) == 2))
      stop("group requires two-level")
      
   if ( any(is.na(x)) )
      stop("Missing values in x") 
      
  
   grouplabel = c(1, 2)[factor(group)]; table(grouplabel)

   y = as.matrix(t(x)) ; dim(y)       # d by n
   n = ncol(y)
   d = nrow(y)
   
   n1 = sum(grouplabel==1)
   n2 = sum(grouplabel==2)
   x1 = matrix(x[grouplabel==1, ], n1, d); dim(x1)
   x2 = matrix(x[grouplabel==2, ], n2, d); dim(x2)
          
     # mean of two group (using y)
     g1 = apply(as.matrix(y[, grouplabel==1]), 1, mean ) # mean vector - length of d
     g2 = apply(as.matrix(y[, grouplabel==2]), 1, mean ) # mean vector - length of d
     
     
     # Euclidean distance (1)
     dcen = sqrt(t(g1 - g2)%*%(g1 - g2)); dcen   
     
     # Educlidean distance (2)
     # xbar = aggregate(cbind(x) ~ grouplabel, FUN=mean)[,-1] ; dim(xbar)      # dcen2 = dist(xbar, method = "euclidean") ; dcen2
 
 return(dcen)
  
   
}

############################################
# calculate Mahalanobis distance between two groups in a high dimensional space
# also work for a single member cluster
# also work in low dimension (n > d)
############################################

dist_mah <- function(x, group, alpha){
   # x: n by d ( d > n)
   # group is a binary group label with the length of n1, n2
  
   if (!(length(unique(group)) == 2))
      stop("group requires two-level")
      
      if ( any(is.na(x)))
      stop("Missing value(s) in x") 
  
   grouplabel = c(1, 2)[factor(group)]; table(grouplabel)

   y = as.matrix(t(x)) ; dim(y)       # d by n
   n = ncol(y)
   d = nrow(y)
   
   n1 = sum(grouplabel==1)
   n2 = sum(grouplabel==2)
   x1 = matrix(x[grouplabel==1, ], n1, d); dim(x1)
   x2 = matrix(x[grouplabel==2, ], n2, d); dim(x2)
   
   # computing sample covariance w/ ridge correction
   if (missing(alpha)){ alpha = sqrt(log(d)/n)  }  # default 
   else{ alpha = alpha }   

   
   
   if (d < n) {  # -------------------- (A)
     
     if (n1 == 1){
         S = cov(x2)   # use the 2nd group
         invS = solve(S)
         } else if (n2 == 1){
         S = cov(x1)  # use the 1st group
         invS = solve(S) 
         } else {     
        S1 = cov(x1)
        S2 = cov(x2)
         S = ((n1-1)*S1 + (n2-1)*S2)/(n1+n2-2)        # d x d sample covariance matrix
      invS = solve(S)
     }
     
     # Mahalabnobis distance of two group (using y)
     g1 = apply(as.matrix(y[, grouplabel==1]), 1, mean ) # mean vector - length of d
     g2 = apply(as.matrix(y[, grouplabel==2]), 1, mean ) # mean vector - length of d
     
     dmah = sqrt(t(g1 - g2)%*%invS%*%(g1 - g2)); dmah 
     
   }else{ # d >= n
     
     # singular value decomposition
     svdy = svd(y)        # input - d by n matrix
     U = svdy$u           # d by n
     D = diag(svdy$d)     # n by n
     V = svdy$v           # n by n
     y1 = D %*% t(V)      # rd by n   ---- new data
     rd = nrow(y1)        # reduced dimension = n instead of d
     
     y11 = matrix(y1[, grouplabel==1], rd, n1); dim(y11)
     y12 = matrix(y1[, grouplabel==2], rd, n2); dim(y12)
     
     # ridge correction parameter
     alpha = alpha 
     
     if (n1 == 1) {
           S = cov(t(y12)) + alpha*diag(rd)  # use the 2nd group only
           invS = solve(S) 
           } else if (n2 == 1) {
           S = cov(t(y11)) + alpha*diag(rd) # use the 1st group only
           invS = solve(S)
           }else{
           S1 = cov(t(y11))
           S2 = cov(t(y12))
            S = ((n1-1)*S1 + (n2-1)*S2)/(n1+n2-2) + alpha*diag(rd)
          invS = solve(S); dim(invS)
       }
     
     # ridge Mahalabnobis distance of two group (using y1)
     g1 = apply(y11, 1, mean)  # mean vector - length of rd
     g2 = apply(y12, 1, mean)  # mean vector - length of rd
     
     dmah = sqrt(t(g1 - g2)%*%invS%*%(g1 - g2)); dmah  # 35.68 (n1=1)
     
     
     # Euclidean distance
     # sqrt(t(g1 - g2) %*%(g1 - g2))  # 25.43 , 28.722(n1=1)
     # dist( rbind(g1, g2) )     # 25.43, 28.722(n1=1)
     
   }  # ------------------------------ (end of A)
   
return(dmah)
   
   
   
   
}

############################################
# calculate MDP distance between two groups in a high dimentional space
# also work for a single member cluster
############################################

dist_mdp <- function(x, group){
   # x: n by d ( d > n)
   # group is a binary group label with the length of n1, n2
  
  if (!(length(unique(group)) == 2))
    stop("group requires two-level")
    
    if ( any(is.na(x)))
      stop("Missing value(s) in x") 
  
   grouplabel = c(1, 2)[factor(group)]
   
   z = as.matrix(t(x)) ; dim(z)       # d by n
   n = ncol(z); n
   d = nrow(z); d
   
   n1 = sum(grouplabel==1)
   n2 = sum(grouplabel==2)
   x1 = matrix(x[grouplabel==1, ], n1, d); dim(x1)
   x2 = matrix(x[grouplabel==2, ], n2, d); dim(x2)
 
   # using  z : d by n
   zbar = apply(z, 1, mean )   # length of d
   zbarmat = matrix(rep(zbar, n), ncol=n, byrow = FALSE) # d by n
   zc = z - zbarmat
   #ell = matrix(c(rep(1, n1), rep(-1, n2)), n, 1)  
   ell = matrix(ifelse(grouplabel == 1, 1, -1), n, 1)
   #library(MASS) # ginv()
   dmdp = 2/norm(ginv(t(zc))%*%ell, type="F"); dmdp  # 24.1
   return(dmdp)
}


############################################
# Outlier detection in high-dimensional data
#
# Uses the package distance functions dist_cen(), dist_mah() and dist_mdp():
# the leave-one-out distance of observation j is the distance between the
# single-member group {j} and the group of all other observations.
#
# HDout(data) is the same as
# HDout(data, method = 1:3, Mout = 4, Nsim = 30, plotY = TRUE)
############################################

HDout <- function(data, method = 1:3, Mout = 4, Nsim = 30, plotY = TRUE) {
  # data   : p by n data matrix (p > n); columns are observations
  # method : 1 = centroid (dist_cen), 2 = ridge Mahalanobis (dist_mah),
  #          3 = maximal data piling (dist_mdp); default 1:3 runs all three
  # Mout   : maximum number of outliers to test (default 4)
  # Nsim   : number of bootstrap simulations (default 30)
  # plotY  : TRUE (default) displays the results in a plot
  #
  # returns a list:
  #   dj_c, dj_rM, dj_MDP : sorted distances at the last iteration (Mout);
  #                         0 for a method that was not run
  #   outliers            : observations flagged at steps 1..Mout, per method

  x <- as.matrix(data)
  p <- nrow(x)
  n <- ncol(x)

  if (anyNA(x)) stop("Missing values in input data")
  if (is.null(method) || length(method) == 0) method <- 1:3
  if (is.null(Mout)  || length(Mout)  == 0) Mout  <- 4
  if (is.null(Nsim)  || length(Nsim)  == 0) Nsim  <- 30
  if (is.null(plotY) || length(plotY) == 0) plotY <- TRUE
  plotY <- isTRUE(as.logical(plotY))
  if (!all(method %in% 1:3)) stop("'method' must be in 1:3")

  METHOD <- c("Centroid", "rM", "MDP")

  dj_c <- 0
  dj_rM <- 0
  dj_MDP <- 0
  outliers <- list()

  ## ---------------------------------------------------------------- helpers

  # Leave-one-out distances of the columns of xx (p by m), using the
  # package distance functions (which take an m by p matrix).
  # Observation j is moved to the last row so that it is the second group;
  # dist_mdp() builds its +1/-1 label vector assuming group 1 comes first.
  loo_dist <- function(xx, imethod) {
    m <- ncol(xx)
    group <- c(rep(1, m - 1), 2)
    # ridge constant of HDout.m: sqrt(log(p) / (number of other observations))
    alpha <- sqrt(log(p) / (m - 1))
    vapply(seq_len(m), function(j) {
      xj <- t(xx[, c(seq_len(m)[-j], j), drop = FALSE])   # m by p
      d <- switch(imethod,
                  dist_cen(xj, group),
                  dist_mah(xj, group, alpha = alpha),
                  dist_mdp(xj, group))
      as.numeric(d)
    }, numeric(1))
  }

  # Bootstrap parameters from x without the flagged observation(s)
  # The covariance is S = Xc %*% t(Xc) / (m - 1) + alpha1 * I (alpha1 = 0
  # when n > p), where Xc is the centered p by m data. It is kept in this
  # low-rank form, so the p by p matrix S is never built.
  boot_par <- function(flag) {
    xoutj <- x[, flag, drop = FALSE]
    m <- ncol(xoutj)
    mu <- rowMeans(xoutj)
    B <- (xoutj - mu) / sqrt(m - 1)              # B %*% t(B) = cov(t(xoutj))
    alpha1 <- if (n > p) 0 else sqrt(log(p) / m)
    list(mu = mu, B = B, alpha1 = alpha1)
  }

  # p by m random sample from N(mu, S)   (MATLAB: mvnrnd(mu', S, m)')
  # (exact: B z1 + sqrt(alpha1) z2 has covariance B t(B) + alpha1 I)
  rmvn <- function(par, m) {
    k <- ncol(par$B)
    R <- par$mu + par$B %*% matrix(stats::rnorm(k * m), k, m)
    if (par$alpha1 > 0) {
      R <- R + sqrt(par$alpha1) * matrix(stats::rnorm(p * m), p, m)
    }
    R
  }

  ## ---------------------------------------------------------------- plotting
  if (plotY) {
    op <- graphics::par(mfrow = c(length(method), Mout),
                        mar = c(2.5, 3.5, 2, 1), mgp = c(2, 0.6, 0),
                        cex.axis = 0.8)
    on.exit(graphics::par(op), add = TRUE)
  }

  ## ---------------------------------------------------------------- methods
  for (imethod in method) {

    flag <- seq_len(n)
    subx <- x
    flagged <- integer(0)

    for (Nout in seq_len(Mout)) {
      subn <- n - Nout + 1

      # Leave-out distance on the Y-axis from the data
      dj <- loo_dist(subx, imethod)
      o <- order(dj)
      dj <- dj[o]
      idx <- flag[o]
      maxid <- idx[subn]             # maybe it is the outlier
      flagged <- c(flagged, maxid)

      # Bootstrap parameters from x without the potential outlier(s)
      flag <- flag[flag != maxid]
      par_b <- boot_par(flag)

      # Simulated leave-out distances on the X-axis
      dstarv <- matrix(0, Nsim, subn)
      for (id in seq_len(Nsim)) {
        R <- rmvn(par_b, subn)
        dstarv[id, ] <- sort(loo_dist(R, imethod))
      }
      dj_ <- colMeans(dstarv)

      # Update data
      subx <- x[, flag, drop = FALSE]

      # Figure
      if (plotY) {
        lim <- range(c(dj_, dj))
        lim <- lim + c(-0.05, 0.08) * diff(lim)
        graphics::plot(dj_[-subn], dj[-subn], pch = 1, cex = 0.9,
                       xlim = lim, ylim = lim, asp = 1,
                       xlab = "", ylab = if (Nout == 1) METHOD[imethod] else "",
                       main = paste0("N=", n + 1 - Nout),
                       cex.main = 1, font.main = 1)
        graphics::points(dj_[subn], dj[subn], pch = 19, cex = 1.4)
        graphics::text(dj_[subn], dj[subn], labels = idx[subn],
                       pos = 3, offset = 0.5, cex = 0.8, xpd = NA)
      }
    }

    if (imethod == 1) dj_c <- dj
    if (imethod == 2) dj_rM <- dj
    if (imethod == 3) dj_MDP <- dj
    outliers[[METHOD[imethod]]] <- flagged
  }

  invisible(list(dj_c = dj_c, dj_rM = dj_rM, dj_MDP = dj_MDP,
                 outliers = outliers))
}
###################################
