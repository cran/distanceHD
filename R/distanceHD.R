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
         S = ((n1-1)*S1 + (n2-2)*S2)/(n1+n2-2)        # d x d sample covariance matrix
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
            S = ((n1-1)*S1 + (n2-2)*S2)/(n1+n2-2) + alpha*diag(rd)
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
   ell = matrix(c(rep(1, n1), rep(-1, n2)), n, 1)  
   # library(MASS) # ginv()
   dmdp = 2/norm(ginv(t(zc))%*%ell, type="F"); dmdp  # 24.1
   return(dmdp)
}
###################################

