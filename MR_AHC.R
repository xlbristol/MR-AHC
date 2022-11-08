

### MR-AHC method for detecting IV clusters and IV estimation using two sample summary statistics in Mendelian randomization (MR) analyses.

### MR_AHC
### Function: MR-AHC method for detecting IV clusters and IV estimation using two sample summary statistics in Mendelian randomization (MR) analyses.

### Input: (1) betaY: A numeric vector of SNP-outcome associations (p_z by 1 vector).
###        (2) betaX: A numeric vector of SNP-exposure associations (p_z by 1 vector).
###        (3) seY: A numeric vector of the standarnd errors of the SNP-outcome associations (p_z by 1 vector).
###        (4) seX: A numeric vector of the standarnd errors of the SNP-exposure associations (p_z by 1 vector).
###        (5) n: A numeric scalar specifiying the sample size.
###        (6) alpha: A numeric scalar between 0 and 1 specifying the significance level for the confidence interval (default = 0.05).
###        (7) tuning: A numeric scalar specifiying the threshold p-value for the Q test (default = 0.1/log(length(betaY))).
###        (8) weak: Logical. If weak = TRUE, the weak IV robust Q test is used (default = TRUE).
###        (9) smallcluster: A numeric scalar specifiying small clusters (default = 4).

### Output: (1) Cluster_number: The number of clusters detected by the algorithm (including the junk cluster).
###         (2) AHC_cluster: All the clustesrs detected by the algorithm. Instruments names (or identities, if variable names are not specified) are listed in each cluster.
###         (3) Null_cluster: Instrumental variables in the null cluster.
###         (4) Junk_cluster: Instrumental variables in the junk cluster.
###         (5) AHC_results: A matrix that summarizes the clustering and estimation results, including:
###           (a) length.IV: The number of IVs in each cluster.
###           (b) Q.IV: The p-value for the Q test of the instruments in each cluster.
###           (c) I.IV: The I statistic of the instruments in each cluster.
###           (d) beta.IV: the point estimate estimated with each cluster.
###           (e) se.IV: the standard error for the causal estimate in each cluster
###           (f) t.IV: the t statistic.
###         (6) F: the F statistic for all the IVs.

### Source code for the internal functions 'dissim', 'getnns' and 'hierclust' is available at
### https://github.com/mda-sw/library/blob/master/hierclus/hcluswtd.r

MR_AHC <- function(betaY, betaX, seY, seX, n, alpha = 0.05, tuning = 0.1/log(n), weak = TRUE, smallcluster = 4){
  
  # Define Constants
  pz <- length(betaY);
  P <- 1; ## this is for single outcome--single exposure case
  
  # F statistic
  F_bar = mean(betaX^2/seX^2)
  
  # Functions for the Q Test
  if (weak){
    # weak-IV robust Q
    Qsq = function(bx,by,sx,sy){
      biv        = by/bx
      k          = length(bx)
      df         = k - 1
      PL2 = function(a){
        b = a[1]
        w = 1/(sy^2/bx^2 + (b^2)*sx^2/bx^2)
        q = sum(w*(biv - b)^2)
      }
      Bhat = optimize(PL2,interval=c(-10,10))$minimum
      W = 1/(sy^2/bx^2 + (Bhat^2)*sx^2/bx^2)
      QE = sum(W*(biv - Bhat)^2)
      QEp = 1-pchisq(QE,df)
      Isq = (QE - df)/QE
      Isq =  max(0,Isq)
      sehat = sqrt(solve(sum(W)))
      return(list(QE, Bhat, sehat, Isq))
    }
  }else{
    Qsq = function(bx,by,sx,sy){
      k          = length(bx)
      df         = k - 1
      y          = matrix(by/bx, nrow = k)
      s          = matrix(sy/bx, nrow = k)
      w          = 1/s^2; sum.w  = sum(w)
      mu.hat     = sum(y*w)/sum.w
      Q          = sum(w*(y-mu.hat)^2)
      se.hat     = sqrt(solve(sum.w))
      Isq        = (Q - df)/Q
      Isq        =  max(0,Isq)
      return(list(Q, mu.hat, se.hat, Isq))
    }
  }
  
  # Q test involving all instruments
  Qsq.all <- Qsq(betaX, betaY, seX, seY);
  Qsq.all.p <- pchisq(Qsq.all[[1]], (pz - P), lower.tail = FALSE);
  
  if (Qsq.all.p > tuning){
    # all the IVs are selected as in the same cluster. No need for clustering.
    clusterIV <- list(c(1:pz));
    Q.IV <- list(Qsq.all.p);
    clusternames <- paste0('cluster', 1);
    beta.IV <- Qsq.all[[2]];
    se.IV <- Qsq.all[[3]];
    I.IV <- Qsq.all[[4]];
  }else{
    # AHC procedure
    
    comb <- combn(pz, P, repeats.allowed=F); Len <- ncol(comb);
    bvec <- matrix(betaY/betaX, nrow = pz);
    sevec <- matrix(abs(seY/betaX), nrow = pz);
    
    dend = hierclust(matrix(bvec, ncol = 1), matrix(1/sevec^2, ncol = 1));
    
    # This is to generate the selection path, i.e. the "tree" (algorithm 1)
    trees = list();
    for (b in 1:(Len - 1)){ 
      dstep <- cutree(dend, k=b); 
      treesb = list();
      for (bb in 1:length(unique(dstep))){
        treesb[[bb]] = comb[, which(dstep == (unique(dstep))[bb])];
      }
      trees[[b]] = treesb;
    }
    
    # Downward testing (algorithm 2)
    start <- 2;
    clusterIV = Q.IV = beta.IV = se.IV = I.IV = clusternames <- list();
    while (start <= length(trees)) {
      branch = trees[[start]];
      index_single = which(lapply(branch, length) == 1);
      if (length(index_single) != 0) branch = branch[-index_single];
      for (k in 1:length(branch)){
        leaf <- sort(unique(as.vector(branch[[k]])));
        check.subset <- function(y){all(leaf %in% y)};
        index = any(unlist(lapply(clusterIV, check.subset)));
        if (length(leaf)< pz && index == FALSE){
          Q.branch <- Qsq(betaX[leaf], betaY[leaf], seX[leaf], seY[leaf])[[1]];
          beta.branch <- Qsq(betaX[leaf], betaY[leaf], seX[leaf], seY[leaf])[[2]];
          se.branch <- Qsq(betaX[leaf], betaY[leaf], seX[leaf], seY[leaf])[[3]];
          p.branch <- pchisq(Q.branch, df = (length(leaf) - P), lower.tail = FALSE);
          I.branch <- Qsq(betaX[leaf], betaY[leaf], seX[leaf], seY[leaf])[[4]];
          if (p.branch > tuning){
            clusterIV[[length(clusterIV)+1]] <- leaf;
            Q.IV[[length(Q.IV)+1]] <- p.branch;
            beta.IV[[length(beta.IV)+1]] <- beta.branch;
            se.IV[[length(se.IV)+1]] <- se.branch;
            I.IV[[length(se.IV)+1]] <- I.branch;
            clusternames[[length(clusternames)+1]] <- paste0('cluster', length(clusternames)+1)
          }
        }
      }
      clusterIV;
      start <- start + 1;
    }
  }
  
  # Estimation results for the selected clusters
  length.IV <- unlist(lapply(clusterIV, length));
  beta.IV <- unlist(beta.IV);
  Q.IV <- unlist(Q.IV);
  I.IV <- unlist(I.IV)
  se.IV <- unlist(se.IV);
  t.IV <- abs(beta.IV/se.IV);
  
  # Summary of the clustering and estimation results
  AHC_results <- cbind(1:length(clusterIV), length.IV, Q.IV, I.IV, beta.IV, se.IV, t.IV);
  colnames(AHC_results) <- c("ID", "length.IV", "Q.IV", "I.IV ", "beta.IV", "se.IV", "t.IV")
  AHC_cluster <- clusterIV
  AHC_cluster_junk <- sort(c(1:pz)[-unlist(AHC_cluster)]) # IVs that do not belong to any clusters are classified as in the junk cluster.
  
  if (length(AHC_cluster_junk) != 0){AHC_cluster[[length(AHC_cluster) + 1]] <- AHC_cluster_junk}else{AHC_cluster_junk <- NULL}
  
  null_index <- which(pchisq((AHC_results[,"t.IV"])^2, df = 1, lower.tail = FALSE) > 0.1/log(n)) # Clusters that do not give significant estimates are classified as the null cluster.
  if (length(null_index) != 0){AHC_cluster_null <- sort(unlist(AHC_cluster[c(null_index)]))}else{AHC_cluster_null <- NULL}
  
  # trimming the small clusters
  small.index = which(lapply(AHC_cluster, length) <= smallcluster)
  if (length(small.index) != 0){
    
    AHC_results <- AHC_results[-small.index, ]
    
    AHC_cluster_junk_new <- sort(unlist(AHC_cluster[small.index]))
    AHC_cluster_junk = sort(union(AHC_cluster_junk_new, AHC_cluster_junk))
    AHC_cluster <- AHC_cluster[-small.index]
    AHC_cluster[[length(AHC_cluster) + 1]] <- AHC_cluster_junk
    
    index.null_junk = which (AHC_cluster_null %in% AHC_cluster_junk == TRUE)
    if(length(index.null_junk) != 0){AHC_cluster_null = AHC_cluster_null[-index.null_junk]}
    if(length(AHC_cluster_null)==0){AHC_cluster_null = NULL}
    
  }
  
  Nr.AHC <- length(AHC_cluster)
  
  # Report results
  results = list( Cluster_number = Nr.AHC,
                  AHC_cluster = AHC_cluster,
                  Null_cluster = AHC_cluster_null,
                  Junk_cluster = AHC_cluster_junk,
                  F = F_bar,
                  AHC_results = AHC_results)
  return(results)
  
}


# Internal Function (1): the dissimilarity matrix with weighting. Not to be called by users.
dissim <- function(a, wt) {
  # Inputs.   a: matrix, for which we want distances on rows,
  #           wt: masses of each row.
  # Returns.  matrix of dims. nrow(a) x nrow(a) with wtd. sqd. Eucl. distances.
  # FM, 2003/11/16
  
  n <- nrow(a)
  m <- ncol(a)
  adiss <- matrix(0, n, n)
  
  for (i1 in 2:n) {
    adiss[i1,i1] <- 0.0
    for (i2 in 1:(i1-1)) {
      adiss[i1,i2] <- 0.0
      for (j in 1:m) {
        # We use the squared Euclidean distance, weighted.
        adiss[i1,i2] <- adiss[i1,i2] + (wt[i1]*wt[i2])/(wt[i1]+wt[i2]) *
          (a[i1,j]-a[i2,j])^2
      }
      adiss[i2,i1] <- adiss[i1,i2]
    }
  }
  adiss
}

# Internal Function (2), not to be called by users.
getnns <- function(diss, flag) {
  # Inputs.  diss: full distance matrix.
  #          flag: "live" rows indicated by 1 are to be processed.
  # Returns. List of: nn, nndiss.
  #          nn:   list of nearest neighbor of each row.
  #          nndiss: nearest neigbbor distance of each row.
  # FM, 2003/11/16
  
  nn <- rep(0, nrow(diss))
  nndiss <- rep(0.0, nrow(diss))
  MAXVAL <- 1.0e12
  if (nrow(diss) != ncol(diss)) stop("Invalid input first parameter.")
  if (nrow(diss) != length(flag)) stop("Invalid inputs 1st/2nd parameters.")
  # if (nrow(diss) != length(nn)) stop("Invalid inputs 1st/3rd parameters.")
  # if (nrow(diss) != length(nndiss)) stop("Invalid inputs 1st/4th parameters.")
  
  for (i1 in 1:nrow(diss)) {
    if (flag[i1] == 1) {
      minobs <- -1
      mindis <- MAXVAL
      for (i2 in 1:ncol(diss)) {
        if ( (diss[i1,i2] < mindis) && (i1 != i2) ) {
          mindis <- diss[i1,i2]
          minobs <- i2
        }
      }
      nn[i1] <- minobs
      nndiss[i1] <- mindis
    }
  }
  list(nn = nn, nndiss = nndiss)
}

# Internal Function (3): generate the clustering path. Not to be called by users.
hierclust <- function(a, wt) {
  
  MAXVAL <- 1.0e12
  
  n <- nrow(a)                               
  diss <- dissim(a, wt)                      # call to function dissim
  flag <- rep(1, n)                          # active/dead indicator
  a <- rep(0, n-1)                           # left subnode on clustering
  b <- rep(0, n-1)                           # right subnode on clustering
  ia <- rep(0, n-1)                          # R-compatible version of a
  ib <- rep(0, n-1)                          # R-compatible version of b
  lev <- rep(0, n-1)                         # level or criterion values
  card <- rep(1, n)                          # cardinalities
  mass <- wt
  order <- rep(0, n)                         # R-compatible order for plotting
  
  nnsnnsdiss <- getnns(diss, flag)           # call to function getnns
  clusmat <- matrix(0, n, n)                 # cluster memberships
  for (i in 1:n) clusmat[i,n] <- i           # init. trivial partition
  
  for (ncl in (n-1):1) {                      # main loop 
    # check for agglomerable pair
    minobs <- -1;  
    mindis <- MAXVAL;
    for (i in 1:n) {
      if (flag[i] == 1) {
        if (nnsnnsdiss$nndiss[i] < mindis) {
          mindis <- nnsnnsdiss$nndiss[i]
          minobs <- i
        }
      }
    }
    # find agglomerands clus1 and clus2, with former < latter
    if (minobs < nnsnnsdiss$nn[minobs]) {
      clus1 <- minobs
      clus2 <- nnsnnsdiss$nn[minobs]
    }
    if (minobs > nnsnnsdiss$nn[minobs]) {
      clus2 <- minobs
      clus1 <- nnsnnsdiss$nn[minobs]
    }
    # So, agglomeration of pair clus1 < clus2 defines cluster ncl
    
    #------------------------------------ Block for subnode labels 
    a[ncl] <- clus1                       # aine, or left child node
    b[ncl] <- clus2                       # benjamin, or right child node
    # Now build up ia, ib as version of a, b which is R-compliant
    if (card[clus1] == 1) ia[ncl] <- (-clus1)     # singleton
    if (card[clus2] == 1) ib[ncl] <- (-clus2)     # singleton
    if (card[clus1] > 1) {                # left child is non-singleton
      lastind <- 0
      for (i2 in (n-1):(ncl+1)) {        # Must have n-1 >= ncl+1 here
        if (a[i2] == clus1) lastind <- i2    # Only concerns a[i2]
      }
      ia[ncl] <- n - lastind             # label of non-singleton
    }
    if (card[clus2] > 1) {                # right child is non-singleton
      lastind <- 0
      for (i2 in (n-1):(ncl+1)) {        # Must have n-1 >= ncl+1 here
        if (a[i2] == clus2) lastind <- i2    # Can only concern a[i2]
      }
      ib[ncl] <- n - lastind             # label of non-singleton
    }
    if (ia[ncl] > 0 || ib[ncl] > 0) {     # Check that left < right
      left <- min(ia[ncl],ib[ncl])
      right <- max(ia[ncl],ib[ncl])
      ia[ncl] <- left                    # Just get left < right
      ib[ncl] <- right
    }
    #--------------------------------------------------------------------
    
    lev[ncl] <- mindis
    for (i in 1:n) {
      clusmat[i,ncl] <- clusmat[i,ncl+1]
      if (clusmat[i,ncl] == clus2) clusmat[i,ncl] <- clus1
    }
    # Next we need to update diss array
    for (i in 1:n) {
      if ( (i != clus1) && (i != clus2) && (flag[i] == 1) ) {
        diss[clus1,i] <- 
          ((mass[clus1]+mass[i])/(mass[clus1]+mass[clus2]+mass[i]))*diss[clus1,i] +
          ((mass[clus2]+mass[i])/(mass[clus1]+mass[clus2]+mass[i]))*diss[clus2,i] -
          (mass[i]/(mass[clus1]+mass[clus2]+mass[i]))*diss[clus1,clus2] 
        diss[i,clus1] <- diss[clus1,i]
      }
    }
    mass[clus1] <- mass[clus1] + mass[clus2]    # Update mass of new cluster
    card[clus1] <- card[clus1] + card[clus2]    # Update card of new cluster
    # Cluster label clus2 is knocked out; following not nec. but no harm
    flag[clus2] <- 0
    nnsnnsdiss$nndiss[clus2] <- MAXVAL
    mass[clus2] <- 0.0
    for (i in 1:n) {
      diss[clus2,i] <- MAXVAL
      diss[i,clus2] <- diss[clus2,i]
    }
    # Finally update nnsnnsdiss$nn and nnsnnsdiss$nndiss
    # i.e. nearest neighbors and the nearest neigh. dissimilarity
    nnsnnsdiss <- getnns(diss, flag)
  }
  
  temp <- cbind(a,b)
  merge2 <- temp[nrow(temp):1, ]
  temp <- cbind(ia,ib)
  merge <- temp[nrow(temp):1,]
  dimnames(merge) <- NULL
  # merge is R-compliant; later suppress merge2
  
  #-------------------------------- Build R-compatible order from ia, ib
  orderlist <- c(merge[n-1,1], merge[n-1,2])
  norderlist <- 2
  for (i in 1:(n-2)) {           # For precisely n-2 further node expansions
    for (i2 in 1:norderlist) {       # Scan orderlist
      if (orderlist[i2] > 0) {     # Non-singleton to be expanded
        tobeexp <- orderlist[i2]
        if (i2 == 1) {
          orderlist <- c(merge[tobeexp,1],merge[tobeexp,2],
                         orderlist[2:norderlist])
        }
        if (i2 == norderlist) {
          orderlist <- c(orderlist[1:(norderlist-1)],
                         merge[tobeexp,1],merge[tobeexp,2])
        }
        if (i2 > 1 && i2 < norderlist) {
          orderlist <- c(orderlist[1:(i2-1)], 
                         merge[tobeexp,1],merge[tobeexp,2],
                         orderlist[(i2+1):norderlist])
        }
        norderlist <- length(orderlist)
      }
    }
  }
  orderlist <- (-orderlist)
  class(orderlist) <- "integer"
  
  xcall <- "hierclust(a,wt)"
  class(xcall) <- "call"
  #clusmat=clusmat
  #labels=as.character(1:n)
  
  retlist <- list(merge=merge,height=as.single(lev[(n-1):1]),order=orderlist,
                  labels=dimnames(a)[[1]],method="minvar",call=xcall,
                  dist.method="euclidean-factor")
  retlist <- list(merge=merge,height=lev[(n-1):1],order=orderlist)
  class(retlist) <- "hclust"
  retlist
}
