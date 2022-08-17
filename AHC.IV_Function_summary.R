

### Agglomerative Hierarchical Clustering method for detecting IV clusters and IV estimation using summary statistics.

### AHC.IV
### Function: Agglomerative Hierarchical Clustering method for detecting the IV clusters and IV estimation using summary data.

### Input: (1) betaY: A numeric vector of SNP-outcome associations (p_z by 1 vector).
###        (2) betaX: A numeric vector of SNP-exposure associations (p_z by 1 vector).
###        (3) seY: A numeric vector of the standarnd errors of the SNP-outcome associations (p_z by 1 vector).
###        (4) seX: A numeric vector of the standarnd errors of the SNP-exposure associations (p_z by 1 vector).
###        (5) alpha: A numeric scalar between 0 and 1 specifying the significance level for the t test (default = 0.05).
###        (6) tuning: A numeric scalar specifiying the threshold p-value for the Q test (default = 0.1/log(length(betaY))).
###        (7) weak: Logical. If weak = TRUE, the weak IV robust Q test is used (default = TRUE).

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

AHC.IV <- function(betaY, betaX, seY, seX, alpha = 0.05, tuning = 0.1/log(length(betaY)), weak = TRUE){

  # Define Constants
  pz <- length(betaY);
  P <- 1; ## this is for single outcome--single exposure case
  
  # F statistic
  F_bar = mean(betaX^2/seX^2)
  
  # Clustering
  comb <- combn(pz, P, repeats.allowed=F); Len <- ncol(comb);
  
  bvec <- matrix(betaY/betaX, nrow = pz);
  sevec <- matrix(seY/betaX, nrow = pz);
  
  dist <- dist(bvec, method = "euclidean");
  dend <- hclust(dist, method = "ward.D2");
  
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
  
  # This is the downward testing step (algorithm 2)
  Qsq.all <- Qsq(betaX, betaY, seX, seY);
  Qsq.all.p <- pchisq(Qsq.all[[1]], (pz - P), lower.tail = FALSE);
  
  if (Qsq.all.p > tuning){
    # all the IVs are selected as in the same cluster
    clusterIV <- list(c(1:pz));
    Q.IV <- list(Qsq.all.p);
    clusternames <- paste0('cluster', 1);
    beta.IV <- Qsq.all[[2]];
    se.IV <- Qsq.all[[3]];
    I.IV <- Qsq.all[[4]];
  }else{
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
  
  null_index <- which(AHC_results[,"t.IV"] < abs(qnorm(alpha/2))) # Clusters that do not give significant estimates are classified as the null cluster.
  if (length(null_index) != 0){AHC_cluster_null <- sort(unlist(AHC_cluster[c(null_index)]))}else{AHC_cluster_null <- NULL}
  
  Nr.AHC <- length(AHC_cluster)
  
  # Report results
  results = list( Cluster_number = Nr.AHC,
                  AHC_cluster = clusterIV,
                  Null_cluster = AHC_cluster_null,
                  Junk_cluster = AHC_cluster_junk,
                  F = F_bar,
                  AHC_results = AHC_results)
  return(results)
  
}



######## Running example: Scenario 4 of the simulations in the MR-Clust paper ############
######## Running example: Scenario 4 of the simulations in the MR-Clust paper ############
######## Running example: Scenario 4 of the simulations in the MR-Clust paper ############

L = 90  # number of IVs
K = 5  # number of IV clusters (including the junk cluster)
N = 5000  # sample size
beta = c(0.4, -0.4, 0.8, 0) #  true causal effects of the non-junk clusters
tau = 2

SNPs_group = mySNPs = matrix(c(rep(1, 10), rep(0, L-10),
                               rep(0,10), rep(1, 20), rep(0, L-30),
                               rep(0, 10), rep(0, 20), rep(1, 40), rep(0, L-70),
                               rep(0, 10), rep(0, 20), rep(0, 40), rep(1,10), rep(0, L-80),
                               rep(0, 10), rep(0, 20), rep(0, 40), rep(0,10), rep(1,10)), 
                               nrow = K, byrow = T)

set.seed(14609)

# generate betaX and seX
effect_x_mean = rnorm(L, 0, 1)
MAF = runif(L, 0.05, 0.5)
betaX = rnorm(L, effect_x_mean, 1/(N*MAF*(1-MAF)))
seX = sqrt(1/(N*MAF*(1-MAF)))

# generate betaY
betaY = rep(NA, L)
for (k in 1:(K-1)){
  SNPs_k = which(SNPs_group[k,]==1)
  betaY[SNPs_k] = rnorm(length(SNPs_k), beta[k]*betaX[SNPs_k], tau/(N*MAF[SNPs_k]*(1-MAF[SNPs_k])))
}

# generate betaY for the junk cluster
SNP_junk = which(SNPs_group[K,] == 1)
junk_mean = rnorm(10, 0, 1)
betaY[SNP_junk] = rnorm(length(SNP_junk), junk_mean*betaX[SNP_junk], tau/(N*MAF[SNP_junk]*(1-MAF[SNP_junk])))

# generate seY
seY = sqrt(tau/(N*MAF*(1-MAF)))

# run AHC
AHC_results <- AHC.IV(betaY, betaX, seY, seX)

