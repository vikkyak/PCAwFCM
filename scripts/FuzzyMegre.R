FMerge <- function(dat, cl_id, cent, U, L){
  
  P <- mat.or.vec(L, L) # log values
  a <- combn(seq(1:L), 2)
  v <- c()
  S <- c()
  K <- mat.or.vec(dim(a)[1], dim(a)[2])
  for (l in 1:(0.5 * (L * L - L))) {
    #for (l in 1:21){
    ind1 <- a[1, l]
    ind2 <- a[2, l]
    x1 <- dat[which(cl_id == ind1),]
    cl <- paste("Cluster", ind1, sep = " ")
    u1 <- U[, cl][names(U[, cl]) %in% rownames(x1)]
    us1 <- sum(U[, cl][names(U[, cl]) %in% rownames(x1)])
    rm(cl)
    x2 <- dat[which(cl_id == ind2),]
    cl <- paste("Cluster", ind2, sep = " ")
    u2 <- U[, cl][names(U[, cl]) %in% rownames(x2)]
    us2 <- sum(U[, cl][names(U[, cl]) %in% rownames(x2)])
    rm(cl)
    # normalized membership degree corresponding to the clusters
    nu1 <- us1 / (us1 + us2)
    nu2 <- us2 / (us1 + us2)
    S  = c(u1,u2)
    K[,l] <- c(nu1, nu2)
  }
  to_merge <- c(a[which(K == max(K), arr.ind = T)], a[which(K == min(K), arr.ind = T)])
  Out <- list(to_merge, c(max(K), min(K)))
  Out
}
