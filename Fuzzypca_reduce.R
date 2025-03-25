Fpca_reduce <- function (D_t, nbt, q)
{
  Y <- prep(D_t, scale = "none", center = TRUE)
  pca_out <- pca(Y,
                 method = "svd",
                 center = FALSE,
                 nPcs = q)
  x <- pca_out@scores

    Cl_history <- list()
    for (t in 1:nbt) {
      dat <- x
      q <- ncol(dat)
      K <- q + 1
      FCM <- fcmM(dat, centers = K)
      du <- rowSums(FCM$u[, duplicated(t(FCM$u))])
      ind <-max(which.max(table(attr(uniquecombs(t(FCM$u)), "index"))))
      FCM$u <- unique(FCM$u[, !duplicated(t(FCM$u))])
      col <- paste("Cluster", ind, sep = " ")
      FCM$u[, col] <- FCM$u[, col] + du
      FCM$v <- FCM$v[!duplicated(FCM$v), ]
      FCM$d <- unique(FCM$d[, !duplicated(t(FCM$d))])
      # If number of cluster is not equal to the dimension of FCM$u
      if (length(unique(FCM$cluster)) == dim(FCM$u)[2])   {
        colnames(FCM$u) <-
          paste("Cluster", 1:length(unique(FCM$cluster)), sep = " ")
        rownames(FCM$v) <-
          paste("Cluster", 1:length(unique(FCM$cluster)), sep = " ")
        aa <- sort(unique(FCM$cluster))
        for (i in 1:length(unique((FCM$cluster)))) {
          FCM$cluster[FCM$cluster == aa[i]] <- K + i
        }
        rm(aa)
        aa <- sort(unique(FCM$cluster))
        for (i in 1:length(unique((FCM$cluster)))) {
          FCM$cluster[FCM$cluster == aa[i]] <- i
        }
      } else{
        FCM$u <-
          FCM$u[, paste("Cluster", sort(unique(FCM$cluster)), sep = " ")]
        colnames(FCM$u) <-
          paste("Cluster", 1:length(unique(FCM$cluster)), sep = " ")
        FCM$v <-
          FCM$v[paste("Cluster", sort(unique(FCM$cluster)), sep = " "),]
        rownames(FCM$v) <-
          paste("Cluster", 1:length(unique(FCM$cluster)), sep = " ")
        aa <- sort(unique(FCM$cluster))
        for (i in 1:length(unique((FCM$cluster)))) {
          FCM$cluster[FCM$cluster == aa[i]] <- K + i
        }
        rm(aa)
        aa <- sort(unique(FCM$cluster))
        for (i in 1:length(unique((FCM$cluster)))) {
          FCM$cluster[FCM$cluster == aa[i]] <- i
        }
      }
      cl_id <- FCM$cluster
      cent <- FCM$v
      U <- FCM$u
      Cl_mat <- c(cl_id)
      L <- length(unique(cl_id))
      Q <- L - 2
      # Start auto encoder part
      for (i in 1:Q) {
        mrg <- FMerge(dat, cl_id, cent, U, L)
        cl_id[which(cl_id == mrg[[1]][2])] <- mrg[[1]][1]
        cent[mrg[[1]][1], ] <-
          mrg[[2]][1] * cent[mrg[[1]][1], ] + mrg[[2]][2] * cent[mrg[[1]][2], ]
        U[, mrg[[1]][1]] <-
          mrg[[2]][1] * U[, mrg[[1]][1]] + mrg[[2]][2] * U[, mrg[[1]][2]]
        cent <- cent[-mrg[[1]][2], -ncol(cent)]
        U <- U[,-mrg[[1]][2]]
        a <- unique(cl_id)
        Omega <- seq(1, max(a), 1)
        b <- setdiff(Omega, a)
        N <- length(b)
        if (N > 0) {
          for (ii in 1:N) {
            cl_id[cl_id > b[N + 1 - ii]] <- cl_id[cl_id > b[N + 1 - ii]] - 1
          }
        }
        colnames(U) <-
          paste("Cluster", 1:length(unique(cl_id)), sep = " ")
        L <-
          length(unique(cl_id))               # reduce number of cluster (K= K-1)
        dat <-
          dat[, -ncol(dat)]                 # reduce PCA dimension (q= q-1)
        Cl_mat <- cbind(Cl_mat, cl_id)
      }                                          # end auto encoder part
      Cl_history[[t]] <- Cl_mat
    }
    return(Cl_history)
}
