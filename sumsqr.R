sumsqr <- function (x, v, clusters) 
{
  sumsqr <- function(x) sum(scale(x, scale = FALSE)^2)
  bwss <- sumsqr(v[clusters, ])
  wss <- sapply(split(as.data.frame(x), clusters), sumsqr)
  twss <- sum(wss)
  tss <- bwss + twss
  ss <- list(bwss, wss, twss, tss)
  names(ss) <- c("between.ss", "within.ss", "tot.within.ss", 
                 "tot.ss")
  return(ss)
}
