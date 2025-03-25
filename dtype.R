dtype <- function (dmetric = "euclidean", pw = 2) 
{
  dmetrics <- c("euclidean", "sqeuclidean", "manhattan", "minkowski", 
                "sqchord", "chebyshev", "pearsonchi", "neymanchi", "sqchi", 
                "divergence", "addsymchi", "prosymchi", "clark", "canberra", 
                "sorensen", "lorentzian", "cosine", "correlation")
  dmetric <- match.arg(dmetric, dmetrics)
  if (dmetric == "euclidean") 
    dtype <- 2
  else if (dmetric == "sqeuclidean") 
    dtype <- 2
  else if (dmetric == "manhattan") 
    dtype <- 2
  else if (dmetric == "minkowski") 
    dtype <- pw
  else if (dmetric == "chebyshev") 
    dtype <- 2
  else if (dmetric == "pearsonchi") 
    dtype <- 2
  else if (dmetric == "neymanchi") 
    dtype <- 2
  else if (dmetric == "sqchi") 
    dtype <- 2
  else if (dmetric == "addsymchi") 
    dtype <- 2
  else if (dmetric == "prosymchi") 
    dtype <- 2
  else if (dmetric == "divergence") 
    dtype <- 2
  else if (dmetric == "clark") 
    dtype <- 2
  else if (dmetric == "sqchord") 
    dtype <- 2
  else if (dmetric == "canberra") 
    dtype <- 1
  else if (dmetric == "sorensen") 
    dtype <- 1
  else if (dmetric == "lorentzian") 
    dtype <- 1
  else if (dmetric == "cosine") 
    dtype <- 2
  else if (dmetric == "correlation") 
    dtype <- 2
  return(dtype)
}