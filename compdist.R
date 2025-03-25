compdist <- function (a, b, dmetric = "euclidean", pw = 2) 
{
  if (missing(a)) 
    stop("Missing data arguments to compute the distances")
  if (missing(b)) 
    b <- a
  if (!is.numeric(a) || !is.numeric(b)) 
    stop("Input data arguments must be numeric to compute the distances")
  dmetrics <- c("euclidean", "sqeuclidean", "manhattan", "minkowski", 
                "chebyshev", "pearsonchi", "neymanchi", "sqchi", "divergence", 
                "addsymchi", "prosymchi", "clark", "canberra", "sorensen", 
                "lorentzian", "sqchord", "cosine", "correlation")
  dmetric <- match.arg(dmetric, dmetrics)
  if (dmetric == "euclidean") 
    distance <- sqrt(sum(t(a - b) * (a - b)))
  else if (dmetric == "sqeuclidean") 
    distance <- sum(t(a - b) * (a - b))
  else if (dmetric == "manhattan") 
    distance <- sum(abs(a - b))
  else if (dmetric == "minkowski") 
    distance <- sum(abs(a - b)^pw)^(1/pw)
  else if (dmetric == "chebyshev") 
    distance <- max(abs(a - b))
  else if (dmetric == "pearsonchi") 
    distance <- sum((t(a - b) * (a - b)/b))
  else if (dmetric == "neymanchi") 
    distance <- sum(t(a - b) * (a - b)/a)
  else if (dmetric == "sqchi") 
    distance <- sum(t(a - b) * (a - b)/(a + b))
  else if (dmetric == "addsymchi") 
    distance <- sum((t(a - b) * (a - b)) * ((a + b)/(a * 
                                                       b)))
  else if (dmetric == "prosymchi") 
    distance <- 2 * sum((t(a - b) * (a - b))/(a + b))
  else if (dmetric == "divergence") 
    distance <- 2 * sum((t(a - b) * (a - b))/(t(a + b) * 
                                                (a + b)))
  else if (dmetric == "clark") 
    distance <- sqrt(sum(abs((a - b)/(a + b))^2))
  else if (dmetric == "sqchord") 
    distance <- sum(sqrt(a) - sqrt(b))^2
  else if (dmetric == "canberra") 
    distance <- sum(abs(a - b)/(a + b))
  else if (dmetric == "sorensen") 
    distance <- sum(abs(a - b))/sum(a + b)
  else if (dmetric == "lorentzian") 
    distance <- sum(log(1 + abs(a - b), exp(1)))
  else if (dmetric == "cosine") 
    distance <- sum(a %*% b)/(sqrt(sum(a * a)) * sqrt(sum(b * 
                                                            b)))
  else if (dmetric == "correlation") 
    distance <- (1 - cor(a, b))/2
  return(distance)
}
