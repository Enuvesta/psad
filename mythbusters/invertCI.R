is.between <- function(x, a, b) {
  return(x > a & x < b)
}

isCovered <- function(CIfun, CL, theta0){
  CI <- CIfun(CL)
  return(is.between(theta0, CI[1], CI[2]))
}

invertCI <- function(CIfun, precision, theta0){
  i <- 1
  start <- .Machine$double.eps
  end   <- 1 - .Machine$double.eps

  while (i<=precision){
    candidates  <- seq(from=start, to=end, length.out=11)
    currentCL   <- candidates[head(which(sapply(X = candidates, FUN = function(x) isCovered(CIfun, x, theta0))), 1)]
    if (length(currentCL)==0){
      currentCL <- 1-.Machine$double.eps
      break
    }else{
      if (currentCL == start){
        end   <- currentCL + 10^(-i)
      }else{
        start <- currentCL - 10^(-i)
        end   <- currentCL
      }
      i <- i+1
    }
  }

  p <- 1-currentCL
  return(p)
}
