scaled.covariance.matrix <- function (alpha, parsed.tree) {
  T <- parsed.tree$T;
  bt <- parsed.tree$branch.times;			 
  if (alpha == 0) {
    V <- bt;
  } else {
    a <- 2*alpha;
    V <- exp(-a*T) * expm1(a*bt) / a;
  }
}

