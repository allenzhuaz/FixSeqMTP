#' Bisection algorithm (FDR)
#'
#' Bisection algorithm to find the solution for the adjusted p-value for FDR controlling g-FSMTPs.
#'
#'@usage
#' bisection.FDR(f, a=0, b=1, p, k, j, n = 1000, tol)
#'@param f the objective function to be optimized for the solution.
#'@param a mininum of the interval which contains the solution from bisection algorithm.
#'@param b maxinum of the interval which contains the solution from bisection algorithm.
#'@param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#'@param k pre-specified number of acceptances allowed in the testing procedure (cannot exceed the length of \code{p})
#'@param j the index of the hypothesis.
#'@param n the number of sections that the interval which from bisection algorithm.
#'@param tol the desired accuracy.
#'@return a solution of the objective function which is between the interval from \code{a} to \code{b}.
#'@author Yalin Zhu
#'@seealso \code{\link{bisection.FWER}}
#'@export
bisection.FDR <- function(f, a=0, b=1, p, k, j, n = 1000, tol) {
  for (i in 1:n) {
    c <- (a + b) / 2 # Calculate midpoint

    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the
    # function and return the root.
    if ((f(c,p,k)[j] == 0) || ((b - a) / 2) < tol) {
      return(c)
    }

    # If another iteration is required,
    # check the signs of the function at the points c and a and reassign
    # a or b accordingly as the midpoint to be used in the next iteration.
    else{
      ifelse(sign(f(c,p,k)[j]) == sign(f(a,p,k)[j]),
             a <- c,
             b <- c)
    }
  }
  # If the max number of iterations is reached and no root has been found,
  # return message and end function.
  print('Too many iterations')
}

#' Objective function to be optimized (4)
#'
#' Objective function to be optimized for the adjusted p-values for FDR controlling g-FSMTP under arbitrary dependence. (See Theorem 3.1 and Theorem 4.1 in Lynch et al. (2016))
#'
#'@usage
#'  optim.arbidept.adjp(alpha, p, k)
#'@param alpha the parameter we need to solve for the adjusted p-values.
#'@param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#'@param k pre-specified number of acceptances allowed in the testing procedure (cannot exceed the length of \code{p})
#'@return difference between adjusted p-value and significant level alpha.
#'@author Yalin Zhu
#'@references
#'  Lynch, G., Guo, W., Sarkar, S. K., & Finner, H. (2016).
#'  The Control of the False Discovery Rate in Fixed Sequence Multiple Testing.
#'  \emph{arXiv preprint} arXiv:1611.03146.
#'@export
optim.arbidept.adjp <- function(alpha, p, k){
  m <- length(p)
  invp <- adjp <- numeric(m); s <- c()
  s[1] <- 0 # r  number of rejection, s number of acceptance
  adjp[1] <- invp[1] <- k*p[1]
  for ( i in 2:m ){
    s[i] <- sum(adjp[1:i-1] > alpha) # count number of acceptances for firtst (i-1) hypotheses.
    if (i <= k){
      invp[i] <- k*p[i]
    }
    else if(i>k){
      invp[i] <- (m-i+1)*k*p[i]/(m-k+1)
    }
    adjp[i] <- invp[i]*(s[i]<k)+1*(s[i]>=k)
  }
  return(adjp-alpha)
}

#' Objective function to be optimized (5)
#'
#' Objective function to be optimized for the adjusted p-values for FDR controlling g-FSMTP under independence. (See Theorem 3.2 and Theorem 4.2 in Lynch et al. (2016))
#'
#'@usage
#'  optim.indept.adjp(alpha, p, k)
#'@param alpha the parameter we need to solve for the adjusted p-values.
#'@param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#'@param k pre-specified number of acceptances allowed in the testing procedure (cannot exceed the length of \code{p})
#'@return difference between adjusted p-value and significant level alpha.
#'@author Yalin Zhu
#'@references
#'  Lynch, G., Guo, W., Sarkar, S. K., & Finner, H. (2016).
#'  The Control of the False Discovery Rate in Fixed Sequence Multiple Testing.
#'  \emph{arXiv preprint} arXiv:1611.03146.
#'@export
optim.indept.adjp <- function(alpha, p, k){
  m <- length(p)
  invp <- adjp <- numeric(m); s <- c()
  s[1] <- 0 # r  number of rejection, s number of acceptance
  adjp[1] <- invp[1] <- k*p[1]/(1-s[1]-(1-k)*p[1])
  for ( i in 2:m ){
    s[i] <- sum(adjp[1:i-1] > alpha) # count number of acceptances for firtst (i-1) hypotheses.
    invp[i] <- k*p[i]/(i-s[i]-(i-k)*p[i])
    adjp[i] <- invp[i]*(s[i]<k)+1*(s[i]>=k)
  }
  return(adjp-alpha)
}

#' Adjusted P-values for Fixed Sequence FDR Controlling Procedure under Arbitrary Dependence
#'
#'  Given a set of pre-ordered p-values and accuracy for the result, returns adjusted p-values using the generalized fixed sequence multiple testing procedures under arbitrary dependence (See Theorem 3.1 and 4.1 in Lynch et al. (2016)). The function also provides an option to make decisions given a pre-specified significant level \eqn{\alpha}.
#'
#'@usage
#'  FSFDR.arbidept.p.adjust(p, alpha=0.05, k=1, tol = 1e-6, make.decision = TRUE)
#'@param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#'@param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#'@param k pre-specified number of acceptances allowed in the testing procedure (cannot exceed the length of \code{p})
#'@param tol desired accuracy. The default value is \code{1e-6 }.
#'@param make.decision logical; if  \code{TRUE} (default), then the output include the decision rules compared adjusted p-values with significant level \eqn{alpha}
#'@return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}) if  \code{make.decision = FALSE}, or  a data frame including original p-values, adjusted p-values and decision rules if \code{make.decision = TRUE}.
#'@details
#' The generalized fixed sequence FDR controlling procedure stops on the \eqn{k}-th acceptances and automatically accepts the rest of hypotheses, where \eqn{k} is a pre-specified positive integer. When \eqn{k=1}, the generalized procedure becomes conventional one (Theorem 3.1 in Lynch et al. (2016)), which stops testing once one acceptance appears.
#'  This method strongly controls FDR under arbitrary dependence.
#'@seealso    \code{\link{FSFWER.arbidept.p.adjust}} for fixed sequence FWER controlling procedures.
#'@author Yalin Zhu
#'@references
#'  Lynch, G., Guo, W., Sarkar, S. K., & Finner, H. (2016).
#'  The Control of the False Discovery Rate in Fixed Sequence Multiple Testing.
#'  \emph{arXiv preprint} arXiv:1611.03146.
#'@examples
#' ## generate a pre-ordered pvalue vector for 50 hypotheses, where 80% are true nulls
#' set.seed(1234); m <- 50; pi0 <- 0.8; m0 <- m*pi0; m1 <- m-m0
#' mu <- c(4*0.9^(1:m1), rep(0,m0))
#' Zstat <- rnorm(n = m, mean = mu)
#' Pval <- 1-pnorm(Zstat)
#' ## conventional fixed sequence procedure
#' FSFDR.arbidept.p.adjust(p = Pval, alpha = 0.05)
#' ## generalized fixed sequence procedure allowing stop at 5th acceptance
#' FSFDR.arbidept.p.adjust(p = Pval, alpha = 0.05, k=5)
#'@export
FSFDR.arbidept.p.adjust <- function(p, alpha=0.05, k=1, tol = 1e-6, make.decision = TRUE){
  m <- length(p); opt.adjp <- numeric(m);
  if (alpha < 0 | alpha > 1) stop("significant level 'alpha' should be between 0 and 1")
  if (k > m | k < 1) stop("number of acceptances should be positive, and cannot exceed total number of hypotheses")
  for (j in 1:m){
    opt.adjp[j] <- bisection.FDR(optim.arbidept.adjp, p=p, k=k, j=j, tol=tol, a = 0, b = 1)
  }
  if (make.decision==TRUE){
    return(data.frame(raw.p = p, adjust.p=opt.adjp, decision=ifelse(opt.adjp<=alpha, "reject","accept")))
  } else{return(opt.adjp)}
}


#' Critical Values for Fixed Sequence FDR Controlling Procedure under Arbitrary Dependence
#'
#'  Given a set of pre-ordered p-values and accuracy for the result, return the corresponding critical values using the generalized fixed sequence FDR controlling procedure under arbitrary dependence (See Theorem 3.1 and 4.1 in Lynch et al. (2016)). The function also provides an option to make decisions given a pre-specified significant level \eqn{\alpha}.
#'
#'@usage
#'  FSFDR.arbidept.cv(p, k=1, alpha = 0.05, make.decision = TRUE)
#'@param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#'@param alpha significant level used to calculate the critical values to make decisions, the default value is 0.05.
#'@param k pre-specified number of acceptances allowed in the testing procedure (cannot exceed the length of \code{p})
#'@param make.decision logical; if  \code{TRUE} (default), then the output include the decision rules compared adjusted p-values with significant level \eqn{alpha}
#'@return
#' A numeric vector of the critical values (of the same length as \code{p}) if  \code{make.decision = FALSE}, or  a data frame including original p-values, critical values and decision rules if \code{make.decision = TRUE}.
#'@seealso  \code{\link{FSFWER.arbidept.cv}}  for fixed sequence FWER controlling procedures.
#'@author Yalin Zhu
#'@references
#'  Lynch, G., Guo, W., Sarkar, S. K., & Finner, H. (2016).
#'  The Control of the False Discovery Rate in Fixed Sequence Multiple Testing.
#'  \emph{arXiv preprint} arXiv:1611.03146.
#'@examples
#' ## generate a pre-ordered pvalue vector for 50 hypotheses, where 80% are true nulls
#' set.seed(1234); m <- 50; pi0 <- 0.8; m0 <- m*pi0; m1 <- m-m0
#' mu <- c(4*0.9^(1:m1), rep(0,m0))
#' Zstat <- rnorm(n = m, mean = mu)
#' Pval <- 1-pnorm(Zstat)
#' ## conventional fixed sequence procedure
#' FSFDR.arbidept.cv(p = Pval, alpha = 0.05)
#' ## generalized fixed sequence procedure allowing stop at 5th acceptance
#' FSFDR.arbidept.cv(p = Pval, alpha = 0.05, k=5)
#'@export
FSFDR.arbidept.cv <- function(p, k=1, alpha = 0.05, make.decision = TRUE){
  m <- length(p); cv <- numeric(m)
  if (alpha < 0 | alpha > 1) stop("significant level 'alpha' should be between 0 and 1")
  if (k > m | k < 1) stop("number of acceptances should be positive, and cannot exceed total number of hypotheses")
  for (i in 1:m){
    if (i <= k){cv[i] <- alpha/k}
    else{cv[i] <- min((m-k+1)/(m-i+1)*alpha/k,1)}
  }
  decision <- c()
  acc <- 0
  for ( i in 1:m ){
    decision[i] <- ifelse(p[i] <= cv[i], "reject", "accept")
    acc <- acc+(p[i] > cv[i])
    if (acc == k) break
  }
  if (i < m){
    decision[i:m] <- "accept"
  }
  if (make.decision==TRUE){
    return(data.frame(raw.p = p, critical.value=cv, decision))
  } else{return(cv)}
}

#' Adjusted P-values for Fixed Sequence FDR Controlling Procedure under Independence
#'
#'  Given a set of pre-ordered p-values and accuracy for the result, returns adjusted p-values using the generalized fixed sequence multiple testing procedures under independence for true nulls (See Theorem 3.2 and 4.2 in Lynch et al. (2016)). The function also provides an option to make decisions given a pre-specified significant level \eqn{\alpha}.
#'
#'@usage
#'  FSFDR.indept.p.adjust(p, alpha=0.05, k=1, tol = 1e-6, make.decision = TRUE)
#'@param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#'@param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#'@param k pre-specified number of acceptances allowed in the testing procedure (cannot exceed the length of \code{p})
#'@param tol desired accuracy. The default value is \code{1e-6 }.
#'@param make.decision logical; if  \code{TRUE} (default), then the output include the decision rules compared adjusted p-values with significant level \eqn{alpha}
#'@return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}) if  \code{make.decision = FALSE}, or  a data frame including original p-values, adjusted p-values and decision rules if \code{make.decision = TRUE}.
#'@details
#' The generalized fixed sequence FDR controlling procedure stops on the \eqn{k}-th acceptances and automatically accepts the rest of hypotheses, where \eqn{k} is a pre-specified positive integer. When \eqn{k=1}, the generalized procedure becomes conventional one (Theorem 3.2 in Lynch et al. (2016)), which stops testing once one acceptance appears.
#'  This method strongly controls FDR if the true null p-values are mutually independent and are independent of the false null p-values. When k=1, the conventional procedure strongly controls FDR if  the p-values are negatively associated on the true null p-values.
#'@seealso    \code{\link{FSFWER.arbidept.p.adjust}} for fixed sequence FWER controlling procedures.
#'@author Yalin Zhu
#'@references
#'  Lynch, G., Guo, W., Sarkar, S. K., & Finner, H. (2016).
#'  The Control of the False Discovery Rate in Fixed Sequence Multiple Testing.
#'  \emph{arXiv preprint} arXiv:1611.03146.
#'@examples
#' ## generate a pre-ordered pvalue vector for 50 hypotheses, where 80% are true nulls
#' set.seed(1234); m <- 50; pi0 <- 0.8; m0 <- m*pi0; m1 <- m-m0
#' mu <- c(4*0.9^(1:m1), rep(0,m0))
#' Zstat <- rnorm(n = m, mean = mu)
#' Pval <- 1-pnorm(Zstat)
#' ## conventional fixed sequence procedure
#' FSFDR.indept.p.adjust(p = Pval, alpha = 0.05)
#' ## generalized fixed sequence procedure allowing stop at 5th acceptance
#' FSFDR.indept.p.adjust(p = Pval, alpha = 0.05, k=5)
#'@export
FSFDR.indept.p.adjust <- function(p, alpha=0.05, k=1, tol = 1e-6, make.decision = TRUE){
  m <- length(p); opt.adjp <- numeric(m);
  if (alpha < 0 | alpha > 1) stop("significant level 'alpha' should be between 0 and 1")
  if (k > m | k < 1) stop("number of acceptances should be positive, and cannot exceed total number of hypotheses")
  for (j in 1:m){
    opt.adjp[j] <- bisection.FDR(optim.indept.adjp, p=p, k=k, j=j, tol=tol, a = 0, b = 1)
  }
  if (make.decision==TRUE){
    return(data.frame(raw.p = p, adjust.p=opt.adjp, decision=ifelse(opt.adjp<=alpha, "reject","accept")))
  } else{return(opt.adjp)}
}

#' Critical Values for Fixed Sequence FDR Controlling Procedure under Independence
#'
#'  Given a set of pre-ordered p-values and accuracy for the result, return the corresponding critical values using the generalized fixed sequence FDR controlling procedure under independence for true nulls (See Theorem 3.2 and 4.2 in Lynch et al. (2016)). The function also provides an option to make decisions given a pre-specified significant level \eqn{\alpha}.
#'
#'@usage
#'  FSFDR.indept.cv(p, k=1, alpha = 0.05, tol = 1e-6, make.decision = TRUE)
#'@param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#'@param alpha significant level used to calculate the critical values to make decisions, the default value is 0.05.
#'@param k pre-specified number of acceptances allowed in the testing procedure (cannot exceed the length of \code{p})
#'@param tol desired accuracy. The default value is \code{1e-6 }.
#'@param make.decision logical; if  \code{TRUE} (default), then the output include the decision rules compared adjusted p-values with significant level \eqn{alpha}
#'@return
#' A numeric vector of the critical values (of the same length as \code{p}) if  \code{make.decision = FALSE}, or  a data frame including original p-values, critical values and decision rules if \code{make.decision = TRUE}.
#'@seealso  \code{\link{FSFWER.arbidept.cv}}  for fixed sequence FWER controlling procedures.
#'@author Yalin Zhu
#'@references
#'  Lynch, G., Guo, W., Sarkar, S. K., & Finner, H. (2016).
#'  The Control of the False Discovery Rate in Fixed Sequence Multiple Testing.
#'  \emph{arXiv preprint} arXiv:1611.03146.
#'@examples
#' ## generate a pre-ordered pvalue vector for 50 hypotheses, where 80% are true nulls
#' set.seed(1234); m <- 50; pi0 <- 0.8; m0 <- m*pi0; m1 <- m-m0
#' mu <- c(4*0.9^(1:m1), rep(0,m0))
#' Zstat <- rnorm(n = m, mean = mu)
#' Pval <- 1-pnorm(Zstat)
#' ## conventional fixed sequence procedure
#' FSFDR.indept.cv(p = Pval, alpha = 0.05)
#' ## generalized fixed sequence procedure allowing stop at 5th acceptance
#' FSFDR.indept.cv(p = Pval, alpha = 0.05, k=5)
#'@export
FSFDR.indept.cv <- function(p, k=1, alpha = 0.05, tol = 1e-6, make.decision = TRUE){
  m <- length(p); cv <- numeric(m)
  if (alpha < 0 | alpha > 1) stop("significant level 'alpha' should be between 0 and 1")
  if (k > m | k < 1) stop("number of acceptances should be positive, and cannot exceed total number of hypotheses")
  m <- length(p)
  r <- 0
  for ( i in 1:m ){
    cv[i] <- (r+1)*alpha/(k+(i-k)*alpha)
    r <- r + (bisection.FDR(optim.indept.adjp, p=p, k=k, j=i, tol=tol, a = 0, b = 1) <= alpha)
  }

  decision <- c()
  acc <- 0
  for ( i in 1:m ){
    decision[i] <- ifelse(p[i] <= cv[i], "reject", "accept")
    acc <- acc+(p[i] > cv[i])
    if (acc == k) break
  }
  if (i < m){
    decision[i:m] <- "accept"
  }
  if (make.decision==TRUE){
    return(data.frame(raw.p = p, critical.value=cv, decision))
  } else{return(cv)}
}



