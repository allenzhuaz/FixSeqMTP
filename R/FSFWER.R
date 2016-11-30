#' Bisection algorithm (FWER)
#'
#' Bisection algorithm to find the solution for the adjusted p-value for FWER controlling g-FSMTPs.
#'
#'@usage
#' bisection.FWER(f, a=0, b=1, p, beta, j, n = 1000, tol)
#'@param f the objective function which to be optimized for the solution.
#'@param a mininum of the interval which cantains the solution from bisection algorithm.
#'@param b maxinum of the interval which cantains the solution from bisection algorithm.
#'@param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#'@param beta  pre-specified constant satisfying \eqn{0 \le \beta <1}, only for \code{method="accept"}.
#'@param j index of the hypothesis.
#'@param n  number of sections that the interval which from bisection algorithm.
#'@param tol desired accuracy.
#'@return a solution of the objective function which between the interwal from \code{a} to \code{b}.
#'@author Yalin Zhu
#'@seealso \code{\link{bisection.FDR}}
#'@export
bisection.FWER <- function(f, a=0, b=1, p, beta, j, n = 1000, tol) {

  for (i in 1:n) {
    c <- (a + b) / 2 # Calculate midpoint

    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the
    # function and return the root.
    if ((f(c, p, beta)[j] == 0) || ((b - a) / 2) < tol) {
      return(c)
    }

    # If another iteration is required,
    # check the signs of the function at the points c and a and reassign
    # a or b accordingly as the midpoint to be used in the next iteration.
    else{
    ifelse(sign(f(c, p, beta)[j]) == sign(f(a, p, beta)[j]),
           a <- c,
           b <- c)
    }
  }
  # If the max number of iterations is reached and no root has been found,
  # return message and end function.
  print('Too many iterations')
}

#' Objective function to be optimized (1)
#'
#' Objective function to be optimized for the adjusted p-values for FWER controlling g-FSMTP based on the numbers of rejections only. (See Procedure A1 in Qiu et al. (2015))
#'
#'@usage
#'  optim.reject.adjp(alpha, p, beta)
#'@param alpha the parameter we need to solve for the adjusted p-values.
#'@param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#'@param beta  pre-specified constant satisfying \eqn{0 \le \beta <1}, only for \code{method="accept"}.
#'@return difference between adjusted p-value and significant level alpha.
#'@author Yalin Zhu
#'@references
#'   Qiu, Z., Guo, W., & Lynch, G. (2015).
#'   On generalized fixed sequence procedures for controlling the FWER.
#'   \emph{Statistics in medicine}, 34(30), 3968-3983.
#'@export
optim.reject.adjp <- function(alpha, p, beta){
  m <- length(p)
 adjp <- numeric(m); s <- c()
  s[1] <- 0 # r  number of rejection, s number of acceptance
  adjp[1] <- (m+s[1])*p[1]
  for ( i in 2:m ){
    s[i] <- sum(adjp[1:i-1] > alpha) # count number of acceptances for firtst (i-1) hypotheses.
    # rejections r[i]<-i-1-s[i]
    adjp[i] <- (m-(i-1-s[i]))*p[i]
  }
  return(adjp-alpha)
}

#' Objective function to be optimized (2)
#'
#' Objective function to be optimized for the adjusted p-values for FWER controlling g-FSMTP based on the numbers of acceptances only. (See Procedure A2 in Qiu et al. (2015))
#'
#'@usage
#'  optim.accept.adjp(alpha, p, beta)
#'@param alpha the parameter we need to solve for the adjusted p-values.
#'@param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#'@param beta  pre-specified constant satisfying \eqn{0 \le \beta <1}, only for \code{method="accept"}.
#'@return difference between adjusted p-value and significant level alpha.
#'@author Yalin Zhu
#'@references
#'   Qiu, Z., Guo, W., & Lynch, G. (2015).
#'   On generalized fixed sequence procedures for controlling the FWER.
#'   \emph{Statistics in medicine}, 34(30), 3968-3983.
#'@export
optim.accept.adjp <- function(alpha, p, beta){
  m <- length(p)
 adjp <- numeric(m); s <- c()
  s[1] <- 0 # r  number of rejection, s number of acceptance
  adjp[1] <- (1-beta^m)/((1-beta)*beta^s[1])*p[1]
  for ( i in 2:m ){
    s[i] <- sum(adjp[1:i-1] > alpha) # count number of acceptances for firtst (i-1) hypotheses.
    # rejections r[i]<-i-1-s[i]
    adjp[i] <- (1-beta^m)/((1-beta)*beta^s[i])*p[i]
  }
  return(adjp-alpha)
  }

#' Objective function to be optimized (3)
#'
#' Objective function to be optimized for the adjusted p-values for FWER controlling g-FSMTP based on the numbers of both rejections and acceptances. (See Procedure A3 in Qiu et al. (2015))
#'
#'@usage
#'  optim.both.adjp(alpha, p, beta)
#'@param alpha the parameter we need to solve for the adjusted p-values.
#'@param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#'@param beta  pre-specified constant satisfying \eqn{0 \le \beta <1}, only for \code{method="accept"}.
#'@return difference between adjusted p-value and significant level alpha.
#'@author Yalin Zhu
#'@references
#'   Qiu, Z., Guo, W., & Lynch, G. (2015).
#'   On generalized fixed sequence procedures for controlling the FWER.
#'   \emph{Statistics in medicine}, 34(30), 3968-3983.
#'@export
optim.both.adjp <- function(alpha, p, beta){
  m <- length(p)
  adjp <- numeric(m); s <- c()
  s[1] <- 0 # r  number of rejection, s number of acceptance
  adjp[1] <- p[1]/(1/(m+s[1])+(m+s[1]-1-2*s[1])/m^2)
  for ( i in 2:m ){
    s[i] <- sum(adjp[1:i-1] > alpha) # count number of acceptances for firtst (i-1) hypotheses.
    # rejections r[i]<-i-1-s[i]
    adjp[i] <- p[i]/(1/(m-(i-1-s[i]))+(m-(i-1-s[i])-1-2*s[i])/m^2)
  }
  return(adjp-alpha)
  }

#' Adjusted P-values for Fixed Sequence FWER Controlling Procedures under Arbitrary Dependence
#'
#'  Given a set of pre-ordered p-values and accuracy for the result, returns adjusted p-values using one of several generalized fixed sequence multiple testing procedures. The function also provides an option to make decisions given a pre-specified significant level \eqn{\alpha}.
#'
#'@usage
#'  FSFWER.arbidept.p.adjust(p, alpha=0.05, beta=0.5, tol = 1e-6,
#'   method = c("reject","accept","both"), make.decision = TRUE)
#'@param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#'@param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#'@param beta  pre-specified constant satisfying \eqn{0 \le \beta <1}, only for \code{method="accept"}. The default value is \code{0.5}.
#'@param tol desired accuracy. The default value is \code{1e-6 }.
#'@param method adjustment method. See details.
#'@param make.decision logical; if  \code{TRUE} (default), then the output include the decision rules compared adjusted p-values with significant level \eqn{alpha}
#'@return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}) if  \code{make.decision = FALSE}, or  a data frame including original p-values, adjusted p-values and decision rules if \code{make.decision = TRUE}.
#'@details
#' The adjustment methods for Fixed Sequence multiple testing include the Procedure A1 only using numbers of rejections ("reject"), Procedure A2 only using numbers of acceptances ("accept") and Procedure A3 using both numbers of rejections and numbers of acceptances ("both").
#'  The three methods strongly control FWER under arbitrary dependence.
#'  The constant \code{beta} need to be specified only for the Procedure A2 (\code{"accept"}),  one can ignore this argument when using other methods.
#'@seealso \code{\link{FSFDR.arbidept.p.adjust}} and   \code{\link{FSFDR.arbidept.p.adjust}} for fixed sequence FDR controlling procedures.
#'@author Yalin Zhu
#'@references
#'   Qiu, Z., Guo, W., & Lynch, G. (2015).
#'   On generalized fixed sequence procedures for controlling the FWER.
#'   \emph{Statistics in medicine}, 34(30), 3968-3983.
#'@examples
#'   ## Clinical trial example in Qiu et al. (2015)
#' Pval <- c(0.0008, 0.0135, 0.0197, 0.7237, 0.0003, 0.2779, 0.0054, 0.8473)
#' FSFWER.arbidept.p.adjust(p=Pval, alpha=0.05, method = "reject")
#' FSFWER.arbidept.p.adjust(p=Pval, alpha=0.05, beta=0.1, method = "accept")
#' FSFWER.arbidept.p.adjust(p=Pval, alpha=0.05, beta=0.5, method = "accept")
#' FSFWER.arbidept.p.adjust(p=Pval, alpha=0.05, beta=0.9, method = "accept")
#' FSFWER.arbidept.p.adjust(p=Pval, alpha=0.05, method = "both")
#'@export
FSFWER.arbidept.p.adjust <- function(p, alpha=0.05, beta=0.5, tol = 1e-6, method = c("reject","accept","both"), make.decision = TRUE)
  {
    m <- length(p); opt.adjp <- numeric(m);
     method <- match.arg(method)
     if (beta < 0 | beta >= 1) stop("pre-sepecified constant 'beta' should be between 0 and 1")
     if (alpha < 0 | alpha > 1) stop("significant level 'alpha' should be between 0 and 1")
  for (j in 1:m){
    opt.adjp[j] <- switch(method,
                          reject = bisection.FWER(optim.reject.adjp, p=p, j=j, beta=beta, tol=tol, a = 0, b = 1),
                          accept = bisection.FWER(optim.accept.adjp, p=p, j=j, beta=beta, tol=tol, a = 0, b = 1),
                          both = bisection.FWER(optim.both.adjp, p=p, j=j, beta=beta, tol=tol, a = 0, b = 1))
  }
  if (make.decision==TRUE){
    return(data.frame(p.value = p, adjust.p.value=opt.adjp, decision=ifelse(opt.adjp<=alpha, "reject","accept")))
  } else{return(opt.adjp)}
}

#' Critical Values for Fixed Sequence FWER Controlling Procedures under Arbitrary Dependence
#'
#'  Given a set of pre-ordered p-values and accuracy for the result, returns critical values using one of several generalized fixed sequence multiple testing procedures. The function also provides an option to make decisions given a pre-specified significant level \eqn{\alpha}.
#'
#'@usage
#'  FSFWER.arbidept.cv(p, alpha=0.05, beta=0.5, tol = 1e-6,
#'   method = c("reject","accept","both"), make.decision = TRUE)
#'@param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#'@param alpha significant level used to calculate the critical values to make decisions, the default value is 0.05.
#'@param beta  pre-specified constant satisfying \eqn{0 \le \beta <1}, only for \code{method="accept"}. The default value is \code{0.5}.
#'@param tol desired accuracy. The default value is \code{1e-6 }.
#'@param method critical value calculation method. See details.
#'@param make.decision logical; if  \code{TRUE} (default), then the output include the decision rules compared adjusted p-values with significant level \eqn{alpha}
#'@return
#' A numeric vector of the critical values (of the same length as \code{p}) if  \code{make.decision = FALSE}, or  a data frame including original p-values, critical values and decision rules if \code{make.decision = TRUE}.
#'@details
#' The critical value calculation methods for Fixed Sequence multiple testing include the Procedure A1 only using numbers of rejections ("reject"), Procedure A2 only using numbers of acceptances ("accept") and Procedure A3 using both numbers of rejections and numbers of acceptances ("both").
#'  The three methods strongly control FWER under arbitrary dependence.
#'  The constant \code{beta} need to be specified only for the Procedure A2 (\code{"accept"}),  one can ignore this argument when using other methods.
#'@seealso \code{\link{FSFDR.arbidept.cv}} and   \code{\link{FSFDR.indept.cv}} for fixed sequence FDR controlling procedures.
#'@author Yalin Zhu
#'@references
#'   Qiu, Z., Guo, W., & Lynch, G. (2015).
#'   On generalized fixed sequence procedures for controlling the FWER.
#'   \emph{Statistics in medicine}, 34(30), 3968-3983.
#'@examples
#'   ## Clinical trial example in Qiu et al. (2015)
#' Pval <- c(0.0008, 0.0135, 0.0197, 0.7237, 0.0003, 0.2779, 0.0054, 0.8473)
#' FSFWER.arbidept.cv(p=Pval, alpha=0.05, method = "reject")
#' FSFWER.arbidept.cv(p=Pval, alpha=0.05, beta=0.1, method = "accept")
#' FSFWER.arbidept.cv(p=Pval, alpha=0.05, beta=0.5, method = "accept")
#' FSFWER.arbidept.cv(p=Pval, alpha=0.05, beta=0.9, method = "accept")
#' FSFWER.arbidept.cv(p=Pval, alpha=0.05, method = "both")
#'@export
FSFWER.arbidept.cv <- function(p, alpha=0.05, beta=0.5, tol = 1e-6, method = c("reject","accept","both"), make.decision = TRUE){
    if (beta < 0 | beta >= 1) stop("pre-sepecified constant 'beta' should be between 0 and 1")
    if (alpha < 0 | alpha > 1) stop("significant level 'alpha' should be between 0 and 1")
    m <- length(p); cv <- numeric(m)
    method <- match.arg(method)
  r <- s <- 0 # r  number of rejection, s number of acceptance
  for (i in 1:m){
  cv[i] <- switch(method,
                  reject = alpha/(m-r),
                  accept = (1-beta)*beta^s*alpha/(1-beta^m),
                  both = (1/(m-r)+(m-r-1)/m^2-2*s/m^2)*alpha)
  if(p[i] <= cv[i]){r <- r+1}
  else{s <- s+1}
  }
  if (make.decision==TRUE){
    return(data.frame(p.value = p, critical.value=cv, decision=ifelse(p<=cv, "reject","accept")))
  } else{return(cv)}
}

