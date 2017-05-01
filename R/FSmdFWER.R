#' Critical values for Fixed Sequence mdFWER Controlling Procedure under Arbitrary Dependence Along with Directional Decisions Regarding Parameters of Interest
#'
#' Given a set of pre-ordered test statistics and the corresponding p-values, returns critical values using the directional fixed sequence multiple testing procedures under arbitrary dependence (See Procedure 1 and Theorem 1 in Grandhi et al. (2016)). The function also provides an option to make decisions and determine the sign given a pre-specified significant level \eqn{\alpha} and the test statistics.
#'
#' @usage
#' FSmdFWER.arbidept.cv(p, test.stat, alpha=0.05, make.decision = TRUE)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param test.stat numeric vector of test statistics, which are used to determine the direction of decisions, with the same length of \code{p}.
#' @param alpha significant level used to compare with Critical values to make decisions, the default value is 0.05.
#' @param make.decision logical; if \code{TRUE} (default), then the output include the decision rules compared original p-values with the critical values, and directions of the decision based on the sign of test statistics.
#' @return
#' A numeric vector of the critical values (of the same length as \code{p}) if  \code{make.decision = FALSEALSE}, or a data frame including original p-values, critical values, test statistics and directional decision rules if \code{make.decision = TRUE}.
#' @seealso \code{\link{FSmdFWER.indept.cv}} for fixed sequence mdFWER controlling procedures under independence.
#' @author Yalin Zhu
#' @references
#'  Grandhi, A., Guo, W., & Romano, J. P. (2016).
#'  Control of Directional Errors in Fixed Sequence Multiple Testing.
#'  \emph{arXiv preprint} arXiv:1602.02345.
#' @examples
#' ## Clinical trial example in Grandhi et al. (2016)
#' Pval <- c(0.0008, 0.0135, 0.0197, 0.7237, 0.0003, 0.2779, 0.0054, 0.8473)
#' Zstat <- c(3.4434, 2.5085, 2.3642, -0.3543, 3.7651, 1.0900, 2.8340, 0.1930)
#' FSmdFWER.arbidept.cv(p = Pval, test.stat = Zstat, make.decision = TRUE)
#' @export
FSmdFWER.arbidept.cv <- function(p, test.stat, alpha=0.05, make.decision = TRUE){
  m <- length(p); cv <- numeric(m)
  for (i in 1:m){
    cv[i]<- alpha/(2^(i-1))
  }
  # find the index of Stop
  if (all(p<=cv)){
    decision <- rep("reject", m)
  }
  else{
    st <- which.max(p > cv)
    decision <- c(rep("reject", st-1), rep("accept", m-st+1))
  }
  if (make.decision==TRUE){
    return(data.frame(raw.p = p, critical.value=cv, decision=decision,
                      test.stat = test.stat, direction = ifelse(decision == "reject", ifelse(test.stat > 0, "+", "-"), NA)))
  } else{return(cv)}
}


#' Adjusted P-values for Fixed Sequence mdFWER Controlling Procedure under Arbitrary Dependence Along with Directional Decisions Regarding Parameters of Interest
#'
#' Given a set of pre-ordered test statistics and the corresponding p-values, returns adjusted p-values using the directional fixed sequence multiple testing procedures under arbitrary dependence (See Procedure 1 and Theorem 1 in Grandhi et al. (2016)). The function also provides an option to make decisions and determine the sign given a pre-specified significant level \eqn{\alpha} and the test statistics.
#'
#' @usage
#' FSmdFWER.arbidept.p.adjust(p, test.stat, alpha=0.05, make.decision = TRUE)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param test.stat numeric vector of test statistics, which are used to determine the direction of decisions, with the same length of \code{p}.
#' @param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#' @param make.decision logical; if \code{TRUE} (default), then the output include the decision rules compared adjusted p-values with significant level \eqn{alpha}, and directions of the decision based on the sign of test statistics.
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}) if  \code{make.decision = FALSEALSE}, or a data frame including original p-values, adjusted p-values, test statistics and directional decision rules if \code{make.decision = TRUE}.
#' @seealso \code{\link{FSmdFWER.indept.p.adjust}} for fixed sequence mdFWER controlling procedures under independence.
#' @author Yalin Zhu
#' @references
#'  Grandhi, A., Guo, W., & Romano, J. P. (2016).
#'  Control of Directional Errors in Fixed Sequence Multiple Testing.
#'  \emph{arXiv preprint} arXiv:1602.02345.
#' @examples
#' ## Clinical trial example in Grandhi et al. (2016)
#' Pval <- c(0.0008, 0.0135, 0.0197, 0.7237, 0.0003, 0.2779, 0.0054, 0.8473)
#' Zstat <- c(3.4434, 2.5085, 2.3642, -0.3543, 3.7651, 1.0900, 2.8340, 0.1930)
#' FSmdFWER.arbidept.p.adjust(p = Pval, test.stat = Zstat, make.decision = TRUE)
#' @export
FSmdFWER.arbidept.p.adjust <- function(p, test.stat, alpha=0.05, make.decision = TRUE){
  m <- length(p); adjp <- numeric(m)
  for (i in 1:m){
    adjp[i]<- min(max(2^(i-1)*p[1:i]), 1)
  }
  if (make.decision==TRUE){
    return(data.frame(raw.p = p, adjust.p=adjp, decision=ifelse(adjp<=alpha, "reject","accept"),
                      test.stat = test.stat, direction = ifelse(adjp <= alpha, ifelse(test.stat > 0, "+", "-"), NA)))
  } else{return(adjp)}
}




#' Critical values for Fixed Sequence mdFWER Controlling Procedure under Independence Along with Directional Decisions Regarding Parameters of Interest
#'
#' Given a set of pre-ordered test statistics and the corresponding p-values, returns critical values using the directional fixed sequence multiple testing procedures under independence (See Procedure 2 and Theorem 2 in Grandhi et al. (2016)). The function also provides an option to make decisions and determine the sign given a pre-specified significant level \eqn{\alpha} and the test statistics.
#'
#' @usage
#' FSmdFWER.indept.cv(p, test.stat, alpha=0.05, make.decision = TRUE)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param test.stat numeric vector of test statistics, which are used to determine the direction of decisions, with the same length of \code{p}.
#' @param alpha significant level used to compare with Critical values to make decisions, the default value is 0.05.
#' @param make.decision logical; if \code{TRUE} (default), then the output include the decision rules compared original p-values with the critical values, and directions of the decision based on the sign of test statistics.
#' @return
#' A numeric vector of the critical values (of the same length as \code{p}) if  \code{make.decision = FALSEALSE}, or a data frame including original p-values, critical values, test statistics and directional decision rules if \code{make.decision = TRUE}.
#' @seealso \code{\link{FSmdFWER.arbidept.cv}} for fixed sequence mdFWER controlling procedures under arbitrary dependence.
#' @author Yalin Zhu
#' @references
#'  Grandhi, A., Guo, W., & Romano, J. P. (2016).
#'  Control of Directional Errors in Fixed Sequence Multiple Testing.
#'  \emph{arXiv preprint} arXiv:1602.02345.
#' @examples
#' ## Clinical trial example in Grandhi et al. (2016)
#' Pval <- c(0.0008, 0.0135, 0.0197, 0.7237, 0.0003, 0.2779, 0.0054, 0.8473)
#' Zstat <- c(3.4434, 2.5085, 2.3642, -0.3543, 3.7651, 1.0900, 2.8340, 0.1930)
#' FSmdFWER.indept.cv(p = Pval, test.stat = Zstat, make.decision = TRUE)
#' @export
FSmdFWER.indept.cv <- function(p, test.stat, alpha=0.05, make.decision = TRUE){
  m <- length(p); cv <- numeric(m)
  for (i in 1:m){
    cv[i]<- alpha
  }
  # find the index of Stop
  if (all(p<=cv)){
    decision <- rep("reject", m)
  }
  else{
    st <- which.max(p > cv)
    decision <- c(rep("reject", st-1), rep("accept", m-st+1))
  }
  if (make.decision==TRUE){
    return(data.frame(raw.p = p, critical.value=cv, decision=decision,
                      test.stat=test.stat, direction = ifelse(decision == "reject", ifelse(test.stat > 0, "+", "-"), NA)))
  } else{return(cv)}
}

#' Adjusted P-values for Fixed Sequence mdFWER Controlling Procedure under Independence Along with Directional Decisions Regarding Parameters of Interest
#'
#' Given a set of pre-ordered test statistics and the corresponding p-values, returns adjusted p-values using the directional fixed sequence multiple testing procedures under independence (See Procedure 2 and Theorem 2 in Grandhi et al. (2016)). The function also provides an option to make decisions and determine the sign given a pre-specified significant level \eqn{\alpha} and the test statistics.
#'
#' @usage
#' FSmdFWER.indept.p.adjust(p, test.stat, alpha=0.05, make.decision = TRUE)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param test.stat numeric vector of test statistics, which are used to determine the direction of decisions, with the same length of \code{p}.
#' @param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#' @param make.decision logical; if \code{TRUE} (default), then the output include the decision rules compared adjusted p-values with significant level \eqn{alpha}, and directions of the decision based on the sign of test statistics.
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}) if  \code{make.decision = FALSEALSE}, or a data frame including original p-values, adjusted p-values, test statistics and directional decision rules if \code{make.decision = TRUE}.
#' @seealso \code{\link{FSmdFWER.indept.p.adjust}} for fixed sequence mdFWER controlling procedures under independence.
#' @author Yalin Zhu
#' @references
#'  Grandhi, A., Guo, W., & Romano, J. P. (2016).
#'  Control of Directional Errors in Fixed Sequence Multiple Testing.
#'  \emph{arXiv preprint} arXiv:1602.02345.
#' @examples
#' ## Clinical trial example in Grandhi et al. (2015)
#' Pval <- c(0.0008, 0.0135, 0.0197, 0.7237, 0.0003, 0.2779, 0.0054, 0.8473)
#' Zstat <- c(3.4434, 2.5085, 2.3642, -0.3543, 3.7651, 1.0900, 2.8340, 0.1930)
#' FSmdFWER.indept.p.adjust(p = Pval, test.stat = Zstat, make.decision = TRUE)
#' @export
FSmdFWER.indept.p.adjust <- function(p, test.stat, alpha=0.05, make.decision = TRUE){
  m <- length(p); adjp <- numeric(m)
  for (i in 1:m){
    adjp[i]<- max(p[1:i])
  }
  if (make.decision==TRUE){
    return(data.frame(raw.p = p, adjust.p=adjp, decision=ifelse(adjp <= alpha, "reject","accept"),
                      test.stat = test.stat, direction = ifelse(adjp <= alpha, ifelse(test.stat > 0, "+", "-"), NA)))
  } else{return(adjp)}
}

