% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FSFWER.R
\name{optim.reject.adjp}
\alias{optim.reject.adjp}
\title{Objective function to be optimized (1)}
\usage{
optim.reject.adjp(alpha, p, beta)
}
\arguments{
\item{alpha}{the parameter we need to solve for the adjusted p-values.}

\item{p}{numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.}

\item{beta}{pre-specified constant satisfying \eqn{0 \le \beta <1}, only for \code{method="accept"}.}
}
\value{
difference between adjusted p-value and significant level alpha.
}
\description{
Objective function to be optimized for the adjusted p-values for FWER controlling g-FSMTP based on the numbers of rejections only. (See Procedure A1 in Qiu et al. (2015))
}
\author{
Yalin Zhu
}
\references{
Qiu, Z., Guo, W., & Lynch, G. (2015).
  On generalized fixed sequence procedures for controlling the FWER.
  \emph{Statistics in medicine}, 34(30), 3968-3983.
}

