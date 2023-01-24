\name{summary.brada}
\alias{summary.brada}
\alias{summary}
\title{
summary.brada
}
\description{
Prints the results of a Bayesian response-adaptive design analysis.
}
\usage{
\method{summary}{brada}(object, ...)
}
\arguments{
  \item{object}{
    An Object of class \code{"brada"}.
  }
  \item{...}{
    Additional parameters, see \code{"summary(object, ...)"}.
  }
}
\details{
Summarises the results of a Bayesian response-adaptive design analysis.
}
\value{
  Prints the results onto the console.
}
\author{
  Riko Kelter
}
\examples{
\donttest{pp_design = brada(Nmax = 40, batchsize = 5, nInit = 10, 
               p_true = 0.2 , p0 = 0.2, p1 = 0.2, 
               nsim = 100,
               a0 = 1, b0 = 1, 
               theta_T = 0.90, theta_L = 0.1, theta_U = 1, 
               method = "PP",
               cores = 2)
summary(pp_design)}
}