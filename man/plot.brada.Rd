\name{plot.brada}
\alias{plot.brada}
\alias{plot}
\title{
plot.brada
}
\description{
Plots the results of a Bayesian response-adaptive design analysis.
}
\usage{
\method{plot}{brada}(x, trajectories = 100, ...)
}
\arguments{
  \item{x}{
    An Object of class \code{"brada"}.
  }
  \item{trajectories}{
    Number of trajectories to be plotted. Defaults to 100.
  }
   \item{...}{
    Additional parameters, see \code{"plot(x, ...)"}.
  }
}
\value{
  Returns a plot.
}
\author{
  Riko Kelter
}
\examples{
\donttest{design = brada(Nmax = 40, batchsize = 5, nInit = 10, 
               p_true = 0.2 , p0 = 0.2, p1 = 0.2, 
               nsim = 100,
               a0 = 1, b0 = 1, 
               theta_T = 0.90, theta_L = 0.1, theta_U = 1, 
               method = "PP",
               cores = 2)
plot(design)}    
}
