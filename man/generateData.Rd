\name{generateData}
\alias{generateData}
\title{
generateData
}
\description{
Generates a matrix of trial data.
}
\usage{
generateData(p,Nmax,nsim,seed=420)
}
\arguments{
  \item{p}{
    true response probability
  }
  \item{Nmax}{
    Maximum trial size.
  }
  \item{nsim}{
    Number of Monte Carlo iterations.
  }
  \item{seed}{
    Random number generator seed.
  }
}
\value{
  Returns a matrix with simulated trial data.
}
\author{
Riko Kelter
}
\examples{
generateData(p=0.2,Nmax=40,nsim=100,seed=420)
}