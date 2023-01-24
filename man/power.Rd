\name{power}
\alias{power}
\title{
power
}
\description{
Performs a power analysis for a \code{brada} object.
}
\usage{
power(brada_object, p_true, nsim=100, cores=2)
}
\arguments{
  \item{brada_object}{
    An object of class \code{brada}
  }
  \item{p_true}{
    the true response probability used for the power analysis
  }
  \item{nsim}{
    the number of Monte Carlo simulation, defaults to 100.
  }
  \item{cores}{
    CPU cores used for computation. Defaults to 2.
  }
}
\value{
  Returns an object of class \code{brada}.
}
\author{
  Riko Kelter
}
\examples{
\donttest{design = brada(Nmax = 30, batchsize = 5, nInit = 10, 
               p_true = 0.2 , p0 = 0.2, p1 = 0.2, 
               nsim = 1000,
               a0 = 1, b0 = 1, 
               theta_T = 0.90, theta_L = 0.1, theta_U = 1, 
               method = "PP",
               cores = 1)
    design_power = power(design, p_true = 0.4, nsim = 1000)
    plot(design_power)}
}