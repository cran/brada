\name{monitor}
\alias{monitor}
\title{
monitor
}
\description{
Monitors a running trial with a binary endpoint and calculates the predictive probability or predictive evidence that the trial will result in a success. Reports whether to stop early for futility or efficacy based on a vector of binary observations.
}
\usage{
monitor(brada_object, obs)
}
\arguments{
  \item{brada_object}{
    An object of class \code{brada}.
  }
  \item{obs}{
    A vector of binary observations, where 1 is a success (response) and 0 a failure (no response).
  }
}
\value{
  No return value, prints the result of the monitoring to the console.
}
\author{
  Riko Kelter
}
\examples{
\donttest{design = brada(Nmax = 40, batchsize = 5, nInit = 10, 
                   p_true = 0.2 , p0 = 0.2, p1 = 0.2, 
                   nsim = 100,
                   a0 = 1, b0 = 1, 
                   theta_T = 0.95, theta_L = 0.05, theta_U = 0.975, 
                   method = "PP",
                   cores = 2)
monitor(design, obs = c(0,1,1,0,0,1,0,1,1,1))}
}