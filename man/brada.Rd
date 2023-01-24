\name{brada}
\alias{brada}
\title{
brada
}
\description{
Performs a Bayesian response-adaptive design analysis for trials with a binary endpoint.
}
\usage{
brada(a0=1,b0=1,Nmax=40,batchsize=5,nInit,p_true,p0,p1,
theta_T=0.90,theta_L=0.1,theta_U=1,nsim=100,
seed=42,method="PP",refFunc="flat",nu=0,
shape1=1,shape2=1,truncation=1,cores=2)
}
\arguments{
  \item{a0}{
    shape1 parameter of the beta prior.
  }
  \item{b0}{
    shape2 parameter of the beta prior.
  }
  \item{Nmax}{
    Maximum trial size.
  }
  \item{batchsize}{
    sample size after which an interim analysis is performed.
  }
  \item{nInit}{
    Initial sample size at which the first interim analysis is performed.
  }
  \item{p_true}{
    True binary response probability used for simulation.
  }
  \item{p0}{
    Right boundary of the null hypothesis to be tested.
  }
  \item{p1}{
    Left boundary of the alternative hypothesis to be tested.
  }
  \item{theta_T}{
    Threshold used in the designs for including trajectories as evidential.
  }
  \item{theta_L}{
    Stopping threshold for futility.
  }
  \item{theta_U}{
    Stopping threshold for efficacy.
  }
  \item{nsim}{
    Number of Monte Carlo iterations.
  }
  \item{seed}{
    Random number generator seed.
  }
  \item{cores}{
    Number of CPU cores to be used for computation. Defaults to 2, but 4 or larger is recommended.
  }
  \item{method}{
    Can be either "PP" or "PPe", depending on whether the predictive probability approach or the predictive evidence value design is desired. Note that the former is a special case of the latter.
  }
  \item{refFunc}{
    A string, either "flat", "beta", "binaryStep", "relu", "palu" or "lolu". See vignettes for explanation.
  }
  \item{nu}{
    A numeric value larger or equal to zero, indicating which evidence threshold if used in the predictive evidence value design.
  }
  \item{shape1}{
    shape1 parameter of the beta reference function, if used.
  }
  \item{shape2}{
    shape2 parameter of the beta reference function, if used.
  }
  \item{truncation}{
    Truncation point in case an artificial neural network reference function is used.
  }
}
\value{
  Returns an object of class brada.
}
\author{
  Riko Kelter
}
\examples{
\donttest{pp_design = brada(Nmax = 30, batchsize = 5, nInit = 10, 
               p_true = 0.2 , p0 = 0.2, p1 = 0.2, 
               nsim = 10,
               a0 = 1, b0 = 1, 
               theta_T = 0.90, theta_L = 0.1, theta_U = 1, 
               method = "PP",
               cores = 2)
summary(pp_design)}
}