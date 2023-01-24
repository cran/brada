\name{calibrate}
\alias{calibrate}
\title{
calibrate
}
\description{
Calibrates a brada object to achieve specified false-positive and false-negative rates.
}
\usage{
calibrate(brada_object, nsim = 100, cores = 2, seq, 
alpha=NULL, beta=NULL, calibration = "nu")
}
\arguments{
  \item{brada_object}{
    An object of class \code{brada}
  }
  \item{nsim}{
    Number of Monte Carlo iterations
  }
  \item{cores}{
    Number of cores used for computation
  }
  \item{seq}{
    Sequence of values for the evidence threshold \code{nu}, or sequence of the futility threshold \code{theta_L}, depending on which value is passed to the function in the \code{calibration} argument. For example, \code{seq=seq(0,1,0.1)} in combination with \code{calibration="nu"} and \code{alpha=0.1} implies that the function tries to calibrate the false-positive rate to decrease below \code{alpha=0.1} by iteratively increasing \code{nu} from zero to one in steps of size 0.1.
  }
  \item{alpha}{
    Upper bound for false-positive rate. Note that it is only possible to specify either \code{alpha} or \code{beta}. When \code{alpha} is passed as an argument, \code{calibration} should take the value \code{"nu"}, while when \code{beta} is passed as an argument, \code{calibration} should take the value \code{"theta_L"}.
  }
  \item{beta}{
    Upper bound for false-negative rate
  }
  \item{calibration}{
    String which specifies which parameter to calibrate. Can take the values \code{"nu"} or \code{"theta_L"}.
  }
}
\value{
  Prints the output to the console and returns the false-positive rate or false-negative rate of the calibrated design, depending on which value the \code{calibration} argument takes.
}
\author{
Riko Kelter
}
