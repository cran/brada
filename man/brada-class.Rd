\name{brada-class}
\Rdversion{1.0}
\docType{class}
\alias{brada-class}
\title{Class \code{"brada-class"}}
\description{Class for modelling the results of a Bayesian response-adaptive design analysis}
\section{Objects from the Class}{Store the results of a Bayesian response-adaptive design analysis}
\section{Slots}{
\describe{\item{\code{data}:}{Object of class \code{"list"} holding the results of the Bayesian response-adaptive design analysis. \code{a0} and \code{b0} store the beta prior shape parameters, \code{Nmax} and \code{batchsize} store the maximum trial size and the batchsize used for interim analyses. \code{nInit} is the minimum sample size at which the first interim analysis is conducted. \code{p_true} is the true response probability used for simulation, \code{p0} is the right boundary of the null hypothesis and \code{p1} the left boundary of the alternative hypothesis. ...}}}
\keyword{classes}

