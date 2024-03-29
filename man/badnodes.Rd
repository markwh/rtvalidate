% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate.R
\name{badnodes}
\alias{badnodes}
\title{Returns a vector of worst-performing nodes (by error)}
\usage{
badnodes(valdata, variable = "width", n = 4, which = c("abs", "min",
  "max"), standardize = FALSE)
}
\arguments{
\item{valdata}{As returned by \code{rt_valdata()}}

\item{variable}{Which variable's errors define "bad" nodes?}

\item{n}{Number of bad nodes to return}

\item{which}{"abs" for worst absolute errors, "min" for worst
negative errors, "max" for worst positive errors}

\item{standardize}{scale by estimated uncertainty? Default is \code{FALSE}.}
}
\description{
Returns a vector of worst-performing nodes (by error)
}
