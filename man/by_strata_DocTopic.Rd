% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis_covariates.R
\name{by_strata_DocTopic}
\alias{by_strata_DocTopic}
\title{Estimate document-topic distribution by strata (for covariate models)}
\usage{
by_strata_DocTopic(x, by_var, labels, by_values = NULL, ...)
}
\arguments{
\item{x}{the output from the covariate keyATM model (see \code{\link[=keyATM]{keyATM()}}).}

\item{by_var}{character. The name of the variable to use.}

\item{labels}{character. The labels for the values specified in \code{by_var} (ascending order).}

\item{by_values}{numeric. Specific values for \code{by_var}, ordered from small to large. If it is not specified, all values in \code{by_var} will be used.}

\item{...}{other arguments passed on to the \code{\link[=predict.keyATM_output]{predict.keyATM_output()}} function.}
}
\value{
strata_topicword object (a list).
}
\description{
Estimate document-topic distribution by strata (for covariate models)
}
