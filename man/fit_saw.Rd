% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post_saw.R
\name{fit_saw}
\alias{fit_saw}
\title{Post-SAW Estimation Procedure.}
\usage{
fit_saw(
  y,
  X,
  Z = NULL,
  time_effect = TRUE,
  id_effect = TRUE,
  s_thresh = "residual",
  return_info = FALSE
)
}
\arguments{
\item{y}{Matrix of labels. Has dimension T x N.}

\item{X}{List of feature matrices. The pth entry corresponds to the design
matrix of the pth covariate and has dimension T x N.}

\item{Z}{Instruments corresponding to the argument X. If NULL all X variables
are their own instrument.}

\item{time_effect}{Boolean indicating if a time effect is to be estimated.}

\item{id_effect}{Boolean indicating if an individual effect is to be
estimated and returned.}

\item{s_thresh}{Tuning parameter for the threshold lambda. Default to
"residual". Can be either numeric or "residual" or "smalln".}

\item{return_info}{Return additional info on model fit.}
}
\description{
Exports main function of the package, which performs the SAW procedure and
afterwards the Post-SAW estimation.
}
