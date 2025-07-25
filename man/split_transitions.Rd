% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_occupancies.R
\name{split_transitions}
\alias{split_transitions}
\title{Extract origin and destination states from transition names}
\usage{
split_transitions(p_list, delim = NULL)
}
\arguments{
\item{p_list}{A named list of transitions where names follow the \code{"from->to"} convention.}

\item{delim}{Character delimiter between \code{from} and \code{to} in transition names. If not specified we try to detect it from names of \code{p_list}}
}
\value{
A named list with two elements:
\describe{
\item{\code{from}}{Character vector of unique origin states}
\item{\code{to}}{Character vector of unique destination states}
}
}
\description{
Parses the names of transitions in \code{p_list} to identify the unique sets of origin (\code{from}) and destination (\code{to}) states.
}
