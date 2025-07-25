% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_occupancies.R
\name{fill_missing_transitions}
\alias{fill_missing_transitions}
\title{Fill in missing transitions with zero probabilities}
\usage{
fill_missing_transitions(p_list, from_states, to_states, delim = NULL)
}
\arguments{
\item{p_list}{A named list of numeric vectors. Each element corresponds to a transition
between states, e.g., \code{"P->W"}, \code{"W->D"}.}

\item{from_states}{Character vector of origin (transient) states.}

\item{to_states}{Character vector of destination states (can include absorbing).}

\item{delim}{Character string separating origin and destination states in the names of transitions. If not specified we try to detect it from names of \code{p_list}}
}
\value{
A completed \code{p_list} list with any missing transitions added as zero vectors.
}
\description{
Ensures that all possible transitions from a set of \code{from_states} to a set of \code{to_states} are present in the \code{p_list}, filling any missing transitions with zero probability vectors. This omits self-transitions.
}
