% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GridOperations.R
\name{gridSearchEvaluation}
\alias{gridSearchEvaluation}
\title{Run Grid Search}
\usage{
gridSearchEvaluation(
  opt_object,
  output_filename = NULL,
  n_cores = 1,
  redo_run = F
)
}
\arguments{
\item{opt_object}{An optimization task object that includes a valid simulation grid. The grid should be created using \code{drawSimulationGrid()}.}

\item{output_filename}{Optional character string specifying the name of the CSV file to save results. Default is \code{"GridSearchResult.csv"}.}

\item{n_cores}{Integer. Number of CPU cores to use for parallel simulation. If \code{n_cores > 1}, uses \code{mclapply()} for parallel execution. Default is 1 (sequential).}

\item{redo_run}{Logical. If TRUE, re-evaluates all grid points even if results already exist. Default is FALSE.}
}
\value{
The updated \code{opt_object} with updated \code{grid_run_status}. The results are appended to a CSV file in the working directory.
}
\description{
Runs PKPD simulations over a grid of treatment parameter combinations stored in the optimization object.
Function \code{drawSimulationGrid} need to first be run before running this function to generate the grid table.
The results are evaluated using the objective function, and penalties are written to a result file.
Previously evaluated points are skipped unless \code{redo_run = TRUE}.
}
