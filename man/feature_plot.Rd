% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{feature_plot}
\alias{feature_plot}
\title{Feature plot}
\usage{
feature_plot(
  expr_mat,
  meta,
  markers,
  plot = xy_plot,
  sort = TRUE,
  scale = TRUE,
  ...
)
}
\arguments{
\item{expr_mat}{A sample x gene expression matrix.}

\item{meta}{Metadata with embedding coordinates.}

\item{markers}{A character vector with genes to display.}

\item{plot}{Function to use for generating individual plots.}

\item{sort}{Logical. Whether to sort highest expressing cells up.}

\item{scale}{Logical. Whether to scale the expression values.}

\item{...}{Other arguments passed to patchwork::wrap_plots().}
}
\value{
A panel of feature plots for each gene.
}
\description{
Feature plot
}
