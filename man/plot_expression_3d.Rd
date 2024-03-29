% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_expression_3d}
\alias{plot_expression_3d}
\title{Plot gene expression across the mouse brain in 3D}
\usage{
plot_expression_3d(
  stage,
  gene,
  annotation_level = NULL,
  annotation_colors = many,
  expression_colors = inferno,
  sizes = c(10, 1000),
  both_hemispheres = TRUE,
  ...
)
}
\arguments{
\item{stage}{The developmental stage to plot.}

\item{annotation_level}{The structure annotation level to color code.}

\item{annotation_colors}{Color map for structure annotation.}

\item{expression_colors}{Colors for expression scale.}

\item{sizes}{The size range for points.}

\item{both_hemispheres}{Logical. Whether to plot both hemispheres (TRUE) of only one (FALSE).}

\item{...}{Other arguments passed to plotly::plot_ly().}

\item{genes}{A character vector with the genes to plot.}

\item{view}{String indicating the perspective to show. Valid values are
'sagittal', 'coronal', 'traverse', 'z' , 'x', 'y', 'slice', '3D'.}
}
\value{
A gene expression plot in 3D.
}
\description{
Plot gene expression across the mouse brain in 3D
}
