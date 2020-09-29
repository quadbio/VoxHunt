#' Project single cell gene expression data to spatial expression patterns in the mouse brain
#'
#' Calculates correlation between single-cell gene expression and
#' in situ hybridization data of the developing mouse brain from
#' the Allen Brain Atlas.
#'
#' @rdname voxel_map
#' @export voxel_map
#'
voxel_map <- function(object, ...){
    UseMethod(generic = 'voxel_map', object = object)
}

#' Project single cell gene expression data to BrainSpan
#'
#' Calculates correlation between single-cell gene expression and reference samples from
#' BrainSpan
#'
#' @rdname brainspan_map
#' @export brainspan_map
#'
brainspan_map <- function(object, ...){
    UseMethod(generic = 'brainspan_map', object = object)
}

#' Project single cell gene expression data to developing mouse brain data (LaManno & Siletti et al. 2020)
#'
#' Calculates correlation between single-cell gene expression and clusters from LaManno & Siletti et al. 2020
#'
#' @rdname mousebrain_map
#' @export mousebrain_map
#'
mousebrain_map <- function(object, ...){
    UseMethod(generic = 'mousebrain_map', object = object)
}

#' Summarize data over groups
#'
#' @rdname summarize_groups
#' @export summarize_groups
#'
summarize_groups <- function(object, ...){
    UseMethod(generic = 'summarize_groups', object = object)
}

#' Assign cells to structure
#'
#' @rdname assign_cells
#' @export assign_cells
#'
assign_cells <- function(object, ...){
    UseMethod(generic = 'assign_cells', object = object)
}

#' Summarize data over structures
#'
#' @rdname summarize_structures
#' @export summarize_structures
#'
summarize_structures <- function(object, ...){
    UseMethod(generic = 'summarize_structures', object = object)
}

#' Plot brain structure annotation
#'
#' Plot the voxel map of the developing mouse brain colorcoded by structure annotation.
#'
#' @rdname plot_annotation
#' @export plot_annotation
#'
plot_annotation <- function(object = 'E13', ...){
    UseMethod(generic = 'plot_annotation', object = object)
}

#' Plot similarity map of single cells
#'
#' @rdname plot_map
#' @export plot_map
#'
plot_map <- function(object, ...){
    UseMethod(generic = 'plot_map', object = object)
}

#' Plot similarity map of single cells to voxels in 3D
#'
#' @rdname plot_map_3d
#' @export plot_map_3d
#'
plot_map_3d <- function(object, ...){
    UseMethod(generic = 'plot_map_3d', object = object)
}

#' Plot similarity to brain structures
#'
#' @rdname plot_structure_similarity
#' @export plot_structure_similarity
#'
plot_structure_similarity <- function(object, ...){
    UseMethod(generic = 'plot_structure_similarity', object = object)
}

#' Devonvolute bulk RNA-seq data based on reference data
#'
#' @rdname deconvolute
#' @export deconvolute
#'
deconvolute <- function(object, ...){
    UseMethod(generic = 'deconvolute', object = object)
}








