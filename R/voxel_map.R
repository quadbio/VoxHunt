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

#' @param stage The developmental stage in the ABA to map single cells to.
#' @param groups A character or factor vector or for grouping of cells,
#' e.g. clusters, cell types.
#' @param method A character string indicating which correlation coefficient to compute.
#' @param genes_use A character vector with genes to use for computing the correlation.
#' We recommend to use 150 - 500 genes.
#' @param allow_neg Logical. Whether to allow negative correlations or set them to 0.
#'
#' @return A VoxelMap object with a cell x voxel correlation matrix and metadata.
#'
#' @rdname voxel_map
#' @export
#' @method voxel_map default
#'
voxel_map.default <- function(
    object,
    stage = 'E13',
    groups = NULL,
    method = 'pearson',
    genes_use = NULL,
    allow_neg = F
){

    if (!exists('DATA_LIST') | !exists('PATH_LIST')){
        stop('Data has not been loaded. Please run load_aba_data() first.')
    }

    stage <- stage_name(stage)
    if (is.null(DATA_LIST[[stage]])){
        DATA_LIST[[stage]] <<- readRDS(PATH_LIST[[stage]])
    }

    voxel_mat <- DATA_LIST[[stage]]$matrix
    inter_genes <- intersect(colnames(object), colnames(voxel_mat))
    if (!is.null(genes_use)){
        inter_genes <- intersect(inter_genes, genes_use)
    }
    expr_mat <- t(object[, inter_genes])
    voxel_mat[voxel_mat < 1] <- 0
    voxel_mat <- t(voxel_mat[, inter_genes])

    corr_mat <- safe_cor(expr_mat, voxel_mat, method = method, allow_neg = allow_neg)
    if (is.null(groups)){
        cell_meta <- tibble(
            cell = rownames(corr_mat)
        )
    } else {
        cell_meta <- tibble(
            cell = rownames(corr_mat),
            group = groups
        )
    }
    utils::data(voxel_meta, envir = environment())
    voxel_meta = voxel_meta[match(colnames(corr_mat), voxel_meta$voxel), ]

    vox_map <- list(
        corr_mat = corr_mat,
        cell_meta = cell_meta,
        voxel_meta = voxel_meta,
        genes = inter_genes
    )
    class(vox_map) <- 'VoxelMap'

    return(vox_map)
}

#' @param group_name A string indicating the metadata column for grouping the cells,
#' e.g. clusters, cell types.
#'
#' @rdname voxel_map
#' @export
#' @method voxel_map Seurat
#'
voxel_map.Seurat <- function(
    object,
    stage = 'E13',
    group_name = NULL,
    method = 'pearson',
    genes_use = NULL,
    allow_neg = F
){
    expr_mat <- t(Seurat::GetAssayData(object, slot = 'data'))
    if (is.null(group_name)){
        groups <- Seurat::Idents(object)
    } else {
        groups <- object[[group_name]][, 1]
    }
    vox_cor <- voxel_map(
        object = Matrix::Matrix(expr_mat, sparse = T),
        stage = stage,
        groups = groups,
        allow_neg = allow_neg,
        method = method,
        genes_use = genes_use
    )
}


#' Print VoxelMap objects
#'
#' @rdname print
#' @export
#' @method print VoxelMap
#'
print.VoxelMap <- function(object){
    n_cells <- dim(object$corr_mat)[1]
    n_voxels <- dim(object$corr_mat)[2]
    n_genes <- length(object$genes)
    stage  <- unique(object$voxel_meta$stage)
    cat(paste0(
        'A VoxelMap object\n', n_cells, ' cells mapped to\n',
        n_voxels, ' voxels in the ', stage, ' mouse brain\nbased on ',
        n_genes, ' features\n'
    ))
}

#' Summarize VoxelMap object over groups
#'
#' @rdname summarize_groups
#' @export summarize_groups
#'
summarize_groups <- function(object, ...){
    UseMethod(generic = 'summarize_groups', object = object)
}

#' @import Matrix
#'
#' @param groups A metadata column or character vector to group the cells,
#' e.g. clusters, cell types.
#' @param fun Function used to aggregate the groups.
#'
#' @return A tibble with group summaries
#'
#' @rdname summarize_groups
#' @export
#' @method summarize_groups VoxelMap
#'
summarize_groups.VoxelMap <- function(
    object,
    groups = NULL,
    fun = colMeans
){

    if (is.null(groups) & 'group'%in%colnames(object$cell_meta)){
        groups <- object$cell_meta$group
    } else if (is.null(groups) & !'group'%in%colnames(object$cell_meta)){
        groups <- ' '
    }

    cluster_cor <- aggregate_matrix(object$corr_mat, groups=groups, fun=fun)

    plot_df <- cluster_cor %>%
        as.matrix() %>%
        tibble::as_tibble(rownames='voxel') %>%
        tidyr::gather(group, corr, -voxel) %>%
        dplyr::mutate(group=factor(group, levels=levels(factor(groups))))
    plot_df <- suppressMessages(dplyr::left_join(plot_df, object$voxel_meta))

    return(plot_df)
}


#' Assign cells to structure
#'
#' @rdname assign_cells
#' @export assign_cells
#'
assign_cells <- function(object, ...){
    UseMethod(generic = 'assign_cells', object = object)
}


#' @import Matrix
#'
#' @param object A VoxelMap opject
#'
#' @return A tibble with cells assigned to voxels
#'
#' @rdname assign_cells
#' @export
#' @method assign_cells VoxelMap
#'
assign_cells.VoxelMap <- function(object){

    which_max_corr <- colnames(object$corr_mat)[apply(object$corr_mat, 1, which.max)]
    max_corr <- apply(object$corr_mat, 1, max)
    max_corr_df <- tibble(
        voxel = which_max_corr,
        corr = max_corr,
        cluster = rownames(object$corr_mat)
    )

    max_corr_df <- suppressMessages(left_join(max_corr_df, object$voxel_meta))

    return(max_corr_df)
}


#' Summarize VoxelMap object over structures
#'
#' @rdname summarize_structures
#' @export summarize_structures
#'
summarize_structures <- function(object, ...){
    UseMethod(generic = 'summarize_structures', object = object)
}


#' @import Matrix
#'
#' @param annotation_level The structure annotation level to summarize voxels to.
#' @param fun Function to use for summarizing voxels.
#'
#' @return A tibble with structure summaries
#'
#' @rdname summarize_structures
#' @export
#' @method summarize_structures VoxelMap
#'
summarize_structures.VoxelMap <- function(
    object,
    annotation_level = 'custom_3',
    fun = colMeans
){

    corr_mat <- t(object$corr_mat)
    voxel_meta <- dplyr::group_by_(object$voxel_meta, annotation_level) %>%
        dplyr::filter(voxel%in%rownames(corr_mat)) %>%
        dplyr::filter(dplyr::n() > 5)
    cluster_cor <- aggregate_matrix(
        corr_mat[voxel_meta$voxel, ],
        groups = voxel_meta[[annotation_level]],
        fun = fun
    )
    plot_df <- cluster_cor %>%
        as.matrix() %>%
        tibble::as_tibble(rownames='cell') %>%
        tidyr::gather(struct, corr, -cell) %>%
        dplyr::mutate(struct=factor(struct, levels=levels(factor(voxel_meta[[annotation_level]]))))
    plot_df <- suppressWarnings(suppressMessages(
        dplyr::left_join(plot_df, object$voxel_meta, by=c(struct=annotation_level)
    )))

    return(plot_df)
}
