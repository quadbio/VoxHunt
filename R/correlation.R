#' Correlate single-cell gene expression with in situ hybridization data
#'
#' @rdname voxel_map
#' @export voxel_map
#'
voxel_map <- function(object, ...){
    UseMethod(generic = 'voxel_map', object = object)
}

#' Correlate single-cell gene expression with in situ hybridization data
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
    allow_neg = F,
    genes_use = NULL
){

    if (!exists('DATA_LIST') | !exists('PATH_LIST')){
        stop('Data has not been loaded. Please run load_aba_data() first.')
    }

    stage <- stage_name(stage)
    if (is.null(DATA_LIST[[stage]])){
        DATA_LIST[[stage]] <<- read_loom(PATH_LIST[[stage]])
    }

    voxel_mat <- DATA_LIST[[stage]]$matrix
    inter_genes <- intersect(colnames(object), colnames(voxel_mat))
    if (!is.null(genes_use)){
        inter_genes <- intersect(inter_genes, genes_use)
    }
    expr_mat <- as.matrix(t(object[, inter_genes]))
    voxel_mat[voxel_mat < 1] <- 0
    voxel_mat <- as.matrix(t(voxel_mat[, inter_genes]))

    vox_corr <- safe_cor(expr_mat, voxel_mat)

    return(vox_corr)
}

#' Correlate single-cell gene expression with in situ hybridization data
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
    allow_neg = F,
    genes_use = NULL
){
    expr_mat <- t(GetAssayData(object, slot = 'data'))
    if (is.null(group_name)){
        groups <- Idents(org_srt)
    } else {
        groups <- object[[group_name]][, 1]
    }
    vox_cor <- voxel_map(
        object = expr_mat,
        stage = stage,
        groups = groups,
        allow_neg = allow_neg,
        method = method,
        genes_use = genes_use
    )
}

#' Safe correlation function which returns a sparse matrix without missing values
#'
safe_cor <- function(
    x,
    y,
    method = 'pearson',
    allow_neg = F
){
    corr_mat <- stats::cor(x, y, method = method)
    corr_mat[is.na(corr_mat)] <- 0
    corr_mat <- Matrix::Matrix(corr_mat, sparse=T)
    if (!allow_neg){
        corr_mat[corr_mat < 0] <- 0
    }
    return(corr_mat)
}


