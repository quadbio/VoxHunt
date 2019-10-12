#' Correlate single-cell gene expression with in situ hybridization data
#'
#' @rdname voxel_correlation
#' @export voxel_correlation
#'
voxel_correlation <- function(object, ...){
    UseMethod(generic = 'voxel_correlation', object = object)
}

#' Correlate single-cell gene expression with in situ hybridization data
#'
#' @rdname voxel_correlation
#' @export
#' @method voxel_correlation default
#'
voxel_correlation.default <- function(
    object,
    stage = 'E13',
    groups = NULL,
    method = 'pearson',
    allow_neg = F,
    genes_use = NULL
){
    inter_genes <- intersect(colnames(object), colnames(voxel_mat))
    if (!is.null(genes_use)){
        inter_genes <- intersect(inter_genes, genes_use)
    }
    expr_mat <- as.matrix(t(object[, inter_genes]))
    voxel_mat[voxel_mat < 1] <- 0
    voxel_mat <- as.matrix(t(voxel_mat[, inter_genes]))
    vox_cor <- stats::cor(expr_mat, voxel_mat, method = method)
    vox_cor[is.na(vox_cor)] <- 0
    vox_cor <- Matrix::Matrix(vox_cor, sparse=T)
    if (!allow_neg){
        vox_cor[vox_cor < 0] <- 0
    }
    return(vox_cor)
}

#' Correlate single-cell gene expression with in situ hybridization data
#'
#' @rdname voxel_correlation
#' @export
#' @method voxel_correlation Seurat
#'
voxel_correlation.Seurat <- function(
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
    vox_cor <- voxel_correlation(
        object = expr_mat,
        voxel_mat = voxel_mat,
        groups = groups,
        allow_neg = allow_neg,
        method = method,
        genes_use = genes_use
    )
}



