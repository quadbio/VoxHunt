#' Correlate single-cell gene expression with in situ hybridization data
#'
#' @export
#'
voxel_correlation <- function(expr_mat, voxel_mat,
                              allow_neg = F,
                              method = 'pearson',
                              genes_use = NULL
                              ){
    inter_genes <- intersect(colnames(expr_mat), colnames(voxel_mat))
    if (!is.null(genes_use)){
        inter_genes <- intersect(inter_genes, genes_use)
    }
    expr_mat <- as.matrix(t(expr_mat[, inter_genes]))
    voxel_mat <- as.matrix(t(voxel_mat[, inter_genes]))
    vox_cor <- stats::cor(expr_mat, voxel_mat, method = method)
    vox_cor[is.na(vox_cor)] <- 0
    if (!allow_neg){
        vox_cor[vox_cor < 0] <- 0
    }
    return(vox_cor)
}



