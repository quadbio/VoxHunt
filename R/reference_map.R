#' @param object A matrix with query expression data.
#' @param reference_object A matrix with reference expression data.
#' @param reference_meta A data frame with reference metadata.
#' @param groups A character or factor vector or for grouping of cells,
#' e.g. clusters, cell types.
#' @param method A character string indicating which correlation coefficient to compute.
#' @param genes_use A character vector with genes to use for computing the correlation.
#' We recommend to use 150 - 500 genes.
#' @param allow_neg Logical. Whether to allow negative correlations or set them to 0.
#' @param pseudobulk_groups Logical. Whether to summarize the group expression before computing the correlation.
#'
#' @return A ReferenceMap object with a cell x cell correlation matrix and metadata.
#'
#' @rdname reference_map
#' @export
#' @method reference_map default
#'
reference_map.default <- function(
    object,
    reference_object,
    query_meta = NULL,
    reference_meta = NULL,
    groups = NULL,
    reduction = NULL,
    method = 'pearson',
    genes_use = NULL,
    allow_neg = FALSE,
    pseudobulk_groups = TRUE
){

    inter_genes <- intersect(colnames(object), colnames(reference_object))
    if (!is.null(genes_use)){
        inter_genes <- intersect(inter_genes, genes_use)
    }
    expr_mat <- t(object[, inter_genes])
    ref_mat <- t(reference_object[, inter_genes])

    expr_mat <- object[, inter_genes]
    if (pseudobulk_groups){
        expr_mat <- aggregate_matrix(expr_mat, groups = groups, fun = Matrix::colMeans)
    } else {
        expr_mat <- t(expr_mat)
    }

    corr_mat <- safe_cor(expr_mat, ref_mat, method = method, allow_neg = allow_neg)

    if (is.null(query_meta) | pseudobulk_groups){
        query_meta <- tibble::tibble(
            cell = rownames(corr_mat)
        )
    }

    if (is.null(reference_meta)){
        reference_meta <- tibble::tibble(
            cell = colnames(corr_mat)
        )
    }

    if (!is.null(groups)){
        if (pseudobulk_groups){
            groups <- factor(levels(factor(groups)), levels=levels(factor(groups)))
        }
        query_meta$group <- groups
    }

    ref_map <- list(
        corr_mat = corr_mat,
        query_meta = query_meta,
        ref_meta = reference_meta,
        genes = inter_genes,
        reduction = reduction,
        single_cell = !pseudobulk_groups
    )
    class(ref_map) <- 'ReferenceMap'

    return(ref_map)
}


#' @param group_name A string indicating the metadata column for grouping the cells,
#' e.g. clusters, cell types.
#'
#' @rdname reference_map
#' @export
#' @method reference_map Seurat
#'
reference_map.Seurat <- function(
    object,
    reference_object,
    slot = 'data',
    assay = 'RNA',
    group_name = NULL,
    reduction = 'umap',
    method = 'pearson',
    genes_use = NULL,
    allow_neg = FALSE,
    pseudobulk_groups = TRUE
){
    expr_mat <- t(SeuratObject::LayerData(object, layer = slot, assay = assay))
    query_meta <- as_tibble(object@meta.data, rownames='cell')
    ref_mat <- t(SeuratObject::LayerData(reference_object, layer = slot, assay = assay))
    ref_meta <- as_tibble(reference_object@meta.data, rownames='cell')
    reduction <- as_tibble(reference_object[['umap']]@cell.embeddings, rownames='cell')[,1:3]
    colnames(reduction)[2:3] <- c('x', 'y')

    if (is.null(group_name)){
        groups <- Seurat::Idents(object)
    } else {
        groups <- object[[group_name]][, 1]
    }

    ref_cor <- reference_map(
        object = Matrix::Matrix(expr_mat, sparse = T),
        reference_object = Matrix::Matrix(ref_mat, sparse = T),
        reference_meta = ref_meta,
        query_meta = query_meta,
        groups = groups,
        reduction = reduction,
        allow_neg = allow_neg,
        method = method,
        genes_use = genes_use,
        pseudobulk_groups = pseudobulk_groups
    )

    return(ref_cor)
}


#' Print ReferenceMap objects
#'
#' @rdname print
#' @export
#' @method print ReferenceMap
#'
print.ReferenceMap <- function(object){
    n_cells <- nrow(object$corr_mat)
    n_ref <- ncol(object$corr_mat)
    n_genes <- length(object$genes)
    query_str <- ifelse(object$single_cell, 'cells', 'groups')
    cat(paste0(
        'A ReferenceMap object\n', n_cells, ' query ',
        query_str, ' mapped to\n',
        n_ref, ' reference cells \nbased on ',
        n_genes, ' features\n'
    ))
}


#' @import Matrix
#'
#' @param groups A metadata column or character vector to group the cells,
#' e.g. clusters, cell types.
#' @param summarize A character vector indicating whether to summarize the
#' query (`'query'`) or the reference (`'reference'`) cells.
#' @param fun Function used to aggregate the groups.
#'
#' @return A tibble with group summaries
#'
#' @rdname summarize_groups
#' @export
#' @method summarize_groups ReferenceMap
#'
summarize_groups.ReferenceMap <- function(
    object,
    summarize = 'query',
    groups = NULL,
    fun = colMeans
){
    if (summarize == 'query'){
        corr_mat <- object$corr_mat
        meta <- object$query_meta
    } else if (summarize == 'reference'){
        corr_mat <- t(object$corr_mat)
        meta <- object$ref_meta
    }

    if (is.null(groups)){
        group_vec <- meta$group
    } else if (length(groups)==ncol(corr_mat)){
        group_vec <- groups
    } else {
        group_vec <- meta[[groups]]
    }

    if (object$single_cell & summarize=='query') {
        cnames <- c('query_group', 'ref_group')
    } else if (object$single_cell & summarize=='reference') {
        cnames <- c('ref_group', 'cell')
    } else if (!object$single_cell & summarize=='query') {
        cnames <- c('query_group', 'ref_group')
    } else if (!object$single_cell & summarize=='reference') {
        cnames <- c('ref_group', 'query_group')
    }

    cluster_cor <- t(aggregate_matrix(corr_mat, groups=group_vec, fun=fun))

    plot_df <- cluster_cor %>%
        as.matrix() %>%
        tibble::as_tibble(rownames=cnames[1]) %>%
        tidyr::pivot_longer(!c(cnames[1]), names_to=cnames[2], values_to='corr')

    return(plot_df)
}

