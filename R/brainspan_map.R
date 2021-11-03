#' @param stage A integer vector indicating the stages (pcw) to map to..
#' @param groups A character or factor vector or for grouping of cells,
#' e.g. clusters, cell types.
#' @param method A character string indicating which correlation coefficient to compute.
#' @param genes_use A character vector with genes to use for computing the correlation.
#' We recommend to use 150 - 500 genes.
#' @param allow_neg Logical. Whether to allow negative correlations or set them to 0.
#' @param pseudobulk_groups Logical. Whether to summarizse the group expression before computing the correlation.
#'
#' @return A BrainSpanMap object with a cell x ref correlation matrix and metadata.
#'
#' @rdname brainspan_map
#' @export
#' @method brainspan_map default
#'
brainspan_map.default <- function(
    object,
    stages = c(8:16),
    groups = NULL,
    method = 'pearson',
    genes_use = NULL,
    allow_neg = FALSE,
    pseudobulk_groups = TRUE
){

    utils::data(brainspan, envir = environment())

    ref_meta <- filter(brainspan$row_meta, age_num%in%stages)
    ref_mat <- brainspan$matrix[ref_meta$ref, ]

    inter_genes <- intersect(colnames(object), colnames(ref_mat))
    if (!is.null(genes_use)){
        inter_genes <- intersect(inter_genes, genes_use)
    }
    expr_mat <- t(object[, inter_genes])
    ref_mat[ref_mat < 1] <- 0
    ref_mat <- t(ref_mat[, inter_genes])

    expr_mat <- object[, inter_genes]
    if (pseudobulk_groups){
        expr_mat <- aggregate_matrix(expr_mat, groups = groups, fun = Matrix::colMeans)
    } else {
        expr_mat <- t(expr_mat)
    }

    corr_mat <- safe_cor(expr_mat, ref_mat, method = method, allow_neg = allow_neg)
    if (is.null(groups)){
        cell_meta <- tibble(
            cell = rownames(corr_mat)
        )
    } else {
        if (pseudobulk_groups){
            groups <- factor(levels(factor(groups)), levels=levels(factor(groups)))
        }
        cell_meta <- tibble(
            cell = rownames(corr_mat),
            group = groups
        )
    }

    ref_map <- list(
        corr_mat = corr_mat,
        cell_meta = cell_meta,
        ref_meta = ref_meta,
        genes = inter_genes
    )
    class(ref_map) <- 'BrainSpanMap'
    return(ref_map)
}


#' @param group_name A string indicating the metadata column for grouping the cells,
#' e.g. clusters, cell types.
#'
#' @rdname brainspan_map
#' @export
#' @method brainspan_map Seurat
#'
brainspan_map.Seurat <- function(
    object,
    stages = c(8:16),
    group_name = NULL,
    method = 'pearson',
    genes_use = NULL,
    allow_neg = FALSE,
    pseudobulk_groups = TRUE
){
    expr_mat <- t(Seurat::GetAssayData(object, slot = 'data'))
    if (is.null(group_name)){
        groups <- Seurat::Idents(object)
    } else {
        groups <- object[[group_name]][, 1]
    }
    ref_cor <- brainspan_map(
        object = Matrix::Matrix(expr_mat, sparse = T),
        stages = stages,
        groups = groups,
        allow_neg = allow_neg,
        method = method,
        genes_use = genes_use,
        pseudobulk_groups = pseudobulk_groups
    )
    return(ref_cor)
}


#' Print BrainSpanMap objects
#'
#' @rdname print
#' @export
#' @method print BrainSpanMap
#'
print.BrainSpanMap <- function(object){
    n_cells <- dim(object$corr_mat)[1]
    n_ref <- dim(object$corr_mat)[2]
    n_genes <- length(object$genes)
    cat(paste0(
        'A BrainSpanMap object\n', n_cells, ' cells mapped to\n',
        n_ref, ' reference samples \nbased on ',
        n_genes, ' features\n'
    ))
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
#' @method summarize_groups BrainSpanMap
#'
summarize_groups.BrainSpanMap <- function(
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
        tibble::as_tibble(rownames='ref') %>%
        tidyr::gather(group, corr, -ref) %>%
        dplyr::mutate(group=factor(group, levels=levels(factor(groups))))
    plot_df <- suppressMessages(dplyr::left_join(plot_df, object$ref_meta))

    return(plot_df)
}


#' @import Matrix
#'
#' @param annotation_level The structure annotation level to summarize samples to.
#' @param fun Function to use for summarizing samples.
#'
#' @return A tibble with structure summaries
#'
#' @rdname summarize_structures
#' @export
#' @method summarize_structures BrainSpanMap
#'
summarize_structures.BrainSpanMap <- function(
    object,
    annotation_level = c('structure_group', 'structure_name', 'structure_acronym'),
    fun = colMeans
){
    annotation_level <- match.arg(annotation_level)
    if (annotation_level == 'structure_name'){
        annotation_level <- 'structure_acronym'
    }
    corr_mat <- t(object$corr_mat)
    ref_meta <- dplyr::group_by_at(object$ref_meta, annotation_level) %>%
        dplyr::filter(ref%in%rownames(corr_mat)) %>%
        dplyr::filter(dplyr::n() > 5)
    cluster_cor <- aggregate_matrix(
        corr_mat[ref_meta$ref, ],
        groups = ref_meta[[annotation_level]],
        fun = fun
    )
    plot_df <- cluster_cor %>%
        as.matrix() %>%
        tibble::as_tibble(rownames='cell') %>%
        tidyr::gather(struct, corr, -cell) %>%
        dplyr::mutate(struct=factor(struct, levels=levels(factor(ref_meta[[annotation_level]]))))
    plot_df <- suppressWarnings(suppressMessages(dplyr::left_join(plot_df, object$cell_meta)))

    return(plot_df)
}
