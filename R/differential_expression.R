#' Finds differentially expressed features for all structures on a given annotation level
#'
#' @export
#'
structure_markers <- function(expr_mat,
                                   annotation_level = NULL,
                                   structure_groups = NULL
){
    utils::data(voxel_meta, envir = environment())

    if (is.null(annotation_level)){
        groups <- structure_groups
    } else {
        meta <- filter(voxel_meta, voxel%in%rownames(expr_mat))
        groups <- meta[[annotation_level]]
    }

    expr_mat[expr_mat < 1] <- 0
    expr_mat <- log2(expr_mat + 1)
    de_df <- de(expr_mat, groups)
    de_df <- arrange(de_df, group, desc(auc))

    return(de_df)
}

#' Performs differential expression analysis on a sample x feature matrix
#'
de <- function(expr_mat, groups){
    de_df <- presto::wilcoxauc(t(expr_mat), groups)
    de_df <- dplyr::select(de_df, 'gene'=feature, group, 'avg_exp'=avgExpr, 'fc'=logFC, auc,
                    pval, padj, 'prcex_self'=pct_in, 'prcex_other'=pct_out)
    return(tibble::as_tibble(de_df))
}
