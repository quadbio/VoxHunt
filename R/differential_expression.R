#' Find differentially expressed features for all structures on a given annotation level
#'
#' Performs DE analysis between structure annotation levels and returns a tibble.
#'
#' @param stage The developmental stage from the ABA.
#' @param annotation_level The structure annotation level to perform DE an analysis between.
#' @param structure_groups A character or factor vector to provide custom grouping of voxels.
#'
#' @return A tibble with AUROC and p-values.
#'
#' @export
#'
structure_markers <- function(
    stage = 'E13',
    annotation_level = 'custom_3',
    structure_groups = NULL
){
    utils::data(voxel_meta, envir = environment())

    if (!exists('DATA_LIST') | !exists('PATH_LIST')){
        stop('Data has not been loaded. Please run load_aba_data() first.')
    }

    stage <- stage_name(stage)
    if (is.null(DATA_LIST[[stage]])){
        DATA_LIST[[stage]] <<- readRDS(PATH_LIST[[stage]])
    }

    voxel_mat <- DATA_LIST[[stage]]$matrix

    if (!is.null(structure_groups)){
        groups <- structure_groups
    } else {
        meta <- filter(voxel_meta, voxel%in%rownames(voxel_mat))
        groups <- meta[[annotation_level]]
    }

    voxel_mat[voxel_mat < 1] <- 0
    voxel_mat <- log2(voxel_mat + 1)
    de_df <- de(voxel_mat, groups)
    de_df <- arrange(de_df, group, desc(auc))

    return(de_df)
}

#' Performs differential expression analysis on a sample x feature matrix
#'
#' @return A tibble with p-values and AUROC values.
#'
de <- function(expr_mat, groups){
    de_df <- presto::wilcoxauc(t(expr_mat), groups)
    de_df <- dplyr::select(de_df, 'gene'=feature, group, 'avg_exp'=avgExpr, 'fc'=logFC, auc,
                    pval, padj, 'prcex_self'=pct_in, 'prcex_other'=pct_out)
    return(tibble::as_tibble(de_df))
}

#' Finds variable features for all structures on a given annotation level
#'
#' @param stage A string indicating the developmental stage in the ABA.
#' @param nfeatures The number of features to return.
#'
#' @return A tibble with variable features.
#'
#' @export
#'
variable_genes <- function(
    stage = 'E13',
    nfeatures = 2000
){

    if (!exists('DATA_LIST') | !exists('PATH_LIST')){
        stop('Data has not been loaded. Please run load_aba_data() first.')
    }

    stage <- stage_name(stage)
    if (is.null(DATA_LIST[[stage]])){
        DATA_LIST[[stage]] <<- readRDS(PATH_LIST[[stage]])
    }

    voxel_mat <- DATA_LIST[[stage]]$matrix

    var_feat <- Seurat::FindVariableFeatures(t(voxel_mat), verbose=F) %>%
        dplyr::as_tibble(rownames='gene') %>%
        dplyr::arrange(dplyr::desc(vst.variance.standardized)) %>%
        dplyr::filter(dplyr::row_number() <= nfeatures) %>%
        {colnames(.) <- stringr::str_replace_all(colnames(.), '\\.', '_'); .}

    return(var_feat)
}
