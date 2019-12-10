#' Finds differentially expressed features for all structures on a given annotation level
#'
#' @export
#'
structure_markers <- function(
    stage = 'E13pt5',
    annotation_level = 'custom_3',
    structure_groups = NULL
){
    utils::data(voxel_meta, envir = environment())

    if (!exists('DATA_LIST') | !exists('PATH_LIST')){
        stop('Data has not been loaded. Please run load_aba_data() first.')
    }

    stage <- stage_name(stage)
    if (is.null(DATA_LIST[[stage]])){
        DATA_LIST[[stage]] <<- read_loom(PATH_LIST[[stage]])
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
de <- function(expr_mat, groups){
    de_df <- presto::wilcoxauc(t(expr_mat), groups)
    de_df <- dplyr::select(de_df, 'gene'=feature, group, 'avg_exp'=avgExpr, 'fc'=logFC, auc,
                    pval, padj, 'prcex_self'=pct_in, 'prcex_other'=pct_out)
    return(tibble::as_tibble(de_df))
}

#' Finds differentially expressed features for all structures on a given annotation level
#'
#' @export
#'
variable_genes <- function(
    stage = 'E13pt5',
    method = 'vst',
    nfeatures = 2000
){

    if (!exists('DATA_LIST') | !exists('PATH_LIST')){
        stop('Data has not been loaded. Please run load_aba_data() first.')
    }

    stage <- stage_name(stage)
    if (is.null(DATA_LIST[[stage]])){
        DATA_LIST[[stage]] <<- read_loom(PATH_LIST[[stage]])
    }

    voxel_mat <- DATA_LIST[[stage]]$matrix

    var_feat <- Seurat::FindVariableFeatures(t(voxel_mat), verbose=F) %>%
        dplyr::as_tibble(rownames='gene') %>%
        dplyr::arrange(dplyr::desc(vst.variance.standardized)) %>%
        dplyr::filter(dplyr::row_number() <= nfeatures) %>%
        {colnames(.) <- stringr::str_replace_all(colnames(.), '\\.', '_'); .}

    return(var_feat)
}
