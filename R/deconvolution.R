#' @param expr_mat Expression matrix to deconvolute.
#' @param genes_use Marker genes to use for deconvolution.
#' @param reference Which reference to use.
#' @param stage Developmental stage to use as reference.
#' @param other_types Whether to infer proportions of type 'other'.
#'
#' @return A tibble with deconvoluted proportions for each reference structure.
#'
deconvolute.default <- function(
    expr_mat, genes_use,
    reference = 'ABA',
    stage = 'E13',
    annotation_level = 'custom_2',
    pseudo_tpm = FALSE,
    other_types = FALSE
){
    ref <- get_aba_ref(
        stage = stage,
        annotation_level = annotation_level,
        pseudo_tpm = pseudo_tpm
    )

    decon_df <- run_epic(
        expr_mat, ref$expr, genes_use,
        ref_sd = ref$sd,
        other_types = other_types
    )

    return(decon_df)
}

#' Prepare ABA ISH reference for deconvolution
#'
#' @param stage Developmental stage to use as reference.
#' @param annotation_level Annotation level to summarize to.
#'
#' @return A tibble with deconvoluted proportions for each reference structure.
#'
get_aba_ref <- function(
    stage = 'E13',
    annotation_level = 'custom_2',
    pseudo_tpm = FALSE
){
    if (!exists('DATA_LIST') | !exists('PATH_LIST')){
        stop('Data has not been loaded. Please run load_aba_data() first.')
    }

    stage <- stage_name(stage)
    if (is.null(DATA_LIST[[stage]])){
        DATA_LIST[[stage]] <<- readRDS(PATH_LIST[[stage]])
    }

    voxel_mat <- DATA_LIST[[stage]]$matrix
    voxel_annot <- DATA_LIST[[stage]]$row_meta[[annotation_level]]

    if (pseudo_tpm){
        voxel_mat <- (voxel_mat / Matrix::colSums(voxel_mat)) * 1e6
    }

    ref_expr <- aggregate_matrix(voxel_mat, groups=voxel_annot)
    ref_sd <- aggregate_matrix(as.matrix(voxel_mat), groups=voxel_annot, fun=matrixStats::colSds)
    rownames(ref_sd) <- rownames(ref_expr)

    ref <- list(
        expr = ref_expr,
        sd = ref_sd
    )
    return(ref)
}

#' Run EPIC deconvolution
#'
#' @param expr_mat Input matrix to deconvolute.
#' @param ref_expr Reference expression matrix.
#' @param genes_use Marker genes to use.
#' @param ref_sd Variability of reference expression.
#' @param other_types Whether to infer proportions of type 'other'.
#'
#' @return A tibble with deconvoluted proportions for each reference structure.
#'
run_epic <- function(
    expr_mat, ref_expr, genes_use,
    ref_sd = NULL,
    other_types = FALSE
){

    if (!is.null(ref_sd)){
        ref_sd <- as.matrix(ref_sd)
    }

    epic_decon <- suppressWarnings(
        EPIC::EPIC(
            bulk = as.matrix(expr_mat),
            reference = list(
                refProfiles = as.matrix(ref_expr),
                sigGenes = genes_use,
                refProfiles.var = ref_sd
            ),
            withOtherCells = other_types
        )
    )
    prop_df <- epic_decon$mRNAProportions %>%
        as_tibble(rownames='sample') %>%
        gather(struct, prop, -sample)

    return(prop_df)

}
