deconvolute.default <- function(
    expr_mat, marker_genes,
    reference = 'ABA',
    stage = 'E13',
    annotation_level = 'custom_2',
    other_types = FALSE
){
    ref <- get_aba_ref(
        stage = stage,
        annotation_level = annotation_level
    )

    decon_df <- run_epic(
        expr_mat, ref$expr, marker_genes,
        ref_sd = ref$sd,
        other_types = other_types
    )

    return(decon_df)

}


get_aba_ref <- function(
    stage = 'E13',
    annotation_level = 'custom_2'
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

    ref_expr <- aggregate_matrix(voxel_mat, groups=voxel_annot)
    ref_sd <- aggregate_matrix(as.matrix(voxel_mat), groups=voxel_annot, fun=colSds)
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
#' @param marker_genes Marker genes to use.
#' @param ref_sd Variability of reference expression.
#' @param other_types Whether to infer proportions of type 'other'.
#'
#' @return A tibble with deconvoluted proportions for each reference structure.
#'
run_epic <- function(
    expr_mat, ref_expr, marker_genes,
    ref_sd = NULL,
    other_types = FALSE
){

    if (!is.null(ref_sd)){
        ref_sd <- as.matrix(ref_sd)
    }

    epic_decon <- EPIC::EPIC(
        bulk = as.matrix(expr_mat),
        reference = list(
            refProfiles = as.matrix(ref_expr),
            sigGenes = marker_genes,
            refProfiles.var = ref_sd
        ),
        withOtherCells = other_types
    )

    prop_df <- epic_decon$mRNAProportions %>%
        as_tibble(rownames='sample') %>%
        gather(struct, prop, -sample)

    return(prop_df)

}
