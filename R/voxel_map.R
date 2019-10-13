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

    corr_mat <- safe_cor(expr_mat, voxel_mat)
    if (is.null(groups)){
        cell_meta <- tibble(
            cell = rownames(corr_mat)
        )
    } else {
        cell_meta <- tibble(
            cell = rownames(corr_mat),
            groups = groups
        )
    }
    utils::data(voxel_meta, envir = environment())
    voxel_meta = dplyr::filter(voxel_meta, sage==stage)

    vox_map <- list(
        corr_mat = corr_mat,
        cell_meta = cell_meta,
        voxel_meta = voxel_meta,
        genes = inter_genes
    )
    class(vox_map) <- 'VoxelMap'

    return(vox_map)
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


#' Print VoxelMap objects
#'
#' @rdname print
#' @export
#' @method print VoxelMap
#'
print.VoxelMap <- function(object){
    n_cells <- dim(object$corr_mat)[1]
    n_voxels <- dim(object$corr_mat)[2]
    n_genes <- length(object$genes)
    stage  <- unique(object$voxel_meta$sage)
    cat(paste0(
        'A VoxelMap object\n', n_cells, ' cells mapped to ',
        n_voxels, ' voxels in the ', stage, ' mouse brain\nbased on ',
        n_genes, ' features.\n'
    ))
}


