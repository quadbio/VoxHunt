#' Define the path to voxel data and load lazily if required
#'
#' @export
#'
load_aba_data <- function(
    dir,
    lazy = T
){
    e11_path <- file.path(dir, 'grid_data_E11pt5.loom')
    e13_path <- file.path(dir, 'grid_data_E13pt5.loom')
    e15_path <- file.path(dir, 'grid_data_E15pt5.loom')
    e18_path <- file.path(dir, 'grid_data_E18pt5.loom')
    p4_path <- file.path(dir, 'grid_data_P4.loom')
    p14_path <- file.path(dir, 'grid_data_P14.loom')
    p28_path <- file.path(dir, 'grid_data_P28.loom')
    p56_path <- file.path(dir, 'grid_data_P56.loom')
    PATH_LIST <<- list(
        E11 = e11_path,
        E13 = e13_path,
        E15 = e15_path,
        E18 = e18_path,
        P4 = p4_path,
        P14 = p14_path,
        P28 = p28_path,
        P56 = p56_path
    )
    DATA_LIST <<- list(
        E11 = NULL,
        E13 = NULL,
        E15 = NULL,
        E18 = NULL,
        P4 = NULL,
        P14 = NULL,
        P28 = NULL,
        P56 = NULL
    )
    if (!lazy){
        DATA_LIST <<- list(
            E11 = read_loom(e11_path),
            E13 = read_loom(e13_path),
            E15 = read_loom(e15_path),
            E18 = read_loom(e18_path),
            P4 = read_loom(p4_path),
            P14 = read_loom(p14_path),
            P28 = read_loom(p28_path),
            P56 = read_loom(p56_path)
        )
    }
}


#' Aggregate matrix over groups
#'
#' @import loomR
#'
read_loom <- function(
    filename,
    layer = 'matrix',
    row_id = 'voxel',
    col_id = 'gene'
){
    loom_file <- connect(filename)
    row_meta <- tibble::as_tibble(loom_file$get.attribute.df(2))
    colnames(row_meta)[1] <- row_id
    col_meta <- tibble::as_tibble(loom_file$get.attribute.df(1))
    colnames(col_meta)[1] <- col_id
    mat <- t(loom_file$get.sparse(layer))
    colnames(mat) <- col_meta[[col_id]]
    rownames(mat) <- row_meta[[row_id]]
    loom_file$close_all()

    out_list <- list(
        matrix = mat,
        row_meta = row_meta,
        col_meta = col_meta
    )

    return(out_list)
}


#' Aggregate matrix over groups
#'
#' @import Matrix
#' @import dplyr
#'
aggregate_matrix <- function(
    x,
    groups = NULL,
    fun = colMeans
){
    if (!is.null(groups)){
        agg_mat <- sapply(levels(factor(groups)), function(g){
            fun(x[which(groups==g), ])
        })
    } else {
        agg_mat <- as.matrix(fun(x))
        colnames(agg_mat) <- ' '
    }
    agg_mat <- Matrix::Matrix(agg_mat, sparse=T)
    return(agg_mat)
}


#' Get gene expression from matrix and return data frame
#'
#' @import dplyr
#'
get_markers <- function(
    expr_mat,
    markers,
    scale = T
){
    if (str_detect(class(expr_mat), 'atrix')){
        data <- expr_mat
    } else {
        data <- t(expr_mat@assays$RNA@data)
    }
    expr_mat <- as.matrix(data[, colnames(data)%in%markers])
    if (scale){
        expr_mat <- scale(expr_mat)
    }
    expr_mat <- as_tibble(expr_mat, rownames="cell") %>%
        gather(gene, expr, -cell)
    return(expr_mat)
}


#' Pick proper stage name
#'
stage_name <- function(
    stage
){
    stage_names <- list(
        E11 = c('E11pt5', 'E11.5', 'E11', 11, '11'),
        E13 = c('E13pt5', 'E13.5', 'E13', 13, '13'),
        E15 = c('E15pt5', 'E15.5', 'E15', 15, '15'),
        E18 = c('E18pt5', 'E18.5', 'E18', 18, '18'),
        P4 = c('P4', 4, '4'),
        P14 = c('P14', 14, '14'),
        P28 = c('P28', 28, '28'),
        P56 = c('P56', 56, '56')
    )
    for (n in names(stage_names)){
        if (stage %in% stage_names[[n]]){
            return(n)
        }
    }
}


#' Safe correlation function which returns a sparse matrix without missing values
#'
safe_cor <- function(
    x,
    y,
    method = 'pearson',
    allow_neg = F
){
    corr_mat <- stats::cor(x, y, method = method)
    corr_mat[is.na(corr_mat)] <- 0
    corr_mat <- Matrix::Matrix(corr_mat, sparse=T)
    if (!allow_neg){
        corr_mat[corr_mat < 0] <- 0
    }
    return(corr_mat)
}
