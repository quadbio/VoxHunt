#' Define the path to voxel data and load lazily if required
#'
#' @param dir The directory with loom files for each stage.
#' @param lazy Logical. If TRUE, the data is only loaded into memory when required,
#' if FALSE, the data is loaded right away.
#'
#' @export
#'
load_aba_data <- function(
    dir,
    lazy = T
){
    e11_path <- file.path(dir, 'E11.loom')
    e13_path <- file.path(dir, 'E13.loom')
    e15_path <- file.path(dir, 'E15.loom')
    e18_path <- file.path(dir, 'E18.loom')
    p4_path <- file.path(dir, 'P4.loom')
    p14_path <- file.path(dir, 'P14.loom')
    p28_path <- file.path(dir, 'P28.loom')
    p56_path <- file.path(dir, 'P56.loom')
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


#' Read loom files to list
#'
#' @import loomR
#'
#' @param filename The name of the loom file
#' @param layer The name of the layer to load as matrix
#' @param row_id A name for the colum in the row_meta indicating the row names.
#' @param col_id A name for the colum in the col_meta indicating the column names.
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
#'
#' @param groups A character vector with the groups to aggregate over.
#' @param fun The aggregation function to be applied to each chunk of the matrix.
#'
aggregate_matrix <- function(
    x,
    groups = NULL,
    fun = Matrix::colMeans
){
    if (length(groups) == nrow(x)){
        agg_mat <- sapply(levels(factor(groups)), function(g){
            fun(x[which(groups==g), ])
        })
        agg_mat <- Matrix::Matrix(agg_mat, sparse=T)
    } else if (length(groups) <= 1){
        agg_mat <- fun(x)
        agg_mat <- Matrix::Matrix(agg_mat, sparse=T)
        colnames(agg_mat) <- groups
        rownames(agg_mat) <- colnames(x)
    } else {
        stop('Length of groups must be either nrow(x) or 1.')
    }
    return(agg_mat)
}


#' Get gene expression from matrix and return data frame
#'
#' @import dplyr
#'
#' @param expr_mat A samples x genes gene expression matrix.
#' @param markers A character vector with the names of the genes to obtain.
#' @param scale Logicel. Whether to scale the expression values.
#'
get_markers <- function(
    expr_mat,
    markers,
    scale = T
){
    genes_use <- colnames(expr_mat)[colnames(expr_mat)%in%markers]
    expr_mat <- as.matrix(expr_mat[, genes_use])
    colnames(expr_mat) <- genes_use
    if (scale){
        expr_mat <- scale(expr_mat)
    }
    expr_mat <- as_tibble(expr_mat, rownames="voxel") %>%
        gather(gene, expr, -voxel)
    return(expr_mat)
}


#' Pick proper stage name
#'
stage_name <- function(stage){
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


#' Fast correlation and covariance calcualtion for sparse matrices
#'
#' @param x Sparse matrix or character vector
#' @param y Sparse matrix or character vector
#'
sparse_covcor <- function(x, y=NULL) {
    if (!is(x, "dgCMatrix")) stop("x should be a dgCMatrix")
    if (is.null(y)) {
        xtx <- crossprod(x)
        n <- nrow(x)
        cMeans <- colMeans(x)
        covmat <- (as.matrix(xtx) - n*tcrossprod(cMeans))/(n-1)
        sdvec <- sqrt(diag(covmat))
        cormat <- covmat/crossprod(t(sdvec))
        return(list(cov=covmat, cor=cormat))
    } else {
        if (!is(y, "dgCMatrix")) stop("y should be a dgCMatrix")
        if (nrow(x) != nrow(y)) stop("x and y should have the same number of rows")
        n <- nrow(x)
        cMeansX <- colMeans(x)
        cMeansY <- colMeans(y)
        covmat <- (as.matrix(crossprod(x, y)) - n * tcrossprod(cMeansX, cMeansY))/(n-1)
        sdvecX <- sqrt(diag((as.matrix(crossprod(x)) - n*tcrossprod(cMeansX))/(n-1)))
        sdvecY <- sqrt(diag((as.matrix(crossprod(y)) - n*tcrossprod(cMeansY))/(n-1)))
        cormat <- covmat / outer(sdvecX, sdvecY)
        return(list(cov=covmat, cor=cormat))
    }
}


#' Fast correlation for sparse matrices
#'
sparse_cor <- function(x, y=NULL) {
    sparse_covcor(x, y)$cor
}


#' Safe correlation function which returns a sparse matrix without missing values
#'
#' @param x Sparse matrix or character vector
#' @param y Sparse matrix or character vector
#' @param method Method to use for calculating the correlation coefficient.
#' @param allow_neg Logical. Whether to allow negative values or set them to 0.
#'
safe_cor <- function(
    x,
    y,
    method = 'pearson',
    allow_neg = F
){
    if (method == 'pearson'){
        x <- Matrix::Matrix(x, sparse = T)
        y <- Matrix::Matrix(y, sparse = T)
        corr_mat <- sparse_cor(x, y)
    } else {
        x <- as.matrix(x)
        y <- as.matrix(y)
        corr_mat <- stats::cor(x, y, method = method)
    }
    corr_mat[is.na(corr_mat)] <- 0
    corr_mat <- Matrix::Matrix(corr_mat, sparse=T)
    if (!allow_neg){
        corr_mat[corr_mat < 0] <- 0
    }
    return(corr_mat)
}
