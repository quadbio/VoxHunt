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
    lazy = TRUE
){
    e11_path <- file.path(dir, 'E11.rds')
    e13_path <- file.path(dir, 'E13.rds')
    e15_path <- file.path(dir, 'E15.rds')
    e18_path <- file.path(dir, 'E18.rds')
    p4_path <- file.path(dir, 'P4.rds')
    p14_path <- file.path(dir, 'P14.rds')
    p28_path <- file.path(dir, 'P28.rds')
    p56_path <- file.path(dir, 'P56.rds')
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
    is_file <- purrr::map_lgl(PATH_LIST, file.exists)
    if (any(!is_file)){
        warning('Warning: The provided directory is missing some expected files.')
    }
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
            E11 = readRDS(e11_path),
            E13 = readRDS(e13_path),
            E15 = readRDS(e15_path),
            E18 = readRDS(e18_path),
            P4 = readRDS(p4_path),
            P14 = readRDS(p14_path),
            P28 = readRDS(p28_path),
            P56 = readRDS(p56_path)
        )
    }
}

#' Load the mousebrain data (LaManno & Siletti et al. 2020)
#'
#' @param dir The directory with loom files for each stage.
#'
#' @export
#'
load_mousebrain_data <- function(
    dir
){
    all_path <- file.path(dir, 'dev_all.loom')
    all_loom <- hdf5r::H5File$new(all_path)
    agg_path <- file.path(dir, 'dev_all.agg.loom')
    agg_loom <- hdf5r::H5File$new(agg_path)

    all_tsne <- all_loom[['col_attrs/TSNE']][,]
    all_region <- all_loom[['col_attrs/Region']][]
    all_class <- all_loom[['col_attrs/Class']][]
    all_cluster <- all_loom[['col_attrs/Clusters']][]

    all_meta <- tibble::tibble(
        tSNE1 = all_tsne[1, ],
        tSNE2 = all_tsne[2, ],
        region = all_region,
        class = all_class,
        cluster = as.character(all_cluster)
    )

    utils::data(human2mouse, envir = environment())

    gene_names_mouse <- agg_loom[['row_attrs/Gene']][]
    gene_names_human <- human2mouse$human_symbol[match(gene_names_mouse, human2mouse$mouse_symbol)]
    gene_names <- gene_names_human[!is.na(gene_names_human)]
    gene_idx <- match(gene_names, gene_names_human)

    agg_expression <- agg_loom[['matrix']][, gene_idx]
    colnames(agg_expression) <- make.names(gene_names)
    rownames(agg_expression) <- agg_loom[['col_attrs/Clusters']][]

    MOUSEBRAIN_DATA <<- list(
        matrix = Matrix(agg_expression, sparse=T),
        meta = all_meta
    )

    all_loom$close_all()
    agg_loom$close_all()
}


#' Aggregate matrix over groups
#'
#' @import Matrix
#'
#' @param groups A character vector with the groups to aggregate over.
#' @param fun The aggregation function to be applied to each chunk of the matrix.
#'
#' @return A summary matrix.
#'
aggregate_matrix <- function(
    x,
    groups = NULL,
    fun = Matrix::colMeans
){
    if (length(groups) == nrow(x)){
        agg_mat <- sapply(levels(factor(groups)), function(g){
            chunk <- x[which(groups==g), ]
            if (is.null(dim(chunk))){
                return(chunk)
            } else {
                return(fun(chunk))
            }
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


#' Get gene expression values from matrix
#'
#' @import dplyr
#'
#' @param expr_mat A samples x genes gene expression matrix.
#' @param markers A character vector with the names of the genes to obtain.
#' @param scale Logical. Whether to scale the expression values.
#'
#' @return A tibble with expression of selected genes.
#'
get_markers <- function(
    expr_mat,
    markers,
    scale = TRUE
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
#' @param stage An alias for the developmental stage.
#'
#' @return A valid stage name.
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
#' @param x Sparse matrix or character vector.
#' @param y Sparse matrix or character vector.
#'
#' @return A list containing a covariance and correlation matrix.
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
#' @param x Sparse matrix or character vector.
#' @param y Sparse matrix or character vector.
#' @param method Method to use for calculating the correlation coefficient.
#' @param allow_neg Logical. Whether to allow negative values or set them to 0.
#'
#' @return A correlation matrix.
#'
safe_cor <- function(
    x,
    y,
    method = 'pearson',
    allow_neg = F
){
    if (method == 'pearson'){
        x <- Matrix::Matrix(x, sparse = TRUE)
        y <- Matrix::Matrix(y, sparse = TRUE)
        corr_mat <- sparse_cor(x, y)
    } else {
        x <- as.matrix(x)
        y <- as.matrix(y)
        corr_mat <- stats::cor(x, y, method = method)
    }
    corr_mat[is.na(corr_mat)] <- 0
    corr_mat <- Matrix::Matrix(corr_mat, sparse=TRUE)
    if (!allow_neg){
        corr_mat[corr_mat < 0] <- 0
    }
    return(corr_mat)
}

#' Fast correlation for sparse matrices
#'
sparse_cor <- function(x, y=NULL) {
    sparse_covcor(x, y)$cor
}


#' Scale function that works for vectors and grouped dataframes
#'
#' @param x A numeric vector.
#'
#' @return A scales numeric vector.
#'
zscale <- function(x){
    return((x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE))
}

#' Convert dataframe to matrix
#'
#' @param df A dataframe.
#'
#' @return A matrix.
#'
as_matrix <- function(df){
    mat <- as.matrix(df[2:ncol(df)])
    rownames(mat) <- df[[1]]
    return(mat)
}

#' Print diagnostic message
#'
#' @param text Text to print
#' @param verbose Whether to print
#'
log_message <- function(text, verbose=TRUE){
    if (verbose){
        message(text)
    }
}








