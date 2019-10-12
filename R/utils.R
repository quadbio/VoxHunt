#' Define the path to voxel data and load lazily if required
#'
#' @export
#'
load_aba_data <- function(
    directory,
    lazy = T
){
    assign('DATA_DIR', directory, envir = .GlobalEnv)
    if (lazy){
        return()
    } else {
        e11 <- read_loom(file.path(DATA_DIR, 'grid_data_E11pt5.loom'))
        assign('E11', e11, envir = .GlobalEnv)
        e13 <- read_loom(file.path(DATA_DIR, 'grid_data_E13pt5.loom'))
        assign('E13', e13, envir = .GlobalEnv)
        e15 <- read_loom(file.path(DATA_DIR, 'grid_data_E15pt5.loom'))
        assign('E15', e15, envir = .GlobalEnv)
        e18 <- read_loom(file.path(DATA_DIR, 'grid_data_E18pt5.loom'))
        assign('E18', e18, envir = .GlobalEnv)
        p4 <- read_loom(file.path(DATA_DIR, 'grid_data_P4.loom'))
        assign('P4', p4, envir = .GlobalEnv)
        p14 <- read_loom(file.path(DATA_DIR, 'grid_data_P14.loom'))
        assign('P14', p14, envir = .GlobalEnv)
        p28 <- read_loom(file.path(DATA_DIR, 'grid_data_P28.loom'))
        assign('P28', p28, envir = .GlobalEnv)
        p56 <- read_loom(file.path(DATA_DIR, 'grid_data_P56.loom'))
        assign('P56', p56, envir = .GlobalEnv)
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
    row_meta <- tibble::as_tibble(loom_file$get.attribute.df(2)) %>%
        dplyr::rename(row_id='CellID')
    col_meta <- tibble::as_tibble(loom_file$get.attribute.df(1)) %>%
        dplyr::rename(col_id='Gene')
    mat <- t(loom_file$get.sparse(layer))
    colnames(mat) <- col_meta[[row_id]]
    rownames(mat) <- row_meta[[col_id]]
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

# extract_legend <- function(x){
#     tmp <- ggplot::ggplot_gtable(ggplot::ggplot_build(x))
#     leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#     legend <- tmp$grobs[[leg]]
#     return(legend)}
