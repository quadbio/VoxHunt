#' Aggregate matrix over groups
#'
#' @import Matrix
#' @import dplyr
#'
aggregate_matrix <- function(x,
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
get_markers <- function(expr_mat, markers, scale=T){
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
